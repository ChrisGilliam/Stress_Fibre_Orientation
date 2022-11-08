% clears data stored in matlab:
clear; 
% clears command window
clc;

addpath([cd,filesep,'CircStat']);
addpath([cd,filesep,'export_fig']);

%%
%--------------------------------------------------------------------------
% Algorithm Parameters: 
% Note that these parameter names coinside with the parameters outlined in:
% W. Karlon et al, 'Measurement of Orientation and Distribution of Cellular 
% Alignment and Cytoskeletal Organization', Annals of Biomedical 
% Engineering, Vol. 27, pp. 712–720, 1999
%--------------------------------------------------------------------------
W = 20;     % Local window size (local window is W by W pixels in size)
s = 6;      % Size of 2D Gaussian filters (filters are s by s pixel in size)
sigma = 3;  % Sigma value of the 2D Gaussians (assuming symmetrical Gaussians)

p_local_threshold = 0.01;   % signifcance value for local Watson U2 test


%%
%--------------------------------------------------------------------------
% Load Image from file: 
%--------------------------------------------------------------------------
% Select image from folder:
[filename, pathname] = uigetfile({'*.*'},'Select image to process');

% Read image into matlab using pathname:
im = imread([pathname,filename]);
filename = filename(logical(cumsum(ismember(filename,'.'))==0));

% Create Folder for Results:
if exist(filename,'dir') == 0
    mkdir(filename);
end

% convert to double and use only one channel:

if size(im,3) > 1
    im1 = double(im(:,:,3));
    im1 = fliplr(im1.');
end
im = double(im(:,:,1));
im = fliplr(im.');

% determine size of image:
[M,N] = size(im);

% Determine x and y coordinates for the image
% Note current in pixels:
% [x,y] = meshgrid(1:N, 1:M);

%%
%--------------------------------------------------------------------------
% Determine Edges of cells: 
%--------------------------------------------------------------------------
[mask_edges,canny_threshold] = Determining_corticalFibers(im);
pos_canny = mask2position(mask_edges,[M,N]);


%%
%--------------------------------------------------------------------------
% Create Gaussian Filters: 
%--------------------------------------------------------------------------
disp('Calculating Orientation of Stress Fibres....');
% Create i and j coordinates for the filters:
[i,j] = meshgrid(-s:s,-s:s);
% Create filters:
hy = (2.*i./sigma^2).*exp( -(i.^2 + j.^2)./sigma^2);    % filter in x
hx = (2.*j./sigma^2).*exp( -(i.^2 + j.^2)./sigma^2);    % filter in y


%%
%--------------------------------------------------------------------------
% Calculate image gradients and determine per pixel estimate of orientation 
%--------------------------------------------------------------------------

% Compute gradients using filtering:
Gx = imfilter(im,hx,'symmetric');   % Gradient in x
Gy = imfilter(im,hy,'symmetric');   % Gradient in y

% calculate magnitude and angle for gradient assuming it is a 2D vector:
G = sqrt(abs(Gx).^2 + abs(Gy).^2);      % Magnitude
Phi = atan(Gy./Gx);                     % angle (in radians)
Phi1 = atan2(Gy,Gx);                     % angle (in radians)

%%
%--------------------------------------------------------------------------
% Determine dominant orientations in local window: 
%--------------------------------------------------------------------------

% determine a vector of possible orientations:
angle_vector = 0:179;

% Initialise 3D Matrix to hold data (efficient in matlab)
A = zeros(M,N,numel(angle_vector));

% for loop to calculate the individual images for each value of theta:
for index = 1:numel(angle_vector)
    
    % obtain local value of theta and conver to radians:
    theta = angle_vector(index)*pi/180;
    
    % calculate 
    holder = G.*exp(2.*cos(2.*(theta - Phi)))./exp(2);
    A(:,:,index) = local_sum(holder,W);
end

[~,i2] = max(A,[],3);

% Obtain dominant angle:
Dom_angle = angle_vector(i2);
disp('Finished Calculating Orientations of Stress Fibres');

%%
%--------------------------------------------------------------------------
% Performing local Watson U2 test: 
%--------------------------------------------------------------------------
disp('Starting local Watson U2 testing...');
h = waitbar(0,'Calculating Local Watson U2 test...');
p1 = zeros(size(Phi));
for i1 = 1:M
    waitbar(i1/M) ; 
    for i2 = 1:N
        index11 = max(i1-W./2,1);
        index12 = min(i1+W./2,M);
        index1 = index11:index12;
        index21 = max(i2-W./2,1);
        index22 = min(i2+W./2,N);
        index2 = index21:index22;
        
        Phi_holder = Phi(index1,index2);
        G_holder = G(index1,index2);

        Phi_holder1 = Phi1(index1,index2);
%         [p11(i1,i2) U UC] = circ_raotest(Phi_holder1(:));
        [p1(i1,i2),~]=WatsonU2Test(Phi_holder1(:));
        
    end
end
close(h);
disp('Finished local Watson U2 tests');
%%
%--------------------------------------------------------------------------
% Determining mask of valid pixels to analyse: 
%--------------------------------------------------------------------------

% calculate mask based on local Rao's Spacing test:
% mask_RS = 1 -> orientation in local region departs from uniform circular
% distribution with significance set by p_local_threshold;
mask_RS = logical(p1<=p_local_threshold);
        
% calculate mask based on comparing local mean to threshold:
% mask_int1 = 1 -> Mean intensity of local region greater than 0.4*mean of whole image
% mask_int2 = 1 -> Mean intensity of local region less than 2 standard
% derivations from mean of whole image.
max_I = max(im(:));
threshold1 = 0.4*mean(im(:)./max_I);
mask_int1 = logical(local_mean(im./max_I,W)>threshold1);

threshold2 = mean(im(:)./max_I)+std(im(:)./max_I).*2;
mask_int2 = logical(local_mean(im./max_I,W)<threshold2);

% Determine boundary pixels based on W:
mask_boundary = BoundaryMask(im,W./2);
Usable_mask = ~mask_edges.*mask_RS.*(mask_int2.*mask_int1).*~mask_boundary;

%%
%--------------------------------------------------------------------------
% Plotting and saving results: 
%--------------------------------------------------------------------------

% re-map dominant angle to -90 to 90:
Dom_angle(Dom_angle>90) = Dom_angle(Dom_angle>90) - 180;

% convert other angles to degrees:
Phi = Phi.*180/pi;

% Plotting raw results:
figure; imagesc(Phi); axis square; h = colorbar('fontsize',12);
caxis([-90,90]); title('Raw Angle Estimate');
h.Label.String = 'Orientation Angle (\circ)';
colormap hsv;
set(gcf, 'Color', 'w');
export_fig([cd,filesep,filename,filesep,'Raw_Orientation.tiff'],'-r300');

% Plotting dominant Angles:
figure; imagesc(Dom_angle); axis square; h = colorbar('fontsize',12);
caxis([-90,90]); title('Dominant Angle');
h.Label.String = 'Orientation Angle (\circ)';
colormap hsv;
set(gcf, 'Color', 'w');
export_fig([cd,filesep,filename,filesep,'Dominant_Orientation.tiff'],'-r300');

% Plotting Dominant angles on top of image:
plot_OrientationImage(im,Dom_angle,Usable_mask);
export_fig([cd,filesep,filename,filesep,'Dominant_Orientation_overlayImage.tiff'],'-r300');

% Plotting Polar histogram:
plot_PolarHist(Phi1(Usable_mask==1));
export_fig([cd,filesep,filename,filesep,'Polar_Histogram.tiff'],'-r300');

% Plotting histogram of Dominant angles:
plot_histogram(Dom_angle(Usable_mask==1));
export_fig([cd,filesep,filename,filesep,'Linear_Histogram.tiff'],'-r300');

%%
%--------------------------------------------------------------------------
% Calculate Circular Statistics:
%--------------------------------------------------------------------------
% Global Watson U2 test:
[p_WU2,U2]=WatsonU2Test(Phi1(Usable_mask==1));

% Circular mean:
[Circ_mean, upper_Circ_mean, lower_Circ_mean] = circ_mean(Dom_angle(Usable_mask==1)*pi./180);
Circ_mean = Circ_mean.*180./pi;
upper_Circ_mean = upper_Circ_mean.*180./pi;
lower_Circ_mean = lower_Circ_mean.*180./pi;

% Resultant Vector Length:
Circ_r = circ_r(Dom_angle(Usable_mask==1)*pi./180);

% Variance based on Resultant Vector length:
Circ_r_var = circ_var(Dom_angle(Usable_mask==1)*pi./180);

% Standard Deviation based on Resultant Vector length:
Circ_r_std = circ_std(Dom_angle(Usable_mask==1)*pi./180);

% % of orientated fibres in an interval:
factor = 100/sum(Usable_mask(:));
interval = 15:15:90;
Percent_interval = zeros(1,length(interval));
for index = 1:length(interval)
    if index == 1
        Percent_interval(1,index) = sum(logical(abs(Dom_angle(Usable_mask==1))<=interval(index))).*factor;
    else
        Percent_interval(1,index) = sum(logical(abs(Dom_angle(Usable_mask==1))<=interval(index) & abs(Dom_angle(Usable_mask==1))>interval(index-1)))./sum(Usable_mask(:))*100;  
    end
end

save_to_Excel(p_WU2,U2,Circ_mean,upper_Circ_mean,lower_Circ_mean,Circ_r,Circ_r_var,Circ_r_std,Percent_interval,filename);

clear('A','index','index1','index11','index12','index2','index21','index22','im1',...
    'Phi_holder','Phi_holder1','pos_canny','holder','G_holder','h','j','i','i1','i2');
save([cd,filesep,filename,filesep,'Results.mat']);

%% Embedded functions 
function J=local_sum(im,W)
%% Performs summation over local region
J=imfilter(im,ones(W,W),'symmetric');
end

function J=local_mean(im,W)
%% Performs mean over local region
J=imfilter(im,ones(W,W)./(W*W),'symmetric');
end

function pos = mask2position(mask,s_im)
%% Converts mask to a set of positions in x and y
wd=floor(s_im./2);
[X,Y]=meshgrid(-wd(2):s_im(2)-1-wd(2),-wd(1):s_im(1)-1-wd(1));
[y_cord,x_cord]=find(mask==1);
pos=[Y(y_cord+(x_cord-1)*s_im(1))+wd(1)+1,X(y_cord+(x_cord-1)*s_im(1))+wd(2)+1];
end

function mask=BoundaryMask(theta,W)
%% Determines boundary samples
[M,N]=size(theta);
m=(1:M)'*ones(1,N);
n=ones(M,1)*(1:N);
mask=logical(~(m>=W+1&m<=M-W &n>=W+1&n<=N-W));
end

function [p,U2]=WatsonU2Test(theta)
%% Performs Watson US test on the angles theta
% Note that theta is assumed to be in radians

% Determine number of angles
N = numel(theta);

% map theta so that it is between 0 to 2pi:
theta = mod( mod( theta(:), 2*pi ) + 2*pi , 2*pi );

% Order and map to 0 to 1:
U = sort(theta) ./ (2*pi);

% Calculate mean of U:
mu_U = mean(U);

% Rearranged version of the Watson U2 formula from Mardia & Jupp,
% Directional Statistics 2000
%     U2 = sum( U.^2 )- N.*mu_U.^2 - 2.*sum( (1:N).*U )./N + (N+1).*mu_U + N/12;
U2 = U.'*U - N.*mu_U.^2 - 2.*((1:N)*U)./N + (N+1).*mu_U + N/12;

%modified U2 according to Stephens, 1970
U2 = (U2 - 0.1./N + 0.1./(N.^2) ).*(1+0.8./N);

% Approximate p value - valid for large N
p = 2.*exp(-2.*U2.*(pi.^2));
end

function plot_OrientationImage(im,theta,Usable_mask)
%% Overlays the orientations on an image
h_hsv = hsv(180);
index = round((theta+90)/(180)*180);
RGB = ind2rgb(index,h_hsv);
mask1 = ones(size(im));
mask1(~Usable_mask==1) = nan;
for m = 1:3
    RGB(:,:,m) = RGB(:,:,m).*mask1;
end
thres = quantile(im(Usable_mask==1),0.95);
figure; image(RGB.*im./thres); caxis([-90,90]); h = colorbar('fontsize',12); 
colormap(h_hsv); 
axis square; axis off
h.Label.String = 'Orientation Angle (\circ)';
set(gcf, 'Color', 'w');
end

function plot_PolarHist(theta)
%% function plots a circular histogram of the data - requires the CircStat toolbox
figure;
obj1 = CircHist(theta./pi.*180, 90, 'areAxialData', true,'drawR',false);
% remove offset between bars and plot-center
rl = rlim;
obj1.setRLim([0, rl(2)]);
obj1.scaleBarSide = 'right';
obj1.polarAxs.ThetaZeroLocation = 'right';
obj1.setThetaLabel('Direction', 'bottomleft');
end

function plot_histogram(theta)
%% function plots a normal histogram of the data
figure; histogram(theta,30,'Normalization','probability');
xlabel('Orientation Angle (\circ)');
ylabel('Frequency Stress Fibre Orientation Angle (%)');
y_range = get(gca,'ytick').*100;
set(gca,'yticklabel', y_range);
line([-15-2.6, -15-2.6],[0,max(y_range)./100],'color','k','linewidth',1.5,'linestyle','--');
line([15+2.6, 15+2.6],[0,max(y_range)./100],'color','k','linewidth',1.5,'linestyle','--');
set(gcf, 'Color', 'w');
end

function [mask_edges,thres_array] = Determining_corticalFibers(im)
%% function is a semi-automated process to determine the cortical fibers of the cells
[M,N] = size(im);

[mask{2},t] = edge(im,'canny');
thres_array = [max(t).*0.5, (max(t)+1)./4, 0.5];
mask{1} = edge(im,'canny',thres_array(1));
mask{3} = edge(im,'canny',thres_array(3));
holder = {'Round 1:', 'Round 2:'};
for rot = 1:2
    for l = 1:3
        pos_canny = mask2position(mask{l},[M,N]);
        figure(l); imagesc(abs(im),[0,255]),axis image,colormap(gray)%
        hold on
        plot(pos_canny(:,2),pos_canny(:,1),'g*','MarkerSize',2,'LineWidth',1);
    end
    disp(holder(rot));
    disp('Please choose which pair of figures best represents the cortical fibers of the cells')
    disp('The options are:');
    disp('A: 1 and 2');
    disp('B: 2 and 3');
    loop_trig = 1;
    while loop_trig
        prompt = 'Please type either A or B: ';
        result = input(prompt, 's');

        switch result
            case {'A','a'}
                loop_trig = 0;
                thres_array = [thres_array(1), (thres_array(1)+thres_array(2))./2, thres_array(2)];
                mask{1} = edge(im,'canny',thres_array(1));
                mask{2} = edge(im,'canny',thres_array(2));
                mask{3} = edge(im,'canny',thres_array(3));

            case {'B','b'}
                loop_trig = 0;
                thres_array = [thres_array(2), (thres_array(2)+thres_array(3))./2, thres_array(3)];
                mask{1} = edge(im,'canny',thres_array(1));
                mask{2} = edge(im,'canny',thres_array(2));
                mask{3} = edge(im,'canny',thres_array(3));
            otherwise 
                disp('Please try again.');
        end
    end
end

for l = 1:3
    pos_canny = mask2position(mask{l},[M,N]);
    figure(l); imagesc(abs(im),[0,255]),axis image,colormap(gray)%
    hold on
    plot(pos_canny(:,2),pos_canny(:,1),'g*','MarkerSize',2,'LineWidth',1);
end

% last choice:
disp('Final Choice');
disp('Please choose which figure best represents the cortical fibers of the cells')
disp('The options are:');
disp('A: 1');
disp('B: 2');
disp('C: 3');
loop_trig = 1;
while loop_trig
    prompt = 'Please type either A, B or C: ';
    result = input(prompt, 's');

    switch result
        case {'A','a'}
            loop_trig = 0;
            thres_array = thres_array(1);
            mask_edges = mask{1};

        case {'B','b'}
            loop_trig = 0;
            thres_array = thres_array(2);
            mask_edges = mask{2};

        case {'C','c'}
            loop_trig = 0;
            thres_array = thres_array(3);
            mask_edges = mask{3};

        otherwise 
            disp('Please try again.');
    end
end

% dilate edges to remove influcence when analysing stress fibre
% orientation:
mask_edges = logical(imfilter(double(mask_edges),ones(5,5),'symmetric')>0);
close all;
end

function save_to_Excel(p_WU2,U2,Circ_mean,upper_Circ_mean,lower_Circ_mean,Circ_r,Circ_r_var,Circ_r_std,Percent_interval,filename)
% save data to excel file
Statistics = {'Circular Mean (degrees)', 'Upper Bound on Circular Mean (95% confidence interval, in degrees)',  'Lower Bound on Circular Mean (95% confidence interval, in degrees)',...
    'Resultant Vector Length', 'Variance of Resultant Vector Length', 'Standard Deviation of Resultant Vector Length',...
    'Percent of Orientated fibres in +/-15 degrees','Percent of Orientated fibres in +/- 15-30 degrees',...
    'Percent of Orientated fibres in +/- 30-45 degrees','Percent of Orientated fibres in +/- 45-60 degrees',...
    'Percent of Orientated fibres in +/- 60-75 degrees','Percent of Orientated fibres in +/- 75-90 degrees',...
    'p-value for Watson U2','U2 value'};
Values = [Circ_mean,upper_Circ_mean,lower_Circ_mean,Circ_r,Circ_r_var,Circ_r_std,Percent_interval,p_WU2,U2];

T = table(Statistics.',Values.');
writetable(table({'Statistics'},{'Values'}),[cd,filesep,filename,filesep,'Statistics.xlsx'],'Sheet',1,'Range','A1','WriteVariableNames',false);
writetable(T,[cd,filesep,filename,filesep,'Statistics.xlsx'],'Sheet',1,'Range','A2','WriteVariableNames',false);

end