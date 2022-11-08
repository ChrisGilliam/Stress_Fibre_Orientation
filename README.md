# Stress_Fibre_Orientation
Image processing code used to generate the results in "[A microfluidic system for studying the effect of disturbed flow on endothelial cells](https://www.frontiersin.org/articles/10.3389/fbioe.2019.00081/full)", Frontiers in Bioengineering and Biotechnology, Vol. 7, apr 2019.
Authors: F. Tovar-Lopez, P. Thurgood,  C. Gilliam, N. Nguyen, E. Pirogova, K. Khoshmanesh, and S. Baratchi
Code written by C. Gilliam (c.gilliam.1@bham.ac.uk)

The code estimates and quantifies the orientation of actin stress fibres within a cell. The code is based on the algorithm outlined in W. Karlon et al, 'Measurement of Orientation and Distribution of Cellular Alignment and Cytoskeletal Organization', Annals of Biomedical Engineering, Vol. 27, pp. 712–720, 1999. Note that the code identifies the cortical fibres of the cell to remove them from the processing.

## Code Operations ##
The code is run using the script Main_File.m. The script operates as follows:
1) Asks the user to select an image for processing. A reference image "laminar shear-2.png" is included.
2) Runs a semi-automated process to determine the cortical fibers of the cells in the image. This process involves estimating the edges of the cell using a Canny filter with different parameters and asking the user to select the best results (i.e. user tunes the parameters of the Canny filter).
3) Computes a per pixel estimate of the dominant orientation angle of the stress fibres based on the algorithm outlined in W. Karlon et al, 'Measurement of Orientation and Distribution of Cellular Alignment and Cytoskeletal Organization', Annals of Biomedical Engineering, Vol. 27, pp. 712–720, 1999.
4) Performs a local Watson U2 on the dominant orientation angles. The code uses the Watson U2 formula outlined in Mardia & Jupp, Directional Statistics 2000, with the modification outlined according to Stephens, 1970.
5) Determines mask of valid pixels to be analyzed. The mask removes pixels that have been identified as cortical fibres; pixels where the orientation angles that fail the Watson U2 test; edge pixels; and pixels that have been saturated.
6) Plots the results and saves the graphs and accompanying statistics.

## Parameters ##
Algorithm parameters:

| Symbol in Paper | Value | Comment |
|:---:|---|---|
| _W_ | 20 | Local window size (local window is W by W pixels in size) |
| _s_ | 6 | Size of 2D Gaussian filters (filters are s by s pixel in size)|
| _sigma_ | 3 | Sigma value of the 2D Gaussians (assuming symmetrical Gaussians)|
| p_local_threshold | 0.01 | significance value for local Watson U2 test |

Note that these parameter names coincide with the parameters outlined in:
W. Karlon et al, 'Measurement of Orientation and Distribution of Cellular Alignment and Cytoskeletal Organization', Annals of Biomedical Engineering, Vol. 27, pp. 712–720, 1999

## Outputs ##
The code creates a new folder with the same name as the image being processed. In this folder the code saves 5 graphs
1) Dominant_Orientation.tiff - Estimated dominant orientation angle produced by the code
2) Dominant_Orientation_overlayImage.tiff - Estimated dominant orientation angle overlayed on the input image. The valid pixel mask is applied to the dominant orientation angle so that only valid dominant orientation angles are shown.
3) Linear_Histogram.tiff - Histogram of the frequency of the orientation angle of stress fibres.
4) Polar_Histogram.tiff - Polar histogram of the frequency of the orientation angle of stress fibres.
5) Raw_Orientation.tiff - Raw orientation angle estimated initial before computing the dominant orientation angle.

The code also outputs the following statistics in the excel file "Statistics.xlsx":
1) Circular Mean (degrees)
2) Upper Bound on Circular Mean (95% confidence interval, in degrees)
3) Lower Bound on Circular Mean (95% confidence interval, in degrees)
4) Resultant Vector Length
5) Variance of Resultant Vector Length
6) Standard Deviation of Resultant Vector Length
7) Percent of Orientated fibres in +/-15 degrees
8) Percent of Orientated fibres in +/- 15-30 degrees
9) Percent of Orientated fibres in +/- 30-45 degrees
10) Percent of Orientated fibres in +/- 45-60 degrees
11) Percent of Orientated fibres in +/- 60-75 degrees
12) Percent of Orientated fibres in +/- 75-90 degrees
13) p-value for Watson U2
14) U2 value

Finally the code saves the Matlab results in "Results.mat"

## Dependencies ##
The circular statistics are computed using the "CircStat" Matlab toolbox. The toolbox is included wtih this code but should be referenced as follows:
P. Berens, CircStat: A MATLAB Toolbox for Circular Statistics. 2009 31 (2009) 21

All of the graphs and images are exported from Matlab using the "export_fig" Matlab toolbox.
