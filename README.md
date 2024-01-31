[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10581517.svg)](https://doi.org/10.5281/zenodo.10581517)

# infected_count

## Description

This script allows to segment cells using a [Cellpose](https://doi.org/10.1038/s41592-020-01018-x) environment, segment bacteria using an [Otsu threshold](https://doi.org/10.1109/TSMC.1979.4310076) and count the number of positive cells for infection.

## Requirements

The script is made in Python and requires [Fiji](https://doi.org/10.1038/nmeth.2019) to be ran. On top of this, multiple update sites need to be activated following [this guide](https://imagej.net/update-sites/#following-an-update-site): 
* PTBIOP

The script also read images from OMERO which needs to be configured. An [OMERO-insight](https://github.com/ome/omero-insight) Fiji JAR also needs to be copied in the OMERO plugins folder.

Once activated, just drag and drop the script in the main Fiji window and click on the RUN button.

As Fiji is operating system independant, this can be run on Windows, Mac and Linux. No specific hardware is necessary and script should take a couple of minutes to finish.

## Run script

### Input

The script will ask for multiple inputs :
* OMERO credentials
* Average size of the cells to find (in microns)
* Threshold for cells to be counted as positive
* Path to the Cellpose environment
* Path for the results to be saved

### Runtime 

Once started, the script will fetch the images from OMERO. It will segment the cells using the cytoplasm and nuclei channels using Cellpose. It will also segment the bacteria signal after background subtraction using an Otsu threshold.
The mean intensity of the binary image will be measured for each detected cells and an empirical threshold is used to consider the cells as being infected or not. 

### Output

The script will save the ROIs of each cells, classified as infected or non infected, allowing for easier verification. If selected at first, a segmented image of the bacteria as well as a tiff version of the original image will be saved in the same folder.
On top of this, a CSV will also be saved reporting the number of cells for each image and how many were infected.
