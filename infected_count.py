# ─── SCRIPT PARAMETERS ──────────────────────────────────────────────────────────

#@ String(label="Username", description="please enter your username") USERNAME
#@ String(label="Password", description="please enter your password", style="password") PASSWORD
#@ String(label="Info about file(s)", description="Link got from OMERO, or image IDs separated by commas") OMERO_link
#@ Double(label="Approximate size of the cells", description="Average size of #the cells", value=40) cell_diam
#@ Double(label="Threshold for positive bacteria infection", value=0.2) bact_thresh
#@ File(label="Path to cellpose environment", style="directory", description="Path to conda env") cellpose_env
#@ File(label="Path for storage", style="directory", description="Where to store the results") destination
#@ Boolean(label="Save results ?", value=True) save_results

#@ CommandService command

# ─── REQUIREMENTS ───────────────────────────────────────────────────────────────

# List of update sites needed for the code
# * PTBIOP
# * OMERO jar
# * IMCF

# ─── Imports ──────────────────────────────────────────────────────────────────

import os
import sys
import csv
import time

from itertools import izip

from ij import IJ, Prefs
from ij.measure import ResultsTable
from ij.plugin import Duplicator, ImageCalculator
from ij.plugin.frame import RoiManager

# Bioformats imports
from loci.plugins import BF, LociExporter
from loci.plugins.in import ImporterOptions
from loci.plugins.out import Exporter

# Omero Dependencies
from omero.gateway import Gateway
from omero.gateway import LoginCredentials
from omero.log import SimpleLogger

from ch.epfl.biop.ij2command import Labels2Rois
from ch.epfl.biop.wrappers.cellpose.ij2commands import Cellpose_SegmentImgPlusAdvanced, CellposePrefsSet

# ─── Functions ────────────────────────────────────────────────────────────────

def parse_url(omero_str):
    """Parse an OMERO URL with one or multiple images selected

    Parameters
    ----------
    omero_str : str
        String which is either the link gotten from OMERO or image IDs separated by commas

    Returns
    -------
    str[]
        List of all the images IDs parsed from the string
    """
    if omero_str.startswith("https"):
        image_ids = omero_str.split("image-")
        image_ids.pop(0)
        image_ids = [s.split('%')[0].replace("|", "") for s in image_ids]
    else:
        image_ids = omero_str.split(",")
    return image_ids

def omero_connect():
    """Connect to OMERO using the credentials entered

    Returns
    -------
    gateway
        OMERO gateway
    """
    # Omero Connect with credentials and simpleLogger
    cred = LoginCredentials()
    cred.getServer().setHostname(HOST)
    cred.getServer().setPort(PORT)
    cred.getUser().setUsername(USERNAME.strip())
    cred.getUser().setPassword(PASSWORD.strip())
    simpleLogger = SimpleLogger()
    gateway = Gateway(simpleLogger)
    gateway.connect(cred)
    return gateway

def open_image_plus(HOST, USERNAME, PASSWORD, groupId, imageId):
    """Open an ImagePlus from an OMERO server

    Parameters
    ----------
    HOST : str
        Adress of your OMERO server
    USERNAME : str
        Username to use in OMERO
    PASSWORD : str
        Password
    groupId : double
        OMERO group ID
    imageId : int
        ID of the image to open
    """
    stackview = "viewhyperstack=true stackorder=XYCZT "
    datasetorg = "groupfiles=false swapdimensions=false openallseries=false concatenate=false stitchtiles=false"
    coloropt = "colormode=Default autoscale=true"
    metadataview = "showmetadata=false showomexml=false showrois=true setroismode=roimanager"
    memorymanage = "virtual=false specifyranges=false setcrop=false"
    split = " splitchannels=false splitfocalplanes=false splittimepoints=false"
    other = "windowless=true"
    options = ("location=[OMERO] open=[omero:server=%s\nuser=%s\npass=%s\ngroupID=%s\niid=%s] %s %s %s %s %s %s %s " %
               (HOST, USERNAME, PASSWORD, groupId, imageId, stackview, datasetorg, coloropt, metadataview, memorymanage, split, other))
    IJ.runPlugIn("loci.plugins.LociImporter", options)

def progress_bar(progress, total, line_number, prefix=''):
    """Progress bar for the IJ log window

    Parameters
    ----------
    progress : int
        Current step of the loop
    total : int
        Total number of steps for the loop
    line_number : int
        Number of the line to be updated
    prefix : str, optional
        Text to use before the progress bar, by default ''
    """

    size = 30
    x    = int(size*progress/total)
    IJ.log("\\Update%i:%s\t[%s%s] %i/%i\r" % (line_number, prefix, "#"*x, "."*(size-x), progress, total))

def timed_log(message):
    """Print a log message with a timestamp added

    Parameters
    ----------
    message : str
        Message to print
    """
    IJ.log(time.strftime("%H:%M:%S", time.localtime()) + ": " + message)

def BFExport(implus, savepath):
    """Export using BioFormats

    Parameters
    ----------
    implus : ImagePlus
        ImagePlus of the file to save
    savepath : str
        Path where to save the image

    """
    paramstring = "outfile=[" + savepath + "] windowless=true compression=Uncompressed saveROI=false"


    print('Savepath: ', savepath)
    plugin     = LociExporter()
    plugin.arg = paramstring
    exporter   = Exporter(plugin, implus)
    exporter.run()

# ─── Variables ────────────────────────────────────────────────────────────────

# OMERO server info
HOST    = "xxxx.xxxx.xxxx.xx"
PORT    = 4064
groupId = "-1"

# Channel order
nuclei_chnl   = 1
bacteria_chnl = 2
membrane_chnl = 3

# ─── Main Code ────────────────────────────────────────────────────────────────

IJ.log("\\Clear")
timed_log("Script starting")

command.run(CellposePrefsSet, False,
            "envType", "conda",
            "cellposeEnvDirectory", cellpose_env,
            "version", "2.0")

gateway = omero_connect()

image_ids_array = parse_url(str(OMERO_link))
image_ids_array.sort()

filenames_list = []
total_cell_list = []
infected_cell_list_otsu = []
# infected_cell_list_yen = []

rm_cells = RoiManager(False)

out_csv = os.path.join(destination.getPath(), "Results_infected.csv")

IJ.run("Set Measurements...", "mean")

try:

    for image_index, image_id in enumerate(image_ids_array):

        progress_bar(image_index + 1, len(image_ids_array), 2, "Processing : ")

        rt_stats = ResultsTable.getResultsTable("Results")
        if rt_stats:
            rt_stats.reset()
        rm_cells.reset()

        IJ.log("\\Update3:Loading the image from OMERO")
        open_image_plus(HOST, USERNAME, PASSWORD, groupId, image_id)

        imp = IJ.getImage()
        imp.hide()
        cal = imp.getCalibration()
        basename = imp.getTitle().replace(" ","_")

        filenames_list.append(imp.getTitle())

        pxl_cell_diameter = int(round(cell_diam / cal.pixelWidth))

        IJ.log("\\Update3:Running Cellpose")
        res = command.run(Cellpose_SegmentImgPlusAdvanced, False,
            "imp", imp,
            "diameter", pxl_cell_diameter,
            "model", "cyto2",
            "nuclei_channel", nuclei_chnl,
            "cyto_channel", membrane_chnl,
            "dimensionMode", "2D"
        ).get()
        imp_cells = res.getOutput("cellpose_imp")
        imp_cells.setTitle("cell_labels")
        imp_cells.setCalibration(cal)
        command.run(Labels2Rois, True, "rm", rm_cells, "imp", imp_cells).get()

        imp_bact = Duplicator().run(
            imp,
            bacteria_chnl, bacteria_chnl,
            1, imp.getNSlices(),
            1, 1
        )
        imp_bact_bg = imp_bact.duplicate()
        IJ.run(imp_bact_bg, "Gaussian Blur...", "sigma=20")
        imp_bact_otsu = ImageCalculator.run(imp_bact, imp_bact_bg, "Subtract")
        IJ.run(imp_bact_otsu, "Median...", "radius=2")

        # imp_bact_yen = imp_bact_otsu.duplicate()
        IJ.setAutoThreshold(imp_bact_otsu, "Otsu dark no-reset")
        Prefs.blackBackground = True
        IJ.run(imp_bact_otsu, "Convert to Mask", "")

        # IJ.setAutoThreshold(imp_bact_yen, "Yen dark no-reset")
        # Prefs.blackBackground = True
        # IJ.run(imp_bact_yen, "Convert to Mask", "")

        rm_cells.runCommand(imp_bact_otsu, "Measure")
        rt_stats = ResultsTable.getResultsTable()
        mean_int_bact_otsu = rt_stats.getColumn("Mean")
        rt_stats.reset()

        # rm_cells.runCommand(imp_bact_yen, "Measure")
        # rt_stats = ResultsTable.getResultsTable()
        # mean_int_bact_yen = rt_stats.getColumn("Mean")
        # rt_stats.reset()

        total_cell_list.append(rm_cells.getCount())

        count = sum(i > bact_thresh for i in mean_int_bact_otsu)
        infected_cell_list_otsu.append(count)

        for index, mean_value in enumerate(mean_int_bact_otsu):
            if mean_value > bact_thresh:
                rm_cells.getRoi(index).setGroup(1)
            else:
                rm_cells.getRoi(index).setGroup(2)

        rm_cells.getRoi(index).setGroupName(1, "Pair")
        rm_cells.getRoi(index).setGroupName(2, "Impair")


        # count = sum(i > bact_thresh for i in mean_int_bact_yen)
        # infected_cell_list_yen.append(count)

        if save_results:
            rm_cells.runCommand("Save", os.path.join(destination.getPath(), basename + "_ROIs.zip"))
            BFExport(imp_bact_otsu, os.path.join(destination.getPath(), basename + "_bact_otsu.tif"))
            # BFExport(imp_bact_yen, os.path.join(destination.getPath(), basename + "_bact_yen.tif"))
            BFExport(imp, os.path.join(destination.getPath(), basename + "_original_image.tif"))

    with open(out_csv, 'wb') as f:
        writer = csv.writer(f, delimiter=";")
        writer.writerow(
            ["Filename","Total number of cells", "Number of infected cells Otsu"]
        )
        writer.writerows(izip(filenames_list, total_cell_list, infected_cell_list_otsu))

    imp.close()
    imp_bact_otsu.close()
    imp_cells.close()
    # imp_bact_yen.close()

finally:
    gateway.disconnect()

timed_log("FINISHED")
