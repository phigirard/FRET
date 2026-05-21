#*******************************************************************************
#
#   Philippe GIRARD
#   Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
#   FRET_LSM_Timelapse.py
#   Release v14.0
#
#   Script for FRET analysis of spectral confocal LSM images in Fiji/ImageJ.
#   It loads spectral datasets or separate donor/acceptor files, performs
#   optional bleaching correction and background subtraction, and computes
#   pixel-wise FRET metrics (FRET index, A/D or D/A ratio).
#
#   Copyright 2026 - BSD-3-Clause license
#
#******************************************************************************/


# ---------------------------------------------------------------------------
# Fiji script parameters (injected via SciJava @ parameters)
# ---------------------------------------------------------------------------
#@ String msg1 (visibility=MESSAGE, value="----------------------------------------------------------- Input File(s): -----------------------------------------------------------", required=False)
#@ String fileType (label="File type: ",choices={ "Spectral Confocal LSM/CZI   ", "Separate TIF files Donor/Acceptor"}, style="radioButtonHorizontal", persist=True)
#@ String msg2 (visibility=MESSAGE, value="                                                                          ", required=False)
#@ String msg3 (visibility=MESSAGE, value="------------------------------------- Correction of the intensity decay due to photobleaching: ---------------------------------------", required=False)
#@ Boolean bleachCorr (label="Correction ?", description="Bleaching correction",value=False, persist=True)
#@ String CorrectionMethod (label="Correction method",choices={ "Simple Ratio", "Exponential Fit", "Histogram Matching" }, style="listBox", persist=True)
#@ String BleachCorrChoice (label="Only for Simple Ratio", choices={"Manual (Background values below)", "Manual (ROI selection)"}, style="radioButtonHorizontal")
#@ String msg4 (visibility=MESSAGE, value="                                                                         ", required=False)
#@ String msg5 (visibility=MESSAGE, value="------------------------------------------ Manual or Automatic background subtraction: --------------------------------------------", required=False)
#@ String ChoiceSub (label=" ", choices={"Manual (values below)", "Manual (ROI selection)","Automatic (ROI from threshold)", "Automatic (Rolling ball)"}, style="radioButtonHorizontal")
#@ Integer BGValueDonor (label="Background value (Donor):", value=100, persist=True)
#@ Integer BGValueAcceptor (label="Background value (Acceptor):", value=100, persist=True)
#@ Integer rollingBall (label="Rolling ball radius", description="Subtract Background",value=50, persist=True)
#@ String msg6 (visibility=MESSAGE, value="                                                                          ", required=False)
#@ String msg7 (visibility=MESSAGE, value="------------------------------------------------ Manual or Automatic threshold: ---------------------------------------------------", required=False)
#@ Boolean manualThreshold (label="Apply manual threshold value:", description="Manual threshold",value=False, persist=True)
#@ Integer thresholdValue (label="Threshold value:", value=100, persist=True)
#@ String msg8 (visibility=MESSAGE, value="                                                                          ", required=False)
#@ String msg9 (visibility=MESSAGE, value="------------------------------------------------------------- Other: --------------------------------------------------------------", required=False)
#@ Double timeLapse (label="Timelapse (min)", description="What is the timelapse between 2 images?",value=5, persist=False)
#@ String FRETchoice (label="FRET metric:",choices={"FRET index = 100 x A/(A+D)   ", "FRET ratio = A/D   ", "FRET ratio = D/A" }, style="radioButtonHorizontal", persist=True)
#@ Boolean calibrationBar (label="Display Calibration Bar ?", description="Calibration Bar",value=True, persist=True)
#@ String msg10 (visibility=MESSAGE, value="                                                                          ", required=False)

#@ UIService uiService
#@ LogService log
#@ PrefService prefs

# ---------------------------------------------------------------------------
# Librairies
# ---------------------------------------------------------------------------

# Python standard library
import sys
import os
import csv
import math

# Java AWT & charting
from java.io import File
from java.lang import Float
from java.awt import Font
from java.awt import Color
from java.awt.image import BufferedImage
from java.awt.geom import Rectangle2D
from java.awt.geom import Ellipse2D

from org.jfree.chart import ChartPanel
from org.jfree.chart import JFreeChart
from org.jfree.chart.axis import NumberAxis
from org.jfree.chart.plot import XYPlot
from org.jfree.chart.renderer.xy import XYLineAndShapeRenderer
from org.jfree.data.xy import XYDataset
from org.jfree.data.xy import XYSeries
from org.jfree.data.xy import XYSeriesCollection


# ImageJ/Fiji core classes
from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import Prefs
from ij import WindowManager
from ij.io import Opener
from ij.io import OpenDialog
from ij.process import ImageProcessor
from ij.process import ColorProcessor
from ij.process import ImageConverter
from ij.process import AutoThresholder
from ij.process import ImageStatistics
from ij.gui import GenericDialog
from ij.gui import WaitForUserDialog
from ij.gui import YesNoCancelDialog
from ij.gui import PlotWindow
from ij.gui import Roi
from ij.plugin import ZAxisProfiler
from ij.plugin import ZProjector
from ij.plugin import RGBStackMerge
from ij.plugin import RGBStackConverter
from ij.plugin import HyperStackReducer
from ij.plugin import Duplicator
from ij.plugin import RoiEnlarger
from ij.plugin import ImageCalculator as IC
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.measure import Calibration
from ij.plugin.filter import ThresholdToSelection
from ij.plugin.filter import Analyzer
from ij.plugin.frame import ThresholdAdjuster
from fiji.util.gui import GenericDialogPlus

# Bio-Formats / Loci
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.in import ImporterOptions
from loci.plugins import BF

# CMCI Library (Kota Miura) for Bleach Correction
from emblcmci import BleachCorrection_SimpleRatio
from emblcmci import BleachCorrection_ExpoFit
from emblcmci import BleachCorrection_MH


# ---------------------------------------------------------------------------
# Module constants and internal configuration
# ---------------------------------------------------------------------------

# ImageJ preferences
Prefs.blackBackground = True

# Bio-Formats import settings
SET_CONCAT = False
SHOW_OME_XML = False
AUTO_SCALE = False

# Channel selection behavior
SELECT_CHANNELS_INTERACTIVELY = True

# Display settings
DEFAULT_LUT = "Fire"

# Available photobleaching correction methods (must match UI choices)
CORRECTION_METHODS = (
	"Simple Ratio",
	"Exponential Fit",
	"Histogram Matching"
)

# Available FRET metrics (must match UI choices)
FRET_METRICS = (
	"FRET index = 100 x A/(A+D)   ",
	"FRET ratio = A/D   ",
	"FRET ratio = D/A"
)

# Available background subtraction methods (must match UI choices)
BACKGROUND_SUBTRACTION_METHODS = (
	"Manual (values below)",
	"Manual (ROI selection)",
	"Automatic (ROI from threshold)",
	"Automatic (Rolling ball)"
)

# CSV headers for background/threshold log file
CSV_FIELDNAMES = (
	"Frame #",
	"Subtraction method",
	"Donor Background",
	"Acceptor Background",
	"Min threshold",
	"Max Threshold"
)

# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

def log_info(message):
    """Log an informational message to both the Fiji Log window and stdout."""
    msg = "[INFO] " + str(message)
    IJ.log(msg)
    print(msg)


def log_step(message):
    """Log a major pipeline step with visual separation."""
    IJ.log("")
    IJ.log("=== " + str(message) + " ===")
    print("\n=== " + str(message) + " ===")


def log_warning(message):
    """Log a warning message."""
    msg = "[WARNING] " + str(message)
    IJ.log(msg)
    print(msg)


def log_error(message):
    """Log an error message."""
    msg = "[ERROR] " + str(message)
    IJ.log(msg)
    print(msg)


# ---------------------------------------------------------------------------
# Helper functions 
# ---------------------------------------------------------------------------

#### Fonctions for PART 1: data preparation & spectral channel selection

def getImpIndexes(imagefile_, sizeC_, sizeZ_, sizeT_, idxseries_, selectIdx_, imageDir_):
	"""Determine donor and acceptor channel indices from a spectral LSM image.

	A Z-projection is generated for the selected series, and the user is asked
	to draw a ROI containing the fluorescent signal. The mean intensity across
	all spectral channels is then plotted to display the emission spectrum.
	The user selects donor and acceptor channel indices based on this plot.

	Parameters
	----------
	imagefile_ : str
		Path to the spectral LSM image.
	sizeC_ : int
		Number of channels in the spectral stack.
	sizeZ_ : int
		Number of Z-slices.
	sizeT_ : int
		Number of time frames.
	idxseries_ : int
		Index of the series to analyze.
	selectIdx_ : bool
		If True, interactively select donor/acceptor channels.
	imageDir_ : str
		Directory used to save  outputs.

	Returns
	-------
	idxDonor_ : int
		Selected donor channel index.
	idxAcceptor_ : int
		Selected acceptor channel index.
	"""
	log_step("Channel selection from spectral LSM")
	log_info("Loading series %d for channel profile extraction" % idxseries_)

	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxseries_, True)
	options.setZBegin(idxseries_, 0)
	options.setZEnd(idxseries_, sizeZ_ - 1)
	options.setTBegin(idxseries_, int(sizeT_ / 2))
	options.setTEnd(idxseries_, int(sizeT_ / 2))
	impsLSM = BF.openImagePlus(options)

	IJ.run(impsLSM[0], "Grays", "stack")
	impProj_ = ZProjector.run(impsLSM[0], "max")
	IJ.run(impProj_, "Enhance Contrast", "saturated=0.35")

	idxDonor_ = 3
	idxAcceptor_ = 7

	if selectIdx_:
		IJ.setTool("rectangle")
		impProj_.show()
		IJ.run("Brightness/Contrast...")
		waitDialog = WaitForUserDialog("ROI", "Select a ROI with the rectangle tool")
		waitDialog.show()
		roi_ = impProj_.getRoi()
		
		if roi_ is None:
			IJ.run(impProj_, "Select All", "")
		else:
			impsLSM[0].setRoi(roi_)

		impsLSM[0].setCalibration(Calibration())
		impsLSM[0].setDimensions(1, 1, sizeC_)
		plot = ZAxisProfiler.getPlot(impsLSM[0])
		plot.setXYLabels("Channel", "Mean")
		xvalues = plot.getXValues()
		yvalues = plot.getYValues()

		series = XYSeries("Mean intensity channel")
		for i in range(len(xvalues)):
			series.add(xvalues[i], yvalues[i])

		dataset = XYSeriesCollection(series)
		yaxis = NumberAxis("Mean")
		xaxis = NumberAxis("# Channel")
		offset = 0.5
		xaxis.setRange(1 - offset, sizeC_ + offset)

		r = XYLineAndShapeRenderer()
		r.setSeriesPaint(0, Color.RED)
		r.setSeriesShape(0, Ellipse2D.Double(-3.0, -3.0, 6.0, 6.0))

		xyplot = XYPlot(dataset, xaxis, yaxis, r)
		xyplot.setBackgroundPaint(Color.white)
		chart = JFreeChart(xyplot)
		chart.removeLegend()

		impPlot = IJ.createImage("Mean intensity channel", "RGB", 512, 512, 1)
		imagePlot = impPlot.getBufferedImage()
		chart.draw(imagePlot.createGraphics(),
			Rectangle2D.Float(0, 0, impPlot.getWidth(), impPlot.getHeight()))
		impPlot.setImage(imagePlot)

		# Save spectrum image alongside analysis folder
		IJ.saveAs(impPlot, "TIFF", os.path.join(imageDir_, "channel_spectrum.tif"))
		impPlot.show()

		gui = GenericDialog("Select FRET donor/acceptor channels")
		gui.addSlider("Donor channel: ", 1, sizeC_, 3)
		gui.addSlider("Acceptor channel: ", 1, sizeC_, 7)
		gui.showDialog()

		impPlot.close()
		impProj_.killRoi()
		impProj_.close()

		idxDonor_ = int(gui.getNextNumber())
		idxAcceptor_ = int(gui.getNextNumber())
		log_info("Selected donor channel: %d" % idxDonor_)
		log_info("Selected acceptor channel: %d" % idxAcceptor_)

	for i in range(len(impsLSM)):
		impsLSM[i].flush()

	return idxDonor_, idxAcceptor_


def extractImpFromIndex(imagefile_, idxChannel_, idxSeries_):
	"""Extract a single-channel ZCT hyperstack from a spectral LSM file."""
	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxSeries_, True)
	options.setAutoscale(AUTO_SCALE)
	options.setShowOMEXML(SHOW_OME_XML)
	options.setConcatenate(SET_CONCAT)
	options.setCBegin(idxSeries_, idxChannel_ - 1)
	options.setCEnd(idxSeries_, idxChannel_ - 1)
	return BF.openImagePlus(options)[0]


def adjustSizeNum(S_, length_):
	"""Pad a numeric string with leading zeros until it reaches the given length."""
	newS_ = str(S_)
	while len(newS_) < length_:
		newS_ = "0" + newS_
	return newS_


def createFolder(imagePath_, basename_):
	"""Create an analysis folder next to the input image and return its path."""
	srcDir_ = os.path.dirname(imagePath_)
	anDir_ = os.path.join(srcDir_, basename_)
	if not os.path.exists(anDir_):
		os.makedirs(anDir_)
	return anDir_



#### Fonctions for PART 2: background selection, bleaching correction, etc.

def getBackgroundROI(imp_):
	"""Interactively select a ROI to measure background intensity."""
	imp_.show()
	IJ.setTool("freehand")
	myWait = WaitForUserDialog("Select a ROI",
		"Select a ROI for background subtraction")
	myWait.show()
	roi_ = imp_.getRoi()
	imp_.killRoi()
	imp_.hide()
	return roi_


def subtractBG(imp_, roi_):
	"""Subtract the mean gray intensity measured inside a ROI from an image."""
	imp_.setRoi(roi_)
	mean_ = imp_.getStatistics(Measurements.MEAN).mean
	IJ.run(imp_, "Select None", "")
	IJ.run(imp_, "Subtract...", "value=" + str(mean_) + " slice")
	return mean_


def applyROI2NAN(imp_, roi_):
	"""Fill a ROI with NaN values in a 32-bit image."""
	ip_ = imp_.getProcessor()
	ip_.setColor(float("nan"))
	ip_.fill(roi_)
	IJ.run(imp_, "Select None", "")
	return


def bleachCorrection(imp_, CorrectionMethodIdx_, backROI_):
	"""Apply photobleaching correction to a time-lapse stack."""
	impcorrected = imp_.duplicate()

	if CorrectionMethodIdx_ == 0:
		imp_.setRoi(backROI_)
		mean_ = imp_.getStatistics(Measurements.MEAN).mean
		imp_.killRoi()
		BCSR = BleachCorrection_SimpleRatio(impcorrected, mean_)
		BCSR.correctBleach()
	elif CorrectionMethodIdx_ == 1:
		BCEF = BleachCorrection_ExpoFit(impcorrected)
		BCEF.core()
		impfit = WindowManager.getImage("y = a*exp(-bx) + c")
		if impfit is not None:
			impfit.hide()
			impfit.close()
	elif CorrectionMethodIdx_ == 2:
		BCMH = BleachCorrection_MH(impcorrected)
		BCMH.doCorrection()

	return impcorrected


def applyThreshold(imp_, minthres_, maxthres_):
	"""Apply a raw threshold on an image and set background pixels to NaN."""
	IJ.setRawThreshold(imp_, minthres_, maxthres_, None)
	IJ.run(imp_, "NaN Background", "")
	IJ.run(imp_, "Despeckle", "")
	return

	
#### PART 3 :  FRET metric computation functions

def CalculationFRETmetric(impD_, impA_, FRETmetric_):
    """Compute the selected FRET metric from donor and acceptor stacks.

    Depending on the user choice, this function computes:
      - FRET index = 100 × A / (A + D)
      - FRET ratio A/D
      - FRET ratio D/A

    Zero-valued and saturated pixels are excluded to avoid infinities and
    extreme ratios.

    Parameters
    ----------
    impD_ : ImagePlus
        Donor image stack (32-bit).
    impA_ : ImagePlus
        Acceptor image stack (32-bit).
    FRETmetric_ : str
        Name of the FRET metric as selected in the UI.

    Returns
    -------
    ImagePlus
        FRET image stack corresponding to the chosen metric.
    """
    if FRETmetric_ == FRET_METRICS[0]:
        # 1) Compute (Donor + Acceptor) image (denominator)
        impD_ = IC.run(impD_, impA_, "Add create 32-bit stack")
        # Remove zero-valued pixels to avoid infinities
        IJ.setRawThreshold(impD_, 1, Float.MAX_VALUE, None)
        IJ.run(impD_, "NaN Background", "stack")
        # 2) Compute Acceptor / (Donor + Acceptor)
        impA_ = IC.run(impA_, impD_, "Divide create 32-bit stack")
        IJ.setRawThreshold(impA_, 0, 1, None)
        IJ.run(impA_, "NaN Background", "stack")
        IJ.run(impA_, "Multiply...", "value=100 stack")
        IJ.run(impA_, "Enhance Contrast", "saturated=0.35 stack")
    elif FRETmetric_ == FRET_METRICS[1]:
        # A/D ratio
        impA_ = CalculationFRETratio(impD_, impA_)
    else:
        # D/A ratio
        impA_ = CalculationFRETratio(impA_, impD_)
    return impA_


def CalculationFRETratio(imp1_, imp2_):
    """Compute a FRET ratio image as imp2_ / imp1_.

    Both images are thresholded to remove zero-valued pixels and extreme
    values, then the ratio is computed and contrast enhanced.

    Parameters
    ----------
    imp1_ : ImagePlus
        Denominator image (e.g., donor or acceptor).
    imp2_ : ImagePlus
        Numerator image (e.g., acceptor or donor).

    Returns
    -------
    ImagePlus
        Ratio image (32-bit) with invalid pixels set to NaN.
    """
    # Remove zero-valued pixels in denominator
    IJ.setRawThreshold(imp1_, 1, Float.MAX_VALUE, None)
    IJ.run(imp1_, "NaN Background", "stack")
    # Compute ratio imp2_ / imp1_
    impRatio_ = IC.run(imp2_, imp1_, "Divide create 32-bit stack")
    # Remove extreme ratio values created by division
    IJ.setRawThreshold(impRatio_, 0, Float.MAX_VALUE, None)
    IJ.run(impRatio_, "NaN Background", "stack")
    IJ.run(impRatio_, "Enhance Contrast", "saturated=0.35 stack")
    return impRatio_


def drawCalibrationBar(statsMin_, statsMax_):
    """Create a calibration bar image for FRET values.

    The bar is a 32-bit image showing the continuous range between the
    minimum and maximum FRET values, labeled at both ends.

    Parameters
    ----------
    statsMin_ : float
        Minimum value for the calibration bar.
    statsMax_ : float
        Maximum value for the calibration bar.

    Returns
    -------
    ImagePlus
        ImagePlus object containing the calibration bar.
    """
    impBar_ = IJ.createImage("Calibration Bar", "32-bit black", 276, 50, 1)
    ipBar_ = impBar_.getProcessor()
    step = (statsMax_ - statsMin_) / 255.0
    for i in range(256):
        for j in range(30):
            ipBar_.setf(11 + i, j, statsMin_ + i * step)
    ipBar_.setColor(256)
    ipBar_.setFont(Font("SansSerif", Font.BOLD, 12))
    ipBar_.drawString("{:.1f}".format(statsMin_), 0, 48)
    ipBar_.drawString("{:.1f}".format(statsMax_), 246, 48)
    impBar_.updateAndDraw()
    impBar_.setDisplayRange(statsMin_, statsMax_)
    return impBar_


# ---------------------------------------------------------------------------
# MAIN SCRIPT
# ---------------------------------------------------------------------------

#### PART 1 : Preparation of data & analysis - Preprocessing
log_step("PART 1 : Preparation of data & analysis - Preprocessing")

# Clear the Fiji Log window
IJ.log("\\\\Clear")


# Clear the Fiji console if not in headless mode
try:
	uiService.getDefaultUI().getConsolePane().clear()
except:
	pass
	

# Close Results Table if opened
if IJ.isResultsWindow():
	IJ.run("Clear Results", "")
	tw = ResultsTable().getResultsWindow()
	if tw is not None:
		tw.close()

## Input handling: spectral LSM vs. separate donor/acceptor images

if fileType == "Spectral Confocal LSM/CZI   " :
	
	od = OpenDialog("Select a spectral LSM image", None)
	lsmPath = od.getPath()
	if lsmPath is None:
		log_warning("No LSM file selected. Aborting.")
		sys.exit(0)	
	log_info("Input LSM file: %s" % lsmPath)
	
	# initialite the reader and get the OME metadata
	reader = ImageReader()
	omeMeta = MetadataTools.createOMEXMLMetadata()
	reader.setMetadataStore(omeMeta)
	reader.setId(lsmPath)
	sizeC = reader.getSizeC()
	sizeT = reader.getSizeT()
	sizeZ = reader.getSizeZ()
	seriesCount = reader.getSeriesCount()
	reader.close()
	log_info("Detected dimensions: C=%d, Z=%d, T=%d, series=%d" %
		(sizeC, sizeZ, sizeT, seriesCount))
	
	
	#select the serie if several
	idxSerie=0
	if seriesCount>1 :
		gui = GenericDialog("Select image serie")
		gui.addSlider("Image series: ", 1, seriesCount, 3)
		gui.showDialog()
		idxSerie = int(gui.getNextNumber()-1)
		log_info("Selected series index: %d" % idxSerie)
	
	#Create Folder for Donor and Acceptor images
	basename = os.path.basename(os.path.splitext(lsmPath)[0]).replace(' ', '_').lower()
	basename += "_S" + adjustSizeNum(str(idxSerie), 2)
	imageDir = createFolder(lsmPath , basename)

	#Open a menu for selecting Donor and Acceptor image indexes
	log_info("Select donor and acceptor channels...")
	idxDonor, idxAcceptor = getImpIndexes(lsmPath,sizeC, sizeZ, sizeT, idxSerie, SELECT_CHANNELS_INTERACTIVELY, imageDir)
 	

	#Extract Donor and Acceptor images 
	log_step("Extract donor and acceptor image stacks")
	impDonor = extractImpFromIndex(lsmPath, idxDonor, idxSerie)
	impAcceptor = extractImpFromIndex(lsmPath, idxAcceptor, idxSerie)

	
	IJ.run(impDonor, "Grays", "stack")
	IJ.run(impAcceptor, "Grays", "stack")
	#Save Donor and Acceptor raw images
	IJ.saveAs(impDonor, "TIFF",os.path.join(imageDir, basename+"_c1.tif")) 
	IJ.saveAs(impAcceptor, "TIFF",os.path.join(imageDir, basename+"_c2.tif")) 
	log_info("Saved raw donor and acceptor images.")

else : 
	log_step("Input from separate donor/acceptor TIF files")
	gui = GenericDialogPlus("Select the Donor/Acceptor images")
	gui.addFileField("Select the donor file", prefs.get(None, "dir.donor", "DefaultImage"))
	gui.addFileField("Select the acceptor file", prefs.get(None, "dir.acceptor", "DefaultImage"))
	
	gui.showDialog()
	if gui.wasOKed():
	    donorPath   = gui.getNextString()
	    acceptorPath = gui.getNextString()
	    prefs.put(None, "dir.donor", donorPath) 
	    prefs.put(None, "dir.acceptor", acceptorPath) 
	else :
		log_warning("User cancelled donor/acceptor file selection. Aborting.")
		sys.exit(0)
	
	
	#create folder for analysis 
	basename = os.path.basename(os.path.splitext(donorPath)[0]).replace(' ', '_').lower()
	imageDir = createFolder(donorPath , basename)

	#open the files 
	log_info("Donor file: %s" % donorPath)
	log_info("Acceptor file: %s" % acceptorPath)
	impDonor = Opener().openImage(donorPath)
	impAcceptor = Opener().openImage(acceptorPath)

## Calibration and stack properties

cal = impAcceptor.getCalibration()
unit = cal.getUnit()
if unit =="micron" :
	unit ="um"
pix2phys = cal.getX(1)	


# Get width, height, frame, bit depth 
nbSlice = impAcceptor.getStackSize()
width = impAcceptor.width
height = impAcceptor.height
depth = impAcceptor.getBitDepth()
log_info("Stack size: %d slices, width=%d, height=%d, bit-depth=%d"
	% (nbSlice, width, height, depth))

# Create  Image stack 
stackDonor = ImageStack(width, height)
stackAcceptor = stackDonor.duplicate()


IJ.run(impDonor, "Enhance Contrast", "saturated=0.35")	
IJ.run(impAcceptor, "Enhance Contrast", "saturated=0.35")	


#### PART 2 :  Bleaching correction and substract background 
log_step("PART 2 : Bleaching correction and background subtraction")

doBleachROI = True
doThreshold = True
# Create Array of Dictionnary for Background/Threshold  CSV File
infoImg = []
for slic in range(nbSlice): 
	if (nbSlice > 1) :
		log_info("Process image %d/%d" % (slic + 1, nbSlice))
		
	# Duplicate the frame number 'slic+1' and convert the image in 32-bit
	impDonor.setSlice(slic+1)
	impDonor_slice = impDonor.crop("whole-slice")
	ImageConverter(impDonor_slice).convertToGray32()
	impAcceptor.setSlice(slic+1)
	impAcceptor_slice = impAcceptor.crop("whole-slice")
	ImageConverter(impAcceptor_slice).convertToGray32()
	
	#Check the bit depth of the images and remove saturated and null pixels (saturated pixel are above  2^depth )
	if slic == 0 and depth > 8 :
		maxPix = impAcceptor_slice.getStatistics(Measurements.MIN_MAX).max
		if maxPix < 4096 :
			depth = 12 # the camera is 12-bits dynamical range =[0,4095]
	maxVal = math.pow(2,depth)-1
	
	applyThreshold(impDonor_slice, 1,  maxVal - 1)
	applyThreshold(impAcceptor_slice, 1,  maxVal - 1)

	
	#Background subtraction
	if (bleachCorr or ChoiceSub == BACKGROUND_SUBTRACTION_METHODS[1]) and doBleachROI :
		backROI = getBackgroundROI(impAcceptor_slice)
		if (nbSlice>1) :
			dial = YesNoCancelDialog(IJ.getInstance(), "Same ROI ?",
				"Do you want to use the same ROI for all the images",
				"  Yes  ", "  No  ")
			doBleachROI = not dial.yesPressed() 
	
	if bleachCorr:
		if nbSlice > 1:
			log_info("Correction of the photobleaching")
			CorrectionMethodIdx = CORRECTION_METHODS.index(CorrectionMethod)
			impDonor = bleachCorrection(impDonor_slice, CorrectionMethodIdx, backROI)
			impAcceptor = bleachCorrection(impAcceptor_slice, CorrectionMethodIdx, backROI)
		else:
			log_info("No photobleaching correction because raw data is not a stack")
	
	impT_slice = impAcceptor_slice.duplicate()
	if doThreshold and not manualThreshold:
		IJ.run(impT_slice, "Enhance Contrast", "saturated=0.35")	
		impT_slice.show()
		ta = ThresholdAdjuster()
		ta.show()
		ta.update()
		IJ.setAutoThreshold(impT_slice, "Moments dark")
		waitDialog = WaitForUserDialog("Manual threshold",
			"Please, adjust the threshold as desired, then press 'OK' (do not press 'Apply')")
		waitDialog.show()
		thres_min = impT_slice.getProcessor().getMinThreshold()
		thres_max = impT_slice.getProcessor().getMaxThreshold()
		ta.close()
		impT_slice.hide()
		if nbSlice > 1:
			d = YesNoCancelDialog(IJ.getInstance(), "Same threshold values ?",
				"Do you want to proceed automatically with the threshold values for all the images",
				"  Yes  ", "  No  ")
			doThreshold = not d.yesPressed()
		log_info("Threshold values: min = %.1f, max = %.1f" % (thres_min, thres_max))
		
	if slic == 0 and manualThreshold :
		thres_min = thresholdValue
		thres_max = maxVal
		
	IJ.setThreshold(impT_slice, thres_min, thres_max)
	roiCells = ThresholdToSelection.run(impT_slice) # Threshold selection in the image = signal
	roiNan = roiCells.getInverse(impT_slice) #get the inverse roi of  roiCells = Background
	impT_slice.close()
	
	if ChoiceSub == BACKGROUND_SUBTRACTION_METHODS[0]:
		IJ.run(impDonor_slice, "Subtract...", "value=" + str(BGValueDonor) + " slice")
		IJ.run(impAcceptor_slice, "Subtract...", "value=" + str(BGValueAcceptor) + " slice")
	elif ChoiceSub == BACKGROUND_SUBTRACTION_METHODS[1]:
		BGValueDonor = subtractBG(impDonor_slice, backROI)
		BGValueAcceptor = subtractBG(impAcceptor_slice, backROI)
		log_info("Background donor = %.1f, acceptor = %.1f" %
			(BGValueDonor, BGValueAcceptor))
	elif ChoiceSub == BACKGROUND_SUBTRACTION_METHODS[2]:
		BGValueDonor = subtractBG(impDonor_slice, roiNan)
		BGValueAcceptor = subtractBG(impAcceptor_slice, roiNan)
		log_info("Background donor = %.1f, acceptor = %.1f" %
			(BGValueDonor, BGValueAcceptor))
	elif ChoiceSub == BACKGROUND_SUBTRACTION_METHODS[3]:
		IJ.run(impDonor_slice, "Subtract Background...", "rolling=" + str(rollingBall) + " stack")
		IJ.run(impAcceptor_slice, "Subtract Background...", "rolling=" + str(rollingBall) + " stack")
		BGValueDonor = rollingBall
		BGValueAcceptor = rollingBall
		log_info("Background subtraction by rolling ball (radius = %d)" % rollingBall)

	# Convert roiNan in NaN values in all images
	applyROI2NAN(impDonor_slice,roiNan)
	applyROI2NAN(impAcceptor_slice,roiNan)

	# Create dictionnary values for the frame 'slic+1'
	BackThres = [str(slic + 1), ChoiceSub, BGValueDonor, BGValueAcceptor,
		thres_min, thres_max]
	infoImg.append(dict(zip(CSV_FIELDNAMES, BackThres)))

	stackDonor.addSlice(impDonor_slice.getProcessor())
	stackAcceptor.addSlice(impAcceptor_slice.getProcessor())

with open(os.path.join(imageDir, "infoFile.csv"), "wb") as csvfile:
	writer = csv.DictWriter(csvfile, fieldnames=CSV_FIELDNAMES)
	writer.writeheader()
	writer.writerows(infoImg)

log_info("Saved background/threshold log: infoFile.csv")

 
#save thresholded Donor image
impDonor_OUT=ImagePlus(impDonor.getTitle(), stackDonor)
IJ.run(impDonor_OUT, "Enhance Contrast", "saturated=0.35 stack")
IJ.saveAs(impDonor_OUT, "TIFF",os.path.join(imageDir, basename+"_c1thres.tif")) 

#save thresholded Acceptor image
impAcceptor_OUT=ImagePlus(impAcceptor.getTitle(), stackAcceptor)
IJ.run(impAcceptor_OUT, "Enhance Contrast", "saturated=0.35 stack")
IJ.saveAs(impAcceptor_OUT, "TIFF",os.path.join(imageDir, basename+"_c2thres.tif"))
	
#### PART 3 :  FRET metric images	
log_step("PART 3 : Measurement of " + FRETchoice)

metric = "ratio"
if FRETchoice == FRET_METRICS[0]:
	metric = "index"
elif FRETchoice == FRET_METRICS[1]:
	metric += "A_D"
else:
	metric += "D_A"
	

impFRET = CalculationFRETmetric(impDonor_OUT, impAcceptor_OUT, FRETchoice)
FRETTitle = "FRET_" + metric + "_" + os.path.basename(basename)
impFRET.setTitle(FRETTitle)

stats = impFRET.getStatistics(Measurements.MIN_MAX)
statsMax = stats.max
statsMin = stats.min
if FRETchoice == FRET_METRICS[0]:
	statsMax = math.ceil(stats.max)
	statsMin = math.floor(stats.min)

impFRET.setDisplayRange(statsMin, statsMax)
impFRET.setCalibration(cal)
IJ.run(impFRET, DEFAULT_LUT, "stack")
IJ.saveAs(impFRET, "TIFF", os.path.join(imageDir, FRETTitle))
impFRET.show()
log_info("Saved FRET stack: %s.tif" % FRETTitle)

Analyzer.setMeasurements(Measurements.AREA + Measurements.MEAN + Measurements.STD_DEV)
Analyzer.setPrecision(5)
rt = ResultsTable()
for slic in range(nbSlice):
	impFRET.setSlice(slic + 1)
	impFRET_slice = impFRET.crop("whole-slice")
	analyzer = Analyzer(impFRET_slice, rt)
	analyzer.measure()

rt.show("Mean FRET index (%)")
IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column")
rt.saveAs(os.path.join(imageDir, "MeanFRETindex.csv"))
log_info("Saved FRET measurements: MeanFRETindex.csv")

if calibrationBar:
	impBar = drawCalibrationBar(statsMin, statsMax)
	IJ.run(impBar, DEFAULT_LUT, "")
	impBar.show()
	IJ.saveAs(impBar, "TIFF", os.path.join(imageDir, "FRET_CalibrationBar"))
	log_info("Saved FRET calibration bar.")

log_step("End of analysis")
log_info("FRET_LSM_Timelapse script completed successfully.")
