#*******************************************************************************
#
#  Philippe GIRARD
#  Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
#  FRET_Wound_Healing.py
#  Release v3.0
#
#  Script for wound healing analysis using ROI meshing and FRET measurement.
#  It loads a FRET image and corresponding wound ROI, divides the ROI into
#  regular bands, measures FRET intensity in each band, and generates a
#  meshing visualization.
#
#  Copyright 2026 - BSD-3-Clause license
#
#******************************************************************************/

# ---------------------------------------------------------------------------
# Fiji / SciJava script parameters (UI fields)
# ---------------------------------------------------------------------------
#@ File roiFile (label="Select the ROI of wound healing:", style="file")
#@ File impFile (label="Select the FRET image:", style="file")
#@ Integer widthBand (label="Region width (in pixels):", value=6, persist=True)
#@ Integer heightHBand (label="Region height (in pixels):", value=6, persist=True)


#@ UIService uiService
#@ LogService log
#@ CommandService command


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------

# ImageJ core classes
from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import WindowManager
from ij import Prefs
from ij.io import Opener
from ij.io import RoiDecoder
from ij.gui import GenericDialog
from ij.gui import Roi
from ij.gui import ShapeRoi
from ij.gui import WaitForUserDialog
from ij.plugin import RoiEnlarger
from ij.plugin import ImageCalculator as IC
from ij.plugin.frame import RoiManager, ThresholdAdjuster
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.measure import Calibration
from ij.plugin.filter import Analyzer
from ij.plugin.filter import MaximumFinder
from ij.plugin.filter import ThresholdToSelection
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.process import ImageProcessor
from ij.process import ImageConverter
from ij.process import FloatPolygon
from ij.process import ImageStatistics


# Standard Python library
import os
import sys
import csv
import math


# Java utilities
from java.io import File
from java.lang import Double
from java.lang import Float
from java.awt import Color


# ---------------------------------------------------------------------------
# Module constants and configuration
# ---------------------------------------------------------------------------

Prefs.blackBackground = True

# Conversion / CSV output settings
IJ.run("Conversions...", "scale")
IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column save_row")


# Message prefixes
LOG_INFO  = "[INFO]  "
LOG_WARN  = "[WARNING] "
LOG_ERROR = "[ERROR] "

# Set upper area limit to infinity so that PA does not discard regions based on size.
MAXSIZE = Double.POSITIVE_INFINITY


# ---------------------------------------------------------------------------
# Logging helpers
# ---------------------------------------------------------------------------

def log_info(message):
	"""Log an informational message to the Fiji log window and stdout."""
	msg = LOG_INFO + str(message)
	IJ.log(msg)
	print(msg)


def log_warning(message):
	"""Log a warning message to the Fiji log window and stdout."""
	msg = LOG_WARN + str(message)
	IJ.log(msg)
	print(msg)


def log_error(message):
	"""Log an error message to the Fiji log window and stdout."""
	msg = LOG_ERROR + str(message)
	IJ.log(msg)
	print(msg)


def log_step(message):
	"""Log a major pipeline step with visual separation."""
	IJ.log("")
	IJ.log("=== " + str(message) + " ===")
	print("\n=== " + str(message) + " ===")


# ---------------------------------------------------------------------------
# All functions for the analysis
# ---------------------------------------------------------------------------

def changeValues(imp_, valIN_, ValOUT_):
	"""Replace all pixels equal to valIN_ in the image by ValOUT_."""
	ip_ = imp_.getProcessor()
	pixels_ = ip_.getPixels()
	for i in range(len(pixels_)):
		if pixels_[i] == valIN_:
			pixels_[i] = ValOUT_
	imp_.updateAndDraw()
	return


def changeValue2NAN(imp_, val_):
	"""Replace all pixels equal to val_ in the image by NaN."""
	changeValues(imp_, val_, Double.NaN)
	return


def measure(imp_, roi_, areaBand_):
	"""Measure the mean gray value inside the given ROI, but skip too‑small bands.

	This function is designed for band‑based FRET analysis to exclude small broken
	ROI segments that do not represent a reliable wound front region.
	"""
	imp_.setRoi(roi_)
	stats_ = imp_.getStatistics(Measurements.MEAN + Measurements.AREA)
	# If the actual ROI area is less than 1/4 of the nominal band area, treat it as invalid.
	if stats_.area < areaBand_ / 4:
		return Double.NaN
	else:
		return stats_.mean


def getRoiBandArray(roi_, nbHBand_, heightHBand_, width_, height_):
	"""Split the input ROI into horizontal bands of height = heightHBand_."""
	roiHBand_ = []
	shaperoi_ = ShapeRoi(roi_)
	for i in range(nbHBand_):
		bandRoi = ShapeRoi(Roi(0, i * heightHBand_, width_, heightHBand_))
		roiAND = shaperoi_.clone().and(bandRoi)
		roiHBand_.append(roiAND)
	return roiHBand_


# ---------------------------------------------------------------------------
# Start of main workflow
# ---------------------------------------------------------------------------

log_step("Wound Healing Analysis - Start")


# Close all images and reset the UI
IJ.run("Close All", "")

# Clear the Fiji Log window
IJ.log("\\\\Clear")

# Clear the Fiji console if not in headless mode
try:
	uiService.getDefaultUI().getConsolePane().clear()
except:
	pass

# Close existing Results window
if IJ.isResultsWindow():
	IJ.run("Clear Results", "")
	tw = ResultsTable().getResultsWindow()
	if tw is not None:
		tw.close()


# Optional global ResultsTable for logging
rt = ResultsTable()


# Set up RoiManager
rm = RoiManager.getInstance()
if rm is None:
	log_info("RoiManager not found; creating a new instance.")
	rm = RoiManager()
rm.reset()
log_info("RoiManager initialized.")


# Convert input File parameters to strings and paths
try:
	roiPath = roiFile.getCanonicalPath()
	impPath = impFile.getCanonicalPath()
	fretname = os.path.basename(os.path.splitext(impPath)[0])
	imageDir = os.path.dirname(impPath)
	log_info("Input ROI file: %s" % roiPath)
	log_info("Input FRET image: %s" % impPath)
	log_info("Output directory: %s" % imageDir)
except Exception as e:
	log_error("Failed to parse file paths: %s" % str(e))
	sys.exit(1)



# Open FRET image with error handling
impFRET = None
try:
	impFRET = Opener().openImage(impPath)
	if impFRET is None:
		raise Exception("OpenImage returned None")
	IJ.run(impFRET, "Grays", "")
	log_info("FRET image loaded: %s" % fretname)
except Exception as e:
	log_error("Could not open FRET image: %s" % str(e))
	sys.exit(1)


# Open wound ROI with error handling
roiCells = None
try:
	roiCells = RoiDecoder.open(roiPath)
	if roiCells is None:
		raise Exception("RoiDecoder.open returned None")
	log_info("Wound healing ROI loaded.")
except Exception as e:
	log_error("Could not open ROI file: %s" % str(e))
	sys.exit(1)

# Extract calibration and image metadata
cal = impFRET.getCalibration()
pix2phys = cal.getX(1)
unit = cal.getUnit()
if unit == "micron":
	unit = "µm"
#log_info("Pixel size: %.3f %s" % (pix2phys, unit))


width  = impFRET.width
height = impFRET.height
depth  = impFRET.getBitDepth()
log_info("Image size: width=%d px, height=%d px, bit-depth=%d" % (width, height, depth))


# Remove small regions inside the FRET image
IJ.run(impFRET, "Select None", "")
stats = impFRET.getStatistics(Measurements.MIN_MAX)
#statsMax = stats.max
statsMin = stats.min
impFRETMask = impFRET.duplicate()
IJ.setThreshold(impFRETMask, statsMin, 1000000000000000000000000000000.0000)
IJ.run(impFRETMask, "Convert to Mask", "")
stats = impFRETMask.getStatistics(Measurements.AREA)
options = PA.ADD_TO_MANAGER
p = PA(options, 0, None, 25, MAXSIZE)
p.analyze(impFRETMask)
rm.setSelectedIndexes(range(rm.getCount()))
rm.runCommand(impFRET,"Combine")
IJ.setBackgroundColor(0, 0, 0)
IJ.run(impFRET, "Clear Outside", "")
changeValue2NAN(impFRET, 0)
rm.reset()

# Number of horizontal bands (height‑wise)
nbHBand = int(height / heightHBand)
lastBand = 0 if height == nbHBand * heightHBand else 1
nbHBand += lastBand



# Create result image (meshing display)
impResult = IJ.createImage("Wound_Meshing", "32-bit black", width, height, 1)
changeValue2NAN(impResult, 0)
ipResult = impResult.getProcessor()

log_info("Meshing result image created.")


# Create horizontal bands of ROI
roiHBand = getRoiBandArray(roiCells, nbHBand, heightHBand, width, height)



# Compute number of bands per line (for later FRETvalue indexing)
nbRoiperLine = []
for iHBand, roiH in enumerate(roiHBand):
	widthROIH = roiH.getBounds().width
	roiXstart = widthROIH
	nbBand = int(widthROIH / widthBand)
	lastBand = 0 if widthROIH == nbBand * widthBand else 1
	nbBand += lastBand
	nbRoiperLine.append(nbBand)
maxband = max(nbRoiperLine)

log_info("Max band count per line: %d" % maxband)



# Mesh and measure FRET values in each band
log_step("Meshing and FRET measurement")

FRETvalue = []
reallengthWH = []
imp_slice = impFRET.duplicate()
# Nominal band area in physical units
bandArea = widthBand * heightHBand * pix2phys * pix2phys
for iHBand, roiH in enumerate(roiHBand):
	widthROIH = roiH.getBounds().width
	roiXstart = widthROIH
	nbBand = int(widthROIH / widthBand)
	lastBand = 0 if widthROIH == nbBand * widthBand else 1
	nbBand += lastBand
	reallengthWH.append(widthROIH)

	for iBand in range(maxband):
		roiXstart -= widthBand
		roiV = ShapeRoi(Roi(roiXstart, 0, widthBand, height))
		roi = roiH.clone().and(roiV)

		if roi.getLength() == 0 or roiXstart + widthBand < 0:
			FRETvalue.append(Double.NaN)
		else:
			rm.addRoi(roi)
			# Measure the mean value, but return NaN if the actual ROI area is too small
			meanVal = measure(imp_slice, roi, bandArea)
			FRETvalue.append(meanVal)
			# Color the corresponding pixels in the meshing image
			ipResult.setColor(meanVal)
			ipResult.fill(roi)


# Save band ROI set and visualize meshing
roisFile = os.path.join(imageDir, "RoiSet_Meshing.zip")
rm.runCommand("Save", roisFile)
log_info("Saved ROI band set: %s" % roisFile)


impResult.setCalibration(cal)
IJ.run(impResult, "Select None", "")
IJ.run(impResult, "Enhance Contrast", "saturated=0.35")
IJ.run(impResult, "Fire", "")

meshPath = os.path.join(imageDir, fretname + "_Meshing.tif")
IJ.saveAs(impResult, "TIFF", meshPath)
log_info("Saved meshing image: %s" % meshPath)
impResult.show()



# Create and show FRET channel table (band by band)
log_step("Building FRET channel table")

channelTable = ResultsTable()

for iX in range(maxband):
	valueX = '{:.2f}'.format(iX * widthBand * pix2phys)
	channelTable.setValue('Width from WH (microns)', iX, valueX)

	for iY in range(nbHBand):
		# Index is: iY * maxband + iX
		val = FRETvalue[iY * maxband + iX]
		channelTable.setValue('FRET-line_' + str(iY + 1), iX, val)

channelTable.show("Channel Analysis Results")

log_info("Analysis table shown.")


# Export measurements to CSV
csvPath = os.path.join(imageDir, "WH_Measurements.csv")
IJ.run("Input/Output...", "jpeg=85 gif=-1 file=.csv save_column")
channelTable.saveAs(csvPath)
log_info("Saved measurements CSV: %s" % csvPath)




# Final logging
log_step("Wound Healing Analysis - End")
log_info("Script completed successfully.")


# --------------------------------------------------------------------------
# If you want to log all computed FRET values in console
# (optional, only for debugging)
# --------------------------------------------------------------------------
VERBOSE = False

if VERBOSE:
	for iY in range(nbHBand):
		line = ""
		for iX in range(maxband):
			idx = iY * maxband + iX
			if math.isnan(FRETvalue[idx]):
				line += " nan   "
			else:
				line += " %5.2f " % FRETvalue[idx]
		log_info("FRET line %d: %s" % (iY, line))
