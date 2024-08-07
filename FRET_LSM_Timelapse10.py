#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	FRET_LSM_Timelapse.py
#	Release v10.0
#
#	Copyright 2022 - BSD-3-Clause license
#                                                                             
#******************************************************************************/

#@ String msg1 (visibility=MESSAGE, value="--------------------------------------- Input File(s): ---------------------------------------", required=False) 
#@ String fileType (label="File type",choices={ "Spectral Confocal LSM   ", "Separate TIF files Donor/Acceptor"}, style="radioButtonHorizontal", persist=True) 
#@ String msg2 (visibility=MESSAGE, value="                                                                          ", required=False) 
#@ String msg3 (visibility=MESSAGE, value="--------------------- Correction of the intensity decay due to photobleaching: ----------------------", required=False) 
#@ Boolean bleachCorr (label="Correction ?", description="Bleaching correction",value=False, persist=True) 
#@ String CorrectionMethod (label="Correction method",choices={ "Simple Ratio", "Exponential Fit", "Histogram Matching" }, style="listBox", persist=True) 
#@ String msg4 (visibility=MESSAGE, value="                                                                         ", required=False) 
#@ String msg5 (visibility=MESSAGE, value="------------------------- Manual or Automatic background subtraction: --------------------------", required=False) 
#@ Boolean manualSub (label="Manual subtraction ?", description="Manual background subtraction",value=False, persist=True) 
#@ Integer rollingBall (label="Rolling ball radius", description="Subtract Background",value=50, persist=True)  
#@ String msg6 (visibility=MESSAGE, value="                                                                          ", required=False) 
#@ String msg7 (visibility=MESSAGE, value="----------------------------------------- Other: -----------------------------------------", required=False) 
#@ Double timeLapse (label="Timelapse (min)", description="What is the timelapse between 2 images?",value=5, persist=False)  
#@ String FRETchoice (label="FRET metric:",choices={"FRET index = 100 x A/(A+D)   ", "FRET ratio = A/D   ", "FRET ratio = D/A" }, style="radioButtonHorizontal", persist=True) 
#@ Boolean calibrationBar (label="Display Calibration Bar ?", description="Calibration Bar",value=True, persist=True) 
#@ String msg8 (visibility=MESSAGE, value="                                                                          ", required=False) 


#@ UIService uiService
#@ LogService log
#@ PrefService prefs 


#---------------------------------------------------------------
#---------------          LIBRARIES            -----------------
#---------------------------------------------------------------

# Python Library -------------------------------------------------------------------------------------------
import sys
import os, csv
import math 
#from random import *

from java.lang import Float
from java.awt import Font

# ImageJ Library ----------------------------------------------------------------------------
from ij import IJ, ImagePlus, Prefs
from ij.io import Opener, OpenDialog
from ij.process import ImageProcessor, ColorProcessor, ImageConverter, AutoThresholder, ImageStatistics
from ij.gui import GenericDialog, WaitForUserDialog, PlotWindow, ProfilePlot, Overlay, Line
from ij.gui import Roi
from ij.plugin import ZAxisProfiler, ZProjector, RGBStackMerge, RGBStackConverter, HyperStackReducer, Duplicator
from ij.plugin import ImageCalculator as IC
from ij.measure import ResultsTable , Measurements, CurveFitter, Calibration
from ij.plugin.filter import Analyzer, MaximumFinder
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin.frame import RoiManager,ThresholdAdjuster

from fiji.util.gui import GenericDialogPlus

# Loci Library ------------------------------------------------------------------------------------------
from loci.formats import ImageReader
from loci.formats import MetadataTools
from loci.plugins.in import ImporterOptions
from loci.plugins import BF


# CMCI Library (Kota Miura) for Bleach Correction --------------------------------------------------------
from emblcmci import BleachCorrection_SimpleRatio,BleachCorrection_ExpoFit,BleachCorrection_MH




#---------------------------------------------------------------
#---------------         CONSTANTS            -----------------
#---------------------------------------------------------------
title = "FRET Analysis for timelapse LSM images"

CorrectionMethods = list([ "Simple Ratio", "Exponential Fit", "Histogram Matching" ])
FRETmetrics = list([ "FRET index = 100 x A/(A+D)   ", "FRET ratio = A/D   ", "FRET ratio = D/A" ])


setconcat = False
showomexml = False
autoscale = False
Prefs.blackBackground = True

percThreshold = 0.05 # percentage of pixels for acceptor threshold 
selectIdx =True

lut = 'Fire'
#---------------------------------------------------------------


#---------------------------------------------------------------
#----------------- All Functions for analysis  -----------------
#---------------------------------------------------------------

#### Fonctions for PART 1 : Preparation of data & analysis - Preprocessing

# Get indexes of Donor and Acceptor images
def getImpIndexes(imagefile_ , sizeC_ , sizeZ_ , sizeT_ , idxseries_ , selectIdx_):
	#read in ImagePlus with arguments	
	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxseries_,True)	
	options.setZBegin(idxseries_, 0)
	options.setZEnd(idxseries_,sizeZ_-1)
	options.setTBegin(idxseries_, int(sizeT_/2))
	options.setTEnd(idxseries_,int(sizeT_/2))  
	impsLSM = BF.openImagePlus(options) # should be one imagestack 
	
	
	#projection in Z dimension if multiple-Z :impProj_
	IJ.run(impsLSM[0], "Grays", "stack")
	impProj_ = ZProjector.run(impsLSM[0],"max") if sizeZ_ >1 else impsLSM[0]
	impProj_.setCalibration(Calibration())
	impProj_.setDimensions(1, sizeC_, 1) #setDimensions(NChannels, NSlices, NFrames) for ZAxisProfiler function
	impProj_.setSlice(int(sizeC_/2))
	IJ.run(impProj_, "Enhance Contrast", "saturated=0.35")
			
	if selectIdx_ :
		IJ.setTool("rectangle")
		impProj_.setSlice(1)
		impProj_.show()	
		waitDialog = WaitForUserDialog("ROI","Select a ROI with the rectangle tool")
		waitDialog.show()
		plot_ = ZAxisProfiler.getPlot(impProj_)
		impPlot_ = plot_.getImagePlus()
		impPlot_.show()
		
		gui = GenericDialog("Select FRET Donor/Acceptor images")
		gui.addSlider("Donor image: ", 1, sizeC_, 3)
		gui.addSlider("Acceptor image: ", 1, sizeC_, 7)
		gui.showDialog()
		impPlot_.close()
		impProj_.killRoi()
		impProj_.close()
		idxDonor_ = int(gui.getNextNumber())
		idxAcceptor_ = int(gui.getNextNumber())
	else :
		idxDonor_, idxAcceptor_ = 3, 7
	for i in range(len(impsLSM)) :
		impsLSM[i].flush()
		
	return idxDonor_, idxAcceptor_

	
#Open ZCT hyperstack with only one channel index 
def extractImpFromIndex(imagefile_ , idxChannel_ , idxSerie_ , idxFrame_ = 0):
	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxSerie_,True)	
	options.setAutoscale(autoscale)
	options.setShowOMEXML(showomexml)
	options.setConcatenate(setconcat)
	if idxFrame_>0 :
		options.setTBegin(idxSerie_, idxFrame_)
		options.setTEnd(idxSerie_,idxFrame_)  
	options.setCBegin(idxSerie_, idxChannel_ -1)
	options.setCEnd(idxSerie_,idxChannel_ -1)
	return BF.openImagePlus(options)[0]


# Adjust the string size with 0
def adjustSizeNum(S_ , length_):
	newS_ = str(S_)
	while len(newS_)<length_ :
		newS_ = "0"+newS_
	return newS_
	
#create folder
def createFolder(imagePath_ , basename_):
	srcDir_ = os.path.dirname(imagePath_)
	anDir_ = os.path.join(srcDir_, basename_) 
	if not os.path.exists(anDir_):
		os.makedirs(anDir_)
	return anDir_



#### Fonctions for PART 2 : Analysis & Results 
# select the background ROI from an image 
def getBackgroundROI(imp_) :
	imp_.setSlice(imp_.getStackSize())
	imp_.show()
	IJ.setTool("rectangle")
	myWait = WaitForUserDialog("Select a ROI", "Select a ROI with the Rectangle tool for background substraction")
	myWait.show()
	roi_  = imp_.getRoi()
	imp_.killRoi
	imp_.hide()
	return roi_

# subtract the mean intensity of the ROI obtained from the function getBackgroundROI
def subtractBackground(imp_, roi_):
	Analyzer.setMeasurements(Measurements.MEAN)
	Analyzer.setPrecision(5)
	rt= ResultsTable.getResultsTable()
	analyzer = Analyzer(imp_,rt)
	for i in range(imp_.getStackSize()):
		rt.reset() 
		imp_.setSlice(i+1)
		imp_.setRoi(roi_)
		analyzer.measure()
		imp_.deleteRoi()
		IJ.run(imp_, "Subtract...", "value="+str(rt.getValue ('Mean', 0))+" slice")
	return

def bleachCorrection(imp_, CorrectionMethodIdx_, backROI_):
	impcorrected = imp_.duplicate()
	if(CorrectionMethodIdx == 0) : #Simple Ratio Method
		rt= ResultsTable.getResultsTable() 
		Analyzer.setMeasurements(Measurements.MEAN)
		Analyzer.setPrecision(5)
		analyzer = Analyzer(imp_,rt)
		imp_.setSlice(1)
		imp_.setRoi(backROI)
		analyzer.measure()
		imp_.deleteRoi()
		BCSR = BleachCorrection_SimpleRatio(impcorrected,rt.getValue('Mean', 0))
		BCSR.correctBleach()
	elif (CorrectionMethodIdx == 1): #Exponential Fitting Method
		BCEF = BleachCorrection_ExpoFit(impcorrected)
		WindowManager.getImage("y = a*exp(-bx) + c").hide()
		BCEF.core()
		impfit = WindowManager.getImage("y = a*exp(-bx) + c")
		impfit.hide()
		impfit.close()
	elif (CorrectionMethodIdx == 2):	#Histogram Matching Method
		BCMH  = BleachCorrection_MH(imp_)
		BCMH.doCorrection()
	return impcorrected
	
def bleachCorrection(imp_): #bleach correction by Exponential Fitting Method
	impcorrected = imp_.duplicate()
	BCEF = BleachCorrection_ExpoFit(impcorrected)
	BCEF.core()
	impfit = WindowManager.getImage("y = a*exp(-bx) + c")
	impfit.hide()
	impfit.close()
	return impcorrected
	
def applyThreshold(imp_, minthres_, maxthres_):
	IJ.setRawThreshold(imp_, minthres_, maxthres_, None)
	stackOK = "stack" if imp_.getStackSize() > 1 else ""
	IJ.run(imp_, "NaN Background", stackOK)	
	IJ.run(imp_, "Despeckle", stackOK)
	return 

# Computation of the fret index
def CalculationFRETmetric(impD_,impA_, FRETmetric_):
	if FRETmetric_ == FRETmetrics[0]:
		impD_ = IC.run(impD_,impA_,  "Add create 32-bit stack")#------- 1) image Donor+Acceptor -> Denominator
		IJ.setRawThreshold(impD_,1,Float.MAX_VALUE, None) # remove 0-value pixels to exclude infinity value in the  divide calculation
		IJ.run(impD_, "NaN Background", "stack")	
		impA_ = IC.run( impA_, impD_, "Divide create 32-bit stack")#------- 2) image Acceptor/(Donor+Acceptor)
		IJ.setRawThreshold(impA_,0,1, None)
		IJ.run(impA_, "NaN Background", "stack")	
		IJ.run(impA_, "Multiply...", "value=100 stack")
		IJ.run(impA_, "Enhance Contrast", "saturated=0.35 stack")
	elif FRETmetric_ == FRETmetrics[1]:
		impA_ = CalculationFRETratio(impD_,impA_)
	else :
		impA_ = CalculationFRETratio(impA_,impD_)
	return impA_

def CalculationFRETratio(imp1_,imp2_) : #image imp2_/imp1_
	IJ.setRawThreshold(imp1_,1,Float.MAX_VALUE, None) # remove 0-value pixels to exclude infinity value in the  divide calculation
	IJ.run(imp1_, "NaN Background", "stack")	
	impRatio_ = IC.run( imp2_, imp1_, "Divide create 32-bit stack")
	IJ.setRawThreshold(impRatio_,0,Float.MAX_VALUE, None) # remove extreme-value pixels (due to the division)
	IJ.run(impRatio_, "NaN Background", "stack")	
	IJ.run(impRatio_, "Enhance Contrast", "saturated=0.35 stack")
	return impRatio_
	

#draw calibration bar
def drawCalibrationBar(statsMin_, statsMax_):
	impBar_ = IJ.createImage("Calibration Bar", "32-bit black", 276, 50, 1)
	ipBar_ = impBar_.getProcessor()
	step= (statsMax_ - statsMin_)/255
	for i in range(256) :
		for j in range(30) :
			ipBar_.setf(11+i,j, statsMin_ + i*step)
	ipBar_.setColor(256)
	ipBar_.setFont(Font("SansSerif", Font.BOLD, 12))
	ipBar_.drawString("{:.1f}".format(statsMin_),0, 48)
	ipBar_.drawString("{:.1f}".format(statsMax_),246, 48)
	impBar_.updateAndDraw()
	impBar_.setDisplayRange(statsMin_, statsMax_)
	return impBar_

#---------------------------------------------------------------
#-----------------       End of Functions      -----------------
#---------------------------------------------------------------

#---------------------------------------------------------------
#     -----------------       Start       -----------------
#---------------------------------------------------------------


#### PART 1 : Preparation of data & analysis - Preprocessing

print "PART 1 : Preparation of data & analysis - Preprocessing"
# clear the console automatically when not in headless mode
uiService.getDefaultUI().getConsolePane().clear()



#close Result Table if opened
if IJ.isResultsWindow() :
	IJ.run("Clear Results", "")
	tw = ResultsTable().getResultsWindow()
	tw.close()

rm = RoiManager.getInstance()
if not rm:
	rm = RoiManager(False)
rm.reset()


if fileType == "Spectral Confocal LSM   " :
	
	od = OpenDialog("Select a spectral LSM image", None)
	lsmPath = od.getPath()

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
	
	
	#select the serie if several
	idxSerie=0
	if seriesCount>1 :
		gui = GenericDialog("Select image serie")
		gui.addSlider("Image series: ", 1, seriesCount, 3)
		gui.showDialog()
		idxSerie = int(gui.getNextNumber()-1)


	#Open a menu for selecting Donor and Acceptor image indexes
	print "Select Donor, Acceptor images and the serie index to analyze"
	idxDonor, idxAcceptor = getImpIndexes(lsmPath,sizeC, sizeZ, sizeT, idxSerie, selectIdx)
 

	#Extract Donor and Acceptor images 
	print "Process Donor and Acceptor images" 
	impDonor = extractImpFromIndex(lsmPath, idxDonor, idxSerie)
	impAcceptor = extractImpFromIndex(lsmPath, idxAcceptor, idxSerie)

	#Create Folder for Donor and Acceptor images
	basename = os.path.basename(os.path.splitext(lsmPath)[0]).replace(' ', '_').lower()
	basename += "_S" + adjustSizeNum(str(idxSerie), 2)
	imageDir = createFolder(lsmPath , basename)
	
	#Save Donor and Acceptor raw images
	IJ.saveAs(impDonor, "TIFF",os.path.join(imageDir, basename+"_c1.tif")) #save Donor raw image
	IJ.saveAs(impAcceptor, "TIFF",os.path.join(imageDir, basename+"_c2.tif")) #save Acceptor raw image

else : 

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
		sys.exit(0)
	
	
	#create folder for analysis 
	basename = os.path.basename(os.path.splitext(donorPath)[0]).replace(' ', '_').lower()
	imageDir = createFolder(donorPath , basename)

	#open the files 
	impDonor = Opener().openImage(donorPath)
	impAcceptor = Opener().openImage(acceptorPath)

cal = impDonor.getCalibration()
unit = cal.getUnit()
if unit =="micron" :
	unit ="um"


#### PART 2 :  Bleaching correction and substract background 
print 'PART 2 : Bleaching correction and substract background'

#Background subtraction
print 'Select a ROI in the background'
if bleachCorr or manualSub :
	backROI = getBackgroundROI(impAcceptor)

#Check the bit depth of the images
depth = impAcceptor.getBitDepth()
impAcceptor.setSlice(1)
impT_slice = impAcceptor.crop("whole-slice")

ImageConverter(impAcceptor).convertToGray32()
ImageConverter(impDonor).convertToGray32()

maxPix = impT_slice.getStatistics(Measurements.MIN_MAX).max
if maxPix < 4096 and depth > 8:
	depth = 12 # the camera is 12-bits dynamical range =[0,4095]
maxVal = math.pow(2,depth)-1
applyThreshold(impDonor, 0,  maxVal - 1)
applyThreshold(impAcceptor, 0,  maxVal - 1)


if bleachCorr :
	print 'Correction of the photobleaching '
	CorrectionMethodIdx = CorrectionMethods.index(CorrectionMethod)
	impDonor = bleachCorrection(impDonor, CorrectionMethodIdx, backROI)
	impAcceptor = bleachCorrection(impAcceptor, CorrectionMethodIdx, backROI)

if manualSub :
	subtractBackground(impDonor,backROI)
	subtractBackground(impAcceptor,backROI)
else : 
	IJ.run(impDonor, "Subtract Background...", "rolling="+str(rollingBall)+" stack")
	IJ.run(impAcceptor, "Subtract Background...", "rolling="+str(rollingBall)+" stack")


impT_slice.show()
ta = ThresholdAdjuster()
ta.show()
ta.update()
IJ.setAutoThreshold(impT_slice, "Default dark")
waitDialog = WaitForUserDialog("Manual threshold", "Please, adjust the threshold as desired, then press 'OK' (do not press 'Apply')") # human thresholding
waitDialog.show()
thres_min = impT_slice.getProcessor().getMinThreshold()
thres_max = impT_slice.getProcessor().getMaxThreshold()
ta.close()
impT_slice.hide()
impT_slice.close()

thres_min = max(thres_min, 1)
applyThreshold(impDonor, thres_min,  min(thres_max, maxVal - 1))
applyThreshold(impAcceptor, thres_min,  min(thres_max, maxVal - 1))

print "The threshold value is = ", thres_min


#IJ.saveAs(impDonor, "TIFF",os.path.join(imageDir, basename+"_c1thres.tif")) #save Donor raw image
#IJ.saveAs(impAcceptor, "TIFF",os.path.join(imageDir, basename+"_c2thres.tif")) #save Acceptor raw image



	
#### PART 3 :  Calculate FRET metric stack	
print 'PART 3 : Measurement of '+ FRETchoice

metric = "ratio" 
if FRETchoice == FRETmetrics[0]:
	metric = "index" 
elif FRETchoice == FRETmetrics[1]: 
	metric+="A_D"
else :
	metric+="D_A"
impFRET = CalculationFRETmetric(impDonor, impAcceptor, FRETchoice)
FRETTitle = "FRET_"+metric+"_"+os.path.basename(basename)
impFRET.setTitle(FRETTitle)
stats = impFRET.getStatistics(Measurements.MIN_MAX)
statsMax = stats.max
statsMin = stats.min
if FRETchoice == FRETmetrics[0]:
	statsMax = math.ceil(stats.max)
	statsMin = math.floor(stats.min)
impFRET.setDisplayRange(statsMin, statsMax)
impFRET.setCalibration(cal)
#IJ.run(impFRET, "Enhance Contrast", "saturated=0.35 stack")
IJ.run(impFRET, lut, "stack")
IJ.saveAs(impFRET, "TIFF",os.path.join(imageDir, FRETTitle))
impFRET.show()


if calibrationBar :
	impBar = drawCalibrationBar(statsMin, statsMax)
	IJ.run(impBar, lut , "")
	impBar.show()
	IJ.saveAs(impBar, "TIFF",os.path.join(imageDir, "FRET_CalibrationBar"))

print 'End of analysis'
