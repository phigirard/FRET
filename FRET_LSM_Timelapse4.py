#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	FRET_LSM_Timelapse.py
#	Release v1.0
#
#	Copyright 2022 - BSD-3-Clause license
#                                                                             
#******************************************************************************/

#@ String msg1 (visibility=MESSAGE, value="--------------------------------------- Input File(s): ---------------------------------------", required=False) 
#@ String fileType (label="File type",choices={ "Spectral Confocal LSM", "Separate TIF files Donor/Acceptor"}, style="radioButtonHorizontal", persist=True) 
#@ File lsmFile (label="LSM File", description="Select a file", style="file") 
#@ File donorFile (label="Donor File", description="Select the donor file", style="file") 
#@ File acceptorFile (label="Acceptor File", description="Select the acceptor file", style="file") 
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
#@ String FRETchoice (label="FRET metric:",choices={"FRET index (A/A+D)", "FRET ratio (A/D)", "FRET ratio (D/A)"}, style="radioButtonHorizontal", persist=True) 
#@ String msg8 (visibility=MESSAGE, value="                                                                          ", required=False) 


#@ UIService uiService
#@ LogService log


#---------------------------------------------------------------
#---------------          LIBRARIES            -----------------
#---------------------------------------------------------------

# Python Library -------------------------------------------------------------------------------------------
# Python Library -------------------------------------------------------------------------------------------
import sys
import os, csv
#from math import *
#from random import *

from java.lang import Double, Float
from java.awt import Polygon, Color
from java.awt.event import AdjustmentListener  

# ImageJ Library ----------------------------------------------------------------------------
from ij import IJ, ImagePlus, ImageStack,  Prefs
from ij.io import Opener
from ij.process import ImageProcessor, ImageConverter, FloatProcessor, AutoThresholder, ImageStatistics
from ij.gui import GenericDialog, WaitForUserDialog, PlotWindow, ProfilePlot, Overlay, Line
from ij.gui import Roi
from ij.plugin import ZAxisProfiler, ZProjector, RGBStackMerge, RGBStackConverter, HyperStackReducer, Duplicator
from ij.plugin import ImageCalculator as IC
from ij.measure import ResultsTable , Measurements, CurveFitter, Calibration
from ij.plugin.filter import Analyzer, MaximumFinder
from ij.plugin.filter import ParticleAnalyzer as PA
from ij.plugin.frame import RoiManager,ThresholdAdjuster


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
FRETmetrics = list([ "FRET index (A/A+D)", "FRET ratio (A/D)", "FRET ratio (D/A)" ])


setconcat = False
showomexml = False
autoscale = False
Prefs.blackBackground = True

percThreshold = 0.05 # percentage of pixels for acceptor threshold 
selectIdx =True
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
		impA_.setDisplayRange(0, 100)
	elif FRETmetric_ == FRETmetrics[1]:
		impA_ = CalculationFRETratio(impD_,impA_)
	else :
		impA_ = CalculationFRETratio(impA_,impD_)
	return impA_

def CalculationFRETratio(imp1_,imp2_) : #image imp2_/imp1_
	IJ.setRawThreshold(imp1_,1,Float.MAX_VALUE, None) # remove 0-value pixels to exclude infinity value in the  divide calculation
	IJ.run(imp1_, "NaN Background", "stack")	
	impRatio_ = IC.run( imp2_, imp1_, "Divide create 32-bit stack")
	IJ.run(impRatio_, "NaN Background", "stack")	
	IJ.run(impRatio_, "Enhance Contrast", "saturated=0.35 stack")
	return impRatio_

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
	ImageConverter(imp_).convertToGray32()
	IJ.setRawThreshold(imp_, minthres_, maxthres_, "Black & White")
	if imp_.getStackSize() > 1 :
		IJ.run(imp_, "NaN Background", "stack")	
		IJ.run(imp_, "Despeckle", "stack")
	else :
		IJ.run(imp_, "NaN Background", "")	
		IJ.run(imp_, "Despeckle", "")
	return 


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


if fileType == "Spectral Confocal LSM" :
	#convert Files from #@ parameters to String 
	imagefile = lsmFile.getCanonicalPath()
	lsmPath = lsmFile.getCanonicalPath()

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
	impFolder = createFolder(lsmPath , idxSerie)
	basename_ = os.path.basename(os.path.splitext(lsmPath)[0]).replace(' ', '_').lower()
	basename_ += "_S" + adjustSizeNum(str(idxSeries), 2)
	basename = os.path.basename(impFolder)
	
	#Save Donor and Acceptor raw images
	IJ.saveAs(impDonor, "TIFF",os.path.join(impFolder, basename+"_c1.tif")) #save Donor raw image
	IJ.saveAs(impAcceptor, "TIFF",os.path.join(impFolder, basename+"_c2.tif")) #save Acceptor raw image

else : 
	#convert Files from #@ parameters to String and extract the main directory of the data
	donorPath = donorFile.getCanonicalPath()
	acceptorPath = acceptorFile.getCanonicalPath()

	#create folder for analysis 
	basename = os.path.basename(os.path.splitext(donorPath)[0]).replace(' ', '_').lower()
	impFolder = createFolder(donorPath , basename)

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
backROI = getBackgroundROI(impAcceptor)


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

impAcceptor.setSlice(1)
impT_slice = impAcceptor.crop("whole-slice")
impT_slice.show()
ta = ThresholdAdjuster()
ta.update()
ta.show()
IJ.setAutoThreshold(impT_slice, "Default dark")
waitDialog = WaitForUserDialog("Manual threshold", "Please, adjust the threshold as desired, then press 'OK' (do not press 'Apply')") # human thresholding
waitDialog.show()
thres_min = impT_slice.getProcessor().getMinThreshold()
thres_max = impT_slice.getProcessor().getMaxThreshold()
ta.close()
impT_slice.hide()
impT_slice.close()

applyThreshold(impDonor, thres_min, thres_max)
applyThreshold(impAcceptor, thres_min, thres_max)

	
#### PART 3 :  Calculate FRET metric stack
metric = "index" if FRETchoice == FRETmetrics[0] else "ratio"
	
print 'PART 3 : Calculate FRET '+metric+' stack'

impFRET = CalculationFRETmetric(impDonor, impAcceptor, FRETchoice)
FRETTitle = "FRET"+metric+"_"+os.path.basename(basename)+".tif"
impFRET.setTitle(FRETTitle)
impFRET.setCalibration(cal)
IJ.run(impFRET, "Enhance Contrast", "saturated=0.35")
IJ.run(impFRET, "Fire", "stack")
IJ.saveAs(impFRET, "TIFF",os.path.join(impFolder, FRETTitle))
impFRET.show()
print 'End of analysis'