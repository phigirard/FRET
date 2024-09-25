#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	FRET_Multiplex_Mariia.py
#	Release v2.0
#
#	Copyright 2024 - BSD-3-Clause license
#                                                                             
#******************************************************************************/

#@ String msg1 (visibility=MESSAGE, value="------------------------------------------- Input File(s): -------------------------------------------", required=False) 
#@ File FRET1Donnor (label="Select the donnor of the 1st FRET pair:", style="file")
#@ File FRET1Acepptor (label="Select the acceptor of the 1st FRET pair:", style="file") 
#@ File FRET2Donnor (label="Select the donnor of the 2nd FRET pair:", style="file")
#@ File FRET2Acepptor (label="Select the acceptor of the 2ndFRET pair:", style="file") 
#@ String msg2 (visibility=MESSAGE, value="                                                                          ", required=False) 
#@ String msg3 (visibility=MESSAGE, value="------------------------------------- Correction of the intensity decay due to photobleaching: ---------------------------------------", required=False) 
#@ Boolean bleachCorr (label="Correction ?", description="Bleaching correction",value=False, persist=True) 
#@ String CorrectionMethod (label="Correction method",choices={ "Simple Ratio", "Exponential Fit", "Histogram Matching" }, style="listBox", persist=True) 
#@ String msg4 (visibility=MESSAGE, value="                                                                         ", required=False) 
#@ String msg5 (visibility=MESSAGE, value="------------------------------------------ Manual or Automatic background subtraction: --------------------------------------------", required=False) 
#@ String ChoiceSub (label=" ", choices={"Manual (values below)", "Manual (ROI selection)","Automatic (ROI from threshold)", "Automatic (Rolling ball)"}, style="radioButtonHorizontal") 
#@ Integer BGValueDonor1 (label="Background value (Donor 1):", value=100, persist=True) 
#@ Integer BGValueAcceptor1 (label="Background value (Acceptor 1):", value=100, persist=True) 
#@ Integer BGValueDonor2 (label="Background value (Donor 2):", value=100, persist=True) 
#@ Integer BGValueAcceptor2 (label="Background value (Acceptor 2):", value=100, persist=True) 
#@ Integer rollingBall (label="Rolling ball radius", description="Subtract Background",value=50, persist=True)  
#@ String msg6 (visibility=MESSAGE, value="                                                                          ", required=False) 
#@ String msg7 (visibility=MESSAGE, value="------------------------------------------------ Manual or Automatic threshold: ---------------------------------------------------", required=False) 
#@ Boolean manualThreshold (label="Apply manual threshold value:", description="Manual threshold",value=False, persist=True) 
#@ Integer thresholdValue1 (label="Threshold value of the first pair:", value=100, persist=True) 
#@ Integer thresholdValue2 (label="Threshold value of the second pair:", value=100, persist=True) 
#@ String msg8 (visibility=MESSAGE, value="                                                                          ", required=False) 
#@ String msg9 (visibility=MESSAGE, value="------------------------------------------------------------- Other: --------------------------------------------------------------", required=False) 
#@ Double timeLapse (label="Timelapse (min)", description="What is the timelapse between 2 images?",value=5, persist=False)  
#@ String FRETchoice (label="FRET metric:",choices={"FRET index = 100 x A/(A+D)   ", "FRET ratio = A/D   ", "FRET ratio = D/A" }, style="radioButtonHorizontal", persist=True) 
#@ Boolean calibrationBar (label="Display Calibration Bar ?", description="Calibration Bar",value=True, persist=True) 
#@ String msg10 (visibility=MESSAGE, value="                                                                          ", required=False) 


#@ UIService uiService
#@ LogService log
#@ PrefService prefs 



# Python Library -------------------------------------------------------------------------------------------
import sys
import os
import csv
import math 

from java.io import File
from java.lang import Float
from java.awt import Font

# ImageJ Library ----------------------------------------------------------------------------
from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import Prefs
from ij.io import Opener
from ij.process import ImageProcessor
from ij.process import ColorProcessor
from ij.process import ImageConverter
from ij.process import AutoThresholder
from ij.process import ImageStatistics
from ij.gui import WaitForUserDialog
from ij.gui import YesNoCancelDialog
from ij.gui import Roi
from ij.plugin import RoiEnlarger
from ij.plugin import ImageCalculator as IC
from ij.measure import ResultsTable
from ij.measure import Measurements
from ij.measure import Calibration
from ij.plugin.filter import ThresholdToSelection
from ij.plugin.frame import RoiManager
from ij.plugin.frame import ThresholdAdjuster

from fiji.util.gui import GenericDialogPlus


# CMCI Library (Kota Miura) for Bleach Correction --------------------------------------------------------
from emblcmci import BleachCorrection_SimpleRatio,BleachCorrection_ExpoFit,BleachCorrection_MH



#---------------------------------------------------------------
#---------------         CONSTANTS            -----------------
#---------------------------------------------------------------
Prefs.blackBackground = True

CorrectionMethods = list([ "Simple Ratio", "Exponential Fit", "Histogram Matching" ])
FRETmetrics = list([ "FRET index = 100 x A/(A+D)   ", "FRET ratio = A/D   ", "FRET ratio = D/A" ])
Subvalues = ["Manual (values below)", "Manual (ROI selection)","Automatic (ROI from threshold)", "Automatic (Rolling ball)"]

# CSV headers for background/Threshold file
fieldnames = ['Frame #', 'Subtraction method', 'Donor Background', 'Acceptor Background', 'Min threshold', 'Max Threshold']
 
pairs = ["1st","2nd"]

lut = 'Fire'
#---------------------------------------------------------------


#---------------------------------------------------------------
#----------------- All Functions for analysis  -----------------
#---------------------------------------------------------------

#### Preparation of data & analysis - Preprocessing


# Adjust the string size with 0
# Adjust the string size with 0
def adjustSizeNum(S_ , length_):
	newS_ = str(S_)
	while len(newS_)<length_ :
		newS_ = "0"+newS_
	return newS_
	



#### Analysis & Results 

# select the background ROI from an image 
def getBackgroundROI(imp_) :
	print 'Select area for background substraction'
	imp_.show()
	IJ.setTool("rectangle")
	myWait = WaitForUserDialog("ROI selection", "Select a ROI with the Rectangle tool for background substraction")
	myWait.show()
	roi_  = imp_.getRoi()
	imp_.deleteRoi()
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


#### Preparation of data & analysis - Preprocessing

print 'Preparation of data & analysis - Preprocessing'
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


#convert Files from #@ parameters to String and extract the main directory of the data
donorSPath = [FRET1Donnor.getCanonicalPath(), FRET2Donnor.getCanonicalPath()]
acceptorSPath = [FRET1Acepptor.getCanonicalPath(),FRET2Acepptor.getCanonicalPath()]


#create folder for analysis 
srcDir = os.path.dirname(donorSPath[0])
basename = os.path.basename(os.path.splitext(donorSPath[0])[0]).lower()
imageDir = os.path.join(srcDir, basename) 
if not os.path.exists(imageDir):
	os.makedirs(imageDir)


for i in range(2):
	#Extract Donor and Acceptor images 
	print 'Process Donor and Acceptor images of '+pairs[i]+' pair' 
	impDonor = Opener().openImage(donorSPath[i])
	IJ.run(impDonor, "Grays", "")
	IJ.run(impDonor, "Smooth", "stack")
	impAcceptor = Opener().openImage(acceptorSPath[i])
	IJ.run(impAcceptor, "Grays", "")
	IJ.run(impAcceptor, "Smooth", "stack")
	
	# Get calibration and set time & pixel units 
	if i==0 :
		cal = impAcceptor.getCalibration()
		unit = cal.getUnit()
		pix2phys = cal.getX(1)	
		if unit == "micron" :
			unit ="um"

	# Remove physical calibration for processing steps
	impDonor.setCalibration(Calibration())
	impAcceptor.setCalibration(Calibration())
	
	print 'Bleaching correction and substract background'
	#Background subtraction
	print 'Select a ROI in the background'
	depth = impAcceptor.getBitDepth()
	if i == 0 :
		backROI = getBackgroundROI(impAcceptor)
		depth = impAcceptor.getBitDepth()
	
	ImageConverter(impAcceptor).convertToGray32()
	ImageConverter(impDonor).convertToGray32()
	impAcceptor.setSlice(1)
	impT_slice = impAcceptor.crop("whole-slice")
	
	#Check the bit depth of the images and remove saturated and null pixels (saturated pixel are above  2^depth )
	maxPix = impT_slice.getStatistics(Measurements.MIN_MAX).max
	if maxPix < 4096 and depth > 8:
		depth = 12 # the camera is 12-bits dynamical range =[0,4095]
	maxVal = math.pow(2,depth)-1
	applyThreshold(impDonor, 0,  maxVal-1)
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

	applyThreshold(impDonor, max(thres_min, 1),  min(thres_max, maxVal - 1))
	applyThreshold(impAcceptor, max(thres_min, 1),  min(thres_max, maxVal - 1))
	
	
	print 'Measurement of '+ FRETchoice

	metric = "ratio" 
	if FRETchoice == FRETmetrics[0]:
		metric = "index" 
	elif FRETchoice == FRETmetrics[1]: 
		metric+="A_D"
	else :
		metric+="D_A"

	impFRET = CalculationFRETmetric(impDonor, impAcceptor, FRETchoice)
	FRETTitle = "FRET_"+metric+"_"+os.path.basename(basename)+"_"+pairs[i]
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
	IJ.run(impFRET, "Fire", "stack")
	IJ.saveAs(impFRET, "TIFF",os.path.join(imageDir, FRETTitle))
	impFRET.show()
	
	if calibrationBar :
		impBar = drawCalibrationBar(statsMin, statsMax)
		IJ.run(impBar, lut , "")
		impBar.show()
		IJ.saveAs(impBar, "TIFF",os.path.join(imageDir, "FRET_CalibrationBar_"+pairs[i]))
	
print 'End of analysis'
