#*******************************************************************************
#                                                                             
#	Philippe GIRARD 
# 	Université Paris Cité, CNRS, Institut Jacques Monod, F-75013 Paris, France
#
# 	FRET_LSM_Timelapse.py
#	Release v13.0
#
#	Copyright 2022 - BSD-3-Clause license
#                                                                             
#******************************************************************************/

#@ String msg1 (visibility=MESSAGE, value="----------------------------------------------------------- Input File(s): -----------------------------------------------------------", required=False) 
#@ String fileType (label="File type: ",choices={ "Spectral Confocal LSM   ", "Separate TIF files Donor/Acceptor"}, style="radioButtonHorizontal", persist=True) 
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


#---------------------------------------------------------------
#---------------          LIBRARIES            -----------------
#---------------------------------------------------------------

# Python Library -------------------------------------------------------------------------------------------
import sys
import os
import csv
import math 

from java.io import File
from java.lang import Float
from java.awt import Font
from java.awt import Color
from java.awt.image import BufferedImage
from java.awt.geom import Rectangle2D, Ellipse2D

from org.jfree.chart import ChartPanel, JFreeChart
from org.jfree.chart.axis import NumberAxis
from org.jfree.chart.plot import XYPlot
from org.jfree.chart.renderer.xy import XYLineAndShapeRenderer
from org.jfree.data.xy import XYDataset, XYSeries, XYSeriesCollection

# ImageJ Library ----------------------------------------------------------------------------
from ij import IJ
from ij import ImagePlus
from ij import ImageStack
from ij import Prefs
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
Prefs.blackBackground = True


CorrectionMethods = list([ "Simple Ratio", "Exponential Fit", "Histogram Matching" ])
FRETmetrics = [ "FRET index = 100 x A/(A+D)   ", "FRET ratio = A/D   ", "FRET ratio = D/A" ]
Subvalues = ["Manual (values below)", "Manual (ROI selection)","Automatic (ROI from threshold)", "Automatic (Rolling ball)"]

# CSV headers for background/Threshold file
fieldnames = ['Frame #', 'Subtraction method', 'Donor Background', 'Acceptor Background', 'Min threshold', 'Max Threshold']

setconcat = False
showomexml = False
autoscale = False


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
	
		
	#projection in Z dimension :impProj_
	IJ.run(impsLSM[0], "Grays", "stack")
	impProj_ = ZProjector.run(impsLSM[0],"max")
	IJ.run(impProj_, "Enhance Contrast", "saturated=0.35")	
	

	idxDonor_, idxAcceptor_ = 3, 7		
	if selectIdx_ :
		IJ.setTool("rectangle")
		impProj_.show()
		IJ.run("Brightness/Contrast...")
		waitDialog = WaitForUserDialog("ROI","Select a ROI with the rectangle tool")
		waitDialog.show()
		
		impsLSM[0].setCalibration(Calibration())
		impsLSM[0].setDimensions(1, 1, sizeC_)
		plot = ZAxisProfiler.getPlot(impsLSM[0])
		plot.setXYLabels('Channel','Mean')
		xvalues = plot.getXValues()
		yvalues = plot.getYValues()
		series = XYSeries("Mean intensity channel")
		for i in range(len(xvalues)) :
			series.add(xvalues[i], yvalues[i])
		dataset = XYSeriesCollection(series) #Add series to dataset
		yaxis = NumberAxis("Mean")
		xaxis = NumberAxis("# Channel")
		offset = 0.5
		xaxis.setRange(1-offset, sizeC_+offset)
		r = XYLineAndShapeRenderer()
		r.setSeriesPaint(0, Color.RED)
		r.setSeriesShape(0, Ellipse2D.Double(-3.0,-3.0,6.0,6.0))
		xyplot = XYPlot(dataset, xaxis, yaxis, r)
		xyplot.setBackgroundPaint(Color.white)
		chart = JFreeChart(xyplot)
		chart.removeLegend()
		impPlot = IJ.createImage("Mean intensity channel", "RGB", 512, 512, 1);
		imagePlot = impPlot.getBufferedImage()
		chart.draw(imagePlot.createGraphics(), Rectangle2D.Float(0, 0, impPlot.getWidth(), impPlot.getHeight()))
		impPlot.setImage(imagePlot)
		impPlot.show()
		
		gui = GenericDialog("Select FRET Donnor/Acceptor images")
		gui.addSlider("Donnor image: ", 1, sizeC_, 3)
		gui.addSlider("Acceptor image: ", 1, sizeC_, 7)
		gui.showDialog()
		impPlot.close()
		impProj_.killRoi()
		impProj_.close()
		idxDonor_ = int(gui.getNextNumber())
		idxAcceptor_ = int(gui.getNextNumber())
			
	for i in range(len(impsLSM)) :
		impsLSM[i].flush()
		
	return idxDonor_, idxAcceptor_



#Open ZCT hyperstack with only one channel index 
def extractImpFromIndex(imagefile_ , idxChannel_ , idxSeries_ ):
	options = ImporterOptions()
	options.setId(imagefile_)
	options.setSeriesOn(idxSeries_,True)	
	options.setAutoscale(autoscale)
	options.setShowOMEXML(showomexml)
	options.setConcatenate(setconcat)
	options.setCBegin(idxSeries_ , idxChannel_ -1)
	options.setCEnd(idxSeries_ ,idxChannel_ -1)  
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
	IJ.run(imp_, "Enhance Contrast", "saturated=0.35")
	imp_.show()
	IJ.setTool("rectangle")
	myWait = WaitForUserDialog("Select a ROI", "Select a ROI with the Rectangle tool for background subtraction")
	myWait.show()
	roi_  = imp_.getRoi()
	imp_.killRoi
	imp_.hide()
	return roi_


def getBackgroundROI(imp_, roi_) :
	imp_.setRoi(roi_)
	mean_ = imp_.getStatistics(Measurements.MEAN).mean
	imp_.killRoi()
	return mean_


# subtract the mean grey intensity fro a specific ROI 
def subtractBG(imp_, roi_):
	imp_.setRoi(roi_)
	RoiEnlarger.enlarge(imp_, -10) #reduction of 10 pixels to be sure not to take into account the edge of the WH
	mean_ = imp_.getStatistics(Measurements.MEAN).mean
	IJ.run(imp_, "Select None", "")
	IJ.run(imp_, "Subtract...", "value="+str(mean_)+" slice")
	return mean_	
	

def applyROI2NAN(imp_, roi_):
	ip_ = imp_.getProcessor()
	ip_.setColor(float('nan'))
	ip_.fill(roi_)
	IJ.run(imp_, "Select None", "")
	return	

def bleachCorrection(imp_, CorrectionMethodIdx_, backValue_):
	impcorrected = imp_.duplicate()
	if(CorrectionMethodIdx == 0) : #Simple Ratio Method
		BCSR = BleachCorrection_SimpleRatio(impcorrected,backValue_)
		BCSR.correctBleach()
	elif (CorrectionMethodIdx == 1): #Exponential Fitting Method
		BCEF = BleachCorrection_ExpoFit(impcorrected)
		WindowManager.getImage("y = a*exp(-bx) + c").hide()
		BCEF.core()
		impfit = WindowManager.getImage("y = a*exp(-bx) + c")
		impfit.hide()
		impfit.close()
	elif (CorrectionMethodIdx == 2):	#Histogram Matching Method
		BCMH  = BleachCorrection_MH(impcorrected)
		BCMH.doCorrection()
	return impcorrected
	
	
def applyThreshold(imp_, minthres_, maxthres_):
	IJ.setRawThreshold(imp_, minthres_, maxthres_, None)
	IJ.run(imp_, "NaN Background", "")	
	IJ.run(imp_, "Despeckle", "")
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
	
	IJ.run(impDonor, "Grays", "stack")
	IJ.run(impAcceptor, "Grays", "stack")
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

# Get calibration and set time & pixel units 
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

# Create  Image stack 
stackDonor = ImageStack(width, height)
stackAcceptor = stackDonor.duplicate()


IJ.run(impDonor, "Enhance Contrast", "saturated=0.35")	
IJ.run(impAcceptor, "Enhance Contrast", "saturated=0.35")	


#### PART 2 :  Bleaching correction and substract background 
print 'PART 2 : Bleaching correction and substract background'

doBleachROI = True
doThreshold = True
# Create Array of Dictionnary for Background/Threshold  CSV File
infoImg = []
for slic in range(nbSlice): 
	if (nbSlice > 1) :
		print "Process image "+str(slic+1)+"/"+str(nbSlice)

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
	if (bleachCorr or ChoiceSub == Subvalues[1]) and doBleachROI :
		backROI = getBackgroundROI(impAcceptor_slice)
		d = YesNoCancelDialog(IJ.getInstance(), "Same ROI ?","Do you want to use the same ROI for all the images", "  Yes  ", "  No  ")
		doBleachROI = not d.yesPressed() 
	
	if bleachCorr :
		print 'Correction of the photobleaching '
		CorrectionMethodIdx = CorrectionMethods.index(CorrectionMethod)
		if BleachCorrChoice == Subvalues[0] :
			BleachCorrValueDonor = BGValueDonor
			BleachCorrValueAcceptor = BGValueAcceptor
		elif BleachCorrChoice == Subvalues[1] :
			BleachCorrValueDonor = getBackgroundROI(impDonor_slice, backROI)
			BleachCorrValueAcceptor = getBackgroundROI(impAcceptor_slice, backROI)
		impDonor = bleachCorrection(impDonor_slice, CorrectionMethodIdx, BleachCorrValueDonor)
		impAcceptor = bleachCorrection(impAcceptor_slice, CorrectionMethodIdx, BleachCorrValueAcceptor)
	
	impT_slice = impAcceptor_slice.duplicate()
	if doThreshold and not manualThreshold:
		#IJ.resetMinAndMax(impT_slice)
		IJ.run(impT_slice, "Enhance Contrast", "saturated=0.35")	
		impT_slice.show()
		ta = ThresholdAdjuster()
		ta.show()
		ta.update()
		IJ.setAutoThreshold(impT_slice, "Triangle dark")
		waitDialog = WaitForUserDialog("Manual threshold", "Please, adjust the threshold as desired, then press 'OK' (do not press 'Apply')") # human thresholding
		waitDialog.show()
		thres_min = impT_slice.getProcessor().getMinThreshold()
		thres_max = impT_slice.getProcessor().getMaxThreshold()
		ta.close()
		impT_slice.hide()
		if (nbSlice > 1) :
			d = YesNoCancelDialog(IJ.getInstance(), "Same threshold values ?","Do you want to proceed automatically with the threshold values for all the images", "  Yes  ", "  No  ")
			doThreshold = not d.yesPressed()
		print "The min threshold value is = ", thres_min
	if slic == 0 and manualThreshold :
		thres_min = thresholdValue
		thres_max = maxVal
	IJ.setThreshold(impT_slice, thres_min, thres_max)
	roiCells = ThresholdToSelection.run(impT_slice) # Threshold selection in the image = signal
	roiNan = roiCells.getInverse(impT_slice) #get the inverse roi of  roiCells = Background
	impT_slice.close()
	
	if ChoiceSub == Subvalues[0] :
		IJ.run(impDonor_slice, "Subtract...", "value="+str(BGValueDonor)+" slice")
		IJ.run(impAcceptor_slice, "Subtract...", "value="+str(BGValueAcceptor)+" slice")
	elif ChoiceSub == Subvalues[1] :
		BGValueDonor = subtractBG(impDonor_slice,backROI)
		BGValueAcceptor = subtractBG(impAcceptor_slice,backROI)
		print "The background value of the donor image  is = ", BGValueDonor
		print "The background value of the acceptor image  is = ", BGValueAcceptor
	elif ChoiceSub == Subvalues[2] :
		BGValueDonor = subtractBG(impDonor_slice,roiNan)
		BGValueAcceptor = subtractBG(impAcceptor_slice,roiNan)
		print "The background value of the donor image  is = ", BGValueDonor
		print "The background value of the acceptor image  is = ", BGValueAcceptor
	elif ChoiceSub == Subvalues[3] :
		IJ.run(impDonor_slice, "Subtract Background...", "rolling="+str(rollingBall)+" stack")
		IJ.run(impAcceptor_slice, "Subtract Background...", "rolling="+str(rollingBall)+" stack")
		BGValueDonor = rollingBall
		BGValueAcceptor = rollingBall

	# Convert roiNan in NaN values in all images
	applyROI2NAN(impDonor_slice,roiNan)
	applyROI2NAN(impAcceptor_slice,roiNan)

	# Create dictionnary values for the frame 'slic+1'
	BackThres = [str(slic+1),ChoiceSub, BGValueDonor, BGValueAcceptor,  thres_min, thres_max]
	infoImg.append(dict(zip(fieldnames, BackThres)))
	
	stackDonor.addSlice(impDonor_slice.getProcessor())
	stackAcceptor.addSlice(impAcceptor_slice.getProcessor())


with  open(os.path.join(imageDir,"infoFile.csv"), 'wb') as csvfile:
	writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
	writer.writeheader()
	writer.writerows(infoImg)	 

impDonor_OUT=ImagePlus(impDonor.getTitle(), stackDonor)
IJ.run(impDonor_OUT, "Enhance Contrast", "saturated=0.35 stack")
IJ.saveAs(impDonor_OUT, "TIFF",os.path.join(imageDir, basename+"_c1thres.tif")) #save thresholded Donor image
impAcceptor_OUT=ImagePlus(impAcceptor.getTitle(), stackAcceptor)
IJ.run(impAcceptor_OUT, "Enhance Contrast", "saturated=0.35 stack")
IJ.saveAs(impAcceptor_OUT, "TIFF",os.path.join(imageDir, basename+"_c2thres.tif")) #save thresholded Acceptor image
	
#### PART 3 :  Calculate FRET metric stack	
print 'PART 3 : Measurement of '+ FRETchoice

metric = "ratio" 
if FRETchoice == FRETmetrics[0]:
	metric = "index" 
elif FRETchoice == FRETmetrics[1]: 
	metric+="A_D"
else :
	metric+="D_A"
impFRET = CalculationFRETmetric(impDonor_OUT, impAcceptor_OUT, FRETchoice)
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

Analyzer.setMeasurements(Measurements.AREA+ Measurements.MEAN +Measurements.STD_DEV)
Analyzer.setPrecision(5)
rt= ResultsTable() 
for slic in range(nbSlice):
	impFRET.setSlice(slic+1)
	impFRET_slice = impFRET.crop("whole-slice")
	analyzer = Analyzer(impFRET_slice,rt)
	analyzer.measure()
rt.show("Mean FRET index (%)")

if calibrationBar :
	impBar = drawCalibrationBar(statsMin, statsMax)
	IJ.run(impBar, lut , "")
	impBar.show()
	IJ.saveAs(impBar, "TIFF",os.path.join(imageDir, "FRET_CalibrationBar"))

print 'End of analysis'
