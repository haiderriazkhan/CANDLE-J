'''
*************************************************************************
 Haider Khan - haiderriazkhan@hotmail.com                               *
 Ruthazer Lab, Montreal Neurological Institute.                         *
 McGill University                                                      *                                      
                                                                        *
 Copyright (C) 2015 Haider Riaz Khan                                    *
                                                                        *
*************************************************************************
 The CANDLE algorithm is described in:                                  *
                                                                        *  
 P. Coupe, Martin Munz, Jose V.Manjon, Edward Ruthazer, D.Louis Collins.*
 A CANDLE for a deeper in-vivo insight. Medical Image Analysis,         *
 16(4):849-64 (2012).                                                   *
*************************************************************************

'''


from ij import IJ, ImageStack, ImagePlus 
from ij.plugin import Filters3D
from ij.process import StackStatistics
import array 
import time
from ij.gui import GenericDialog
from ij.io import FileSaver
import sys
import NativeCodeJNA
import InverseAnscombe



# function that opens up a dialog box where the user inputs filter parameters 
def getOptions():  
    gd = GenericDialog("Options") 
    gd.addMessage("Filter Parameters")  
    gd.addNumericField("Smoothing Parameter", 0.1, 2)  # show 2 decimals
    gd.addNumericField("Patch radius", 2 , 1)
    gd.addNumericField("Search volume radius", 3 , 1)  
    gd.showDialog()  
    
    if gd.wasCanceled():  
        print "User canceled dialog!"  
        sys.exit()  
    # Read out the options    
    beta = gd.getNextNumber()
    patchradius = gd.getNextNumber()
    searchradius = gd.getNextNumber()     
    return beta, patchradius, searchradius  
 


# get input image 
InputImg = IJ.openImage();




# Get Input Image Statistics  
InStats = StackStatistics(InputImg)
print "mean:", InStats.mean, "minimum:", InStats.min, "maximum:", InStats.max

options = getOptions()  
if options is not None:  
    beta, patchradius, searchradius = options  
    


 
# get the image stack within the ImagePlus 
InputStack = InputImg.getStack()
z = InputImg.getNSlices() 
print "number of slices:", z

    

# Instantiate 3D Median Filter plugin  
f3d = Filters3D()

# Start of 3D Median Filter
start_time = time.time()
print "Preprocessing: 3D Median filter"
 
# Retrieve filtered stack 
medianFilteredStack = f3d.filter(InputStack, f3d.MEDIAN, 2, 2, 2) 



# Construct an ImagePlus from the filtered stack 
medianFilteredImage = ImagePlus("MedianFiltered-Image", medianFilteredStack)

# End of 3D Median Filter
elapsed_time = time.time() - start_time
print "Elapsed time:", elapsed_time

# get image dimensions
x = medianFilteredImage.width
y = medianFilteredImage.height



# Get Image Statistics after Median 3D Filter
medianFilterStats = StackStatistics(medianFilteredImage)
print "mean:", medianFilterStats.mean, "minimum:", medianFilterStats.min, "maximum:", medianFilterStats.max

    


# Anscombe transform to convert Poisson noise into Gaussian noise

start_time = time.time()
print "Stabilization: Anscombe transform"

    
IJ.run(medianFilteredImage, "32-bit", ""); 
IJ.run(medianFilteredImage, "Add...", "value=0.375 stack"); 
IJ.run(medianFilteredImage, "Square Root", "stack"); 
IJ.run(medianFilteredImage, "Multiply...", "value=2 stack");

IJ.run(InputImg, "32-bit", ""); 
IJ.run(InputImg, "Add...", "value=0.375 stack"); 
IJ.run(InputImg, "Square Root", "stack"); 
IJ.run(InputImg, "Multiply...", "value=2 stack");       

              

# End of Anscombe transform
elapsed_time = time.time() - start_time
print "Elapsed time:", elapsed_time


Stats = StackStatistics(medianFilteredImage)
print "mean:", Stats.mean, "minimum:", Stats.min, "maximum:", Stats.max




medianFilteredStack = medianFilteredImage.getStack()
InputStack = InputImg.getStack()


# Get the Input and filtered Images as 1D arrays
medfiltArray = array.array('f')
InputImgArray = array.array('f')


for i in xrange(1 , z + 1):
    ip = medianFilteredStack.getProcessor(i).convertToFloat()
    ip2 = InputStack.getProcessor(i).convertToFloat()
    pixels = ip.getPixels()
    pixels2 = ip2.getPixels()
    medfiltArray.extend(pixels)
    InputImgArray.extend(pixels2)


InputImg.flush()


# Noise Estimation and Non-Local Means Filter

print "Going Native ..."
fimg = NativeCodeJNA.NativeCall(InputImgArray, medfiltArray, int(searchradius), int(patchradius), beta , int(x), int(y), int(z))




  

# Optimal Inverse Anscombe Transform
start_time = time.time()
print "Inverse Anscombe"
fimg = InverseAnscombe.OVST(fimg)
elapsed_time = time.time() - start_time
print "Elapsed time:", elapsed_time




outputstack = ImageStack(x, y, z )  

for i in xrange(0, z):  
    # Get the slice at index i and assign array elements corresponding to it.
    outputstack.setPixels(fimg[int(i*x*y):int((i+1)*x*y)], i+1)

print 'Preparing denoised image for display '
outputImp = ImagePlus("Output Image", outputstack)    
print "OutputImage Stats:"
Stats = StackStatistics(outputImp)
print "mean:", Stats.mean, "minimum:", Stats.min, "maximum:", Stats.max



outputImp.setDisplayRange(Stats.min, Stats.max) 

outputImp.show()  

 
fs =  FileSaver(outputImp)
fs.save() 
