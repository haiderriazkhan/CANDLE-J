'''
Created on Jun 27, 2014

@author: haiderriaz
'''


from ij import IJ, ImageStack, ImagePlus 
from ij.plugin import Filters3D, Duplicator
from ij.plugin.filter import ImageMath
from ij.process import StackStatistics, FloatProcessor, ImageProcessor
from ij.process import *
import array
from array import zeros
import jarray 
import time
import math
from ij.gui import GenericDialog
import LocalNoiseEstimation
import InverseAnscombe
#from ij.io import FileSaver   
#from mikera.vectorz import *



# function that opens up a dialog box where the user inputs filter parameters 
def getOptions():  
    gd = GenericDialog("Options") 
    gd.addStringField("Suffix for the filenames of the filtered images:", "Untitled")
    gd.addMessage("Filter Parameters")  
    gd.addNumericField("Smoothing Parameter", 0.1, 2)  # show 2 decimals
    gd.addNumericField("Patch radius", 2 , 1)
    gd.addNumericField("Search volume radius", 3 , 1)
    gd.addMessage("Background")  
    gd.addCheckbox("Fast processing of the dark background", True)  
    gd.showDialog()  
    
    if gd.wasCanceled():  
        print "User canceled dialog!"  
        return  
    # Read out the options  
    name = gd.getNextString()  
    beta = gd.getNextNumber()
    patchradius = gd.getNextNumber()
    searchradius = gd.getNextNumber()  
    background = gd.getNextBoolean()   
    return name, beta, patchradius, searchradius, background  
  


# get input image 
InputImg = IJ.openImage();





# Get Input Image Statistics  
InStats = StackStatistics(InputImg)
print "mean:", InStats.mean, "minimum:", InStats.min, "maximum:", InStats.max

options = getOptions()  
if options is not None:  
    name, beta, patchradius, searchradius, background = options  
    


 
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

    
    

# Background Detection
mask = [1]*(x*y*z)


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



# Estimation of local noise variance

#LocalNoiseEstimation.estimate(medfiltArray, x, y, z)




medianFilteredImage.show()


bob = [4 , 5, 0 , -1 , -6 , 0 , 7 , 200 , 300  , 40]
cock = [45 , 55, 0 , -14 , -64 , 0 , 74 , 29 , 330  , 470]
asympmap = [i for i, e in enumerate(bob) if e > 100]

print bob

lol =  2.0* math.sqrt(3.0/8.0)
biasedmap = [i for i, e in enumerate(bob) if e < lol] 

print biasedmap




 
#fs =  FileSaver(newImage)
#fs.save() 


# Temporary Reusable code
#Stats = StackStatistics()
#print "mean:", Stats.mean, "minimum:", Stats.min, "maximum:", Stats.max
#medfiltArray = map(lambda x: 2 * math.sqrt(x + const  ) ,  medfiltArray)    
#InputImgArray = map(lambda x: 2 * math.sqrt(x + const  ) ,  InputImgArray)
