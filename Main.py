'''
Created on Jun 27, 2014

@author: haiderriaz
'''


from ij import IJ, ImageStack, ImagePlus 
from ij.plugin import Filters3D, Duplicator
#from ij.process import Filters
from ij.process import StackStatistics, FloatProcessor, ImageProcessor
from array import zeros 
import time
import java.lang.Math as javaMath
from net.imglib2.img import ImagePlusAdapter
from net.imglib2.algorithm.fft2 import FFTConvolution
from ij.gui import GenericDialog
#from ij.io import FileSaver   

 

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
dim = InputImg.getNSlices() 
print "number of slices:", dim

# Begin Quarantined Code Block

fact = 64-searchradius
substack = javaMath.ceil(dim/fact)

if substack == 2:
    fact=fact*2
    substack=1
    
if substack == 1:
    fact=dim
     
# End Quarantined Code Block
    

# Instantiate 3D Median Filter plugin  
f3d = Filters3D()

# Start of 3D Median Filter
start_time = time.time()
print "Preprocessing: 3D Median filter"
 
# Retrieve filtered stack 
medianFilteredStack = f3d.filter(InputStack, f3d.MEDIAN, 2, 2, 2) 



# Construct an ImagePlus from the stack 
medianFilteredImage = ImagePlus("MedianFiltered-Image", medianFilteredStack)

# End of 3D Median Filter
elapsed_time = time.time() - start_time
print "Elapsed time:", elapsed_time

# Get Image Statistics after Median 3D Filter
medianFilterStats = StackStatistics(medianFilteredImage)
print "mean:", medianFilterStats.mean, "minimum:", medianFilterStats.min, "maximum:", medianFilterStats.max

if background:
    newStats = StackStatistics(medianFilteredImage)
    mini = newStats.min
    average_val =  newStats.mean
    maxi = newStats.max
    delta = (maxi-mini)/10000
    k=3
    N = javaMath.pow(k, 3)
    ConvMed = Duplicator().run(medianFilteredImage)
    ConvMedWrap = ImagePlusAdapter.wrap(ConvMed)
    
    
    
    ones = FloatProcessor(3 , 3)
    ones.setValue(1)
    ones.fill()
    onesStack = ImageStack(3, 3)
    for i in xrange(1, 4): 
        onesStack.addSlice(ones)
    onesImage = ImagePlus("onesImage", onesStack)
    onesImageWrap = ImagePlusAdapter.wrap(onesImage)
    
    
    
    #FFTConvolution(ConvMedWrap , onesImageWrap).run()
    Stats = StackStatistics(ConvMed)
    print "mean:", Stats.mean, "minimum:", Stats.min, "maximum:", Stats.max
    
    
      
else:
    fp = FloatProcessor(medianFilteredImage.width , medianFilteredImage.height)
    fp.setValue(1)
    fp.fill()
    maskStack = ImageStack(medianFilteredImage.width, medianFilteredImage.height)   
    for i in xrange(1, dim+1): 
        maskStack.addSlice(fp)
    mask =  ImagePlus("mask", maskStack)
    


# Anscombe transform to convert Poisson noise into Gaussian noise
medFiltStack = medianFilteredImage.getStack()   
newMedFiltStack = ImageStack(medianFilteredImage.width, medianFilteredImage.height)
newInputStack = ImageStack(InputImg.width, InputImg.height)
for i in xrange(1 , medianFilteredImage.getNSlices() + 1):
    ip = medFiltStack.getProcessor(i).convertToFloat()
    ip2 = InputStack.getProcessor(i).convertToFloat()
    pixels = ip.getPixels()
    pixels2 = ip2.getPixels()
    for j in xrange (len(pixels)):
        pixels[j] = 2 * javaMath.sqrt(pixels[j] + (3.0/8.0)  )
        pixels2[j] = 2 * javaMath.sqrt(pixels2[j] + (3.0/8.0)  )
           
    newMedFiltStack.addSlice(ip)
    newInputStack.addSlice(ip2)    

medianFilteredImage = ImagePlus("MedianFiltered-Image", newMedFiltStack)
InputImg = ImagePlus("Input-Image", newInputStack)
    
#newStats = StackStatistics(medianFilteredImage)
#print "mean:", newStats.mean, "minimum:", newStats.min, "maximum:", newStats.max
#print medianFilteredImage.getNSlices() 
# Display result of 3D median filter 
medianFilteredImage.show()
#medianFilterStats = StackStatistics(medianFilteredImage)
#print "mean:", medianFilterStats.mean, "minimum:", medianFilterStats.min, "maximum:", medianFilterStats.max 
#print "number of output slices:", medianFilteredImage.getNSlices()

 
#fs =  FileSaver(newImage)
#fs.save() 


# Temporary Reusable code
#Stats = StackStatistics()
#print "mean:", Stats.mean, "minimum:", Stats.min, "maximum:", Stats.max

