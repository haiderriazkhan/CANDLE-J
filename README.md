CANDLE-J
========

CANDLE-J is a multi-language image denoising software designed as an [ImageJ](http://imagej.nih.gov/ij/) (64-bit) plugin and an open source alternative to the [CANDLE](http://www.bic.mni.mcgill.ca/ServicesSoftwareAdvancedImageProcessingTools/CANDLE/) MATLAB program. CANDLE-J is very adept at processing deep in vivo 3D multiphoton microscopy images where the signal to noise ratio is low (SNR). The CANDLE denoising routines have been [tested](http://www.ncbi.nlm.nih.gov/pubmed/22341767) on synthetic data, images acquired on microspheres and in vivo images.  

The denoising algorithms are written in Jython, Java and C. The C methods are accessed with the aide of the Java Native Access library. The Optimized Non Local Means Filter (ONLM) is the bottleneck of the CANDLE-J software and hence the ONLM algorithm is written as a multithreaded C program to improve performance. The Collaborative Approach for eNhanced Denoising under Low-light Excitation (CANDLE) is described in a 2012 Medical Image Analysis   [paper](http://www.bic.mni.mcgill.ca/uploads/ServicesSoftwareAdvancedImageProcessingTools/candle.pdf) by Pierrick *et al.* .    

### Installing CANDLE-J
CANDLE-J is an ImageJ64 plugin and requires ImageJ to be installed on the computer. If ImageJ is not already installed, click [here](http://imagej.nih.gov/ij/download.html) to download the appropriate version. Depending on your operating system, either download the `MacOSX.zip` or `Linux.zip` file from the listings in this repository. Unzip the downloaded file, and place the resultant `CANDLE-J` folder in the `plugins` folder of your local ImageJ directory. Open ImageJ (restart ImageJ if it is already open) and CANDLE-J should be available to use from the `Plugins` dropdown menu. 

**Note** 
- For a Mac OS X; CANDLE-J should be used with `ImageJ64`(64-bit) instead of the version simply named `ImageJ` (32-bit). `ImageJ64` comes bundled with the downloaded ImageJ zip file.
- In a linux machine, use the version of ImageJ that appears with the default linux executable icon. 

Haider Riaz Khan   
haider.riaz@mail.mcgill.ca  
Montreal Neurological Institute  
McGill University  
3801, Rue University, MP121 (lab)  
Montreal, QC H3A 2B4
