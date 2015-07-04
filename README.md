CANDLE-J
========

CANDLE-J is a multi-language image denoising software designed as an [ImageJ64](http://imagej.nih.gov/ij/) plugin and an open source alternative to the [CANDLE](http://www.bic.mni.mcgill.ca/ServicesSoftwareAdvancedImageProcessingTools/CANDLE/) MATLAB program. CANDLE-J is very robust for processing deep in vivo 3D multiphoton microscopy images where the signal to noise ratio is low (SNR). The CANDLE denoising routines have been [tested](http://www.ncbi.nlm.nih.gov/pubmed/22341767) on synthetic data, images acquired on microspheres and in vivo images.  

The denoising algorithms are written in Jython, Java and C. The C methods are accessed with the aide of the Java Native Access library. The Optimized Non Local Means Filter (ONLM) is the bottleneck of the CANDLE-J software and hence the ONLM algorithm is written as a multithreaded C program to improve performance. 

### Installing CANDLE-J
CANDLE-J is an ImageJ64 plugin and requires ImageJ to be installed on the computer. If ImageJ is not already installed, click [here](http://imagej.nih.gov/ij/download.html) to download the appropriate version. Depending on your operating system, either download the `CANDLE-J_MacOSX.zip` or `CANDLE-J_Windows.zip` file from the listings in this repository. Unzip the downloaded file, and place the resultant `CANDLE-J` folder in the `plugins` folder your local ImageJ directory. Open ImageJ (restart if it is already open) and CANDLE-J should be available to use from the `Plugins` dropdown menu. 



Haider Riaz Khan   
haider.riaz@mail.mcgill.ca  
Montreal Neurological Institute  
McGill University  
3801, Rue University, MP121 (lab)  
Montreal, QC H3A 2B4
