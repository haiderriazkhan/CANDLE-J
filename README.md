CANDLE-J
========

CANDLE-J is a multi-language image denoising software designed as an [ImageJ](http://imagej.nih.gov/ij/) (64-bit) plugin and an open source alternative to the CANDLE MATLAB [program](http://www.bic.mni.mcgill.ca/ServicesSoftwareAdvancedImageProcessingTools/CANDLE/). CANDLE-J is very adept at processing deep in vivo 3D multiphoton microscopy images where the signal to noise ratio is low (SNR). The Collaborative Approach for eNhanced Denoising under Low-light Excitation (CANDLE) has been [tested](http://www.bic.mni.mcgill.ca/uploads/ServicesSoftwareAdvancedImageProcessingTools/candle.pdf) on synthetic data, images acquired on microspheres and in vivo images.  

The denoising algorithms are written in Jython, Java and C. The C methods are accessed with the aide of the Java Native Access library. The Optimized Non Local Means Filter (ONLM) is the bottleneck of the CANDLE-J software and hence the ONLM algorithm is written as a multithreaded C program to improve performance. Additionally, the two C programs were compiled with the GCC optimizing option `-O3` turned on. The `-O3` option makes the compiler strive to improve the performance of the object code.       

### Installing CANDLE-J
CANDLE-J is an ImageJ (64-bit) plugin and requires ImageJ to be installed on the computer. If ImageJ is not already installed, click [here](http://imagej.nih.gov/ij/download.html) to download the appropriate version. In the top header of this repository, click **release**. Depending on your operating system, either download the `MacOSX.zip` or `Linux.tar.gz` file from the listing. Unzip the downloaded file, and place the resultant `CANDLE-J` folder in the `plugins` folder of your local ImageJ directory. Open ImageJ (restart ImageJ if it is already open) and CANDLE-J should be available to use from the `Plugins` dropdown menu. 

**Important** 
- For a Mac OS X; CANDLE-J should be used with `ImageJ64`(64-bit) instead of the version simply named `ImageJ` (32-bit). `ImageJ64` comes bundled with the downloaded ImageJ zip file.
- In a Linux machine, use the version of ImageJ that appears with the default Linux executable icon. <img src="Images/LinuxLogo.png" width="24" height="24" />

### Usage
1. Open ImageJ
2. Go to the Plugins dropdown menu.
3. Hover your mouse on CANDLE-J and click CANDLE. ![alt text](/Images/ImageJ_Menu.png)
4. Yo 



Haider Riaz Khan   
haider.riaz@mail.mcgill.ca  
Montreal Neurological Institute  
McGill University  
3801, Rue University, MP121 (lab)  
Montreal, QC H3A 2B4
