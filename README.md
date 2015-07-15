CANDLE-J
========

CANDLE-J is a multi-language image denoising software designed as an [ImageJ](http://imagej.nih.gov/ij/) (64-bit) plugin and an open source alternative to the CANDLE MATLAB [program](http://www.bic.mni.mcgill.ca/ServicesSoftwareAdvancedImageProcessingTools/CANDLE/). CANDLE-J is very adept at processing deep in vivo 3D multiphoton microscopy images where the signal to noise ratio is low (SNR). The Collaborative Approach for eNhanced Denoising under Low-light Excitation (CANDLE) has been [tested](http://www.bic.mni.mcgill.ca/uploads/ServicesSoftwareAdvancedImageProcessingTools/candle.pdf) on synthetic data, images acquired on microspheres and in vivo images.  

The denoising algorithms are written in Jython, Java and C. The C methods are accessed with the aide of the Java Native Access library. The Optimized Non Local Means Filter (ONLM) is the bottleneck of the CANDLE-J software and hence the ONLM algorithm is written as a multithreaded C program to improve performance. Additionally, the two C programs were compiled with the GCC optimizing option `-O3` turned on. The `-O3` option makes the compiler strive to improve the performance of the object code.       

### Install
CANDLE-J is an ImageJ (64-bit) plugin and requires ImageJ to be installed on the computer. If ImageJ is not already installed, click [here](http://imagej.nih.gov/ij/download.html) to download the appropriate version. In the top header of this repository, click **release**. Depending on your operating system, either download the `MacOSX.zip` or `Linux.tar.gz` file from the listing. Unzip the downloaded file, and place the resultant `CANDLE-J` folder in the `plugins` folder of your local ImageJ directory. Open ImageJ (restart ImageJ if it is already open) and CANDLE-J should be available to use from the `Plugins` dropdown menu. 

**Important** 
- For a Mac OS X; CANDLE-J should be used with the default `ImageJ` (64-bit) instead of the version named `ImageJ32` (32-bit). Both come bundled with the downloaded ImageJ zip file.

- In a Linux machine, use the version of ImageJ that appears with the default Linux executable icon. <img src="Images/LinuxLogo.png" width="24" height="24" /> This version runs CANDLE-J from within the ImageJ directory and allows it to use relative file paths. 

### Usage
1- Open ImageJ

2- Go to the Plugins dropdown menu.

3- Hover your mouse on CANDLE-J and click CANDLE.

<div align="left">
        <img width="45%" src="/Images/ImageJ_Menu.png" </img>
</div>
4- Select an input 3D tiff image by navigating in the window that appears.

5 - A GUI to select filtering parameters will appear. Also a dialog box will show up to keep the user updated regarding the progress of the program. 

<div align="center">
        <img width="65%" height="100%" src="/Images/ParameterBox.png" </img>
</div>

- First, select the smoothing parameter (beta). This parameter controls the amount of denoising to be applied. Generally, values between 0.1 and 0.4 are fine for visualization. Higher values may be useful for segmentation or registration purposes.
- Second, select the patch radius. A radius of 1 voxel (i.e., patch of 3x3x3 voxels) is usually appropriate for a greater view of fine structures. A radius of 2 voxels (i.e., patch of 5x5x5 voxels) should be favoured for a restricted field of view of large structures.
- The radius of the search volume can also be changed but it is not recommended to increase it (higher computational time).


6- Once CANDLE-J has run its course, the ouput image will be displayed (can be manipulated further with ImageJ) along with a save dialog box. 

                                        *         *         *  


### Contact

Haider Riaz Khan   
haider.riaz@mail.mcgill.ca  
Montreal Neurological Institute  
McGill University  
3801, Rue University, MP121 (lab)  
Montreal, QC H3A 2B4
