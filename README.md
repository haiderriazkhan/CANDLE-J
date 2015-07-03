CANDLE-J
========

CANDLE-J is a multi-language image denoising software designed as an [Image-J](http://imagej.nih.gov/ij/) plugin and an open source alternative to the [CANDLE](http://www.bic.mni.mcgill.ca/ServicesSoftwareAdvancedImageProcessingTools/CANDLE/) MATLAB program. CANDLE-J is very [robust](http://www.ncbi.nlm.nih.gov/pubmed/22341767) for processing deep in vivo 3D multiphoton microscopy images where the signal to noise ratio is low (SNR).  

The denoising algorithms are written in Jython, C and Java. The C methods are accessed with the aide of the Java Native Access library. The Optimized Non Local Means Filter (ONLM) is the bottleneck of the CANDLE-J software and hence the ONLM algorithm is written as a multithreaded C program to improve performance. 


Haider Riaz Khan   
haider.riaz@mail.mcgill.ca  
Montreal Neurological Institute  
McGill University  
3801, Rue University, MP121 (lab)  
Montreal, QC H3A 2B4
