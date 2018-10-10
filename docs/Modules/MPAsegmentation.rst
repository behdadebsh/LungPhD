================================
Pulmonary_Artery_Segmentation.py
================================

Automatic segmentation of the main pulmonary artery (MPA) can be challenging due to complexity of the vasculature in the chest area. This script is providing a semi-automatic method for this purpose. This code takes advantage of a graphical interface (consisting of the image slice) for user to input a seed point for segmentation. Hence, user should be familiar with the main cardio-pulmonary vasculature anatomy. The seed point is given as an input to the program by user by clicking on the area of interest. There are some details provided step by step on how this cript works below :

- After choosing from which slice to which you are going to do the segmentation (the value k in the code), you run the code.
- The first slice as defined in the previous step will be shown in a canvas.
- Click on one point inside the MPA area in the slice (The local coordinate of your click will be shown on the console).
- You can do multiple click in the previous step but, it will regard them multiple seed points. This option is put in the script for bifurcating vessels so you can do multiple clicks whenever needed to follow the daughter vessels. Note that this option is still a work in progress, therefore, not robust enough yet!
- Then make sure you close the canvas after you are done with inputing seed points.
- It will start growing the mask using level set algorithm and the process is shown on a window.
- Next slice will be viewed and waiting for another seed point from user and you should do it exactly as previous time.
- This iteration will continue while it gets to finish all the slices for all the k values defined before running the script.

All the segmented masks will be saved in a directory given in the code with jpg format and will be ready to be used for other purposes. The segmentation method and algorithm being used is based on the method of Chan & Vese `[1] <https://ieeexplore.ieee.org/document/902291>`_ and is modified for our purpose.

Note:

- For this algorithm to work you need to have contrast enhanced CT images.
