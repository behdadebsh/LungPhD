==============
Lvl_Set_Alg.py
==============

This code implements the paper: "Active Contours Without Edges" By Chan & Vese. This is a nice way to segment images whose foregrounds and backgrounds are statistically different and homogeneous `[1] <https://ieeexplore.ieee.org/document/902291>`_.

lvlset(I, init_mask, max_its, alpha, thresh, color, display)

 :Inputs: 
	:I:           2D image
	:init_mask:   Initialization for the seed point (1 = foreground, 0 = bg)
	:max_its:     Number of iterations to run segmentation for
	:alpha:       (optional) Weight of smoothing term. higer = smoother.  default = 0.2
	:color:       Color for mask boundary visualisation. example: 'r' for red.
	:display:     (optional) displays intermediate outputs. default = true
 :Outputs:
	:seg:        Final segmentation mask (1=fg, 0=bg)

