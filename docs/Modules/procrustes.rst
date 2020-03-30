=========================
UpperArtery_Procrustes.py
=========================

This script reads the centerline calculated from "Centerline_Xi2.py" and the 1D template from the directory on your hard drive. By using procrustes, the 1D template is mapped on the calculated centerline. Procrustes analysis determines a linear transformation (translation, reflection, orthogonal rotation and scaling) of the nodes in mesh to best conform them to the nodes in referenceMesh, using the sum of squared errors as the goodness of fit criterion.

The calculated centerline has 27 nodes for the whole MPA (main, right and left). For the same region, 1D template of the arterial tree has 11 nodes. Therefore, the mapping correspondance is defined as,

       Node number on 1D template	(is mapped to)	Node number on centerline

                 1 ------------------------------------------>  1

                 2 ------------------------------------------>  3

                 3 ------------------------------------------>  5

                 4 ------------------------------------------>  7

                 5 ------------------------------------------>  27

                 6 ------------------------------------------>  11

                 7 ------------------------------------------>  18

                 8 ------------------------------------------>  21

                 9 ------------------------------------------>  15

                 10 -----------------------------------------> 23

                 11 -----------------------------------------> 26


To map as above, a transformation matrix from 1D template to centerline node coordinates in space is calculated.

procrustes(X, Y, scaling, reflection)
    :Inputs:
    
    referenceMesh, mesh
        meshes (as matrices) of target and input coordinates. they must have equal
        numbers of  nodes (rows), but mesh may have fewer dimensions
        (columns) than referenceMesh.

    scaling
        if False, the scaling component of the transformation is forced
        to 1

    reflection
        if 'best' (default), the transformation solution may or may not
        include a reflection component, depending on which fits the data
        best. setting reflection to True or False forces a solution with
        reflection or no reflection respectively.

    :Outputs:

    :d:
        the residual sum of squared errors, normalized according to a
        measure of the scale of referenceMesh, ((referenceMesh - referenceMesh.mean(0))**2).sum()

    :Z:
        the matrix of transformed Y-values

    :tform:
        a dict specifying the rotation, translation and scaling that
        maps X --> Y

