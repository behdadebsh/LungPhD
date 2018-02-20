#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 15:54:59 2018

This script is for Lung Segmentation from DICOM images. Loading the scans from a directory and
transforms the pixels to Hounsfield Units. The features implemented in these codes are written
in a way to read the series of DICOM images from a folder and convert the voxel values into
hounsfield unit numbers. In order to segment the lung, two functions are made to perform this action.
First, markers are made via thresholding and then labeled. Labels are sorted from smaller areas to
bigger areas. Therefore, the 2 largest areas would be lungs. Second, using some morphological tools
and using the output from the first function to make lung mask and segment slices. Another feature of this
script is to visualise the outputs in 2D and 3D graphics.

@author: Behdad Shaarbaf Ebrahimi & Kaggle
            This code has been modified by Behdad Shaarbaf Ebrahimi (UoA - ABI) for academic purposes.
"""

import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as ndimage
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plotly import figure_factory as FF
from plotly.offline import init_notebook_mode, iplot
from skimage import measure, morphology, segmentation

from load_scan import get_pixels_hu, load_scan

init_notebook_mode(connected=True)


# Load the scans in given folder path


path = '/hpc/bsha219/lung/Data/ST12/Raw/DICOMS'
patient_scans = load_scan(path)
patient_images = get_pixels_hu(patient_scans)
slice = input('Which Slice do you want to mask? Enter it here: ')
print ("Original Slice")
plt.imshow(patient_images[slice], cmap='gray')
plt.show()


# Some of the starting Code is taken from ArnavJain, since it's more readable then my own
def generate_markers(image):
    # Creation of the internal Marker
    marker_internal = image < -500
    marker_internal = segmentation.clear_border(marker_internal)
    marker_internal_labels = measure.label(marker_internal)
    areas = [r.area for r in measure.regionprops(marker_internal_labels)]
    areas.sort()
    if len(areas) > 2:
        for region in measure.regionprops(marker_internal_labels):
            if region.area < areas[-2]:
                for coordinates in region.coords:                
                       marker_internal_labels[coordinates[0], coordinates[1]] = 0
    
    marker_internal_labels = measure.label(marker_internal_labels)
    left_lung = marker_internal_labels == 1
    right_lung = marker_internal_labels == 2
    marker_internal = marker_internal_labels > 0
    
    # Creation of the external Marker
    
    external_a = ndimage.binary_dilation(marker_internal, iterations=30)
    external_b = ndimage.binary_dilation(marker_internal, iterations=60)
    
    ex_left_a = ndimage.binary_dilation(left_lung,iterations=30)
    ex_left_b = ndimage.binary_dilation(left_lung,iterations=60)
    
    ex_right_a = ndimage.binary_dilation(right_lung,iterations=30)
    ex_right_b = ndimage.binary_dilation(right_lung,iterations=60)
    
    marker_external = external_b ^ external_a
    marker_ex_left = ex_left_b ^ ex_left_a
    marker_ex_right = ex_right_b ^ ex_right_a
    
    # Creation of the Watershed Marker matrix
    
    marker_watershed = np.zeros((512, 512), dtype=np.int)
    watershed_left = np.zeros((512,512), dtype=np.int)
    watershed_right = np.zeros((512,512), dtype=np.int)
    
    marker_watershed += marker_internal * 255
    marker_watershed += marker_external * 128
    watershed_left += left_lung * 255
    watershed_left += marker_ex_left * 128
    watershed_right += right_lung * 255
    watershed_right += marker_ex_right * 128
    
    return marker_internal, marker_external, marker_watershed, left_lung , right_lung, marker_ex_left, marker_ex_right, watershed_left, watershed_right

# Show some example markers from the middle

#test_patient_internal, test_patient_external, test_patient_watershed, left_lung, right_lung, ex_left, ex_right, watershed_left, watershed_right = generate_markers(patient_images[slice])
#print ("Internal Marker")
#plt.imshow(test_patient_internal, cmap='gray')
#plt.show()
#print ("Left Lung Marker")
#plt.imshow(left_lung, cmap='gray')
#plt.show()
#print ("Right Lung Marker")
#plt.imshow(right_lung, cmap='gray')
#plt.show()


def lung_segment(image):
    #Creation of the markers as shown above:
    marker_internal, marker_external, marker_watershed, left_lung, right_lung, ex_left, ex_right, watershed_left, watershed_right = generate_markers(image)
    
    #Creation of the Sobel-Gradient
    sobel_filtered_dx = ndimage.sobel(image, 1)
    sobel_filtered_dy = ndimage.sobel(image, 0)
    sobel_gradient = np.hypot(sobel_filtered_dx, sobel_filtered_dy)
    sobel_gradient *= 255.0 / np.max(sobel_gradient)
    
    #Watershed algorithm
    watershed = morphology.watershed(sobel_gradient, marker_watershed)
    watershed_left = morphology.watershed(sobel_gradient, watershed_left)
    watershed_right = morphology.watershed(sobel_gradient, watershed_right)
    
    for i in range(512):
        for j in range(512):
            if watershed_left[i,j] == 128:
                watershed_left[i,j] = 0
            if watershed_left[i,j] == 255:
                watershed_left[i,j] = 1
            if watershed_right[i,j] == 128:
                watershed_right[i,j] = 0
            if watershed_right[i,j] == 255:
                watershed_right[i,j] = 1
    
    
    seg_left = image * watershed_left
    seg_right = image * watershed_right
    segmentation = seg_left + seg_right
    segmentation = np.where(segmentation != 0, image, -2000*np.ones((512, 512)))
    vesselness = segmentation > 30
    vesselness = image * vesselness
    #Reducing the image created by the Watershed algorithm to its outline
#    outline = ndimage.morphological_gradient(watershed, size=(3,3))
#    outline = outline.astype(bool)
    
    #Performing Black-Tophat Morphology for reinclusion
    #Creation of the disk-kernel and increasing its size a bit
#    blackhat_struct = [[0, 0, 1, 1, 1, 0, 0],
#                       [0, 1, 1, 1, 1, 1, 0],
#                       [1, 1, 1, 1, 1, 1, 1],
#                       [1, 1, 1, 1, 1, 1, 1],
#                       [1, 1, 1, 1, 1, 1, 1],
#                       [0, 1, 1, 1, 1, 1, 0],
#                       [0, 0, 1, 1, 1, 0, 0]]
#    blackhat_struct = ndimage.iterate_structure(blackhat_struct, 8)
    #Perform the Black-Hat
#    outline += ndimage.black_tophat(outline, structure=blackhat_struct)
    
    #Use the internal marker and the Outline that was just created to generate the lungfilter
#    lungfilter = np.bitwise_or(marker_internal, outline)
    #Close holes in the lungfilter
    #fill_holes is not used here, since in some slices the heart would be reincluded by accident
#    lungfilter = ndimage.morphology.binary_closing(lungfilter, structure=np.ones((5,5)), iterations=3)
    
    #Apply the lungfilter (note the filtered areas being assigned -2000 HU)
#    segmented = np.where(lungfilter == 1, image, -2000*np.ones((512, 512)))
    
    return segmentation, watershed, marker_internal, marker_external, marker_watershed, watershed_left, watershed_right


def vesselness(image):
    segmented, watershed, marker_internal, marker_external, marker_watershed, watershed_left, watershed_right = lung_segment(image)
    
    vessels = segmented > 0
    vessel_mask = ndimage.binary_opening(vessels,iterations=2)
    vessels = image * vessel_mask
#    print ("Vessel mask")
#    plt.imshow(vessel_mask, cmap='gray')
#    plt.show()
#    print ("Vessel Segmentation")
#    plt.imshow(vessels, cmap='gray')
#    plt.show()
    return vessel_mask, vessels

###############################################################################
############## 3D Visualisation ###############################################
###############################################################################
###############################################################################

def make_mesh(image, threshold=30, step_size=1):
    print "Transposing surface"
    p = image.transpose(2, 1, 0)

    print "Calculating surface"
    verts, faces, norm, val = measure.marching_cubes(p, threshold, step_size=step_size, allow_degenerate=True)
    return verts, faces


def plotly_3d(verts, faces):
    x, y, z = zip(*verts)

    print "Drawing"

    # Make the colormap single color since the axes are positional not intensity.
    #    colormap=['rgb(255,105,180)','rgb(255,255,51)','rgb(0,191,255)']
    colormap = ['rgb(236, 236, 212)', 'rgb(236, 236, 212)']

    fig = FF.create_trisurf(x=x,
                            y=y,
                            z=z,
                            plot_edges=False,
                            colormap=colormap,
                            simplices=faces,
                            backgroundcolor='rgb(64, 64, 64)',
                            title="Interactive Visualization")
    iplot(fig)


def plt_3d(verts, faces):
    print "Drawing"
    x, y, z = zip(*verts)
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Fancy indexing: `verts[faces]` to generate a collection of triangles
    mesh = Poly3DCollection(verts[faces], linewidths=0.05, alpha=1)
    face_color = [1, 1, 0.9]
    mesh.set_facecolor(face_color)
    ax.add_collection3d(mesh)

    ax.set_xlim(0, max(x))
    ax.set_ylim(0, max(y))
    ax.set_zlim(0, max(z))
    ax.set_axis_bgcolor((0.7, 0.7, 0.7))
    plt.show()

#class VolumeSlicer(HasTraits):
#    # The data to plot
#    data = Array()
#
#    # The 4 views displayed
#    scene3d = Instance(MlabSceneModel, ())
#    scene_x = Instance(MlabSceneModel, ())
#    scene_y = Instance(MlabSceneModel, ())
#    scene_z = Instance(MlabSceneModel, ())
#
#    # The data source
#    data_src3d = Instance(Source)
#
#    # The image plane widgets of the 3D scene
#    ipw_3d_x = Instance(PipelineBase)
#    ipw_3d_y = Instance(PipelineBase)
#    ipw_3d_z = Instance(PipelineBase)
#
#    _axis_names = dict(x=0, y=1, z=2)
#
#
#    #---------------------------------------------------------------------------
#    def __init__(self, **traits):
#        super(VolumeSlicer, self).__init__(**traits)
#        # Force the creation of the image_plane_widgets:
#        self.ipw_3d_x
#        self.ipw_3d_y
#        self.ipw_3d_z
#
#
#    #---------------------------------------------------------------------------
#    # Default values
#    #---------------------------------------------------------------------------
#    def _data_src3d_default(self):
#        return mlab.pipeline.scalar_field(self.data,
#                            figure=self.scene3d.mayavi_scene)
#
#    def make_ipw_3d(self, axis_name):
#        ipw = mlab.pipeline.image_plane_widget(self.data_src3d,
#                        figure=self.scene3d.mayavi_scene,
#                        plane_orientation='%s_axes' % axis_name)
#        return ipw
#
#    def _ipw_3d_x_default(self):
#        return self.make_ipw_3d('x')
#
#    def _ipw_3d_y_default(self):
#        return self.make_ipw_3d('y')
#
#    def _ipw_3d_z_default(self):
#        return self.make_ipw_3d('z')
#
#
#    #---------------------------------------------------------------------------
#    # Scene activation callbaks
#    #---------------------------------------------------------------------------
#    @on_trait_change('scene3d.activated')
#    def display_scene3d(self):
#        outline = mlab.pipeline.outline(self.data_src3d,
#                        figure=self.scene3d.mayavi_scene,
#                        )
#        self.scene3d.mlab.view(40, 50)
#        # Interaction properties can only be changed after the scene
#        # has been created, and thus the interactor exists
#        for ipw in (self.ipw_3d_x, self.ipw_3d_y, self.ipw_3d_z):
#            # Turn the interaction off
#            ipw.ipw.interaction = 0
#        self.scene3d.scene.background = (0, 0, 0)
#        # Keep the view always pointing up
#        self.scene3d.scene.interactor.interactor_style = \
#                                 tvtk.InteractorStyleTerrain()
#
#
#    def make_side_view(self, axis_name):
#        scene = getattr(self, 'scene_%s' % axis_name)
#
#        # To avoid copying the data, we take a reference to the
#        # raw VTK dataset, and pass it on to mlab. Mlab will create
#        # a Mayavi source from the VTK without copying it.
#        # We have to specify the figure so that the data gets
#        # added on the figure we are interested in.
#        outline = mlab.pipeline.outline(
#                            self.data_src3d.mlab_source.dataset,
#                            figure=scene.mayavi_scene,
#                            )
#        ipw = mlab.pipeline.image_plane_widget(
#                            outline,
#                            plane_orientation='%s_axes' % axis_name)
#        setattr(self, 'ipw_%s' % axis_name, ipw)
#
#        # Synchronize positions between the corresponding image plane
#        # widgets on different views.
#        ipw.ipw.sync_trait('slice_position',
#                            getattr(self, 'ipw_3d_%s'% axis_name).ipw)
#
#        # Make left-clicking create a crosshair
#        ipw.ipw.left_button_action = 0
#        # Add a callback on the image plane widget interaction to
#        # move the others
#        def move_view(obj, evt):
#            position = obj.GetCurrentCursorPosition()
#            for other_axis, axis_number in self._axis_names.iteritems():
#                if other_axis == axis_name:
#                    continue
#                ipw3d = getattr(self, 'ipw_3d_%s' % other_axis)
#                ipw3d.ipw.slice_position = position[axis_number]
#
#        ipw.ipw.add_observer('InteractionEvent', move_view)
#        ipw.ipw.add_observer('StartInteractionEvent', move_view)
#
#        # Center the image plane widget
#        ipw.ipw.slice_position = 0.5*self.data.shape[
#                    self._axis_names[axis_name]]
#
#        # Position the view for the scene
#        views = dict(x=( 0, 90),
#                     y=(90, 90),
#                     z=( 0,  0),
#                     )
#        scene.mlab.view(*views[axis_name])
#        # 2D interaction: only pan and zoom
#        scene.scene.interactor.interactor_style = \
#                                 tvtk.InteractorStyleImage()
#        scene.scene.background = (0, 0, 0)
#
#
#    @on_trait_change('scene_x.activated')
#    def display_scene_x(self):
#        return self.make_side_view('x')
#
#    @on_trait_change('scene_y.activated')
#    def display_scene_y(self):
#        return self.make_side_view('y')
#
#    @on_trait_change('scene_z.activated')
#    def display_scene_z(self):
#        return self.make_side_view('z')
#
#
#    #---------------------------------------------------------------------------
#    # The layout of the dialog created
#    #---------------------------------------------------------------------------
#    view = View(HGroup(
#                  Group(
#                       Item('scene_y',
#                            editor=SceneEditor(scene_class=Scene),
#                            height=250, width=300),
#                       Item('scene_z',
#                            editor=SceneEditor(scene_class=Scene),
#                            height=250, width=300),
#                       show_labels=False,
#                  ),
#                  Group(
#                       Item('scene_x',
#                            editor=SceneEditor(scene_class=Scene),
#                            height=250, width=300),
#                       Item('scene3d',
#                            editor=SceneEditor(scene_class=MayaviScene),
#                            height=250, width=300),
#                       show_labels=False,
#                  ),
#                ),
#                resizable=True,
#                title='Pulmonary Physiology Predictor | Volume Segmentation',
#                )


######## Script:
#segment, watershed, marker_internal, marker_external, marker_watershed, watershed_left, watershed_right = lung_segment(patient_images[slice])


#print ("Lung Mask")
#plt.imshow(watershed, cmap='gray')
#plt.show()
#print ("Left Lung Mask")
#plt.imshow(watershed_left, cmap='gray')
#plt.show()
#print ("Right Lung Mask")
#plt.imshow(watershed_right, cmap='gray')
#plt.show()
#print ("Segmented Lung")
#plt.imshow(segment, cmap='gray')
#plt.show()

#vessel = vesselness(patient_images[slice])

#v, f = make_mesh(patient_images, 30, 2)
#plotly_3d(v, f)

segment_lung = []
segment_vessel = []
for sl in range(len(patient_scans)):
    seg, watershed, marker_internal, marker_external, marker_watershed, watershed_left, watershed_right = lung_segment(patient_images[sl])
    segment_lung.append(seg)
    ves_msk, ves = vesselness(patient_images[sl])
    segment_vessel.append(ves)
    print "Segmented Vessel Slide ", sl
    plt.imshow(ves, cmap='gray')
    plt.show()


v, f = make_mesh(ves, 30, 2)
plt_3d(v, f)
#~ m = VolumeSlicer(data=segment_vessel)
#~ m.configure_traits()
