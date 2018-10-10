=================
Centerline_Xi2.py
=================

This script is written to calculate a centerline points (including the junction point) and radii respective to those points at each point using the *.exnode and .*exelem files reading them using OpenCMISS.Zinc and respective to Xi_2 direction which in this case is towards the length of the tubes for a single bifurcation template. The centerline points are written as *.exdata file format which is compatible with CMGUI.

getCoordinatesForXi2(el_info, el_ids, xi_2_index)

Get the coordinates for the ring specified by the el_ids at the given xi_2 location using the supplied el_info.
    :param el_info: Dictionary of el_ids and a list of coordinates increasing in xi_1 fastest.
    :param el_ids: A list of element ids that describe the ring.
    :param xi_2_index: The index of the xi_2 location to extract the coordinates for.
    :return: A list of coordinates that form the ring.


After having the coordinates in all the nodes' coordinates on the rings, the point center to that specific ring is calculated as the average coordinates of all four nodes on the ring.
