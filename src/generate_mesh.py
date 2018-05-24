import morphic
reload (morphic)
import numpy
import os, sys
import csv
import matplotlib
# matplotlib.use('wx')
# matplotlib.use('Qt4Agg')

class morphic_generate_mesh:

    def __init__(self):
        self.mesh = morphic.Mesh()
        self.count = 0
        self.elements = None
        self.file_name = None

    def generate_mesh(self, file_name, elements, save=True):

        if self.mesh is not None:
            self.mesh = None
        self.mesh = morphic.Mesh()

        data = {}

        if self.elements is not None:
            self.elements = None

        self.elements = elements

        if self.file_name is not None:
            self.file_name = None
        self.file_name = file_name
        self.count = 0

        node_list = []
        node_ind = []
        with open (self.file_name, 'r') as csvfile:
            data[self.file_name] = csv.reader (csvfile, delimiter=' ', quotechar='|')
            for rowx in data[self.file_name]:
                rowy = data[self.file_name].next ()
                rowz = data[self.file_name].next ()
                node = [
                    [float (rowx[1]), float (rowx[2]), float (rowx[3]), float (rowx[4])],
                    [float (rowy[1]), float (rowy[2]), float (rowy[3]), float (rowy[4])],
                    [float (rowz[1]), float (rowz[2]), float (rowz[3]), float (rowz[4])]
                ]

                node_list.append(node)
                node_ind.append(rowx[0])

        for idx, nd in enumerate(node_list):
            self.mesh.add_stdnode(node_ind[idx], nd, group='_default')

        for ii, elem in enumerate(self.elements):
            self.mesh.add_element(ii + 1, ['H3', 'H3'], elem)

        self.mesh.generate()

        if save:
            mesh_save = os.path.normpath(self.file_name + os.sep + os.pardir)
            self.mesh.save(mesh_save + '/' 'BV.mesh')


def visualise_mesh(mesh, fig, visualise, face_colours):
    mid = mesh.label
    Xn = mesh.get_nodes (group='_default')
    Xnid = mesh.get_node_ids (group='_default')

    if visualise:
        # View embryo heart surface mesh
        Xs, Ts = mesh.get_surfaces (res=45)
        # Xl = mesh.get_lines (res=32)

        fig.plot_surfaces('{0}_Faces'.format(mid), Xs, Ts, color=face_colours, opacity=1.)


if __name__ == "__main__":

    elem = [
        ['1', '2', '22', '23'], ['2', '3', '23', '24'], ['3', '4', '24', '25'], ['4', '1', '25', '22'],
        ['22', '23', '5.2', '6.1'], ['23', '24', '6.6', '7.2'], ['24', '25', '7.1', '8.6'], ['25', '22', '8.5', '5.2'],
        ['6.2', '5.1', '11', '10'], ['9.1', '6.3', '12', '11'], ['8.3', '9.1', '13', '12'], ['5.1', '8.4', '10', '13'],
        ['9.2', '6.4', '14', '15'], ['6.5', '7.2', '15', '16'], ['7.1', '8.1', '16', '17'], ['8.2', '9.2', '17', '14'],
        ['11', '10', '19', '18'], ['12', '11', '20', '19'], ['13', '12', '21', '20'], ['10', '13', '18', '21'],
        ['15', '16', '34', '31'], ['16', '17', '31', '32'], ['17', '14', '32', '33'], ['14', '15', '33', '34']
    ]

    ipnode = 'BV.morphic'
    myMesh = morphic_generate_mesh()
    path = '/people/bsha219/lung/Data/Human_PE_Study_HRC/ST12/TLC/Vessel/SurfaceFEMesh/'
    myMesh.generate_mesh('BV.mesh', elem, save=True)

    # mesh2D = morphic.Mesh ('BV.mesh')
    #
    # """ PLOT MESH """
    # visualise = True
    # offset = 0
    # if visualise:
    #     from morphic import viewer
    #     if "fig" not in locals():
    #         fig = viewer.Figure(bgcolor=(1.,1.,1.))
    #
    #     rc = 233.
    #     gc = 150.
    #     bc = 122.
    #     rgb_color = (rc/255., gc/255., bc/255.)
    #     visualise_mesh(mesh2D, fig, visualise, face_colours=rgb_color)