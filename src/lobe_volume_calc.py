import pyvista as pv

# lobes = ['RUL', 'RML', 'RLL', 'LUL', 'LLL']

# my_mesh = pv.read('/hpc/bsha219/lung/Data/CTEPH/Alfred1/Post/Lobe/Lobe_mesh/RUL_surf.ply')
# volume = my_mesh.volume
# print(volume)


# Load the PLY files
mesh1 = pv.read('/hpc/bsha219/lung/Data/CTEPH/Alfred8/Post/Lobe/Lobe mesh/RUL_surf.ply')
mesh2 = pv.read('/hpc/bsha219/lung/Data/CTEPH/Alfred8/Post/Lobe/Lobe mesh/RML_surf.ply')
mesh3 = pv.read('/hpc/bsha219/lung/Data/CTEPH/Alfred8/Post/Lobe/Lobe mesh/RLL_surf.ply')
mesh4 = pv.read('/hpc/bsha219/lung/Data/CTEPH/Alfred8/Post/Lobe/Lobe mesh/LUL_surf.ply')
mesh5 = pv.read('/hpc/bsha219/lung/Data/CTEPH/Alfred8/Post/Lobe/Lobe mesh/LLL_surf.ply')

# Create a MultiBlock object and add the meshes to it
multi_mesh = pv.MultiBlock()
multi_mesh.append(mesh1)
multi_mesh.append(mesh2)
multi_mesh.append(mesh3)
multi_mesh.append(mesh4)
multi_mesh.append(mesh5)

# Create a plotter and add the MultiBlock object to it
plotter = pv.Plotter()
plotter.add_mesh(multi_mesh)

# Display the plot
plotter.show()
