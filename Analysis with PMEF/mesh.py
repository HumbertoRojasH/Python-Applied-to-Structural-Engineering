import gmsh
import pyvista as pv
gmsh.initialize()
gmsh.open("geometria.geo")
gmsh.write("object.vtk")
gmsh.finalize()

import meshio
mesh = meshio.read("object.vtk")
nodos = mesh.points
connection = mesh.get_cells_type("hexahedron")
cells = {"hexahedron": connection}
meshio.write_points_cells("object.vtk", nodos, cells)
