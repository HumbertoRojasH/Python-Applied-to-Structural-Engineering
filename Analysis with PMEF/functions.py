
def Generate_Model(W,H,t,ta,L,n,d,type):
    import cadquery as cq
    from cadquery import exporters
    (L, H, W, t, ta, n, d) = (L.value, H.value, W.value, t.value, ta.value, n.value, d.value)
    dist = L/n
    ppts = [(-dist/2-i*dist,0) for i in range(int(n))]
    pts = [
            (W/2.0, H/2.0),
            (W/2.0, (H/2.0 - t)),
            (ta/2.0, (H/2.0 - t)),
            (ta/2.0, (t - H/2.0)),
            (W/2.0, (t - H/2.0)),
            (W/2.0, H/-2.0),
            (W/-2.0, H/-2.0),
            (W/-2.0, (t - H/2.0)),
            (ta/-2.0, (t - H/2.0)),	
            (ta/-2.0, (H/2.0 - t)),
            (W/-2.0, (H/2.0 - t)),
            (W/-2.0, H/2.0),

            (W/2.0, H/2.0),]
    if type.value == "Castellated Beam":
        result = cq.Workplane("YZ").polyline(pts).close().extrude(L).faces(">Y[2]").workplane().pushPoints(ppts).polygon(6, d).cutThruAll()
    elif type.value == "Cellular Beam":
        result = cq.Workplane("YZ").polyline(pts).close().extrude(L).faces(">Y[2]").workplane().pushPoints(ppts).hole(d)
    exporters.export(result, 'object.step', exporters.ExportTypes.STEP)

def Generate_geo_file(t,ta,type,n,Min_size,Max_size):
    if type.value == "Castellated Beam":
        s1 = 14+int(n.value)*6+1
        s2 = 14+int(n.value)*6+3
    elif type.value == "Cellular Beam":
        s1 = 14+int(n.value)+1
        s2 = 14+int(n.value)+3
    code = """
SetFactory("OpenCASCADE");
Merge "object.step";

//+
Extrude {0, """+str(ta.value)+""", 0} {
  Surface{13}; Layers {1}; Recombine;
}
//+
Extrude {0, 0, """+str(t.value)+"""} {
  Surface{3}; Surface{"""+str(s2)+"""}; Surface{14}; Layers {1}; Recombine;
}
//+
Extrude {0, 0, """+"""-"""+str(t.value)+"""} {
  Surface{8}; Surface{"""+str(s1)+"""}; Surface{12}; Layers {1}; Recombine;
}
//+
Recursive Delete {
  Volume{1}; 
}
Coherence;
Mesh.MeshSizeMin = """+str(Min_size.value)+""";
Mesh.MeshSizeMax = """+str(Max_size.value)+""";
Mesh.Algorithm = 8;
Mesh.Algorithm3D = 4;
Mesh.RecombinationAlgorithm = 0;
Mesh.RecombineAll = 1;
Mesh.SubdivisionAlgorithm = 2;
Mesh 3;
Coherence Mesh;
"""
    f = open('geometria.geo', 'w')
    f.write(code)
    f.close()

def Analysis_FEM(W,H,L,E,v,de,F,fix, FS):
    W,H,L,Elas,poisson,de,F,FS= (W.value, H.value, L.value, E.value, v.value, de.value, F.value, FS.value)
    import meshio
    mesh = meshio.read("object.vtk")
    nodos = mesh.points
    connection = mesh.get_cells_type("hexahedron")

    from scipy.sparse.linalg import spsolve
    from numpy import array, zeros, round
    from pmef.pro import AssembleMatrix, AssembleVector, ApplyBC
    from pmef.pos import deform
    import pyvista as pv

    grav = array([0., 0., -1.])

    class ProblemData:
        SpaceDim = 3
        pde = "Elasticity"

    class ModelData:
        E = Elas
        v = poisson
        density = de
        selfweight = de*9.80665
        gravity = grav

    class ElementData:
        dof = 3
        nodes = 8
        noInt = 8
        type = "Brick8"

    class Mesh:
        NN = len(nodos)
        NC = len(connection)
        Nodos = nodos
        Conex = connection

    BC_data = []
    Force = W*L*F
    nodosz = nodos[round(nodos[:,2],4)==H/2.0,:]
    f_node = Force/len(nodosz)
    for i in range(Mesh.NN):
        [x,y,z]=Mesh.Nodos[i]
        if fix.value == "Start":
            if round(x,4)==0:
                BC_data.append([i, 1, 1, 0.])
                BC_data.append([i, 1, 2, 0.])
                BC_data.append([i, 1, 3, 0.])
        if fix.value == "End - Start":
            if round(x,4)==0:
                BC_data.append([i, 1, 1, 0.])
                BC_data.append([i, 1, 2, 0.])
                BC_data.append([i, 1, 3, 0.])
            if round(x,4)==L:
                BC_data.append([i, 1, 1, 0.])
                BC_data.append([i, 1, 2, 0.])
                BC_data.append([i, 1, 3, 0.])
        if round(z,4)==H/2:
            BC_data.append([i, 0, 3, f_node])

    BC_data = array(BC_data)

    N = Mesh.NN*ElementData.dof
    u = zeros(N,"float64")

    K = AssembleMatrix(Mesh,ElementData, ProblemData, ModelData, "MatrizK")

    f = AssembleVector(Mesh,ElementData, ProblemData, ModelData, "VectorF")

    [K, f] = ApplyBC(K, f, BC_data, Mesh, ElementData, ProblemData, ModelData)

    u = spsolve(K.tocsr(), f)
    defo = deform(Mesh.Nodos,u,FS)

    cells = {"hexahedron": Mesh.Conex}
    meshio.write_points_cells("result.vtk", defo, cells)

    mesh = pv.read("result.vtk")
    mesh["Ux"] = u[0::3]
    mesh["Uy"] = u[1::3]
    mesh["Uz"] = u[2::3]
    mesh_vtp = mesh.extract_surface()
    mesh.save("result.vtk")
    mesh_vtp.save("result.vtp")
    return u
