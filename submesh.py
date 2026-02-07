from firedrake import *

# create boring mesh ... in parallel presumably need to set distribution parameters?
distribution_parameters = {
    "partition": True,
    "overlap_type": (DistributedMeshOverlapType.RIDGE, 1),
}
mesh = UnitSquareMesh(10, 10, distribution_parameters=distribution_parameters)
mdim = mesh.topological_dimension()
assert mdim == 2

# the DGT0 space
DGT0 = FunctionSpace(mesh, "DGT", 0)
assert DGT0 == FunctionSpace(mesh, "HDiv Trace", 0)   # they are the same

# generate a DGT0 function
x, y = SpatialCoordinate(mesh)
f = Function(DGT0).interpolate(10.0 * (x - y)**2)  # valley along diagonal
#print(f.dat.data)

#VTKFile("bar.pvd").write(f)   # generates "NotImplementedError: Can't interpolate onto traces sorry"

# generate the skeleton submesh, i.e. nodes and edges, but not cells
#   reference:  https://github.com/firedrakeproject/firedrake/blob/5144bf85ee557d62d87d6844bd82d39ccdfc4389/tests/firedrake/submesh/test_submesh_interpolate.py#L195
facet_function = Function(DGT0).interpolate(Constant(1.))
facet_value = 999
mesh = RelabeledMesh(mesh, [facet_function], [facet_value])
subm = Submesh(mesh, mdim - 1, facet_value)

#subV = FunctionSpace(subm, "DGT", 0)  # generates "ValueError: Cannot take the trace of a 1-dim cell."

subDG0 = FunctionSpace(subm, "DG", 0)   # seems o.k.
assert DGT0.dim() == subDG0.dim()   # apparently same dimension but different dof ordering?

#fs = Function(subDG0).interpolate(f)   # generates "TypeError: cannot determine truth value of Relational: Abs(1.0*X + 1.0*Y - 1.0) < 1.0e-10"  ... don't know if cross-mesh interpolation *should* work here?

#fs = Function(subDG0).project(f)  # generates "ValueError: Mismatching cells in source (triangle) and target (interval) meshes"

# this writes an apparently valid .pvd, but with wrong values
fs = Function(subDG0)
fs.dat.data[:] = f.dat.data_ro   # COMPLETELY WRONG as a way of interpolating
#print(fs.dat.data)
VTKFile("fooBAD.pvd").write(fs)   # not what it should be

# but we *can* do cross-mesh interpolation from DG0 on mesh to DG0 on the submesh
DG0 = FunctionSpace(mesh, "DG", 0)
fDG0 = Function(DG0).interpolate(10.0 * (x - y)**2)  # valley along diagonal (regenerate f from scratch)
fs2 = Function(subDG0).interpolate(fDG0)  # cross-mesh apparently works
VTKFile("fooGOOD.pvd").write(fs2)

# get coordinates of dofs for both DGT0 and subDG0
def get_coordinates(msh, el_str):
    v = VectorFunctionSpace(msh, el_str, 0)
    X = assemble(interpolate(msh.coordinates,v))
    return X.dat.data_ro

crds_DGT0 = get_coordinates(mesh, "DGT")
#print(crds_DGT0)
crds_subDG0 = get_coordinates(subm, "DG")
#print(crds_subDG0)

# loop through each coordinate of DGT0 and assign value to correct index of DG0 function
fs3 = Function(subDG0)
for ((xi,yi),f_val) in zip(crds_DGT0,f.dat.data_ro):
    i = np.argmin(np.sqrt((xi - crds_subDG0[:,0])**2 + (yi - crds_subDG0[:,1])**2))  # due to rounding error np.where(xi==crds_subDG0[:,0] ...) doesn't work
    fs3.dat.data[i] = f_val
VTKFile("fooBETTER.pvd").write(fs3)