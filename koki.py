from firedrake import *

def expr(m):
        x, y, z = SpatialCoordinate(m)
        return x + y**2 + z**3

degree = 3
distribution_parameters = {
    "partition": True,
    "overlap_type": (DistributedMeshOverlapType.RIDGE, 1),
}
mesh = UnitCubeMesh(8, 8, 8, distribution_parameters=distribution_parameters)
V = FunctionSpace(mesh, "HDiv Trace", 0)
facet_function = Function(V).interpolate(Constant(1.))
facet_value = 999
mesh = RelabeledMesh(mesh, [facet_function], [facet_value])
subm = Submesh(mesh, mesh.topological_dimension() - 1, facet_value)
HDivT3d = FunctionSpace(mesh, "HDiv Trace", degree)
hdivt3d = Function(HDivT3d).interpolate(expr(mesh))
DG2d = FunctionSpace(subm, "DG", degree)
dg2d = Function(DG2d).interpolate(expr(subm))
value3d_int = assemble(inner(hdivt3d('+'), hdivt3d('-')) * dS(facet_value))
value3d_ext = assemble(inner(hdivt3d, hdivt3d) * ds(facet_value))
value2d = assemble(inner(dg2d, dg2d) * dx)
assert abs(value2d - (value3d_int + value3d_ext)) < 5.e-13
DG3d = FunctionSpace(mesh, "DG", degree)
dg3d = Function(DG3d).interpolate(expr(mesh))
dg2d_ = Function(DG2d).interpolate(dg3d)
error = assemble(inner(dg2d_ - expr(subm), dg2d_ - expr(subm)) * dx)**0.5
assert abs(error) < 1.e-14
