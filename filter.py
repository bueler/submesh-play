from firedrake import *
from firedrake.petsc import PETSc
import firedrake.cython.dmcommon as dmcommon
from firedrake.utils import IntType

import numpy as np
import sys

mesh = UnitSquareMesh(2, 2)

# create section which indicates nodes and edges
tdim = mesh.topological_dimension()
assert tdim == 2
entity_dofs = np.zeros(tdim + 1, dtype=IntType)
entity_dofs[:] = [1, 1, 0]  # we will keep the nodes and edges, but not the cells
indicatorSect, _ = dmcommon.create_section(mesh, entity_dofs)
indicatorSect.view()

# Pull Plex from mesh
dm = mesh.topology_dm

# Create a filter label
dm.createLabel("filter")
adaptLabel = dm.getLabel("filter")
adaptLabel.setDefaultValue(0)

# Set label values with function array
DGT1 = FunctionSpace(mesh, "DGT", 1)
indicator = Function(DGT1).interpolate(Constant(1.0))
print(indicator.dat.data)
print(len(indicator.dat.data))
print(DGT1.dim())
dmcommon.mark_points_with_function_array(
    dm, indicatorSect, 0, indicator.dat.data_with_halos, adaptLabel, 1
)

sys.exit(0)

# Create a DMPlexTransform object to apply the filter
opts = PETSc.Options()
opts["dm_plex_transform_active"] = "filter"
opts["dm_plex_transform_type"] = "transform_filter"
dmTransform = PETSc.DMPlexTransform().create(comm=mesh.comm)
dmTransform.setDM(dm)

# For now the only way to set the active label with petsc4py is with PETSc.Options() (DMPlexTransformSetActive() has no binding)
dmTransform.setFromOptions()
dmTransform.setUp()
dmAdapt = dmTransform.apply(dm)

# Labels are no longer needed we need to call destroy on them.
dmAdapt.removeLabel("filter")
dm.removeLabel("filter")
dmTransform.destroy()

# Remove labels to stop further distribution in mesh()
# dm.distributeSetDefault(False) <- Matt's suggestion
dmAdapt.removeLabel("pyop2_core")
dmAdapt.removeLabel("pyop2_owned")
dmAdapt.removeLabel("pyop2_ghost")
# ^ Koki's suggestion

# Pull distribution parameters from original dm
distParams = mesh._distribution_parameters

# Create a new mesh from the adapted dm
refinedmesh = Mesh(dmAdapt, distribution_parameters=distParams, comm=mesh.comm)

# Set transform type back to regular refinemenet
opts["dm_plex_transform_type"] = "refine_regular"
