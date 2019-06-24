#!/usr/bin/env python
# coding: utf-8

# # Meshing a protein using PyGAMer
#
# This tutorial will walk you through the process of using PyGAMer to generate a tetrahedral mesh of the volume surrounding a protein. The protein of interest, PDBID 2JHO, is that of sperm whale myoglobin at 1.4Å resolution.

# Load in PyGAMer and numpy
import pygamer
print("PyGAMer version:", pygamer.__version__())
import numpy as np

# Mesh the protein of interest
mesh = pygamer.readPDB_molsurf('2jho.pdb')

# When using PyGAMer it's important to initialize the orientation of the simplices. Notably `compute_orientation()` will try to apply a self consistent set of orientations. In many cases, you may wish to ensure that the mesh normals are outward facing. This can be achieved by calling `correctNormals()` which essentially computes the volume of the bounded surface and flips the normals if the volume is negative.

# Compute the normal orientation
components, orientable, manifold = mesh.compute_orientation()
mesh.correctNormals()
print(F"The mesh has {components} components, is"
      F" {'orientable' if orientable else 'non-orientable'}, and is"
      F" {'manifold' if manifold else 'non-manifold'}.")


# Note that this mesh has two components. This is because there is some empty space inside of the protein which however is not solvent acessible. Therefore when the surface is extracted, this void space is encapsulated by a separate mesh. You can split the two surfaces accordingly as follows.

meshes = mesh.splitSurfaces()
for i, m in enumerate(meshes):
    print(F"Mesh {i} is {m.getVolume()} A^3 in volume.")
# Keep only the larger mesh
mesh = meshes[0]

# Set selection of all vertices to True so smooth will operate on them.
for v in mesh.vertexIDs:
    v.data().selected = True
# Apply 5 iterations of smoothing
mesh.smooth(max_iter=5, preserve_ridges=True, verbose=True)


print(F"The mesh currently has {mesh.nVertices} vertices, {mesh.nEdges} edges, and {mesh.nFaces} faces.")


# ### Iterative Coarsening
#
# Depending on your problem, you may wish to coarsen or decimate the mesh. Basically reduce the number of vertices. GAMer offers several functions to help you with this goal.

for i in range(5):
    # Coarsen dense regions of the mesh
    mesh.coarse_dense(rate = 2, numiter = 3)
    # Coarsen flat regions of the mesh
    mesh.coarse_flat(rate = 0.1, numiter = 3)
    mesh.smooth(max_iter = 3, preserve_ridges = True, verbose = False)
    print(F"Iteration {i}: {mesh.nVertices} vertices, {mesh.nEdges} edges, and {mesh.nFaces} faces.")

# Smooth the mesh again.
mesh.smooth(max_iter = 10, preserve_ridges = True, verbose = True)

# Set boundary markers of the mesh to 23
for faceID in mesh.faceIDs:
    faceID.data().marker = 23

# Get the root metadata
gInfo = mesh.getRoot()
gInfo.ishole = True    # Don't mesh the inside of
gInfo.marker = -1

# Center mesh at 0,0,0
center, radius = mesh.getCenterRadius()
mesh.translate(-center)


# ## Now let's construct a bounding box to represent the cytosol around the protein

# Generate a surrounding box
box = pygamer.surfacemesh.cube(5)

# Set box boundary markers to 50
for faceID in box.faceIDs:
    faceID.data().marker = 50

# Get and set box metadata
gInfo = box.getRoot()
gInfo.ishole = False
gInfo.marker = 5

# Scale the box by 2 times the radius of the protein mesh
box.scale(radius*2)

# You can write these meshes to file if you wish.
#
# `pygamer.writeOFF('2jho.off', mesh)`
#
# `pygamer.writeOFF('box.off', box)`

# ## Tetrahedralizing

# Construct a list of meshes to pass into TetGen
meshes = [mesh, box]

# Call tetgen to tetrahedralize. The string represents command line arguments passed directly to TetGen.
tetmesh = pygamer.makeTetMesh(meshes, "q1.3/10O8/7AC")

# Tetrahedral meshes can be written out to several formats. Here's an example to write out in Dolfin XML format: `pygamer.writeDolfin('outputTetMesh.xml', tetmesh)`.