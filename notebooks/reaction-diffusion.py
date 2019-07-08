#!/usr/bin/env python
# coding: utf-8

# # Solving a reaction diffusion PDE on a GAMer mesh using FEniCS  
# 
# This tutorial will walk you through the process of setting up and solving a reaction diffusion partial differential equation (PDE) using a mesh imported from GAMer via FEniCS

# In[1]:


import pygamer
#import dolfin as d
import math


# ## Generating the mesh using GAMer

# In[2]:


# generate the surface mesh for a cube
smesh = pygamer.surfacemesh.sphere(4)


# In[ ]:


# Mark all faces on one side of the cube
c1 = 0
c2 = 0
for fid in smesh.faceIDs:
    isInflux = True # flag to mark faces
    indices = fid.indices()
    for idx in indices:
        v = smesh.getVertex([idx])
        #if not math.isclose(v.data().position[0],-1.0,rel_tol=1e-4):
        if v.data().position[0] > 0:
            isInflux = False
            
    if isInflux:
        fid.data().marker = 2
        c1 += 1
    else:
        fid.data().marker = 3
        c2 += 2

print(c1)
print(c2)

# In[ ]:


# construct the volumetric mesh
vmesh = pygamer.makeTetMesh([smesh], "q1.1/15a0.00006O10/2ACVY")
# qradius-edge/dihedral


# In[ ]:


# write volumetric mesh into dolfin readable format
pygamer.writeDolfin('outputTetMesh.xml', vmesh)

pygamer.writeVTK('outputTetMesh.vtk', vmesh)


# ## Setup the model problem in FEniCS

# ### Model problem: strong formulation

# Find $u(\mathbf{x},t)$ such that
# 
# \begin{align}
# \frac{\partial u}{\partial t} &= D\nabla^2 u + f ~~\text{in}~~ \Omega \\
# D(n \cdot \nabla u) &= j_\text{in} ~~\text{on}~~ \partial\Omega_\text{in} \\
# D(n \cdot \nabla u) &= j_\text{out} ~~\text{on}~~ \partial\Omega_\text{out} \\
# \partial\Omega &= \partial\Omega_\text{in} \cup \partial\Omega_\text{out}
# \end{align}
# 
# where $f=-0.5$ is a constant decay rate, $D=5.0$ is the diffusion coefficient, and $j_\text{in}=3.0$ and $j_\text{out}=-1.0$ are the influx/outflux rates for the Neumann boundary conditions.

# In[3]:

import dolfin as d

# Model parameters
D    = 0.01 # diffusion coefficient
u0   = d.Constant(0.0) # initial condition
dt   = 0.05 # time step size
T    = 5.0 # final time
jin  = 5.0 # influx rate
jout = -5.0 # outflux rate
f    = -0.0 # decay rate
t    = 0.0 # current time


# In[4]:


numsteps = round(T/dt)


# ### Model problem: weak/variational formulation

# We can multiply the governing PDE by a test function, $v$, coming from a function space $V$, integrate over the domain and apply the divergence theorem to obtain the following variational formulation:
# 
# \begin{equation}
# \int_\Omega \frac{\partial u}{\partial t}v \,dx = -\int_\Omega D\nabla u \cdot \nabla v \,dx + \int_{\partial\Omega} D(n\cdot\nabla u)v \,ds + \int_\Omega fv \,dx ~~\text{in}~~ \Omega 
# \end{equation}
# 
# We discretize the time-derivative with a backward Euler scheme for stability purposes
# $$\frac{\partial u}{\partial t} \approx \frac{u - u_n}{\Delta t}$$
# where $u_n$ represents the (known) solution computed at the previous timestep. We further simplify the variational formulation by inserting the Neumann boundary conditions and using the shorthand notation where $\langle \cdot , \cdot \rangle$ represents the inner-product over $\Omega$ and $\langle \cdot , \cdot \rangle_{\partial\Omega_\text{in}}$ represents the inner-product on the boundary $\partial\Omega_\text{in}$; terms are separated by their dependence on $u$:
# 
# \begin{equation}
# \langle u,v \rangle  + D\Delta t \langle \nabla u,\nabla v \rangle = \Delta t\langle f,v \rangle +  \Delta t \langle j_\text{in},v \rangle_{\partial\Omega_\text{in}} + \Delta t \langle j_\text{out},v \rangle_{\partial\Omega_\text{out}} + \langle u_n,v \rangle 
# \end{equation}

# In the abstract bilinear form this is written as,
# \begin{align}
# a(u,v) &= \langle u,v \rangle  + D\Delta t \langle \nabla u,\nabla v \rangle \\
# L(v) &= \Delta t\langle f,v \rangle + \Delta t \langle j_\text{in},v \rangle_{\partial\Omega_\text{in}} + \Delta t \langle j_\text{out},v \rangle_{\partial\Omega_\text{out}} + \langle u_n,v \rangle 
# \end{align}

# In[23]:


# Construct linear Lagrange function space over the mesh
dolfin_mesh = d.Mesh('outputTetMesh.xml')
#dolfin_mesh = d.Mesh('/home/jlaughli/gitrepos/spinemodel/surface_diffusion/simple-block-pm2--noflux10.xml')
#dolfin_mesh = d.UnitCubeMesh(10,10,10)
V = d.FunctionSpace(dolfin_mesh,'P',1)


J = d.Expression("jin*sin(t)",t=t,jin=jin,degree=1)

# In[24]:


# Define trial and test functions
u = d.TrialFunction(V)
v = d.TestFunction(V)
# Known solution corresponding to previous timestep. Initialize with initial condition
un = d.interpolate(u0,V)


# Earlier we marked faces on the influx boundary with 2 and faces on the outflux boundaries as 3

# In[25]:

# define integration domains
meshBoundary = d.MeshFunction("size_t",dolfin_mesh,2,dolfin_mesh.domains())
dx = d.Measure('dx',domain=dolfin_mesh)
ds = d.Measure('ds',domain=dolfin_mesh,subdomain_data=meshBoundary)




# define bilinear and linear forms
F = u*v*dx - un*v*dx + D*dt*d.inner(d.grad(u), d.grad(v))*dx - dt*f*u*v*dx #- dt*jout*u*v*ds(2)
#a = u*v*dx + D*dt*d.inner(d.grad(u),d.grad(v))*dx - dt*f*u*v*dx #- dt*jout*u*v*ds(3)
#L = dt*jin*v*ds(2) + un*v*dx #dt*jout*v*ds(3) 
a,L = d.lhs(F), d.rhs(F)


# ### Solve and store solution

# In[26]:


vtkfile = d.File('solution/dolfinOut.pvd')

# store the initial condition
u = d.Function(V)
#u = un
u.assign(un)
vtkfile << (u,t)


# In[27]:

bc = d.DirichletBC(V,"2",2)
#bc=[]

# Main loop
for idx in range(numsteps):
    # step forward in time
    t += dt 
    J.t = t 
    
    d.solve(a==L,u,bcs=bc)
    
    # save solution to file
    vtkfile << (u,t)
    
    print(u.vector().get_local().max())
    #print(dir(u))
    # assign just computed solution to un
    un.assign(u)
    
    print('Time step %d complete ...' % int(idx+1))


# In[ ]:




