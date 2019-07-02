# import
import dolfin as d

# test mesh
#d.File('testmesh.xml') << d.UnitCubeMesh(10,10,10)

# import the tetrahedral mesh from gamer
gmesh = d.Mesh('testmesh.xml')

# file to store solution to
vtkfile = d.File('dolfinOut.pvd')

# model parameters
D = 5 # diffusion coefficient
u0 = d.Constant(0.0) # initial condition
dt = 0.1 # time step size
T = 5.0 # final time
j = 3.0 # influx rate
f = -0.5 # decay rate


# Construct linear Lagrange function spaces
V = d.FunctionSpace(gmesh,'Lagrange',1)

u = d.TrialFunction(V)
v = d.TestFunction(V)

un = d.interpolate(u0,V)

# funcitonal
a = u*v*dx + D*dt*grad(u)*grad(v)*dx
L = f*v*dx + dt*j*v*ds

a = d.lhs(F) # left hand side (bilinear form)
L = d.rhs(F) # right hand side (linear form)

# store initial condition
t = 0.0 # current time
u = d.Function(V)
u = un
vtkfile << (u,t)

# Main loop
#u = d.Function(V)
#for 
