import numpy as np
from mpi4py import MPI

import dolfinx.fem.petsc
import dolfinx.nls.petsc

from dolfinx import mesh, fem, io
from dolfinx.fem.petsc import NonlinearProblem
from dolfinx.nls.petsc import NewtonSolver
from petsc4py import PETSc
from ufl import (
    TestFunction, TrialFunction, derivative,
    inner, dx, grad, SpatialCoordinate
)

# Parameters
R = 1.0
n_polytrope = 0.1

# 1D mesh on [0, R]
domain = mesh.create_interval(MPI.COMM_WORLD, 100, [0.0, R])

# Function space (Lagrange P1)
V = fem.functionspace(domain, ("CG", 1))

# Initial guess for theta
theta = fem.Function(V)
theta.interpolate(lambda x: np.maximum(1.0 - x[0], 0.0)**1.5)

# Test function
v = TestFunction(V)

# Radial coordinate
r = SpatialCoordinate(domain)[0]

# Variational form
F = (r**2) * inner(grad(theta), grad(v)) * dx - (r**2) * theta**n_polytrope * v * dx

# Boundary condition: theta(R) = 0
def boundary_R(x):
    return np.isclose(x[0], R)

dofs_R = fem.locate_dofs_geometrical(V, boundary_R)
bc_R = fem.dirichletbc(PETSc.ScalarType(0.0), dofs_R, V)


def boundary_0(x):
    return np.isclose(x[0], 0.0)

dofs_0 = fem.locate_dofs_geometrical(V, boundary_0)
bc_0 = fem.dirichletbc(PETSc.ScalarType(1.0), dofs_0, V)

# Jacobian (derivative of F)
J = derivative(F, theta)

# Nonlinear problem and solver
problem = NonlinearProblem(F, theta, bcs=[bc_0, bc_R], J=J)
solver = NewtonSolver(MPI.COMM_WORLD, problem)
solver.rtol = 1e-6
solver.max_it = 100
solver.convergence_criterion = "residual"
solver.report = True

# Solve
niter, converged = solver.solve(theta)
print(f"Newton solver: iterations = {niter}, converged = {converged}")

# Save solution to VTK file
with io.VTKFile(domain.comm, "polytrope_solution.vtu", "w") as vtk:
    vtk.write_mesh(domain)
    theta.name = "theta"
    vtk.write_function(theta)


value = theta.eval(
    np.array([0.0, 0.0, 0.0], dtype=np.float64),
    np.array([dofs_0[0]], dtype=np.int32)
)[0]
print("theta(0) =", value)
