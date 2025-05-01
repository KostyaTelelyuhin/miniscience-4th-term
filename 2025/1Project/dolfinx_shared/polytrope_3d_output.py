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
n_polytrope = 3.0

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

# === Generate 3D star model from 1D profile ===

from scipy.interpolate import interp1d
import gmsh
import dolfinx.io



# Получаем радиусы и значения θ
r_vals = domain.geometry.x[:, 0]
theta_vals = theta.x.array

# Интерполируем θ(r)
theta_interp = interp1d(r_vals, theta_vals, kind="cubic", fill_value=0.0, bounds_error=False)

# Создаём сферическую 3D-сетку радиуса R


gmsh.initialize()
gmsh.model.add("star")

# 1. Создаем геометрию — сферу
sphere_tag = gmsh.model.occ.addSphere(0.0, 0.0, 0.0, R)

# 2. Синхронизируем геометрию
gmsh.model.occ.synchronize()

# 3. Добавляем физический объем (иначе dolfinx не сможет импортировать)
gmsh.model.addPhysicalGroup(3, [sphere_tag], tag=1)

# 4. Генерируем 3D сетку
gmsh.model.mesh.generate(3)

# 5. Преобразуем Gmsh → DOLFINx
mesh3d, _, _ = dolfinx.io.gmshio.model_to_mesh(
    gmsh.model, MPI.COMM_WORLD, 0, gdim=3
)



# Пространство функций в 3D
V3 = fem.functionspace(mesh3d, ("CG", 1))
# Координаты узлов, расстояние до центра
coords = mesh3d.geometry.x
r3 = np.linalg.norm(coords, axis=1)

# theta(r)^n → плотность
theta_r = theta_interp(r3)
rho_vals = theta_r**n_polytrope

# Записываем в функцию
rho = fem.Function(V3)
rho.x.array[:] = rho_vals
rho.name = "density"

# Сохраняем 3D модель
with io.VTKFile(mesh3d.comm, "polytrope_star_3d.vtu", "w") as vtk:
    vtk.write_mesh(mesh3d)
    vtk.write_function(rho)

print("3D model saved to polytrope_star_3d.vtu ✅")


value = theta.eval(
    np.array([0.0, 0.0, 0.0], dtype=np.float64),
    np.array([dofs_0[0]], dtype=np.int32)
)[0]
print("theta(0) =", value)
gmsh.finalize()

