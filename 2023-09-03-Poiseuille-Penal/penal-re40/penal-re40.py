'''
Solving the Incompressible Navier-Stokes Equation with Obstacle using FEniCS

In this script, we solve the incompressible Navier-Stokes equation in a channel with a square obstacle using two different approaches. The first approach utilizes a reference mesh, while the second approach employs a Cartesian mesh with a penalized term in the equations to model the obstacle, The penalized solution converges to the reference solution as parameter eta tends to zero.

References

1. Hans Petter Langtangen, Anders Logg, "Solving PDEs in Python – The FEniCS Tutorial Volume" (2017).
2. Philippe Angot, Bruneau, "A Penalization Method to Take into Account Obstacles in Incompressible Viscous Flows" (1999).
3. Katuhiko Goda, "A Multistep Technique with Implicit Difference Schemes for Calculating Two- or Three-Dimensional Cavity Flows" (1978).

We utilize the first reference to understand the mathematical theorems and convergence proofs. The second reference provides insights into the FEniCS library, while the third reference explains the variation formulation for the Incompressible Navier-Stokes equation using Chorin's scheme for time discretization.
'''

##Begin importing necessary libraries
from dolfin import *
import dolfin
import numpy as np

## Print log messages only from the root process in parallel
dolfin.parameters["std_out_all_processes"] = False;

## We are going to define the mesh of the exact solution that we are going to call mesh_ref because it woulb be the reference mesh for our solution. This mesh was generated in Gmsh using triangular elements. 

# Define the reference mesh 
mesh_ref = Mesh('../mesh/mesh-square.xml')
print("Reference mesh imported")

# Print number of elements an plot
print(f"Elements: {mesh_ref.num_cells()}" )

## Now we define a rectangular mesh with triangular elements. This mesh will be used as a domain for the Brinkman problem. 

# Create mesh and define function space
mesh_penal = dolfin.RectangleMesh(dolfin.Point(0.0, 0.0), dolfin.Point(4.0, 1.0), 320, 80, "right/left")
print("Penalized mesh created")

# Print number of elements an plot
print(f"Elements: {mesh_penal.num_cells()}" )

## Define the finite element spaces that we are going to use for the reference problem.

# Define function spaces (P2-P1)
V_ref = dolfin.VectorFunctionSpace(mesh_ref, "Lagrange", 2)
Q_ref = dolfin.FunctionSpace(mesh_ref, "Lagrange", 1)
W_ref = dolfin.VectorFunctionSpace(mesh_ref,'DG',1)

# Define trial and test functions
u_ref = dolfin.TrialFunction(V_ref)
p_ref = dolfin.TrialFunction(Q_ref)
v_ref = dolfin.TestFunction(V_ref)
q_ref = dolfin.TestFunction(Q_ref)

## Define the finite element spaces that we are going to use for the approximated problem.

# Define function spaces (P2-P1)
V_penal = dolfin.VectorFunctionSpace(mesh_penal, "Lagrange", 2)
Q_penal = dolfin.FunctionSpace(mesh_penal, "Lagrange", 1)
W_penal = dolfin.VectorFunctionSpace(mesh_penal,'DG',1)

# Define trial and test functions
u_penal = dolfin.TrialFunction(V_penal)
p_penal = dolfin.TrialFunction(Q_penal)
v_penal = dolfin.TestFunction(V_penal)
q_penal = dolfin.TestFunction(Q_penal)

## Set the parameters for our problem

# Set the time parameters
dt = 0.001
FINAL_TIME = 0.030
TIME_MULT = 10 

# Set the fluid parameters 
U0 = 1
SQUARE_SIDE_HALF_LENGTH = 0.15
SQUARE_CENTRE_X = 0.5
SQUARE_CENTRE_Y = 0.5
Re = 40
nu = U0 * (2 * SQUARE_SIDE_HALF_LENGTH) / Re

# print the viscocity
print("Viscocity = {}".format(nu))

## Set initial conditions
u0_initial = dolfin.Constant((0.0,0.0))
u1_initial = dolfin.Constant((0.0,0.0))
p1_initial = dolfin.Constant(0.0)

## Define the Dirichlet boundary conditions for the velocity and pressure for the reference problem. 

# No-slip boundary condition for velocity
bcu_noslip_top_wall_ref = dolfin.DirichletBC(V_ref, dolfin.Constant((0, 0)), "on_boundary && x[1] > 1.0 - DOLFIN_EPS")
bcu_noslip_bottom_wall_ref = dolfin.DirichletBC(V_ref, dolfin.Constant((0, 0)), "on_boundary && x[1] <  DOLFIN_EPS")
bcu_noslip_hole_ref = DirichletBC(V_ref, dolfin.Constant((0, 0)), "on_boundary && std::abs(x[0] - 0.5) <= 0.15 + DOLFIN_EPS && std::abs(x[1] - 0.5) <= 0.15 + DOLFIN_EPS")

# Inflow boundary condition for velocity
bcu_inflow_ref = dolfin.DirichletBC(V_ref, dolfin.Constant((U0, 0)), "x[0] <  DOLFIN_EPS")

# Outflow boundary condition for pressure
bcp_outflow_ref = dolfin.DirichletBC(Q_ref, dolfin.Constant(0), "x[0] > 4.0 - DOLFIN_EPS")

# Collect boundary conditions
bcu_ref = [bcu_inflow_ref, bcu_noslip_top_wall_ref, bcu_noslip_bottom_wall_ref, bcu_noslip_hole_ref]
bcp_ref = [bcp_outflow_ref]

## Define the Dirichlet boundary conditions for the velocity and pressure for the Brinkman problem. 

# No-slip boundary condition for velocity
bcu_noslip_top_wall_penal = dolfin.DirichletBC(V_penal, dolfin.Constant((0, 0)), "on_boundary && x[1] > 1.0 - DOLFIN_EPS")
bcu_noslip_bottom_wall_penal = dolfin.DirichletBC(V_penal, dolfin.Constant((0, 0)), "on_boundary && x[1] <  DOLFIN_EPS")

# Inflow boundary condition for velocity
bcu_inflow_penal = dolfin.DirichletBC(V_penal, dolfin.Constant((U0, 0)), "x[0] <  DOLFIN_EPS")

# Outflow boundary condition for pressure
bcp_outflow_penal = dolfin.DirichletBC(Q_penal, dolfin.Constant(0), "x[0] > 4.0 - DOLFIN_EPS")

# Collect boundary conditions
bcu_penal = [bcu_inflow_penal, bcu_noslip_top_wall_penal, bcu_noslip_bottom_wall_penal]
bcp_penal = [bcp_outflow_penal]

## Define auxiliary functions to use in the Chorin's method scheme for the time discretization for the incompressible Navier-Stokes equation (reference problem).

# Create auxiliary functions
u0_ref = dolfin.Function(V_ref) # this is the velocity at time n
u1_ref = dolfin.Function(V_ref) # it can be both, tentative velocity and velocity at n+1
p1_ref = dolfin.Function(Q_ref) # this is the pressure at time n+1
w1_ref = dolfin.Function(W_ref) # Define the vorticity function of the velocity at n+1

## Define auxiliary functions to use in the Chorin's method scheme for the time discretization for the incompressible Navier-Stokes equation (penalized problem).

# Create auxiliary functions
u0_penal = dolfin.Function(V_penal) # this is the velocity at time n
u1_penal = dolfin.Function(V_penal) # it can be both, tentative velocity and velocity at n+1
p1_penal = dolfin.Function(Q_penal) # this is the pressure at time n+1
w1_penal = dolfin.Function(W_penal) # Define the vorticity function of velocity at n+1

## Do a test to obtain an error of projecting the penalized solution into the reference domain.

p_penal_test = dolfin.Function(Q_penal)
p_ref_test = dolfin.Function(Q_ref)
p_penal_test_initial = dolfin.Constant(1.0)
p_ref_test_initial = dolfin.Constant(1.0)
p_penal_test.assign(p_penal_test_initial)
p_ref_test.assign(p_ref_test_initial)
p_penal_test_proj = dolfin.interpolate(p_penal_test, Q_ref)
diff_square_test = (p_ref_test - p_penal_test_proj)**2*dx
error_L2_test = np.sqrt(np.abs(dolfin.assemble(diff_square_test)))
print('Test of projecting penalized solution into the reference domain, error = %.3e' % error_L2_test)

## Define some constants and functions for the penalized problem.

# Define coefficients 
k = dolfin.Constant(dt)
f = dolfin.Constant((0.0,0.0))

# Specific coefficients for the penalized problem
funcstring = f"""(x[0] - {SQUARE_CENTRE_X} >= - {SQUARE_SIDE_HALF_LENGTH}  
                        && x[0] - {SQUARE_CENTRE_X} <= {SQUARE_SIDE_HALF_LENGTH})? 
                         ((x[1] - {SQUARE_CENTRE_Y} >= - {SQUARE_SIDE_HALF_LENGTH} 
                        && x[1] - {SQUARE_CENTRE_Y} <= {SQUARE_SIDE_HALF_LENGTH})? 
                        1.0 : 0.0) : 0.0;"""
chi = dolfin.Expression(funcstring, degree=0) 

## Define the variational problem in Chorin's method for the Reference problem. The fisrt step solves for the tentative non-divergence free velocity (in general). The second step is the variational problem for a Poisson equation that solves for the pressure at time m+1. The third step solves the velocity at time m+1.

def set_variational_ref(u, v, p, q, u0, u1, p1, f, Re, k):
    # Tentative velocity step
    F1 = (1/k)*inner(u - u0, v)*dx \
       +  inner(grad(u0)*u0, v)*dx \
       + (1/Re)*inner(grad(u), grad(v))*dx \
       - inner(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)

    # Pressure update
    a2 = inner(grad(p), grad(q))*dx
    L2 = -(1/k)*div(u1)*q*dx

    # Velocity update
    a3 = inner(u, v)*dx
    L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

    a = [a1, a2, a3]
    L = [L1, L2, L3]
    return a, L

## Define the variational problem in Chorin's method for the Brinkman problem. The fisrt step solves for the tentative non-divergence free velocity (in general), note that we are adding a penalized term with an implicit velocity. The second step is the variational problem for a Poisson equation that solves for the pressure at time m+1. The third step solves the velocity at time m+1. 

def set_variational_penal(u, v, p, q, u0, u1, p1, f, eta, chi, Re, k):
    # Tentative velocity step
    F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
         (1/Re)*inner(grad(u), grad(v))*dx + (chi/eta)*inner(u, v)*dx - inner(f, v)*dx
    a1 = lhs(F1)
    L1 = rhs(F1)

    # Pressure update
    a2 = inner(grad(p), grad(q))*dx
    L2 = -(1/k)*div(u1)*q*dx

    # Velocity update
    a3 = inner(u, v)*dx
    L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

    a = [a1, a2, a3]
    L = [L1, L2, L3]
    return a, L

## Assemble the matrices that we are going to use in the linear system that will solve the EDP for the reference problem.

def assemble_matrices(a1, a2, a3):
    # Assemble matrices
    A1 = assemble(a1)
    A2 = assemble(a2)
    A3 = assemble(a3)
    A = [A1, A2, A3]    
    return A

## Define a preconditioner to improve the time computation in solving the linear system.

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

## Define a progress bar.

# Define progress bar
def print_progress_bar(progress):
    bar_width = 50
    pos = int(bar_width * progress)
    bar = '[' + '=' * pos + '>' + ' ' * (bar_width - pos) + ']'
    percentage = int(progress * 100)
    print(f"{bar} {percentage}%\r", end="")
    if progress == 1.0:
        print()

## Time stepping for the reference solution and for the penalized solution:
def time_stepping(A, L, u0, u1, p1, w1, ufile, pfile, wfile, bcu, bcp):
    # Time-stepping
    t = dt
    time_step = 1
    u1_arr = []
    p1_arr = []
    sample_time_arr = []

    # Save the first step    
    ufile.write(u1, 0)
    pfile.write(p1, 0)
    
    while t < FINAL_TIME + DOLFIN_EPS:  

        # Compute tentative velocity step
        b1 = assemble(L[0])
        [bc.apply(A[0], b1) for bc in bcu]
        solve(A[0], u1.vector(), b1, "bicgstab", "default")

        # Pressure correction
        b2 = assemble(L[1])
        [bc.apply(A[1], b2) for bc in bcp]
        [bc.apply(p1.vector()) for bc in bcp]
        solve(A[1], p1.vector(), b2, "bicgstab", prec)

        # Velocity correction
        b3 = assemble(L[2])
        [bc.apply(A[2], b3) for bc in bcu]
        solve(A[2], u1.vector(), b3, "bicgstab", "default")

        # Calculate vorticity from velocity
        curl_u1 = dolfin.curl(u1)
        magnitude_curl_u1 = dolfin.sqrt(dolfin.inner(curl_u1, curl_u1))
        magnitude_curl_u1_proj = dolfin.project(magnitude_curl_u1)
        w1.assign(magnitude_curl_u1_proj)

        # Save the vorticity first step
        if t == dt:    
            wfile.write(w1, 0)

        # Take time samples
        if time_step%TIME_MULT == 0:

            # Save to file
            label_time = time_step/TIME_MULT
            ufile.write(u1, label_time)
            pfile.write(p1, label_time)
            wfile.write(w1, label_time)
            u1_arr.append(u1)  
            p1_arr.append(p1)
 
            # Take time samples
            sample_time_arr.append(round(t, 3))
            
        # Update u0
        u0.assign(u1)
        
        # Time step
        t += dt
        time_step += 1
        
        # Print progress bar
        print_progress_bar(t/(FINAL_TIME))
    ufile.close()
    pfile.close()
    wfile.close()    
    return u1_arr

## Solve the reference problem an iterate eta for different values, calculate L2 error and plot solutions.

# Create files for storing solution
ufile_ref = dolfin.XDMFFile("results/results_ref/velocity.xdmf")
pfile_ref = dolfin.XDMFFile("results/results_ref/pressure.xdmf")
wfile_ref = dolfin.XDMFFile("results/results_ref/vorticity.xdmf")

# Print message
print('\nComputing the reference solution...')

# Initialize functions
u0_ref.assign(u0_initial)
u1_ref.assign(u1_initial)
p1_ref.assign(p1_initial)

# Set the reference problem
a_ref, L_ref = set_variational_ref(u_ref, v_ref, p_ref, q_ref, u0_ref, u1_ref, p1_ref, f, Re, k)
A_ref = assemble_matrices(a_ref[0], a_ref[1], a_ref[2])
u1_arr_ref = time_stepping(A_ref, L_ref, u0_ref, u1_ref, p1_ref, w1_ref, ufile_ref, pfile_ref, wfile_ref, bcu_ref, bcp_ref)

# Output filename
output_filename = "results/error_data.txt"
with open(output_filename, 'w') as out_txt:
    out_txt.write('eta\t')
    out_txt.write('EL2TOF\t')
    out_txt.write('NL2TOS\t')
    out_txt.write('\n')

eta_values = [1e-2, 1e-8]
count = 1
for eta_brinkman in eta_values:
    
    # Create arrays
    u1_arr_penal_proj = []
    p1_arr_penal_proj = []
    error_L2_omega_f_arr = []
    norm_L2_omega_s_arr = []
    error_data = []

    # Set penalized parameter
    eta = dolfin.Constant(eta_brinkman)

    # Set time step to calculate L2 norm in time
    dt_samples = dt * TIME_MULT
    
    # Create files for storing solution
    ufile_penal = dolfin.XDMFFile("results/results_{}/velocity.xdmf".format(eta_brinkman))
    pfile_penal = dolfin.XDMFFile("results/results_{}/pressure.xdmf".format(eta_brinkman))
    wfile_penal = dolfin.XDMFFile("results/results_{}/vorticity.xdmf".format(eta_brinkman))

    # Print message
    print(f'\nComputing... {count} cycles of {len(eta_values)} for eta = {eta_brinkman}')

    # Initialize functions
    u0_penal.assign(u0_initial)
    u1_penal.assign(u1_initial)
    p1_penal.assign(p1_initial)
    
    # Set the penalized problem
    a_penal, L_penal = set_variational_penal(u_penal, v_penal, p_penal, q_penal, u0_penal, u1_penal, p1_penal, f, eta, chi, Re, k)
    A_penal = assemble_matrices(a_penal[0], a_penal[1], a_penal[2])
    u1_arr_penal = time_stepping(A_penal, L_penal, u0_penal, u1_penal, p1_penal, w1_penal, ufile_penal, pfile_penal, wfile_penal, bcu_penal, bcp_penal)

    #project penalized velocity solution on reference domain   
    for u1_penal_sample in u1_arr_penal:    
        u1_penal_proj = dolfin.interpolate(u1_penal_sample, V_ref)
        u1_arr_penal_proj.append(u1_penal_proj)

    # Calculate L2 error in Ω_f
    for index, u1_penal_proj in enumerate(u1_arr_penal_proj):
        diff_square = (u1_arr_ref[index] - u1_penal_proj)**2*dx
        error_L2_omega_f = np.sqrt(np.abs(dolfin.assemble(diff_square))) 
        error_L2_omega_f_arr.append(error_L2_omega_f)
        
    # Calculate L2 error in time in [0,T] X Ω_f 
    error_L2_omega_f_arr_square = [error ** 2 for error in error_L2_omega_f_arr]
    error_L2_time_omega_f = np.sqrt(dt_samples*np.sum(error_L2_omega_f_arr_square))

    # Calculate L2 norm in Ω_s
    for u1_penal_sample in u1_arr_penal:
        diff_square = (chi*u1_penal_sample)**2*dx
        norm_L2_omega_s = np.sqrt(np.abs(dolfin.assemble(diff_square))) 
        norm_L2_omega_s_arr.append(norm_L2_omega_s)

    # Calculate L2 norm in time in [0,T] X Ω_s
    norm_L2_omega_s_arr_square = [error ** 2 for error in norm_L2_omega_s_arr]
    norm_L2_time_omega_s = np.sqrt(dt_samples*np.sum(norm_L2_omega_s_arr_square))
    
    # Save the values in error_data.txt
    with open(output_filename, 'a') as out_txt:
        out_txt.write('%.0e\t' % eta_brinkman)
        out_txt.write('%.2e\t' % error_L2_time_omega_f)
        out_txt.write('%.2e\t' % norm_L2_time_omega_s)
        out_txt.write('\n')
    count += 1
print(f'\nData saved in {output_filename}')
    