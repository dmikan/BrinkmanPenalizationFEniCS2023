##Begin importing necessary libraries
from dolfin import *
import dolfin
import numpy as np

## Print log messages only from the root process in parallel
dolfin.parameters["std_out_all_processes"] = False;

# Define the reference mesh 
mesh_ref = Mesh('../mesh/mesh-square.xml')
print("Reference mesh imported")

# Print number of elements an plot
print(f"Elements: {mesh_ref.num_cells()}" )

# Create mesh and define function space
mesh_penal = Mesh('../mesh/mesh-square-penal.xml')
print("Penalized mesh imported")

# Print number of elements an plot
print(f"Elements: {mesh_penal.num_cells()}" )

# Define function spaces
V_ref_elem =  dolfin.VectorElement("CG", mesh_ref.ufl_cell(), 2)
Q_ref_elem =  dolfin.FiniteElement("CG", mesh_ref.ufl_cell(), 1)
W_ref_elem = dolfin.MixedElement([V_ref_elem, Q_ref_elem])
W_ref = dolfin.FunctionSpace(mesh_ref, W_ref_elem)
K_ref = dolfin.VectorFunctionSpace(mesh_ref,'DG',1) #vorticity in japanese is Kado

# Define trial and test functions
(v_ref, q_ref) = dolfin.TestFunctions(W_ref)
w_ref = dolfin.Function(W_ref)
(u_ref, p_ref) = dolfin.split(w_ref)
k_ref = dolfin.Function(K_ref)

# Define function spaces
V_penal_elem =  dolfin.VectorElement("CG", mesh_penal.ufl_cell(), 2)
Q_penal_elem =  dolfin.FiniteElement("CG", mesh_penal.ufl_cell(), 1)
W_penal_elem = dolfin.MixedElement([V_penal_elem, Q_penal_elem])
W_penal = dolfin.FunctionSpace(mesh_penal, W_penal_elem)
K_penal = dolfin.VectorFunctionSpace(mesh_ref,'DG',1) #vorticity in japanese is Kado

# Define trial and test functions
(v_penal, q_penal) = dolfin.TestFunctions(W_penal)
w_penal = dolfin.Function(W_penal)
(u_penal, p_penal) = dolfin.split(w_penal)
k_penal = dolfin.Function(K_penal)

# Set the fluid parameters 
U0 = 1
SQUARE_SIDE_HALF_LENGTH = 0.15
SQUARE_CENTRE_X = 0.5
SQUARE_CENTRE_Y = 0.5
Re = 80
nu = U0 * (2 * SQUARE_SIDE_HALF_LENGTH) / Re

# print the viscocity
print("Viscocity = {}".format(nu))

# No-slip boundary condition for velocity
bcu_noslip_top_wall_ref = dolfin.DirichletBC(W_ref.sub(0), dolfin.Constant((0, 0)), "on_boundary && x[1] > 1.0 - DOLFIN_EPS")
bcu_noslip_bottom_wall_ref = dolfin.DirichletBC(W_ref.sub(0), dolfin.Constant((0, 0)), "on_boundary && x[1] <  DOLFIN_EPS")
bcu_noslip_hole_ref = DirichletBC(W_ref.sub(0), dolfin.Constant((0, 0)), "on_boundary && std::abs(x[0] - 0.5) <= 0.15 + DOLFIN_EPS && std::abs(x[1] - 0.5) <= 0.15 + DOLFIN_EPS")

# Inflow boundary condition for velocity
inlet_expression  = dolfin.Expression(('4.0 * x[1] * (1.0 - x[1])', 0), degree = 2)
bcu_inflow_ref = dolfin.DirichletBC(W_ref.sub(0), inlet_expression, "x[0] <  DOLFIN_EPS")

# Outflow boundary condition for pressure
bcp_outflow_ref = dolfin.DirichletBC(W_ref.sub(1), dolfin.Constant(0), "x[0] > 4.0 - DOLFIN_EPS")

# Collect boundary conditions
bc_ref = [bcu_inflow_ref, bcu_noslip_top_wall_ref, bcu_noslip_bottom_wall_ref, bcu_noslip_hole_ref, bcp_outflow_ref]

# No-slip boundary condition for velocity
bcu_noslip_top_wall_penal = dolfin.DirichletBC(W_penal.sub(0), dolfin.Constant((0, 0)), "on_boundary && x[1] > 1.0 - DOLFIN_EPS")
bcu_noslip_bottom_wall_penal = dolfin.DirichletBC(W_penal.sub(0), dolfin.Constant((0, 0)), "on_boundary && x[1] <  DOLFIN_EPS")

# Inflow boundary condition for velocity
inlet_expression  = dolfin.Expression(('4.0 * x[1] * (1.0 - x[1])', 0), degree = 2)
bcu_inflow_penal = dolfin.DirichletBC(W_penal.sub(0), inlet_expression, "x[0] <  DOLFIN_EPS")

# Outflow boundary condition for pressure
bcp_outflow_penal = dolfin.DirichletBC(W_penal.sub(1), dolfin.Constant(0), "x[0] > 4.0 - DOLFIN_EPS")

# Collect boundary conditions
bc_penal = [bcu_inflow_penal, bcu_noslip_top_wall_penal, bcu_noslip_bottom_wall_penal, bcp_outflow_penal]

## Define some constants and functions for the penalized problem.

# Specific coefficients for the penalized problem
funcstring = f"""(x[0] - {SQUARE_CENTRE_X} >= - {SQUARE_SIDE_HALF_LENGTH}  
                        && x[0] - {SQUARE_CENTRE_X} <= {SQUARE_SIDE_HALF_LENGTH})? 
                         ((x[1] - {SQUARE_CENTRE_Y} >= - {SQUARE_SIDE_HALF_LENGTH} 
                        && x[1] - {SQUARE_CENTRE_Y} <= {SQUARE_SIDE_HALF_LENGTH})? 
                        1.0 : 0.0) : 0.0;"""
chi = dolfin.Expression(funcstring, degree=0)

## Define the variational problem for the Reference problem. 

def set_variational_ref(u, v, p, q, Re):
    F = dot(dot(u, nabla_grad(u)), v)*dx \
        + nu*inner(grad(u), grad(v))*dx \
        - div(v)*p*dx \
        - q*div(u)*dx
    return F

## Define the variational problem for the Brinkman problem.

def set_variational_penal(u, v, p, q, eta, chi, Re):
    F = dot(dot(u, nabla_grad(u)), v)*dx \
        + nu*inner(grad(u), grad(v))*dx \
        + (chi/eta)*inner(u, v)*dx \
        - div(v)*p*dx \
        - q*div(u)*dx
    return F

## Define the iterative method to solve the nonlinear problem

def newton_solver(F, w, bc, W):
    dw = dolfin.TrialFunction(W)
    dF = dolfin.derivative(F, w)
    nsproblem = dolfin.NonlinearVariationalProblem(F, w, bc, dF)
    solver = dolfin.NonlinearVariationalSolver(nsproblem)
    return solver

# Calculate vorticity from velocity
def compute_curl(u, k):        
    curl_u = dolfin.curl(u)
    #magnitude_curl_u = dolfin.sqrt(dolfin.inner(curl_u, curl_u)) #type:ufl
    #magnitude_curl_u_proj = dolfin.project(magnitude_curl_u) #typee:function
    curl_u_proj = dolfin.project(curl_u)
    k.assign(curl_u_proj)
    return k

## Iteration 
def iterate(w, solver, k, ufile, pfile, kfile):
    num_simulation = 5
    n = 1
    for i in range(num_simulation):
        if n % 1 == 0:
            if MPI.rank(MPI.comm_world) == 0:
                print('Solving weak form n =', n, flush=True)
        solver.solve()
        # Save the last one (steady solution of interest)
        if i == num_simulation - 1:  
            u, p = w.split()
            vorticity = compute_curl(u, k)
            u.rename("u", "Velocity")
            ufile.write(u, n)
            pfile.write(p, n)
            kfile.write(vorticity, n)
        n = n + 1
    ufile.close()
    pfile.close()
    kfile.close()
    u_proj = dolfin.project(u)
    return u_proj

## Solve the reference problem an iterate eta for different values, calculate L2 error and plot solutions.

# Create directories for the reference solutions
ufile_ref = XDMFFile('results/results_ref/velocity.xdmf')
ufile_ref.parameters['rewrite_function_mesh'] = False
ufile_ref.parameters['flush_output'] = True
pfile_ref = XDMFFile('results/results_ref/pressure.xdmf')
pfile_ref.parameters['rewrite_function_mesh'] = False
pfile_ref.parameters['flush_output'] = True
kfile_ref = XDMFFile('results/results_ref/vorticity.xdmf')
kfile_ref.parameters['rewrite_function_mesh'] = False
kfile_ref.parameters['flush_output'] = True

# Print message
print('\nComputing the reference solution...')

# Set the reference problem
F_ref = set_variational_ref(u_ref, v_ref, p_ref, q_ref, Re)
solver_ref = newton_solver(F_ref, w_ref, bc_ref, W_ref)
u1_ref = iterate(w_ref, solver_ref, k_ref, ufile_ref, pfile_ref, kfile_ref)

# Calculate L2 norm in Ω_f (ref)
diff_square = (u1_ref)**2*dx
norm_L2_omega_f_ref = np.sqrt(np.abs(dolfin.assemble(diff_square))) #norm_L2_omega_f (ref) 

# Output filename
output_filename = "results/error_data.txt"
with open(output_filename, 'w') as out_txt:
    out_txt.write('--------------------------------------------------------\n\n')
    out_txt.write('PARAMETER OF SIMULATION:\n')
    out_txt.write('REYNOLDS NUMBER = {}\n'.format(Re))
    out_txt.write('L2 NORM IN OMEGA F OF THE REFERENCE VELOCITY = %.2e\n' % norm_L2_omega_f_ref)
    out_txt.write('--------------------------------------------------------\n\n')
    out_txt.write('FINITE ELEMENT PARAMETERS:\n')
    out_txt.write('ELEMENTS OF REFERENCE MESH = {}\n'.format(mesh_ref.num_cells()))
    out_txt.write('ELEMENTS OF PENALIZED MESH = {}\n'.format(mesh_penal.num_cells()))
    out_txt.write('H_MIN AND H_MAX OF REFERENCE MESH = %.2e, %.2e\n' % (mesh_ref.hmin(), mesh_ref.hmax()))
    out_txt.write('H_MIN AND H_MAX OF PENALIZED MESH = %.2e, %.2e\n' % (mesh_penal.hmin(), mesh_penal.hmax()))
    out_txt.write('VELOCITY FINITE ELEMENT SPACE = {}, {}\n'.format(W_ref.sub(0).ufl_element().family(), W_ref.sub(0).ufl_element().degree()))
    out_txt.write('PRESSURE FINITE ELEMENT SPACE = {}, {}\n'.format(W_ref.sub(1).ufl_element().family(), W_ref.sub(1).ufl_element().degree()))
    out_txt.write('--------------------------------------------------------\n\n')
    out_txt.write('GLOSARY:\n')
    out_txt.write('ETA = PENALIZATION PARAMETER\n')
    out_txt.write('L2NOS = L2 NORM IN OMEGA_S\n')
    out_txt.write('L2EOF = L2 ERROR IN OMEGA_F\n')
    out_txt.write('L2EOF-REL = RELATIVE L2 ERROR IN OMEGA_F \n')
    out_txt.write('--------------------------------------------------------\n\n')
    out_txt.write('ETA\t\t')
    out_txt.write('L2NOS\t\t')
    out_txt.write('L2EOF\t\t')
    out_txt.write('L2EOF-REL\n')

eta_values = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8]
count = 1
for eta_brinkman in eta_values:

    # Set penalized parameter
    eta = dolfin.Constant(eta_brinkman)

    # Create directories for the reference and penalized solutions
    ufile_penal = XDMFFile('results/results_{}/velocity.xdmf'.format(eta_brinkman))
    ufile_penal.parameters['rewrite_function_mesh'] = False
    ufile_penal.parameters['flush_output'] = True
    pfile_penal = XDMFFile('results/results_{}/pressure.xdmf'.format(eta_brinkman))
    pfile_penal.parameters['rewrite_function_mesh'] = False
    pfile_penal.parameters['flush_output'] = True
    kfile_penal = XDMFFile('results/results_{}/vorticity.xdmf'.format(eta_brinkman))
    kfile_penal.parameters['rewrite_function_mesh'] = False
    kfile_penal.parameters['flush_output'] = True

    # Print message
    print(f'\nComputing... {count} cycles of {len(eta_values)} for eta = {eta_brinkman}')
    
    # Set the penalized problem
    F_penal = set_variational_penal(u_penal, v_penal, p_penal, q_penal, eta, chi, Re)
    solver_penal = newton_solver(F_penal, w_penal, bc_penal, W_penal)
    u1_penal = iterate(w_penal, solver_penal, k_penal, ufile_penal, pfile_penal, kfile_penal)

    #project penalized solution into the reference domain 
    u1_penal_proj = dolfin.interpolate(u1_penal, W_ref.sub(0).collapse())

    # Calculate L2 error in Ω_f
    diff_square = (u1_ref - u1_penal_proj)**2*dx
    error_L2_omega_f = np.sqrt(np.abs(dolfin.assemble(diff_square))) 

    # Calculate L2 norm in Ω_s
    diff_square = (chi*u1_penal)**2*dx
    error_L2_omega_s = np.sqrt(np.abs(dolfin.assemble(diff_square))) 

    # Relative error 
    error_L2_omega_f_relative = error_L2_omega_f/norm_L2_omega_f_ref
    
    # Save the values in error_data.txt
    with open(output_filename, 'a') as out_txt:
        out_txt.write('%.0e\t' % eta_brinkman)
        out_txt.write('%.2e\t' % error_L2_omega_s)
        out_txt.write('%.2e\t' % error_L2_omega_f)
        out_txt.write('%.2e\t' % error_L2_omega_f_relative)
        out_txt.write('\n')
    count += 1
if MPI.rank(MPI.comm_world) == 0:
    print('Done!', flush=True)
