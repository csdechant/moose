# Base input file for eigenvalue example tests for multiple waveguide geometries
# RECTANGULAR (Default)
#     Mesh file rectangular.e based on Mesh block:
#         [Mesh]
#           [gmg]
#             type = GeneratedMeshGenerator
#             dim = 2
#             nx = 50
#             ny = 25
#             xmin = 0
#             xmax = 2
#             ymin = 0
#             ymax = 1
#             elem_type = TRI3
#           []
#         []
#     Expected analytic eigenvalue = 12.337005
#     EM Module calculated eigenvalue = 12.363806
# CIRCULAR (Mesh/file=circle.msh, BCs/active='circle eigen_circle')
#     Mesh generated using gmsh
#       radius = 1
#       center = (0, 0)
#     Expected analytic eigenvalue = 5.784025
#     EM Module calculated eigenvalue = 5.824152
# COAXIAL (Mesh/file=coaxial.msh, BCs/active='coaxial eigen_coaxial')
#     Mesh generated using gmsh with coaxial.geo
#       inner_radius = 0.125
#       outer_radius = 0.5
#       center = (0, 0)
#     Expected analytic eigenvalue = 67.108864
#     EM Module calculated eigenvalue = 68.007802

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = circular_axis_V02.msh
  []
  uniform_refine = 2
[]

[Variables]
  [potential]
    order = FIRST
    family = LAGRANGE
    eigen = true
  []
[]

[ICs]
  [u_ic]
    type = FunctionIC
    variable = 'potential'
    function = parsed_function
  []
[]

[Functions]
  [parsed_function]
    type = ParsedFunction
    expression = 'sin(pi*x)-cos(pi*y/2)'
  []
[]

[AuxVariables]
  [Ex]
    order = CONSTANT
    family = MONOMIAL
  []
  [Ey]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [diff]
    type = ADDiffusion
    variable = potential
  []
  [coeff]
    type = ADCoefReaction
    coefficient = -1
    variable = potential
    extra_vector_tags = 'eigen'
  []
[]

[AuxKernels]
  [Ex_aux]
    type = PotentialToFieldAux
    variable = Ex
    gradient_variable = potential
    sign = negative
    component = x
  []
  [Ey_aux]
    type = PotentialToFieldAux
    variable = Ey
    gradient_variable = potential
    sign = negative
    component = y
  []
[]

[BCs]
  [circle]
    type = DirichletBC
    variable = potential
    boundary = 'right_side left_side'
    value = 0
  []
  # [eigen_circle]
  #   type = EigenDirichletBC
  #   variable = potential
  #   boundary = 'right_side left_side'
  # []
[]

[VectorPostprocessors]
  [eigenvalues]
    type = Eigenvalues
  []
[]

# [Executioner]
#   type = Eigenvalue
# []

[UserObjects]
  # [pin_potential_low]
  #   type = NSPressurePin
  #   variable = potential
  #   pin_type = point-value
  #   point = '0 0.5 0'
  # []

  # NSPressure shifts all the variable solution space by the pin value at the pin local,
  # and does not enforce the value during a solve.
  # [pin_potential_high]
  #   type = NSPressurePin
  #   variable = potential
  #   pin_type = point-value
  #   point = '0 1.5 0'
  # []
[]

# [Executioner]
#   type = Eigenvalue
#   solve_type = JACOBI_DAVIDSON
#   # solve_type = KRYLOVSCHUR
#   which_eigen_pairs = SMALLEST_MAGNITUDE
#   # precond_matrix_includes_eigen = true
#   n_eigen_pairs = 3
#   # n_basis_vectors = 100
#   petsc_options = '-eps_monitor_all -eps_view'
#   #petsc_options_iname = '-st_type -eps_target -st_pc_type -eps_mpd -eps_ncv'
#   #petsc_options_value = 'precond 1e-10 lu 12 900'
#
#   petsc_options_iname = '-st_type -st_pc_type'
#   petsc_options_value = 'precond lu'
#
#   # eigen_tol = 1e-5
# []
[Executioner]
  type = Eigenvalue
  which_eigen_pairs = SLEPC_DEFAULT
  n_eigen_pairs = 5
  solve_type = KRYLOVSCHUR
  petsc_options_iname = '-st_type -eps_target -st_ksp_type -st_pc_type'
  petsc_options_value = 'sinvert 0 preonly lu'
  petsc_options = '-eps_view -eps_monitor_all'
  eigen_tol = 1e-15
[]

[Problem]
  type = EigenProblem
  active_eigen_index = 01
[]

[Outputs]
  file_base = circle_index_test_mesh_01
  csv = true
  exodus = true
  # execute_on = FINAL
[]
