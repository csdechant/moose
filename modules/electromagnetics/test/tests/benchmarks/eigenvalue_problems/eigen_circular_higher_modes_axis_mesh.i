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
[]

[Variables]
  [potential]
    order = FIRST
    family = LAGRANGE
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
  # alternative BCs for circle case
  [circle]
    type = DirichletBC
    variable = potential
    boundary = 'right_side left_side'
    value = 0
  []
  [eigen_circle]
    type = EigenDirichletBC
    variable = potential
    boundary = 'right_side left_side'
  []
[]

[VectorPostprocessors]
  [eigenvalues]
    type = Eigenvalues
  []
[]

# [Executioner]
#   type = Eigenvalue
# []

[Executioner]
  type = Eigenvalue
  solve_type = KRYLOVSCHUR
  # solve_type = ARNOLDI
  which_eigen_pairs = SMALLEST_MAGNITUDE
  # precond_matrix_includes_eigen = true
  n_eigen_pairs = 6
  petsc_options = '-eps_monitor_all -eps_view'
  # st_pc_type drastically changes answers
  # -st_type cayley skips some eigenvalues (including dup.), Why?
  petsc_options_iname = '-st_type -eps_target -st_ksp_type -st_pc_type'
  petsc_options_value = 'sinvert 0 preonly lu'
  eigen_tol = 1e-8
[]

[Problem]
  type = EigenProblem
  active_eigen_index = 0
[]

[Outputs]
  file_base = circle_index_test_mesh_00
  csv = true
  exodus = true
  execute_on = FINAL
[]
