# Test for EM module vector kernels CurlCurlField and VectorFunctionReaction
# Manufactured solution: u = y * x_hat - x * y_hat

[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = circular_axis_V03.msh
  []
  second_order = true
  # uniform_refine = 2
[]

[Variables]
  [u]
    family = NEDELEC_ONE
    # family = LAGRANGE_VEC
    order = FIRST
    eigen = true
  []
[]

[Kernels]
  [curl_curl]
    type = CurlCurlField
    variable = u
  []
  [coeff]
    type = VectorFunctionReaction
    variable = u
    function = -1.0
    extra_vector_tags = 'eigen'
  []
[]

[BCs]
  # [sides]
  #   type = VectorCurlPenaltyDirichletBC
  #   variable = u
  #   penalty = 1e8
  #   boundary = 'wall'
  # []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

# [Executioner]
#   type = Steady
#   solve_type = 'NEWTON'
#   petsc_options_iname = '-pc_type'
#   petsc_options_value = 'lu'
# []
[Executioner]
  type = Eigenvalue
  which_eigen_pairs = SLEPC_DEFAULT
  n_eigen_pairs = 5
  # n_basis_vectors = 200
  solve_type = KRYLOVSCHUR
  petsc_options_iname = '-st_type -eps_target -st_ksp_type -st_pc_type'
  petsc_options_value = 'sinvert 9.3 preonly lu'
  petsc_options = '-eps_view -eps_monitor_all'
  eigen_tol = 1e-8
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
#   petsc_options_iname = '-st_type -st_pc_type'
#   petsc_options_value = 'precond lu'
#   # eigen_tol = 1e-5
# []

[Problem]
  type = EigenProblem
  active_eigen_index = 00
[]

[Outputs]
  csv = true
  exodus = true
[]
