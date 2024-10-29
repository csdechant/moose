# This file uses the Exodus file from transient_function_curl_NaturalBC.i for mesh and 
# initial conditions for the vector variable (First Order NEDELEC_ONE)

# NOTE: The new ICs Object 'VectorSolutionIC' does NOT read individual blocks/Domains.
#       There seems to be a bug where vector variables do not exist on points that are
#       on the boundary of blocks (at least for LAGRANGE_VEC). For that reason, 
#       'VectorSolutionIC' reads the whole mesh and places values of the vector where they exist. 
#       This is done with the 'pointValue' function in 'SolutionUserObject' 
#       while letting 'subdomain_ids = nullptr'

[Mesh]
  [geo]
    type = FileMeshGenerator
    file = 'transient_function_curl_NaturalBC_out_at_10steps.e'
    use_for_exodus_restart = true
  []
[]

[Variables]
  # u = (y*t, -x*t, 0)
  [./u]
    family = NEDELEC_ONE
    order = FIRST
  [../]
[]

[ICs]
  [u]
    type = VectorSolutionIC
    variable = u
    solution_uo = soln
    from_variable_x = u_x
    from_variable_y = u_y
  []
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = transient_function_curl_NaturalBC_out_at_10steps.e
    system_variables = 'u_x u_y'
    timestep = 'LATEST'
  []
[]

[Functions]
  # Simple "clockwise rotating" field in XY plane. curl(u) = (0, 0, -2)
  [./field]
    type = ParsedVectorFunction
    expression_x = 'y*t'
    expression_y = '-x*t'
    curl_z = '-2'
  [../]
  [./ffn_x]
    type = ParsedFunction
    expression = 'y*t'
  [../]
  [./ffn_y]
    type = ParsedFunction
    expression = '-x*t'
  [../]
[]

[Kernels]
  [./time]
    type = VectorTimeDerivative
    variable = u
  [../]
    
  [./diff]
    type = VectorFEWave
    variable = u
    x_forcing_func = ffn_x
    y_forcing_func = ffn_y
  [../]
[]

# [BCs]
#   [./top]
#     type = VectorCurlBC
#     curl_value = field
#     variable = u
#     boundary = 'left right top bottom'
#   [../]
# []

[Executioner]
  type = Transient
  num_steps = 10
  start_time = 0.1
  dt = 0.01
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
