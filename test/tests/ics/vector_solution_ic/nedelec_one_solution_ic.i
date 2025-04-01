[Mesh]
  [geo]
    type = FileMeshGenerator
    file = 'gold/nedelec_one_solution_ic_out_0000_mesh.xda'
  []
[]

[Variables]
  [u]
    family = NEDELEC_ONE
    order = FIRST
  []
[]

[ICs]
  [u_ic]
    type = VectorSolutionIC
    from_variable = u
    solution_uo = soln
    variable = u
  []
[]

[Functions]
  [field]
    type = ParsedVectorFunction
    expression_x = 'y*t'
    expression_y = '-x*t'
    curl_z = '-2'
  []
  [ffn_x]
    type = ParsedFunction
    expression = 'y*t'
  []
  [ffn_y]
    type = ParsedFunction
    expression = '-x*t'
  []
[]

[Kernels]
  [time]
    type = VectorTimeDerivative
    variable = u
  []

  [diff]
    type = VectorFEWave
    variable = u
    x_forcing_func = ffn_x
    y_forcing_func = ffn_y
  []
[]

[BCs]
  [top]
    type = VectorCurlBC
    curl_value = field
    variable = u
    boundary = 'left right top bottom'
  []
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = 'gold/nedelec_one_solution_ic_out_0000_mesh.xda'
    es = 'gold/nedelec_one_solution_ic_out_0000.xda'
    system_variables = 'u'
    timestep = LATEST
  []
[]

[Executioner]
  type = Transient
  start_time = 0.1
  num_steps = 1
  dt = 0.01
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  execute_on = 'INITIAL'
[]
