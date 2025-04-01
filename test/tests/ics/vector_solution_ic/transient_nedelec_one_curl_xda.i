[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  xmin = -1
  ymin = -1
  elem_type = QUAD9
[]

[Variables]
  [u]
    family = NEDELEC_ONE
    order = FIRST
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

[Executioner]
  type = Transient
  num_steps = 10
  dt = 0.01
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  file_base = 'nedelec_one_solution_ic_out'
  execute_on = 'FINAL'
  [xda]
    type = XDA
    execute_on = 'FINAL'
  []
[]
