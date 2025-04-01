[Mesh]
  [geo]
    type = FileMeshGenerator
    file = 'gold/lagrange_vec_solution_ic_out_0000_mesh.xda'
  []
[]

[Variables]
  [u]
    family = LAGRANGE_VEC
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

[Kernels]
  [diff]
    type = VectorDiffusion
    variable = u
  []
  [time]
    type = VectorTimeDerivative
    variable = u
  []
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = 'gold/lagrange_vec_solution_ic_out_0000_mesh.xda'
    es = 'gold/lagrange_vec_solution_ic_out_0000.xda'
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
