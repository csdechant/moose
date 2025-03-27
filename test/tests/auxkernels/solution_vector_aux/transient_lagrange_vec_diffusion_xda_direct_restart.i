[Mesh]
  [geo]
    type = FileMeshGenerator
    file = 'transient_vector_diffusion_out_at_10steps_0000_mesh.xda'
  []
[]

[Variables]
  [u]
    family = LAGRANGE_VEC
  []
[]

[ICs]
  [u]
    type = VectorConstantIC
    variable = u
    x_value = 1
    y_value = 2
    z_value = 3
    block = 2
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

[AuxVariables]
  [aux_xda_u]
    order = FIRST
    family = LAGRANGE_VEC
  []
[]

[AuxKernels]
  [aux_xda_u_kernel]
    type = SolutionVectorAux
    variable = aux_xda_u
    solution = soln
    from_variable = u
    direct = true
    execute_on = 'INITIAL'
  []
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = 'transient_vector_diffusion_out_at_10steps_0000_mesh.xda'
    es = 'transient_vector_diffusion_out_at_10steps_0000.xda'
    system_variables = 'u'
    timestep = 'LATEST'
  []
[]

[Executioner]
  type = Transient
  start_time = 0.1
  num_steps = 10
  dt = 0.01
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
