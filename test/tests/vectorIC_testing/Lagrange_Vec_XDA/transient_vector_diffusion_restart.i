# This file uses the Exodus file from transient_vector_diffusion.i for mesh and 
# initial conditions for the vector variable (First Order LAGRANGE_VEC)

# NOTE: The new ICs Object 'VectorSolutionIC' does NOT read individual blocks/Domains.
#       There seems to be a bug where vector variables do not exist on points that are
#       on the boundary of blocks (at least for LAGRANGE_VEC). For that reason, 
#       'VectorSolutionIC' reads the whole mesh and places values of the vector where they exist. 
#       This is done with the 'pointValue' function in 'SolutionUserObject' 
#       while letting 'subdomain_ids = nullptr'

[Mesh]
  [geo]
    # type = FileMeshGenerator
    # file = 'transient_vector_diffusion_out_at_10steps.e'
    # use_for_exodus_restart = true
    type = FileMeshGenerator
    file = 'transient_vector_diffusion_out_at_10steps_0010_mesh.xda'
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
  # [u]
  #   type = VectorSolutionIC
  #   variable = u
  #   solution_uo = soln
  #   from_variable = u
  # []
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
  [uold_x]
    order = FIRST
    family = LAGRANGE
  []
  [uold_y]
    order = FIRST
    family = LAGRANGE
  []

  [uold_vector]
    order = FIRST
    family = LAGRANGE_VEC
  []
[]

[AuxKernels]
  [uold_x]
    type = SolutionAux
    variable = uold_x
    solution = soln
    from_variable = u
    execute_on = 'INITIAL'
  []
  # [uold_y]
  #   type = SolutionAux
  #   variable = uold_y
  #   solution = soln
  #   from_variable = u_y
  #   execute_on = 'INITIAL'
  # []
 
  # [uold_vector]
  #   type = SolutionVectorAux
  #   variable = uold_vector
  #   solution = soln
  #   from_variable = u
  #   execute_on = 'INITIAL'
  # []
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = 'transient_vector_diffusion_out_at_10steps_0010_mesh.xda'
    es = 'transient_vector_diffusion_out_at_10steps_0010.xda'
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
