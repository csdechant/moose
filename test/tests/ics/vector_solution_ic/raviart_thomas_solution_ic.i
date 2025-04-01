k = asin(1)

[Mesh]
  [geo]
    type = FileMeshGenerator
    file = 'gold/raviart_thomas_solution_ic_out_0000_mesh.xda'
  []
[]

[Variables]
  [u]
    family = RAVIART_THOMAS
    order = FIRST
  []
  [p]
    family = MONOMIAL
    order = CONSTANT
  []
  [lambda]
    family = SCALAR
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
  [f]
    type = ParsedVectorFunction
    expression_x =  ${k}*sin(${k}*x)*sin(${k}*y)*t
    expression_y = -${k}*cos(${k}*x)*cos(${k}*y)*t
    div = 2*${k}^2*cos(${k}*x)*sin(${k}*y)*t
  []
[]

[Kernels]
  [time]
    type = VectorTimeDerivative
    variable = u
  []
  
  [coefficient]
    type = VectorFunctionReaction
    variable = u
    sign = negative
  []
  [gradient]
    type = GradField
    variable = u
    coupled_scalar_variable = p
  []
  [divergence]
    type = DivField
    variable = p
    coupled_vector_variable = u
  []
  [forcing]
    type = BodyForce
    variable = p
    function = ${Functions/f/div}
  []
  [mean_zero_p]
    type = ScalarLagrangeMultiplier
    variable = p
    lambda = lambda
  []
[]

[ScalarKernels]
  [constraint]
    type = AverageValueConstraint
    variable = lambda
    pp_name = PP
    value = 0.0
  []
[]

[BCs]
  [sides]
    type = VectorDivPenaltyDirichletBC
    variable = u
    function = f
    penalty = 1e8
    boundary = 'top bottom left right'
  []
[]

[Postprocessors]
  [PP]
    type = ElementIntegralVariablePostprocessor
    variable = p
    execute_on = linear
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[UserObjects]
  [soln]
    type = SolutionUserObject
    mesh = 'gold/raviart_thomas_solution_ic_out_0000_mesh.xda'
    es = 'gold/raviart_thomas_solution_ic_out_0000.xda'
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
  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
[]

[Outputs]
  exodus = true
  execute_on = 'INITIAL'
[]