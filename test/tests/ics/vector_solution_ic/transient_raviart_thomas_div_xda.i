k = asin(1)

[Mesh]
  [gmg]
    type = GeneratedMeshGenerator
    dim = 2
    nx = 10
    ny = 10
    xmax =  1
    ymax =  1
    xmin = -1
    ymin = -1
    elem_type = QUAD9
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

[Executioner]
  type = Transient
  num_steps = 10
  dt = 0.01
  solve_type = NEWTON
  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
[]

[Outputs]
  exodus = true
  file_base = 'raviart_thomas_solution_ic_out'
  execute_on = 'FINAL'
  [xda]
    type = XDA
    execute_on = 'FINAL'
  []
[]