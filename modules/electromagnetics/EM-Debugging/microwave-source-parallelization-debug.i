[Mesh]
  [fmg]
    type = FileMeshGenerator
    file = NCSU_chamber-edits-WO-lid.msh
  []
  coord_type = RZ
  rz_coord_axis = Y
[]

[Problem]
  type = FEProblem
[]

[Functions]
  # This function is used in the current density source equation.
  # In this case, the frequency is 2.45 GHz (omege = 2*pi*2.45e9)
  # and the relative magnetic permeability of stainless steel is about 1
  # (mu = 4*pi*1e-7)
  [omegaMu_steel]
    type = ParsedFunction
    expression = '2*pi*2.45e9 * 4*pi*1e-7'
  []

  # This fuction is used in the displacement current equation.
  # In this case, the frequency is 2.45 GHz (omege = 2*pi*2.45e9),
  # the relative magnetic permeability of a ceramic is about 1 (mu = 4*pi*1e-7),
  # and the relative permittivity of an aluminum oxide ceramic is 9.1 (epsilon = 8.05714e-11)
  [muEpsilonOmega2_ceramic]
    type = ParsedFunction
    expression = '4*pi*1e-7 * 8.05714e-11 * (2*pi*2.45e9)^2'
  []

  # This fuction is used in the displacement current equation.
  # In this case, the frequency is 2.45 GHz (omege = 2*pi*2.45e9),
  # and the permeability / permittivity of a vacuum is used.
  [muEpsilonOmega2_vacuum]
    type = ParsedFunction
    expression = '4*pi*1e-7 * 8.8542e-12 * (2*pi*2.45e9)^2'
  []

  #beta is the wave number
  [beta]
    type = ParsedFunction
    expression = '2*pi*2.45e9/3e8'
  []
  #the real current is 3.94 A (might need to change)
  [curr_real]
    type = ParsedVectorFunction
    expression_y = '1.5e10 * 1.02e-5'
  []
  [curr_imag] # defaults to '0.0 0.0 0.0'
    type = ParsedVectorFunction
  []
[]

[Variables]
  [E_real]
    family = NEDELEC_ONE
    order = FIRST
  []
  [E_imag]
    family = NEDELEC_ONE
    order = FIRST
  []
[]

[Kernels]
  # Calculation the E-Field developed in the port needle
  # as it acts like a current source.
  # Real component calculations
  [pin_curlCurl_real]
    type = CurlCurlField
    variable = E_real
    block = Resonator_Pin
  []
  [pin_source_real]
    type = VectorCurrentSource
    variable = E_real
    component = real
    source_real = curr_real
    source_imag = curr_imag
    function_coefficient = omegaMu_steel
    block = Resonator_Pin
  []
  # Imaginary component calculations
  [pin_curlCurl_imag]
    type = CurlCurlField
    variable = E_imag
    block = Resonator_Pin
  []
  [pin_source_imaginary]
    type = VectorCurrentSource
    variable = E_imag
    component = imaginary
    source_real = curr_real
    source_imag = curr_imag
    function_coefficient = omegaMu_steel
    block = Resonator_Pin
  []

  # Calculation the E-Field developed in the aluminum oxide ceramic.
  # Real component calculations
  [ceramic_curlCurl_real]
    type = CurlCurlField
    variable = E_real
    block = Ceramic
  []
  [ceramic_coeff_real]
    type = VectorFunctionReaction
    variable = E_real
    function = muEpsilonOmega2_ceramic
    sign = negative
    block = Ceramic
  []
  # Imaginary component calculations
  [ceramic_curlCurl_imag]
    type = CurlCurlField
    variable = E_imag
    block = Ceramic
  []
  [ceramic_coeff_imag]
    type = VectorFunctionReaction
    variable = E_imag
    function = muEpsilonOmega2_ceramic
    sign = negative
    block = Ceramic
  []

  # Calculation the E-Field developed in vacuum
  # Real component calculations
  [vacuum_curlCurl_real]
    type = CurlCurlField
    variable = E_real
    block = Plasma
  []
  [vacuum_coeff_real]
    type = VectorFunctionReaction
    variable = E_real
    function = muEpsilonOmega2_vacuum
    sign = negative
    block = Plasma
  []
  # Imaginary component calculations
  [vacuum_curlCurl_imag]
    type = CurlCurlField
    variable = E_imag
    block = Plasma
  []
  [vacuum_coeff_imag]
    type = VectorFunctionReaction
    variable = E_imag
    function = muEpsilonOmega2_vacuum
    sign = negative
    block = Plasma
  []
[]

[BCs]
  [absorbing_real]
    type = VectorEMRobinBC
    variable = E_real
    component = real
    beta = beta
    coupled_field = E_imag
    mode = absorbing
    boundary = 'metal_wall'
  []
  [absorbing_imag]
    type = VectorEMRobinBC
    variable = E_imag
    component = imaginary
    beta = beta
    coupled_field = E_real
    mode = absorbing
    boundary = 'metal_wall'
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Steady
  automatic_scaling = true
  compute_scaling_once = false
  solve_type = 'NEWTON'
  line_search = none
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
[]

[Outputs]
  exodus = true
  print_linear_residuals = true
[]

[Debug]
  show_var_residual_norms = true
[]
