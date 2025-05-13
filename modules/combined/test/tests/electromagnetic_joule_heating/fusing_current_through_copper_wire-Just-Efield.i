# This test is an attempt to make a pure electric field version of
# the test file 'fusing_current_through_copper_wire.i', where that test
# solves using magnetic vector potential. The PDEs are as follow:
#
#   curl(curl(E)) + j*mu*epsilon*omega^2*E = -j*mu*omega*J
#   rho*C*dT/dt - div(k*grad(T)) = Q
#   Q = 0.5*sigma*mag(E)^2
#
# Where:
#   - E is the electric field
#   - mu is the permeability of free space
#   - omega is the angular frequency of the system (60 Hz)
#   - epsilon is the permittivity of free space
#   - j is the sqrt(-1)
#   - J is the supplied DC current
#   - rho is the density of copper
#   - C is the heat capacity of copper
#   - T is the temperature
#   - k is the thermal conductivity of the wire
#   - Q is the Joule heating
#   - sigma is the electric conductivity
# 
# Below are the following notes:
#
#   - For the frequency of 60 Hz, the VectorEMRobinBC needs to be
#     applied to the whole surface of the wire. Otherwise, the simulations
#     will not converge given the current preconditioner options.
#
#   - Currently, the resulting electric field is too high 
#     (the real component of the electric field is ~1e7 V/m, 
#      while in the magnetic vector potential case, 
#      the real componet of the electric field is ~8 V/m)
# 
# Below are two hypothesis for the discrepancy 
# (assuming the magnetic vector potential input file is correct within assumptions):
# 
#   - There is an issue with the boundary conditions at low frequencies. 
#     In particular, a more robust VectorEMRobinBC is needed 
#     (similar to COMSOL formulation, where the conductivity of the wall is included as
#      an input parameter for the BC). It should be noted that while a more robust
#     impendance BC would be great, I don't believe this is the main source of error.
#
#   - The formulation of the electric field is ill-posed, compared to the magnetic 
#     vector potential formulation. In particular, having a supplied DC current in the 
#     electric field formulation is ill-posed. This is because the current definition is
#     J = sigma*E. By forcing the value of J and having only a temperature dependent sigma,
#     the electric field of J is not the same electric field as "curl(curl(E)) + j*mu*epsilon*omega^2*E".
#     We can supply a DC current in the magnetic vector potential formulation without issue
#     because we first seperature the electric field as E = jwA - grad(V), where the
#     DC current, J, is defined using grad(V). Since the magnetic vector potential 
#     simulation forces the value of J, but doesn't solve for grad(V) directly, 
#     we avoid the ill-posing the electric field.

[Mesh]
  # Mesh of the copper wire
  [fmg]
    type = FileMeshGenerator
    file = copper_wire.msh
  []
[]

[Variables]
  # The real and complex components of the electric field in
  # the frequency domain
  [E_real]
    family = NEDELEC_ONE
    order = FIRST
  []
  [E_imag]
    family = NEDELEC_ONE
    order = FIRST
  []

  # The temperature of the air in the copper wire
  [T]
    initial_condition = 293.0 #in K
  []
[]

[Kernels]
  ### Physics to determine the electric field propagation ###
  # The propagation of the real component
  [curl_curl_real]
    type = CurlCurlField
    variable = E_real
  []
  # Displacement current in the wire
  [coeff_real]
    type = ADMatWaveReaction
    variable = E_real
    E_real = E_real
    E_imag = E_imag
    wave_coef_real = wave_equation_coefficient_real
    wave_coef_imag = wave_equation_coefficient_imaginary
    component = real
  []
  # Current supplied to the wire
  [source_real]
    type = VectorCurrentSource
    variable = E_real
    component = real
    source_real = curr_real
    source_imag = curr_imag
    function_coefficient = omegaMu
  []

  # The propagation of the complex component
  [curl_curl_imag]
    type = CurlCurlField
    variable = E_imag
  []
  # Displacement current in the wire
  [coeff_imag]
    type = ADMatWaveReaction
    variable = E_imag
    E_real = E_real
    E_imag = E_imag
    wave_coef_real = wave_equation_coefficient_real
    wave_coef_imag = wave_equation_coefficient_imaginary
    component = imaginary
  []
  # Current supplied to the wire
  [source_imaginary]
    type = VectorCurrentSource
    variable = E_imag
    component = imaginary
    source_real = curr_real
    source_imag = curr_imag
    function_coefficient = omegaMu
  []

  ### Physics to determine the heat transfer ###
  # Heat transfer in the copper wire
  [HeatTdot_in_copper]
    type = ADHeatConductionTimeDerivative
    variable = T
    specific_heat = specific_heat_copper
    density_name = density_copper
  []
  [HeatDiff_in_copper]
    type = ADHeatConduction
    variable = T
    thermal_conductivity = thermal_conductivity_copper
  []
  # Commenting off the Joule heating to decouple the temperature
  # and electric field solve to just study the electric field solutions
  # NOTE: The temperature without Joule heating is 293K throughout the run time
  # [HeatSrc]
  #   type = ADJouleHeatingSource
  #   variable = T
  #   heating_term = 'electric_field_heating'
  # []
[]

[AuxVariables]
  # Decomposing the magnetic vector potential
  # for the electric field calculations
  [E_x_real]
    family = MONOMIAL
    order = FIRST
  []
  [E_y_real]
    family = MONOMIAL
    order = FIRST
  []

  [E_x_imag]
    family = MONOMIAL
    order = FIRST
  []
  [E_y_imag]
    family = MONOMIAL
    order = FIRST
  []

  # The electrical conductivity for the electric
  # field calculations
  [elec_cond]
    family = MONOMIAL
    order = FIRST
  []

  # The electric field profile determined from
  # the magnetic vector potential
  [E_real_total]
    family = NEDELEC_ONE
    order = FIRST
  []
  [E_imag_total]
    family = NEDELEC_ONE
    order = FIRST
  []
[]

[Functions]
  # The supplied current density to the wire
  [curr_real]
    type = ParsedVectorFunction
    expression_x = 60e6 # Units in A/m^2, equivalent to 2850 A in a 5mm diameter wire
  []
  [curr_imag] # defaults to '0.0 0.0 0.0'
    type = ParsedVectorFunction
  []

  # Permittivity of free space
  [eps_real_func]
    type = ParsedFunction
    expression = '8.8542e-12' # Units in F/m
  []
  # Permeability of free space
  [mu_real_func]
    type = ParsedFunction
    expression = '4*pi*1e-7' # Units in N/A^2
  []
  # The angular drive frequency of the system
  [omega_real_func]
    type = ParsedFunction
    expression = '2*pi*60' # Units in rad/s
  []

  # The angular frequency time permeability of free space
  [omegaMu]
    type = ParsedFunction
    symbol_names = 'omega mu'
    symbol_values = 'omega_real_func mu_real_func'
    expression = 'omega*mu'
  []
  # Beta function for absorbing BC
  [beta]
    type = ParsedFunction
    symbol_names = 'omega'
    symbol_values = 'omega_real_func'
    expression = '2*pi*omega/3e8'
  []
[]

[BCs]
  ### Electric field boundary conditions ###
  # Wave absorbing BC for the real component
  [absorbing_left_real]
    type = VectorEMRobinBC
    variable = E_real
    component = real
    beta = beta
    coupled_field = E_imag
    mode = absorbing
    boundary = 'port exit walls'
  []

  # Wave absorbing BC for the complex component
  [absorbing_left_imag]
    type = VectorEMRobinBC
    variable = E_imag
    component = imaginary
    beta = beta
    coupled_field = E_real
    mode = absorbing
    boundary = 'port exit walls'
  []

  ### Temperature boundary conditions ###
  # Convective heat flux BC with copper wire
  # exposed to air
  [surface]
    type = ADConvectiveHeatFluxBC
    variable = T
    boundary = walls
    T_infinity = 293
    heat_transfer_coefficient = 10
  []
[]

[Materials]
  [k]
    type = ADGenericConstantMaterial
    prop_names = 'thermal_conductivity_copper'
    prop_values = '397.48' #in W/(m K)
  []
  [cp]
    type = ADGenericConstantMaterial
    prop_names = 'specific_heat_copper'
    prop_values = '385.0' #in J/(kg K)
  []
  [rho]
    type = ADGenericConstantMaterial
    prop_names = 'density_copper'
    prop_values = '8920.0' #in kg/(m^3)
  []

  # Electrical conductivity (copper is default material)
  [sigma]
    type = ADElectricalConductivity
    temperature = T
    block = copper
  []
  # Material that supplies the correct Joule heating formulation
  [ElectromagneticMaterial]
    type = ElectromagneticHeatingMaterial
    electric_field = E_real
    complex_electric_field = E_imag
    electric_field_heating_name = electric_field_heating
    electrical_conductivity = electrical_conductivity
    formulation = FREQUENCY
    solver = ELECTROMAGNETIC
    block = copper
  []

  # Coefficient for wave propagation
  [WaveCoeff]
    type = WaveEquationCoefficient
    eps_rel_imag = 0
    eps_rel_real = eps_real
    k_real = omega_real
    mu_rel_imag = 0
    mu_rel_real = mu_real
  []
  [eps_real]
    type = ADGenericFunctionMaterial
    prop_names = eps_real
    prop_values = eps_real_func
  []
  [mu_real]
    type = ADGenericFunctionMaterial
    prop_names = mu_real
    prop_values = mu_real_func
  []
  [omega_real]
    type = ADGenericFunctionMaterial
    prop_names = omega_real
    prop_values = omega_real_func
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
  scheme = bdf2
  solve_type = NEWTON
  line_search = NONE
  petsc_options_iname = '-pc_type -pc_factor_shift_type -pc_factor_shift_amount'
  petsc_options_value = 'lu NONZERO 1.e-10'
  dt = 0.1
  end_time = 10
  automatic_scaling = true
  nl_abs_tol = 1e-8
[]

[Outputs]
  exodus = true
  perf_graph = true
[]
