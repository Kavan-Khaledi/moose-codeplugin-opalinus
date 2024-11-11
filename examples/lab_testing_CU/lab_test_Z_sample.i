## Modeling of triaxial consolidated undrained tests
## Cylindrical Element (diameter=3cm, hight=6cm)

pconf = 8 #MPa --> the initial effective confining pressure applied to the sample before shearing
pw = 9 # the initial pore pressure after consilidation phase
pconf_total = ${pconf}+${pw} #MPa Total confining presure
dip_angle = 45

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  PorousFlowDictator = dictator
  biot_coefficient = 1
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = sample3d.msh
  []
  second_order = true
[]

[Variables]
  [disp_x]
    family = LAGRANGE
    order = SECOND
  []
  [disp_y]
    family = LAGRANGE
    order = SECOND
  []
  [disp_z]
    family = LAGRANGE
    order = SECOND
  []
  [porepressure]
    family = LAGRANGE
    order = SECOND
  []
[]

[ICs]
  [porepressure]
    type = FunctionIC
    variable = porepressure
    function = ${pw}
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        add_variables = true
        incremental = true
        eigenstrain_names = ini_stress
      []
    []
  []
[]

[Kernels]
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    use_displaced_mesh = false
    variable = disp_x
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    use_displaced_mesh = false
    variable = disp_y
    component = 1
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    use_displaced_mesh = false
    variable = disp_z
    component = 2
  []
  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = porepressure
  []
  [flux]
    type = PorousFlowAdvectiveFlux
    use_displaced_mesh = false
    variable = porepressure
    gravity = '0 0 0'
    fluid_component = 0
  []
  [poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    variable = porepressure
    fluid_component = 0
  []
[]

[AuxVariables]
  [effective_mean_pressure]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = -${pconf}
  []
  [deviatoric_stress]
    family = MONOMIAL
    order = CONSTANT
    initial_condition = 0
  []

  [total_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
  [internal_plastic_variable]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [effective_mean_pressure]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = effective_mean_pressure
    scalar_type = hydrostatic
    execute_on = 'TIMESTEP_BEGIN'
  []

  [deviatoric_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = deviatoric_stress
    scalar_type = vonmisesStress
    execute_on = 'TIMESTEP_BEGIN'
  []

  [total_strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  []

  [internal_plastic_variable]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 0
    variable = internal_plastic_variable
  []
[]

[BCs]

  [side1]
    type = Pressure
    boundary = 'confining'
    variable = 'disp_x'
    function = ${pconf_total}
  []
  [side2]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'top'
    function = 'if(t<=5000,0,-2.8e-9*(t-5000))'
  []
  [side3]
    type = Pressure
    boundary = 'confining'
    variable = 'disp_y'
    function = ${pconf_total}
  []
  [no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'bottom top'
    value = 0.0
  []
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.0
  []
  [no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom top'
    value = 0.0
  []
  [undrained]
    type = PorousFlowSink #undrained
    boundary = 'top bottom confining'
    flux_function = 0
    variable = porepressure
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []

  [ucsInitialStress]
    type = CartesianLocalCoordinateSystem
    e1 = '1 0 0'
    e2 = '0 1 0'
  []
  [ucsOpalinusMaterial]
    type = CartesianLocalCoordinateSystem
    dip_direction_degree = '0'
    dip_angle_degree = ${dip_angle}
    dip_option = 'e1_e2_plane_e1_horizontal'
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2E3
    density0 = 1000
    thermal_expansion = 0
    viscosity = 9.0E-10 #MPas----> 2.1e-12exp(1808/T)-T=298
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
  []
  [porosity_bulk]
    type = PorousFlowPorosity
    fluid = true
    mechanical = true
    ensure_positive = true
    porosity_zero = 0.11
    solid_bulk = 1.3333E10
  []

  [permeability_bulk]
    type = OpalinusPermeabilityTensor
    permeability1 = 5e-19
    permeability2 = 5e-19
    permeability3 = 5e-20
    local_coordinate_system = 'ucsOpalinusMaterial'
  []

  [relperm0]
    type = PorousFlowRelativePermeabilityCorey
    n = 1
    phase = 0
  []

  [elasticity_tensor]
    type = OpalinusElasticityTensor
    youngs_modulus_in_plane = 11000
    youngs_modulus_normal = 6000
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.25
    shear_module_normal = 2000
    local_coordinate_system = 'ucsOpalinusMaterial'
  []

  [ini_stress]
    type = ComputeEigenstrainFromGeostaticInitialStress
    eigenstrain_name = 'ini_stress'
    local_coordinate_system = 'ucsInitialStress'
    principal_stress_1 = ${pconf}
    principal_stress_2 = ${pconf}
    principal_stress_3 = ${pconf}
  []

  [opalinus]
    type = OpalinusPerfectPlasticStressUpdate
    local_coordinate_system = 'ucsOpalinusMaterial'
    gama_mean = 0.9
    parameter_omega_1 = 0.15
    parameter_b_1 = 6.7
    p_tensile = 6

    yield_function_tol = 1e-3
    smoothing_tol = 0.0
    tip_smoother = 2
    min_step_size = 0.004
    max_NR_iterations = 40
  []

  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'opalinus'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
  []

  [undrained_density_0]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2500
  []
[]

[Preconditioning]
  [SMP]
    type = SMP
    full = true
    petsc_options = '-ksp_snes_ew'
    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_levels -ksp_gmres_restart'
    petsc_options_value = 'gmres hypre boomeramg lu NONZERO 4 301'
  []
[]

[Executioner]
  type = Transient
  #automatic_scaling = true
  solve_type = 'PJFNK'

  line_search = none

  nl_abs_tol = 5e-6
  nl_rel_tol = 1e-10

  l_max_its = 15
  nl_max_its = 15

  start_time = 0.0
  dt = 10000
  end_time = 150000 # hours

  [Quadrature]
    type = SIMPSON
    order = SECOND
  []
[]

[Outputs]
  time_step_interval = 1
  print_linear_residuals = false
  exodus = true
[]
