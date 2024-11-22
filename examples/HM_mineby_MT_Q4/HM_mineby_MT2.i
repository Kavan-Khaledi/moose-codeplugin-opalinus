
[GlobalParams]
  time_unit = days
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = false
  PorousFlowDictator = dictator
[]

!include MT2_mesh.i

active_block_ids = '1 2 3 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115'
inactive_block_names = 'empty'

[Problem]
  solve = true
  kernel_coverage_check = 'SKIP_LIST'
  kernel_coverage_block_list = '${inactive_block_names}'
  material_coverage_check = 'SKIP_LIST'
  material_coverage_block_list = '${inactive_block_names}'
[]

[Variables]
  [disp_x]
    order = FIRST
    family = LAGRANGE
    block = '${active_block_ids}'
  []
  [disp_y]
    order = FIRST
    family = LAGRANGE
    block = '${active_block_ids}'
  []
  [disp_z]
    order = FIRST
    family = LAGRANGE
    block = '${active_block_ids}'
  []
  [porepressure]
    order = FIRST
    family = LAGRANGE
    #scaling = 1e-5
    block = '${active_block_ids}'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        incremental = true
        add_variables = false
        eigenstrain_names = ini_stress
        block = '${active_block_ids}'
      []
    []
  []
[]

# ===== Kernels: Inactive Domains =====
[Kernels]
  [effective_stress_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    component = 0
    block = '${active_block_ids}'
  []

  [effective_stress_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 1
    block = '${active_block_ids}'
  []

  [effective_stress_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    component = 2
    block = '${active_block_ids}'
  []

  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = porepressure
    block = '${active_block_ids}'
    # multiply_by_density = true
  []

  [flux]
    type = PorousFlowAdvectiveFlux
    variable = porepressure
    gravity = '0 0 0'
    block = '${active_block_ids}'
  []

  [poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    variable = porepressure
    fluid_component = 0
    block = '${active_block_ids}'
    # multiply_by_density = true
  []
[]

[ICs]
  [porepressure]
    type = FunctionIC
    variable = porepressure
    function = '2'
    block = '${active_block_ids}'
  []
[]

[AuxVariables]

  [internal_plastic_variable]
    order = CONSTANT
    family = MONOMIAL
  []
  [p]
    order = CONSTANT
    family = MONOMIAL
  []
  [q]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]

  [internal_plastic_variable]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 0
    variable = internal_plastic_variable
    block = '${active_block_ids}'
  []
  [p]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = p
    scalar_type = Hydrostatic
    block = '${active_block_ids}'
  []
  [q]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = q
    scalar_type = VonMisesStress
    block = '${active_block_ids}'
  []
[]

[BCs]

  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'top bottom left right'
    value = 0.0
  []

  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'front back'
    value = 0.0
  []

  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'top bottom left right'
    value = 0.0
  []

  [fix_pw]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    PT_shift = 2
    pt_vals = '-1E3 1E3'
    multipliers = '-1 1'
    fluid_phase = 0
    flux_function = 1e5
    boundary = 'left right top bottom'
  []
  [y_pw]
    type = PorousFlowSink
    boundary = 'front back'
    variable = porepressure
    flux_function = 0.0
  []
  [fix_pw_inner]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    PT_shift = 0.1
    pt_vals = '-1E3 1E3'
    multipliers = '-1 1'
    fluid_phase = 0
    flux_function = 'if(t<=1,0.0,if(y<=min(-21+t+0.01,-5),1e5,0.0))'
    boundary = 'inner_ring'
  []
[]

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2E3
    density0 = 1000
    thermal_expansion = 0
    viscosity = 9.0E-10 #MPas -2.1e-12 * exp(1808/T) --> T = 298
  []
[]

[Materials]

  [temperature]
    type = PorousFlowTemperature
    block = '${active_block_ids}'
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    block = '${active_block_ids}'
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
    block = '${active_block_ids}'
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
    block = '${active_block_ids}'
  []

  [massfrac]
    type = PorousFlowMassFraction
    block = '${active_block_ids}'
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
    block = '${active_block_ids}'
  []

  [porosity_bulk]
    type = PorousFlowPorosity
    fluid = false
    mechanical = true
    ensure_positive = true
    porosity_zero = 0.18
    solid_bulk = 1e8
    block = '${active_block_ids}'
  []

  [undrained_density_0]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2500
    block = '${active_block_ids}'
  []

  [rock_permeability_bulk]
    type = OpalinusPermeabilityTensor
    permeability1 = 5e-19
    permeability2 = 5e-19
    permeability3 = 5e-20

    local_coordinate_system = 'ucsOpalinusMaterial'
    block = '${active_block_ids}'
  []

  [relperm0]
    type = PorousFlowRelativePermeabilityConst
    phase = 0
    kr = 1
    block = '${active_block_ids}'
  []

  [rock_elasticity_tensor]
    type = OpalinusElasticityTensor
    youngs_modulus_in_plane = 1500
    youngs_modulus_normal = 750
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.35
    shear_module_normal = 1000
    local_coordinate_system = 'ucsOpalinusMaterial'
    block = '${active_block_ids}'
  []

  [opalinus_mont_terri]
    type = OpalinusPerfectPlasticStressUpdate
    gama_mean = 0.85
    parameter_omega_1 = 0.15
    parameter_b_1 = 6.7
    p_tensile = 3
    local_coordinate_system = 'ucsOpalinusMaterial'

    yield_function_tol = 1e-4
    smoothing_tol = 0.0
    tip_smoother = 1.0
    max_NR_iterations = 50
    min_step_size = 0.004
    block = '${active_block_ids}'
  []

  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'opalinus_mont_terri'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
    block = '${active_block_ids}'
  []

  [ini_stress]

    type = ComputeEigenstrainFromGeostaticInitialStress
    eigenstrain_name = 'ini_stress'
    local_coordinate_system = 'ucsInitialStress'
    principal_stress_1 = 2.5
    principal_stress_2 = 0.7
    principal_stress_3 = 4.5
    block = '${active_block_ids}'
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
    e1 = '0.939 0 -0.342'
    e2 = '0 1 0'
  []
  [ucsOpalinusMaterial]
    type = CartesianLocalCoordinateSystem
    dip_direction_degree = 270
    dip_angle_degree = 30
    dip_option = 'e1_e2_plane_e1_horizontal'
  []
[]
[MeshModifiers]
  [GlobalSubdomainModifier]
    type = TimedSubdomainModifier
    times = '2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    blocks_from = ' 101 102 103 104 105 106 107 108 109 110 111 112 113 114 115'
    blocks_to = '500 500 500 500 500 500 500 500 500 500 500 500 500 500 500'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []
[]

[Preconditioning]
  [.\SMP]
    type = SMP
    full = true
    petsc_options = '-ksp_snes_ew'

    petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type -sub_pc_type -sub_pc_factor_shift_type -sub_pc_factor_levels -ksp_gmres_restart'
    petsc_options_value = 'gmres hypre boomeramg lu NONZERO 4 301'
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  automatic_scaling = true
  compute_scaling_once = false

  line_search = NONE

  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-7

  l_max_its = 20
  nl_max_its = 15

  start_time = 0.0
  end_time = 300

  [TimeStepper]
    type = FunctionDT
    function = 'if(t<17, 1, 30)'
  []

  #[Quadrature]
  #  type = SIMPSON
  #  order = FIRST
  #[]
[]

[Outputs]
  perf_graph = true
  print_linear_residuals = false
  exodus = true
[]
