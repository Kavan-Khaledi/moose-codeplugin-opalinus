
[GlobalParams]
  time_unit = days
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = false
  PorousFlowDictator = dictator
[]

inactive_block_names = 'tunnel_inactive'
inactive_block_ids = '100'

[Mesh]
  [file]
    type = FileMeshGenerator
    file = "HM_tunnel_x.msh"
  []

  second_order = true

  add_subdomain_names = '${inactive_block_names}'
  add_subdomain_ids = '${inactive_block_ids}'
[]

!include HM_tunnel_x.groups.i

active_block_names = '${RockVolumes} ${TunnelVolumes}'

[Problem]
  kernel_coverage_check = 'SKIP_LIST'
  kernel_coverage_block_list = '${inactive_block_names}'
  material_coverage_check = 'SKIP_LIST'
  material_coverage_block_list = '${inactive_block_names}'
[]

[Variables]
  [disp_x]
    order = SECOND
    family = LAGRANGE
    block = '${active_block_names}'
  []
  [disp_y]
    order = SECOND
    family = LAGRANGE
    block = '${active_block_names}'
  []
  [disp_z]
    order = SECOND
    family = LAGRANGE
    block = '${active_block_names}'
  []
  [porepressure]
    order = SECOND
    family = LAGRANGE
    scaling = 1e-5
    block = '${active_block_names}'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        incremental = true
        #add_variables = true
        eigenstrain_names = ini_stress
        block = '${active_block_names}'
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
    block = '${active_block_names}'
  []

  [effective_stress_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    component = 1
    block = '${active_block_names}'
  []

  [effective_stress_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    component = 2
    block = '${active_block_names}'
  []

  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = porepressure
    block = '${active_block_names}'
  []

  [flux]
    type = PorousFlowFullySaturatedDarcyFlow
    variable = porepressure
    gravity = '0 0 0'
    fluid_component = 0
    block = '${active_block_names}'
  []

  [poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    variable = porepressure
    fluid_component = 0
    block = '${active_block_names}'
  []
[]

[ICs]
  [porepressure]
    type = FunctionIC
    variable = porepressure
    function = '9'
    block = '${active_block_names}'
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
    block = '${active_block_names}'
  []
  [p]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = p
    scalar_type = Hydrostatic
    block = '${active_block_names}'
  []
  [q]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = q
    scalar_type = VonMisesStress
    block = '${active_block_names}'
  []
[]

[BCs]

  [fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = '${XMinSurfaces} ${XMaxSurfaces}'
    value = 0.0
  []

  [fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = '${YMinSurfaces} ${YMaxSurfaces}'
    value = 0.0
  []

  [fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = '${ZMinSurfaces} ${ZMaxSurfaces}'
    value = 0.0
  []

  [fix_pw]
    type = PorousFlowPiecewiseLinearSink
    variable = porepressure
    PT_shift = 9
    pt_vals = '-1E3 1E3'
    multipliers = '-1 1'
    fluid_phase = 0
    flux_function = 1e5
    boundary = '${ZMinSurfaces} ${ZMaxSurfaces}  ${YMaxSurfaces} ${XMinSurfaces} ${XMaxSurfaces}'
  []
  [y_pw]
    type = PorousFlowSink
    boundary = '${YMinSurfaces}'
    variable = porepressure
    flux_function = 0.0
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
    block = '${active_block_names}'
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    block = '${active_block_names}'
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
    block = '${active_block_names}'
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
    block = '${active_block_names}'
  []

  [massfrac]
    type = PorousFlowMassFraction
    block = '${active_block_names}'
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
    block = '${active_block_names}'
  []

  [porosity_bulk]
    type = PorousFlowPorosity
    fluid = true
    mechanical = true
    ensure_positive = true
    porosity_zero = 0.11
    solid_bulk = 1.3333E10
    block = '${active_block_names}'
  []

  [undrained_density_0]
    type = GenericConstantMaterial
    block = '${active_block_names}'
    prop_names = density
    prop_values = 2500
  []

  [rock_permeability_bulk]
    type = OpalinusPermeabilityTensor
    permeability1 = 5e-19
    permeability2 = 5e-19
    permeability3 = 5e-20

    local_coordinate_system = 'ucsOpalinusMaterial'
    block = '${active_block_names}'
  []

  [rock_elasticity_tensor]
    type = OpalinusElasticityTensor
    youngs_modulus_in_plane = 11000
    youngs_modulus_normal = 6000
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.25
    shear_module_normal = 2000
    local_coordinate_system = 'ucsOpalinusMaterial'
    block = '${active_block_names}'
  []

  [opalinus_mont_terri]
    type = OpalinusPerfectPlasticStressUpdate
    gama_mean = 0.9
    parameter_omega_1 = 0.15
    parameter_b_1 = 6.7
    p_tensile = 6
    local_coordinate_system = 'ucsOpalinusMaterial'

    yield_function_tol = 1e-3
    smoothing_tol = 0.0
    tip_smoother = 2
    max_NR_iterations = 50
    min_step_size = 0.004
    block = '${active_block_names}'
  []

  [stress]
    type = ComputeMultipleInelasticStress
    block = '${active_block_names}'
    inelastic_models = 'opalinus_mont_terri'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
  []

  [ini_stress]
    block = '${active_block_names}'
    type = ComputeEigenstrainFromGeostaticInitialStress
    eigenstrain_name = 'ini_stress'
    local_coordinate_system = 'ucsInitialStress'
    principal_stress_1 = 13.5   # effective stresses
    principal_stress_2 = 17.5   # effective stresses
    principal_stress_3 = 13.5   # effective stresses
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
    dip_direction_degree = 0
    dip_angle_degree = 0
    dip_option = 'e1_e2_plane_e1_horizontal'
  []
[]

[MeshModifiers]
  [GlobalSubdomainModifier]
    type = TimedSubdomainModifier
    times = '2 3 4 5 6 7 8 9 10 11 12 13 14 15 16'
    blocks_from = ' tunnel01 tunnel02 tunnel03 tunnel04 tunnel05 tunnel06 tunnel07 tunnel08 tunnel09 tunnel10 tunnel11 tunnel12 tunnel13 tunnel14 tunnel15'
    blocks_to = 'tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive tunnel_inactive'
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
  #automatic_scaling = true

  line_search = none

  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-7

  l_max_its = 20
  nl_max_its = 15
  dt = 1
  start_time = 0.0
  end_time = 17

  [Quadrature]
    type = SIMPSON
    order = SECOND
  []
[]

[Outputs]
  perf_graph = true
  print_linear_residuals = false
  exodus = true
[]
