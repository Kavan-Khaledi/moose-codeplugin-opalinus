
[GlobalParams]
  time_unit = days
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = false
  PorousFlowDictator = dictator
[]


[Mesh]
  [file]
    type = FileMeshGenerator
    file = "HM_tunnel_x.msh"
  []

  second_order = true

  active_block_names = '${RockVolumes} ${TunnelVolumes}'
  inactive_block_names = 'tunnel_inactive'

  add_subdomain_names = '${inactive_block_names}'
[]

!include HM_tunnel_x.groups.i

[Problem]
  kernel_coverage_check = 'SKIP_LIST'
  kernel_coverage_block_list = '${Mesh/inactive_block_names}'
  material_coverage_check = 'SKIP_LIST'
  material_coverage_block_list = '${Mesh/inactive_block_names}'
[]

[Variables]
  [disp_x]
    family = LAGRANGE
    order = SECOND
    block = '${Mesh/active_block_names}'
  []
  [disp_y]
    family = LAGRANGE
    order = SECOND
    block = '${Mesh/active_block_names}'
  []
  [disp_z]
    family = LAGRANGE
    order = SECOND
    block = '${Mesh/active_block_names}'
  []
  [porepressure]
    family = LAGRANGE
    order = SECOND
    scaling = 1e-5
    block = '${Mesh/active_block_names}'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        strain = SMALL
        add_variables = false
        incremental = true
        eigenstrain_names = ini_stress
        block = '${Mesh/active_block_names}'
      []
    []
  []
[]

# ===== Kernels: PorousFlow =====
[Kernels]
  [effective_stress_x]
    type = PorousFlowEffectiveStressCoupling
    block = '${Mesh/active_block_names}'
    variable = disp_x
    component = 0
  []

  [effective_stress_y]
    type = PorousFlowEffectiveStressCoupling
    block = '${Mesh/active_block_names}'
    variable = disp_y
    component = 1
  []

  [effective_stress_z]
    type = PorousFlowEffectiveStressCoupling
    block = '${Mesh/active_block_names}'
    variable = disp_z
    component = 2
  []

  [mass0]
    type = PorousFlowMassTimeDerivative
    block = '${Mesh/active_block_names}'
    fluid_component = 0
    variable = porepressure
  []

  [flux]
    type = PorousFlowFullySaturatedDarcyFlow
    block = '${Mesh/active_block_names}'
    variable = porepressure
    gravity = '0 0 0'
    fluid_component = 0
  []

  [poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    block = '${Mesh/active_block_names}'
    variable = porepressure
    fluid_component = 0
  []
[]

[ICs]
  [porepressure]
    type = FunctionIC
    variable = 'porepressure'
    block = '${Mesh/active_block_names}'
    function = '9'
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
    block = '${Mesh/active_block_names}'
    property = plastic_internal_parameter
    variable = internal_plastic_variable
    index = 0
  []

  [p]
    type = RankTwoScalarAux
    block = '${Mesh/active_block_names}'
    rank_two_tensor = stress
    variable = p
    scalar_type = Hydrostatic
  []

  [q]
    type = RankTwoScalarAux
    block = '${Mesh/active_block_names}'
    rank_two_tensor = stress
    variable = q
    scalar_type = vonMisesStress
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
    boundary = '${YMinSurfaces_fix}'
    variable = porepressure
    flux_function = 0.0
  []
[]

YMinSurfaces_fix = 'rock_i01_f00a
                    rock_i01_f00b
                    rock_i01_f00c
                    rock_i01_f00d
                    rock_i02_f00a
                    rock_i02_f00b
                    rock_i02_f00c
                    rock_i02_f00d' # tunnel_f00'   # @Kavan-Khaledi: remove 'tunnel_f00' to avoid the error shown below

# We caught a libMesh error in ThreadedElementLoopBase:Assertion `i < _val.size()' failed.
# i = 0
# _val.size() = 0

[FluidProperties]
  [simple_fluid]
    type = SimpleFluidProperties
    bulk_modulus = 2E3
    density0 = 1000
    thermal_expansion = 0
    viscosity = 9.0E-10 #MPas -2.1e-12 * exp(1808/T) --> T = 298
  []
[]

# Material: Volume Elements
[Materials]

  [rock_elasticity_tensor]
    type = OpalinusElasticityTensor
    block = '${Mesh/active_block_names}'
    local_coordinate_system = 'ucsOpalinusMaterial'
    youngs_modulus_in_plane = 11000
    youngs_modulus_normal = 6000
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.25
    shear_module_normal = 2000
  []

  [opalinus_mont_terri]
    type = OpalinusPerfectPlasticStressUpdate
    block = '${Mesh/active_block_names}'
    local_coordinate_system = 'ucsOpalinusMaterial'
    gama_mean = 0.9
    parameter_omega_1 = 0.15
    parameter_b_1 = 6.7
    p_tensile = 6
    yield_function_tol = 1e-3
    smoothing_tol = 0.0
    tip_smoother = 2
    max_NR_iterations = 50
    min_step_size = 0.004
  []

  [stress]
    type = ComputeMultipleInelasticStress
    block = '${Mesh/active_block_names}'
    inelastic_models = 'opalinus_mont_terri'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
  []

  [ini_stress]
    block = '${Mesh/active_block_names}'
    type = ComputeEigenstrainFromGeostaticInitialStress
    eigenstrain_name = 'ini_stress'
    local_coordinate_system = 'ucsInitialStress'
    principal_stress_1 = 13.5   # effective stresses
    principal_stress_2 = 17.5   # effective stresses
    principal_stress_3 = 13.5   # effective stresses
  []

  [temperature]
    type = PorousFlowTemperature
    block = '${Mesh/active_block_names}'
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    block = '${Mesh/active_block_names}'
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
    block = '${Mesh/active_block_names}'
  []
  [ppss]
    type = PorousFlow1PhaseFullySaturated
    porepressure = porepressure
    block = '${Mesh/active_block_names}'
  []

  [massfrac]
    type = PorousFlowMassFraction
    block = '${Mesh/active_block_names}'
  []
  [simple_fluid]
    type = PorousFlowSingleComponentFluid
    fp = simple_fluid
    phase = 0
    block = '${Mesh/active_block_names}'
  []

  [porosity_bulk]
    type = PorousFlowPorosity
    fluid = true
    mechanical = true
    ensure_positive = true
    porosity_zero = 0.11
    solid_bulk = 1.3333E10
    block = '${Mesh/active_block_names}'
  []

  [undrained_density_0]
    type = GenericConstantMaterial
    block = '${Mesh/active_block_names}'
    prop_names = density
    prop_values = 2500
  []

  [rock_permeability_bulk]
    type = OpalinusPermeabilityTensor
    block = '${Mesh/active_block_names}'
    permeability1 = 5e-19
    permeability2 = 5e-19
    permeability3 = 5e-20

    local_coordinate_system = 'ucsOpalinusMaterial'
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
