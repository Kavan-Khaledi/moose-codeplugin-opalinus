PorousFlowDictatorName = 'dictator'

[GlobalParams]
  time_unit = days
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = false
  PorousFlowDictator = '${PorousFlowDictatorName}'
[]

inactive_block_names = 'tunnel01_inactive'
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

[AuxVariables]
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxVariables]
  [normal_axis_x]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_axis_y]
    order = CONSTANT
    family = MONOMIAL
  []
  [normal_axis_z]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
    execute_on = timestep_end
    block = '${active_block_names}'
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
    block = '${active_block_names}'
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
    block = '${active_block_names}'
  []
[]

[AuxKernels]
  [normal_axis_x]
    type = MaterialStdVectorAux
    property = normal_local_vector
    variable = normal_axis_x
    index = 0
    execute_on = TIMESTEP_END
    block = '${active_block_names}'
  []
  [normal_axis_y]
    type = MaterialStdVectorAux
    property = normal_local_vector
    variable = normal_axis_y
    index = 1
    execute_on = TIMESTEP_END
    block = '${active_block_names}'
  []
  [normal_axis_z]
    type = MaterialStdVectorAux
    property = normal_local_vector
    variable = normal_axis_z
    index = 2
    execute_on = TIMESTEP_END
    block = '${active_block_names}'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [active_blocks]
        strain = SMALL
        incremental = true
        eigenstrain_names = ini_stress
        block = '${active_block_names}'
      []
    []
  []
[]

#[PorousFlowFullySaturated]
#  coupling_type = HydroMechanical
#  porepressure = porepressure
#  biot_coefficient = 1
#  fp = simple_fluid
#  stabilization = FULL
#  gravity = '0 0 0 ' #-9.81'
#  add_darcy_aux = false
#  eigenstrain_names = ini_stress
#  dictator_name = '${PorousFlowDictatorName}'
#  block = '${active_block_names}'
#  active = 'porepressure'
#[]

# ===== Kernels: PorousFlowFullySaturated =====
[Kernels]
  [effective_stress_x]
    type = PorousFlowEffectiveStressCoupling
    variable = 'disp_x'
    component = 0
    block = '${active_block_names}'
  []

  [effective_stress_y]
    type = PorousFlowEffectiveStressCoupling
    variable = 'disp_y'
    component = 1
    block = '${active_block_names}'
  []

  [effective_stress_z]
    type = PorousFlowEffectiveStressCoupling
    variable = 'disp_z'
    component = 2
    block = '${active_block_names}'
  []

  [mass0]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = 'porepressure'
    block = '${active_block_names}'
  []

  [flux]
    type = PorousFlowFullySaturatedDarcyFlow
    variable = 'porepressure'
    gravity = '0 0 0'
    fluid_component = 0
    block = '${active_block_names}'
  []

  [poro_vol_exp]
    type = PorousFlowMassVolumetricExpansion
    variable = 'porepressure'
    fluid_component = 0
    block = '${active_block_names}'
  []
[]

# ===== Kernels: Gravity =====
#[Kernels]
#  [./gravity]
#    type = Gravity
#    block = '${active_block_names}'
#    use_displaced_mesh = false
#    variable = disp_z
#    value = -${gravitational_acceleration}
#  []
#[]

[ICs]
  [porepressure]
    type = FunctionIC
    variable = porepressure
    function = '2'
    block = '${active_block_names}'
  []
[]

[AuxVariables]
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
    type = DirichletBC
    variable = porepressure
    boundary = '${ZMinSurfaces} ${ZMaxSurfaces} ${YMinSurfaces} ${YMaxSurfaces} ${XMinSurfaces} ${XMaxSurfaces}'
    value = 2.0
  []

  # Pressure applied on the top surface
  [Pressure]
    [pressure_zmax]
      boundary = '${ZMaxSurfaces}'
      function = 1000 #Newtons
    []
  []

  #  [pressure_face]
  #    type = FunctionDirichletBC
  #    variable = porepressure
  #    boundary = 'tunnel_face_y'
  #    function = 'if(t<=1, 2*(1-t), 0.0)'
  #  []
  #
  #  [H_BC_1]
  #    type = PorousFlowPiecewiseLinearSink
  #    variable = porepressure
  #    PT_shift = 2
  #    pt_vals = '-1E3 1E3'
  #    multipliers = '-1 1'
  #    fluid_phase = 0
  #    flux_function = 1E5
  #    boundary = 'bottom_z top_z left_x right_x'
  #  []
  #
  #  [H_BC_2]
  #    type = PorousFlowSink
  #    boundary = 'front_y back_y'
  #    flux_function = 0
  #    variable = porepressure
  #  []
  #
  #  [pp_in]
  #    type = PorousFlowPiecewiseLinearSink
  #    variable = porepressure
  #    PT_shift = 0.0
  #    pt_vals = '-1E3 1E3'
  #    multipliers = '-1 1'
  #    fluid_phase = 0
  #    flux_function = inner_pressure_method1
  #    boundary = 'inner_boundary_tunnel'
  #  []
[]

[Functions]
  [Tunnel01_elasticity_tensor_prefactor]
    type = "StagedFunction"
  []
  [tunnel01_permeability_prefactor]
    type = "StagedFunction"
  []

  # [excavat_method1]
  #   type = ParsedFunction
  #   symbol_names = 'ymin ymax v slop'
  #   symbol_values = '-2   15  1 2'
  #   # excavation face at ymin+(ymax-ymin)*min(t/end_t,1)
  #   # slope is the distance over which the modulus reduces from maxval to minval
  #   expression = 'if(y<ymin+v*t, 0.02, if(y<=ymin+v*t+slop, 0.02+0.98*(y-ymin-v*t)/slop, 1))'
  # []
  # [inner_pressure_method1]
  #   type = ParsedFunction
  #   symbol_names = 'ymin ymax v'
  #   symbol_values = '-2 15 1'
  #   expression = 'if(y<ymin+v*t,1e5,0)'
  # []
  [permeability_prefactor_method1]
    type = ParsedFunction
    symbol_names = 'ymin ymax v'
    symbol_values = '-2   15   1'
    expression = 'if(y<ymin+v*t, 1e5, 1)'
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

# Local coordinate systems (material and initial stress)
[UserObjects]
  #[ucsInitialStress]
  #  type = CartesianLocalCoordinateSystem
  #  e1 = '1 0 0'
  #  e2 = '0 1 0'
  #[]
  [ucsOpalinusMaterial]
    type = CartesianLocalCoordinateSystem
    dip_direction_degree = '90'
    dip_angle_degree = '35'
    dip_option = 'e1_e2_plane_e1_horizontal'
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
    type = PorousFlow1PhaseP
    porepressure = porepressure
    capillary_pressure = pc
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
    type = PorousFlowPorosityConst
    block = '${active_block_names}'
    porosity = 0.15
    PorousFlowDictator = '${PorousFlowDictatorName}'
  []

  [undrained_density_0]
    type = GenericConstantMaterial
    block = '${active_block_names}'
    prop_names = density
    prop_values = 2500
  []

  [rock_permeability_bulk]
    type = OpalinusPermeabilityTensor
    block = '${RockVolumes}'
    permeability_parallel = 1e-19
    permeability_normal = 1e-21
    dip_direction = 90
    dip = 35
  []

  [rock_elasticity_tensor]
    type = OpalinusElasticityTensor
    block = '${RockVolumes}'
    youngs_modulus_in_plane = 2500
    youngs_modulus_normal = 1150
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.35
    shear_module_normal = 1000
    local_coordinate_system = 'ucsOpalinusMaterial'
  []

  [tunnel01_permeability_bulk]
    type = OpalinusPermeabilityTensor
    block = 'tunnel01'
    permeability_parallel = 1e-19
    permeability_normal = 1e-21
    permeability_tensor_prefactor = tunnel01_permeability_prefactor
    dip_direction = 90
    dip = 35
  []

  [tunnel01_elasticity_tensor]
    type = OpalinusElasticityTensor
    block = 'tunnel01'
    youngs_modulus_in_plane = 2500
    youngs_modulus_normal = 1150
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.35
    shear_module_normal = 1000
    local_coordinate_system = 'ucsOpalinusMaterial'
    elasticity_tensor_prefactor = Tunnel01_elasticity_tensor_prefactor
  []

  [tunnelXX_permeability_bulk]
    type = OpalinusPermeabilityTensor
    block = 'tunnel02 tunnel03 tunnel04 tunnel05 tunnel06 tunnel07 tunnel08 tunnel09 tunnel10 tunnel11 tunnel12 tunnel13 tunnel14 tunnel15 tunnel16 tunnel17 tunnel18 tunnel19 tunnel20'
    permeability_parallel = 1e-19
    permeability_normal = 1e-21
    #permeability_tensor_prefactor = permeability_prefactor_method1
    dip_direction = 90
    dip = 35
  []

  [tunnelXX_elasticity_tensor]
    type = OpalinusElasticityTensor
    block = 'tunnel02 tunnel03 tunnel04 tunnel05 tunnel06 tunnel07 tunnel08 tunnel09 tunnel10 tunnel11 tunnel12 tunnel13 tunnel14 tunnel15 tunnel16 tunnel17 tunnel18 tunnel19 tunnel20'
    youngs_modulus_in_plane = 2500
    youngs_modulus_normal = 1150
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.35
    shear_module_normal = 1000
    local_coordinate_system = 'ucsOpalinusMaterial'
    # elasticity_tensor_prefactor = TunnelXX_elasticity_tensor_prefactor
  []

  [stress]
    type = ComputeMultipleInelasticStress
    block = '${active_block_names}'
    inelastic_models = ''
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
  []

  [ini_stress]
    type = OpalinusEigenstrainFromInitialStress
    block = '${active_block_names}'
    eigenstrain_name = ini_stress
    principal_stress_1 = 4.5
    trend_s1 = 0
    plunge_s1 = 90
    principal_stress_2 = 2.5
    trend_s2 = 90
    plunge_s2 = 0
    principal_stress_3 = 0.7
    trend_s3 = 0
    plunge_s3 = 0
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureVG
    m = 0.38
    alpha = 0.0000000005 # MPa^-1
  []
[]

# move elements between subdomains back and forth
[MeshModifiers]
  [StagedSubdomainModifier]
    type = StagedSubdomainModifier
  []
[]

[Stages]

  [Stage0]
    t = 0.0
    [Stage0_FunctionValueChanges]
      type = 'StagedFunctionValueChange'
      function_names = 'Tunnel01_elasticity_tensor_prefactor tunnel01_permeability_prefactor'
      new_values = '1.0                                  1.0                            '
    []
  []

  [Stage1]
    t = 1.0
    [Stage1_FunctionValueChanges]
      type = 'StagedFunctionValueChange'
      function_names = 'Tunnel01_elasticity_tensor_prefactor tunnel01_permeability_prefactor'
      new_values = '0.02                                 0.02                           '
      start_time = 't - 0.1'
      end_time = 't - 0.0001'
    []
    [Stage1_SubdomainModification]
      type = 'StagedSubdomainModification'
      from = 'tunnel01'
      to = 'tunnel01_inactive'
    []
  []
[]

[Preconditioning]
  [.\SMP]
    type = SMP
    full = true
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  petsc_options = '-snes_converged_reason'

  # best overall
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'

  line_search = none

  nl_abs_tol = 1e-4
  nl_rel_tol = 1e-6

  l_max_its = 20
  nl_max_its = 8

  start_time = 0.0
  end_time = 1.0 #17
  [TimeSteppers]
    # [ConstantDT1]
    #   type = ConstantDT
    #   dt = 0.25
    # []
    [StagedTimeSequenceStepper1]
      type = StagedTimeSequenceStepper
    []
  []

  [Quadrature]
    type = SIMPSON
    order = SECOND
  []
[]

[Outputs]
  perf_graph = true

  #exodus = true
  [out]
    type = Exodus
    execute_on = 'TIMESTEP_END'
  []

  # checkpoint = false
  #[out]
  #  type = Checkpoint
  #  num_files = 2
  #  wall_time_interval = 600
  #[]
[]
