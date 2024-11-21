# This test re-assigns the elements of one subdomain to another subdomain.

[GlobalParams]
  time_unit = days
  displacements = 'disp_x disp_y disp_z'
  use_displaced_mesh = false
  PorousFlowDictator = dictator
[]

Box1_inactive_name = 'Box1_inactive'
inactive_domain_block_names = '${Box1_inactive_name}'

[Problem]
  solve = true
  kernel_coverage_check = SKIP_LIST
  kernel_coverage_block_list = '${inactive_domain_block_names}'
  material_coverage_check = SKIP_LIST
  material_coverage_block_list = '${inactive_domain_block_names}'
[]

[Mesh]
  [BaseMesh]
    type = GeneratedMeshGenerator
    subdomain_name = 'BaseMesh'
    elem_type = 'TET10'
    dim = 3
    nx = 3
    ny = 3
    nz = 2
    xmin = -3
    xmax = +3
    ymin = -3
    ymax = +3
    zmin = -2
    zmax = +2
  []

  [Box1]
    type = SubdomainBoundingBoxGenerator
    block_name = 'Box1'
    input = "BaseMesh"
    block_id = 1
    location = "INSIDE"
    bottom_left = "-1.0 -1.0 -0.0"
    top_right = "+1.0 +1.0 +2.0"
  []

  add_subdomain_names = '${inactive_domain_block_names}'
  active_block_names = 'BaseMesh Box1'
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
        generate_output = 'max_principal_stress mid_principal_stress min_principal_stress stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz'
        block = 'BaseMesh Box1'
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
    function = '0'
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
  [total_strain_zz]
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

  [total_strain_zz]
    type = RankTwoAux
    block = '${Mesh/active_block_names}'
    rank_two_tensor = total_strain
    variable = total_strain_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  []

  [p]
    type = RankTwoScalarAux
    block = '${Mesh/active_block_names}'
    rank_two_tensor = stress
    variable = p
    scalar_type = hydrostatic
  []

  [q]
    type = RankTwoScalarAux
    block = '${Mesh/active_block_names}'
    rank_two_tensor = stress
    variable = q
    scalar_type = vonMisesStress
  []

[]

# fix the lower model boundary in y and z direction
[BCs]
  [back_fix_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'back'
    value = 0.0
  []
  [back_fix_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0.0
  []
[]

# fix the left model boundary in x direction
[BCs]
  [left_fix_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0.0
  []
[]

# put some pressure on the right model boundary
[BCs]
  [right_Dirichlet]
    type = FunctionDirichletBC
    variable = disp_x
    boundary = 'right'
    function = right_pressure_function
  []
[]
[Functions]
  [right_pressure_function]
    type = ParsedFunction
    expression = '-0.001 * t'
  []
[]

# no flow at the outside
[BCs]
  [fix_pw]
    type = DirichletBC
    variable = 'porepressure'
    boundary = 'left right back front top bottom'
    value = '0'
  []
[]

[UserObjects]
  [ucsMaterial]
    type = CartesianLocalCoordinateSystem
    e1 = '1 0 0'
    e2 = '0 1 0'
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'porepressure disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
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

# Material: Volume Elements
[Materials]

  [elasticity_tensor]
    type = OpalinusElasticityTensor
    block = '${Mesh/active_block_names}'
    local_coordinate_system = 'ucsMaterial'
    youngs_modulus_in_plane = 11000
    youngs_modulus_normal = 6000
    poisson_ratio_in_plane = 0.15
    poisson_ratio_normal = 0.25
    shear_module_normal = 2000
  []

  [PerfectPlastic]
    type = OpalinusPerfectPlasticStressUpdate
    block = '${Mesh/active_block_names}'
    local_coordinate_system = 'ucsMaterial'
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

  [ini_stress]
    type = ComputeEigenstrainFromInitialStress
    block = '${Mesh/active_block_names}'
    eigenstrain_name = ini_stress
    initial_stress = '0 0 0  0 0 0  0 0 0'
  []

  [stress]
    type = ComputeMultipleInelasticStress
    block = '${Mesh/active_block_names}'
    inelastic_models = 'PerfectPlastic'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
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
    porepressure = 'porepressure'
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

  [permeability_bulk]
    type = OpalinusPermeabilityTensor
    block = '${Mesh/active_block_names}'
    permeability1 = 5e-19
    permeability2 = 5e-19
    permeability3 = 5e-20

    local_coordinate_system = 'ucsMaterial'
  []
[]

# move elements between subdomains back and forth
[MeshModifiers]
  [GlobalSubdomainModifier]
    type = TimedSubdomainModifier
    times = '0.2'
    blocks_from = 'Box1'
    blocks_to = 'Box1_inactive'
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

  solve_type = 'PJFNK' #'NEWTON'
  petsc_options = '-snes_converged_reason'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'

  line_search = none

  l_tol = 1E-5
  l_max_its = 20

  nl_abs_tol = 1E-3
  nl_rel_tol = 1e-7
  nl_max_its = 20

  end_time = 1.0
  dtmin = 0.001
  [TimeSteppers]
    [BlockEventTimeStepper]
      type = TimeSequenceStepper
      time_sequence = '0.05 0.1 0.2 0.4 1.0'
    []
  []

  [Quadrature]
    type = SIMPSON
    order = SECOND
  []
[]

[Outputs]
  perf_graph = true
  exodus = true
[]
