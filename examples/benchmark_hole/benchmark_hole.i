#Benchmark test
# https://docs.itascacg.com/flac3d700/flac3d/zone/test3d/VerificationProblems/CylinderInMohrCoulomb/salencon.html?node3575

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Problem]
  kernel_coverage_check = 'SKIP_LIST'
  kernel_coverage_block_list = 'tunnel_inactive'
  material_coverage_check = 'SKIP_LIST'
  material_coverage_block_list = 'tunnel_inactive'
[]

[Mesh]
  [file]
    type = FileMeshGenerator
    file = benchmark_hole2.msh
  []
  second_order = true
  add_subdomain_names = tunnel_inactive
[]

[Variables]
  [disp_x]
    family = LAGRANGE
    order = SECOND
    block = 'rock tunnel'
  []
  [disp_y]
    family = LAGRANGE
    order = SECOND
    block = 'rock tunnel'
  []
  [disp_z]
    family = LAGRANGE
    order = SECOND
    block = 'rock tunnel'
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        # add_variables = true
        incremental = true
        eigenstrain_names = ini_stress
        block = 'rock tunnel'
      []
    []
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

  [internal_plastic_strain]
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
    block = 'rock tunnel'
  []

  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
    execute_on = timestep_end
    block = 'rock tunnel'
  []

  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
    block = 'rock tunnel'
  []

  [internal_plastic_strain]
    type = MaterialStdVectorAux
    property = plastic_internal_parameter
    index = 0
    variable = internal_plastic_strain
    block = 'rock tunnel'
  []
[]

[BCs]
  [no_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'left right top bottom'
    value = 0.0
  []

  [no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'left right top bottom'
    value = 0.0
  []
  [no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'front face back'
    value = 0.0
  []
[]

[UserObjects]
  [GlobalSubdomainModifier]
    type = TimedSubdomainModifier
    times = '2'
    blocks_from = 'tunnel'
    blocks_to = 'tunnel_inactive'
    execute_on = 'INITIAL TIMESTEP_BEGIN'
  []

  [ucsInitialStress]
    type = CartesianLocalCoordinateSystem
    e1 = '1 0 0'
    e2 = '0 1 0'
  []
  [ucsOpalinusMaterial]
    type = CartesianLocalCoordinateSystem
    e1 = '1 0 0'
    e2 = '0 1 0'
  []
[]

[Materials]

  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    poissons_ratio = 0.21
    youngs_modulus = 6777.9
    block = 'rock tunnel'
  []

  [ini_stress]
    type = ComputeEigenstrainFromGeostaticInitialStress
    eigenstrain_name = 'ini_stress'
    local_coordinate_system = 'ucsInitialStress'
    principal_stress_1 = 30
    principal_stress_2 = 12.6 #we assume plane strain condition at t=0, nu(sigma1+sigma3)
    principal_stress_3 = 30
    block = 'rock tunnel'
  []

  [opalinus]
    type = OpalinusPerfectPlasticStressUpdate
    gama_mean = 1.2
    p_tensile = 3
    local_coordinate_system = 'ucsInitialStress'
    smoothing_tol = 0.0
    tip_smoother = 1
    psi_to_phi = 1
    yield_function_tol = 1.0E-5
    max_NR_iterations = 30
    min_step_size = 0.04
  []

  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'opalinus'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
    block = 'rock tunnel'
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
  solve_type = 'NEWTON'
  line_search = none
  nl_abs_tol = 1e-6
  nl_rel_tol = 1e-10

  l_max_its = 30
  nl_max_its = 15

  start_time = 0.0
  dt = 1
  end_time = 3
[]

[Outputs]
  time_step_interval = 1
  print_linear_residuals = false
  exodus = true
[]
