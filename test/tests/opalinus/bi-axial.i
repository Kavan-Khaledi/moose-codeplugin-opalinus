# Bi-Axial loading of a cube
# for details on the problem see
# https://communities.bentley.com/cfs-file/__key/communityserver-wikis-components-files/00-00-00-05-58/PlxValidation_2D00_Bi_2D00_axial_5F00_compression_5F00_test_5F00_with_5F00_Mohr_2D00_Coulomb_5F00_model_2D00_2018.pdf

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 1
  nz = 1
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 1.0
  zmin = 0.0
  zmax = 1.0
  elem_type = TET10

  # side sets:
  #   name     id   location
  #   back     0    z = zmin
  #   bottom   1    y = ymin
  #   right    2    x = xmax
  #   top      3    y = ymax
  #   left     4    x = xmin
  #   front    5    z = zmax
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
[]

[Physics]

  [SolidMechanics]

    [QuasiStatic]
      [all]
        add_variables = false
        incremental = true
        eigenstrain_names = ini_stress
        generate_output = 'max_principal_stress mid_principal_stress min_principal_stress stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz'
      []
    []
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
  [total_strain_zz]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]

  [total_strain_zz]
    type = RankTwoAux
    rank_two_tensor = total_strain
    variable = total_strain_zz
    index_i = 2
    index_j = 2
    execute_on = timestep_end
  []

  [p]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = p
    scalar_type = hydrostatic
  []

  [q]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    variable = q
    scalar_type = vonMisesStress
  []
[]

[BCs]

  [bottom_normal]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    value = 0
  []

  [top_normal]
    type = DirichletBC
    variable = disp_y
    boundary = 'top'
    value = 0
  []

  [back_front_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'back front'
    value = 0
  []

  [left_normal]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []

  [back_normal]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []

  [right]
    type = Pressure
    boundary = 'right'
    variable = 'disp_x'
    function = '1'
  []
  [top_displacement]
    type = FunctionDirichletBC
    variable = disp_z
    boundary = 'front'
    function = 'if(t<=0.4,0,-0.025*(t-0.4))'
  []
[]

[UserObjects]
  [ucsInitialStress]
    type = CartesianLocalCoordinateSystem
    e1 = '1 0 0'
    e2 = '0 1 0'
  []
[]

[Materials]

  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1000
    poissons_ratio = 0.25
  []

  [PerfectPlastic]
    type = OpalinusPerfectPlasticStressUpdate
    gama_mean = 1
    p_tensile = 1
    local_coordinate_system = 'ucsInitialStress'
    smoothing_tol = 0.0
    tip_smoother = 0.1

    yield_function_tol = 1.0E-5
    min_step_size = 0.004
    max_NR_iterations = 40
  []

  [ini_stress]
    type = ComputeEigenstrainFromInitialStress
    eigenstrain_name = ini_stress
    initial_stress = '-1 0 0  0 0 0  0 0 -1'
  []

  [stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = 'PerfectPlastic'
    perform_finite_strain_rotations = false
    tangent_operator = 'nonlinear'
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
  #automatic_scaling = true
  #compute_scaling_once=false
  solve_type = 'NEWTON'

  petsc_options = '-snes_converged_reason'
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'

  line_search = none

  nl_abs_tol = 1e-3
  nl_rel_tol = 1e-10
  nl_max_its = 20

  start_time = 0.0
  dt = 0.25
  end_time = 2
[]

[Outputs]
  exodus = true
[]
