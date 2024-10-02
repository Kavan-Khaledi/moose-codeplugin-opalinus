# Apply an initial stress that should be defined
# by ComputeEigenstrainFromGeostaticInitialStress,
# and do a transient step to check that nothing
# happens.
# For this test, the local coordinate system is
# not rotated with respect to the global system.
# Gravity is also active. Therefore the initial
# stress field in the vertical direction is not
# constant.

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 5
  ny = 5
  nz = 5
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

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
[]

[Variables]
  [disp_x]
    order = SECOND
    family = LAGRANGE
  []
  [disp_y]
    order = SECOND
    family = LAGRANGE
  []
  [disp_z]
    order = SECOND
    family = LAGRANGE
  []
[]

[Physics]
  [SolidMechanics]
    [QuasiStatic]
      [all]
        displacements = 'disp_x disp_y disp_z'
        eigenstrain_names = 'ini_stress'
        # generate_output = 'max_principal_stress mid_principal_stress min_principal_stress stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz'
      []
    []
  []
[]

[Kernels]
  [gravity]
    type = Gravity
    use_displaced_mesh = false
    variable = disp_z
    value = -9.81
  []
[]

[AuxVariables]
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_xz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yz]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_zz]
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
  []
  [stress_xy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xy
    index_i = 0
    index_j = 1
  []
  [stress_xz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xz
    index_i = 0
    index_j = 2
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [stress_yz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yz
    index_i = 1
    index_j = 2
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
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

  [left_normal]
    type = DirichletBC
    variable = disp_x
    boundary = 'left'
    value = 0
  []

  [right_normal]
    type = DirichletBC
    variable = disp_x
    boundary = 'right'
    value = 0
  []

  [back_normal]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []

  [front_normal]
    type = DirichletBC
    variable = disp_z
    boundary = 'front'
    value = 0
  []
[]

# Local coordinate systems (material and initial stress)
[UserObjects]
  [ucsInitialStress]
    type = CartesianLocalCoordinateSystem
    origin = '0 0 1'
    e1 = '1 0 0'
    e2 = '0 1 0'
  []
[]

[Materials]

  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1
    poissons_ratio = 0.25
  []

  [strain_from_initial_stress]
    type = ComputeEigenstrainFromGeostaticInitialStress
    eigenstrain_name = 'ini_stress'
    local_coordinate_system = 'ucsInitialStress'
    principal_stress_1 = 4.5
    principal_stress_2 = 2.5
    principal_stress_3 = 0.7
    stress_3_increment_z =  ${fparse 0.025 * -9.81} # density * gravity
  []

  [stress]
    type = ComputeLinearElasticStress
  []

  [density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 0.025
  []

[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'

  # best overall
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'

  line_search = none

  nl_abs_tol = 5e-8
  nl_rel_tol = 1e-8
  nl_max_its = 5

  l_max_its = 30
  l_abs_tol = 1e-8
  l_tol = 1e-08

  start_time = 0.0
  dt = 1
  end_time = 1
[]

[Outputs]
  exodus = true
[]
