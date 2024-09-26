# Bi-Axial loading of a cube
# for details on the problem see
# https://communities.bentley.com/cfs-file/__key/communityserver-wikis-components-files/00-00-00-05-58/PlxValidation_2D00_Bi_2D00_axial_5F00_compression_5F00_test_5F00_with_5F00_Mohr_2D00_Coulomb_5F00_model_2D00_2018.pdf

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

[Physics/SolidMechanics/QuasiStatic]
  [./all]
    add_variables = true
    incremental = true
    generate_output = 'max_principal_stress mid_principal_stress min_principal_stress stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz'
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

  [back_normal]
    type = DirichletBC
    variable = disp_z
    boundary = 'back'
    value = 0
  []

  [Pressure]
    [right]
      boundary = right  #xmax
      function = 1000   #Newtons
    []
    [front]
      boundary = front  #zmax
      function = 1000*t #Newtons
    []
  []
[]

# [AuxVariables]
#   [./f0]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./f1]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./f2]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./iter]
#     order = CONSTANT
#     family = MONOMIAL
#   []
#   [./intnl]
#     order = CONSTANT
#     family = MONOMIAL
#   []
# []

# [AuxKernels]
#   [./f0_auxk]
#     type = MaterialStdVectorAux
#     property = plastic_yield_function
#     index = 0
#     variable = f0
#   []
#   [./f1_auxk]
#     type = MaterialStdVectorAux
#     property = plastic_yield_function
#     index = 1
#     variable = f1
#   []
#   [./f2_auxk]
#     type = MaterialStdVectorAux
#     property = plastic_yield_function
#     index = 2
#     variable = f2
#   []
#   [./iter]
#     type = MaterialRealAux
#     property = plastic_NR_iterations
#     variable = iter
#   []
#   [./intnl_auxk]
#     type = MaterialStdVectorAux
#     property = plastic_internal_parameter
#     index = 1
#     variable = intnl
#   []
# []

# [Postprocessors]
#   [./s_I]
#     type = PointValue
#     point = '0 0 0'
#     variable = max_principal_stress
#   []
#   [./s_II]
#     type = PointValue
#     point = '0 0 0'
#     variable = mid_principal_stress
#   []
#   [./s_III]
#     type = PointValue
#     point = '0 0 0'
#     variable = min_principal_stress
#   []
#   [./s_xx]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_xx
#   []
#   [./s_xy]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_xy
#   []
#   [./s_xz]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_xz
#   []
#   [./s_yy]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_yy
#   []
#   [./s_yz]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_yz
#   []
#   [./s_zz]
#     type = PointValue
#     point = '0 0 0'
#     variable = stress_zz
#   []
#   [./f0]
#     type = PointValue
#     point = '0 0 0'
#     variable = f0
#   []
#   [./f1]
#     type = PointValue
#     point = '0 0 0'
#     variable = f1
#   []
#   [./f2]
#     type = PointValue
#     point = '0 0 0'
#     variable = f2
#   []
#   [./iter]
#     type = PointValue
#     point = '0 0 0'
#     variable = iter
#   []
#   [./intnl]
#     type = PointValue
#     point = '0 0 0'
#     variable = intnl
#   []
# []

[UserObjects]
  [./ts]
    type = SolidMechanicsHardeningConstant
    value = 1E6
  []
  [./cs]
    type = SolidMechanicsHardeningConstant
    value = 1e6 #0.5
  []
  [./coh]
    type = SolidMechanicsHardeningConstant
    value = 1000 #1E6
  []
  [./angphi]
    type = SolidMechanicsHardeningConstant
    value = 30
    convert_to_radians = true
  []
  [./angpsi]
    type = SolidMechanicsHardeningConstant
    value = 1
    convert_to_radians = true
  []
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1000000
    poissons_ratio = 0.25
  []
  [./tensile]
    type = CappedMohrCoulombStressUpdate
    tensile_strength = ts
    compressive_strength = cs
    cohesion = coh
    friction_angle = angphi
    dilation_angle = angpsi
    smoothing_tol = 0.001
    yield_function_tol = 1.0E-12
  []
  [./stress]
    type = ComputeMultipleInelasticStress
    inelastic_models = tensile
    perform_finite_strain_rotations = false
  []
[]


[Executioner]
  type = Transient
  #automatic_scaling = true
  #compute_scaling_once=false
  solve_type = 'PJFNK'

  petsc_options = '-snes_converged_reason'

  # best overall
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = ' lu       mumps'

  line_search = none

  nl_abs_tol = 5e-5
  nl_rel_tol = 1e-6
  nl_max_its = 5

  l_max_its = 30
  l_abs_tol = 1e-6
  l_tol = 1e-05

  start_time = 0.0
  dt = .25
  end_time = 5 #10

[]

[Outputs]
  exodus = true
[]
