
[Mesh]
  [part1]
    type = ConcentricCircleMeshGenerator
    num_sectors = 8
    radii = '2.25 11'
    rings = '1 20 2'
    has_outer_square = on
    pitch = 50
    preserve_volumes = on
    smoothing_max_it = 10
  []

  [rotate]
    type = TransformGenerator
    input = part1
    transform = ROTATE
    vector_value = '0 90 0'
  []

  [part1_3d]
    type = MeshExtruderGenerator
    input = rotate
    extrusion_vector = '0 -20 0'
    num_layers = 40
    bottom_sideset = front
    top_sideset = back
  []

  [inner_ring]
    type = ParsedGenerateSideset
    input = part1_3d
    combinatorial_geometry = '((x*x+z*z)>4.9 & (x*x+z*z)<5.2)'
    new_sideset_name = inner_ring
  []

  # [convert]
  #   type = ElementsToTetrahedronsConverter
  #   input = inner_ring
  # []

  [subdomains1]
    type = ParsedSubdomainMeshGenerator
    input = inner_ring
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -19 & y >= -20'
    block_id = 101
  []
  [subdomains2]
    type = ParsedSubdomainMeshGenerator
    input = subdomains1
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -18 & y >= -19'
    block_id = 102
  []
  [subdomains3]
    type = ParsedSubdomainMeshGenerator
    input = subdomains2
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -17 & y >= -18'
    block_id = 103
  []
  [subdomains4]
    type = ParsedSubdomainMeshGenerator
    input = subdomains3
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -16 & y >= -17'
    block_id = 104
  []
  [subdomains5]
    type = ParsedSubdomainMeshGenerator
    input = subdomains4
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -15 & y >= -16'
    block_id = 105
  []
  [subdomains6]
    type = ParsedSubdomainMeshGenerator
    input = subdomains5
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -14 & y >= -15'
    block_id = 106
  []
  [subdomains7]
    type = ParsedSubdomainMeshGenerator
    input = subdomains6
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -13 & y >= -14'
    block_id = 107
  []
  [subdomains8]
    type = ParsedSubdomainMeshGenerator
    input = subdomains7
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -12 & y >= -13'
    block_id = 108
  []
  [subdomains9]
    type = ParsedSubdomainMeshGenerator
    input = subdomains8
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -11 & y >= -12'
    block_id = 109
  []
  [subdomains10]
    type = ParsedSubdomainMeshGenerator
    input = subdomains9
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -10 & y >= -11'
    block_id = 110
  []
  [subdomains11]
    type = ParsedSubdomainMeshGenerator
    input = subdomains10
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -9 & y >= -10'
    block_id = 111
  []
  [subdomains12]
    type = ParsedSubdomainMeshGenerator
    input = subdomains11
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -8 & y >= -9'
    block_id = 112
  []
  [subdomains13]
    type = ParsedSubdomainMeshGenerator
    input = subdomains12
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -7 & y >= -8'
    block_id = 113
  []
  [subdomains14]
    type = ParsedSubdomainMeshGenerator
    input = subdomains13
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -6 & y >= -7'
    block_id = 114
  []
  [subdomains15]
    type = ParsedSubdomainMeshGenerator
    input = subdomains14
    combinatorial_geometry = '(x*x+z*z)<5.1 & y < -5 & y >= -6'
    block_id = 115
  []
  # second_order = true
  add_subdomain_ids = '500'
  add_subdomain_names = 'empty'
[]
