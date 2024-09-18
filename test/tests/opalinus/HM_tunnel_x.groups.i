# Variables containing a set of block names representing groups


block_anchors = ''
BoundarySurfaces = 'B_xmax_b
                    B_xmax_c
                    B_xmin_a
                    B_xmin_d
                    B_zmin_c
                    B_zmin_d
                    rock_i01_f00a
                    rock_i01_f00b
                    rock_i01_f00c
                    rock_i01_f00d
                    rock_i01_f20a
                    rock_i01_f20b
                    rock_i01_f20c
                    rock_i01_f20d
                    rock_i02_f00a
                    rock_i02_f00b
                    rock_i02_f00c
                    rock_i02_f00d
                    rock_i02_f20a
                    rock_i02_f20b
                    rock_i02_f20c
                    rock_i02_f20d
                    tunnel_f00
                    tunnel_f20'

BoundarySurfaces_XMax = 'B_xmax_b
                         B_xmax_c'

BoundarySurfaces_XMin = 'B_xmin_a
                         B_xmin_d'

BoundarySurfaces_YMax = 'rock_i01_f20a
                         rock_i01_f20b
                         rock_i01_f20c
                         rock_i01_f20d
                         rock_i02_f20a
                         rock_i02_f20b
                         rock_i02_f20c
                         rock_i02_f20d
                         tunnel_f20'

BoundarySurfaces_YMin = 'rock_i01_f00a
                         rock_i01_f00b
                         rock_i01_f00c
                         rock_i01_f00d
                         rock_i02_f00a
                         rock_i02_f00b
                         rock_i02_f00c
                         rock_i02_f00d
                         tunnel_f00'

BoundarySurfaces_ZMin = 'B_zmin_c
                         B_zmin_d'

RockVolumes = 'rock01a
               rock01b
               rock01c
               rock01d
               rock02a
               rock02b
               rock02c
               rock02d'

TunnelVolumes = 'tunnel01
                 tunnel02
                 tunnel03
                 tunnel04
                 tunnel05
                 tunnel06
                 tunnel07
                 tunnel08
                 tunnel09
                 tunnel10
                 tunnel11
                 tunnel12
                 tunnel13
                 tunnel14
                 tunnel15
                 tunnel16
                 tunnel17
                 tunnel18
                 tunnel19
                 tunnel20'

XMaxSurfaces = 'B_xmax_b
                B_xmax_c'

XMinSurfaces = 'B_xmin_a
                B_xmin_d'

YMaxSurfaces = 'rock_i01_f20a
                rock_i01_f20b
                rock_i01_f20c
                rock_i01_f20d
                rock_i02_f20a
                rock_i02_f20b
                rock_i02_f20c
                rock_i02_f20d
                tunnel_f20'

YMinSurfaces = 'rock_i01_f00a
                rock_i01_f00b
                rock_i01_f00c
                rock_i01_f00d
                rock_i02_f00a
                rock_i02_f00b
                rock_i02_f00c
                rock_i02_f00d
                tunnel_f00'

ZMaxSurfaces = 'B_zmax_a
                B_zmax_d'

ZMinSurfaces = 'B_zmin_c
                B_zmin_d'

# Fake users of the variables containing a set of block names representing groups

[Functions]
	[FakeUser_block_anchors]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${block_anchors}
	[]
	[FakeUser_BoundarySurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${BoundarySurfaces}
	[]
	[FakeUser_BoundarySurfaces_XMax]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${BoundarySurfaces_XMax}
	[]
	[FakeUser_BoundarySurfaces_XMin]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${BoundarySurfaces_XMin}
	[]
	[FakeUser_BoundarySurfaces_YMax]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${BoundarySurfaces_YMax}
	[]
	[FakeUser_BoundarySurfaces_YMin]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${BoundarySurfaces_YMin}
	[]
	[FakeUser_BoundarySurfaces_ZMin]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${BoundarySurfaces_ZMin}
	[]
	[FakeUser_RockVolumes]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${RockVolumes}
	[]
	[FakeUser_TunnelVolumes]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${TunnelVolumes}
	[]
	[FakeUser_XMaxSurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${XMaxSurfaces}
	[]
	[FakeUser_XMinSurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${XMinSurfaces}
	[]
	[FakeUser_YMaxSurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${YMaxSurfaces}
	[]
	[FakeUser_YMinSurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${YMinSurfaces}
	[]
	[FakeUser_ZMaxSurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${ZMaxSurfaces}
	[]
	[FakeUser_ZMinSurfaces]
		type = ParsedFunction
		expression = 'a'
		symbol_names = 'a'
		symbol_values = '1'
		control_tags = ${ZMinSurfaces}
	[]
[]
