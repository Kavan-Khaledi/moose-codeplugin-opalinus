//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OpalinusPermeabilityTensor.h"
#include "Function.h"

registerMooseObject(MOOSEAPPNAME, OpalinusPermeabilityTensor);

InputParameters
OpalinusPermeabilityTensor::validParams()
{
  InputParameters params = Material::validParams();

  params.addClassDescription(
      "This Material calculates the permeability tensor assuming it is constant");
  params.addParam<FunctionName>(
      "permeability_tensor_prefactor",
      "1.0",
      "Optional function to use as a scalar prefactor on the permeability tensor.");
  params.addParam<Real>("dip_direction",
                        0.0,
                        "clock_wise rotation around z-axis. The y-axis is assumed to point NORTH");
  params.addParam<Real>(
      "dip",
      0.0,
      "counter clock-wise rotation around x-axis. The y-axis is assumed to point NORTH");
  params.addRequiredParam<Real>("permeability_parallel",
                                "Permeability of P-sample (flow parallel to bedding), unit m2");
  params.addRequiredParam<Real>("permeability_normal",
                                "Permeability S-sample (loading normal to bedding), unit m2");
  return params;
}

OpalinusPermeabilityTensor::OpalinusPermeabilityTensor(const InputParameters & parameters)
  : Material(parameters),
    _permeability_qp(declareProperty<RealTensorValue>("PorousFlow_permeability_qp")),
    _dpermeability_qp_dvar(
        declareProperty<std::vector<RealTensorValue>>("dPorousFlow_permeability_qp_dvar")),
    _dpermeability_qp_dgradvar(declareProperty<std::vector<std::vector<RealTensorValue>>>(
        "dPorousFlow_permeability_qp_dgradvar")),
    _prefactor_f(getFunction("permeability_tensor_prefactor")),
    _k_s(parameters.get<Real>("permeability_normal")),
    _k_p(parameters.get<Real>("permeability_parallel")),
    _geological_angles(getParam<Real>("dip_direction"), getParam<Real>("dip"))
{

  _unrotated_permeability = RankTwoTensor(_k_p, 0, 0, 0, _k_p, 0, 0, 0, _k_s);
  const Real phi_1 = -1.0 * _geological_angles(0) * (libMesh::pi / 180.0);
  const Real Phi = _geological_angles(1) * (libMesh::pi / 180.0);
  const Real c1 = std::cos(phi_1);
  const Real c2 = std::cos(Phi);
  const Real s1 = std::sin(phi_1);
  const Real s2 = std::sin(Phi);

  _geological_rotation(0, 0) = c1;
  _geological_rotation(1, 0) = s1;
  _geological_rotation(2, 0) = 0.0;
  _geological_rotation(0, 1) = -c2 * s1;
  _geological_rotation(1, 1) = c1 * c2;
  _geological_rotation(2, 1) = s2;
  _geological_rotation(0, 2) = s1 * s2;
  _geological_rotation(1, 2) = -c1 * s2;
  _geological_rotation(2, 2) = c2;

  _unrotated_permeability.rotate(_geological_rotation);
  _input_permeability = RealTensorValue(_unrotated_permeability(0, 0),
                                        _unrotated_permeability(0, 1),
                                        _unrotated_permeability(0, 2),
                                        _unrotated_permeability(1, 0),
                                        _unrotated_permeability(1, 1),
                                        _unrotated_permeability(1, 2),
                                        _unrotated_permeability(2, 0),
                                        _unrotated_permeability(2, 1),
                                        _unrotated_permeability(2, 2));
}

void
OpalinusPermeabilityTensor::computeQpProperties()
{

  _permeability_qp[_qp] = _input_permeability;

  _permeability_qp[_qp] *= _prefactor_f.value(_t, _q_point[_qp]);
}
