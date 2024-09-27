//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeEigenstrainFromGeostaticInitialStress.h"
#include "RankTwoTensor.h"
#include "Function.h"
#include "Conversion.h" // for stringify
#include "CartesianLocalCoordinateSystem.h"

registerMooseObject(MOOSEAPPNAME, ComputeEigenstrainFromGeostaticInitialStress);

InputParameters
ComputeEigenstrainFromGeostaticInitialStress::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes an eigenstrain from an initial stress");
  params.addParam<std::string>("base_name",
                               "The base_name for the elasticity tensor that will be "
                               "used to compute strain from stress.  Do not provide "
                               "any base_name if your elasticity tensor does not use "
                               "one.");

  params.addRequiredParam<UserObjectName>(
      "local_coordinate_system", "The UserObject that defines the local coordinate system. ");

  params.addParam<Real>(
      "principal_stress_1",
      0.0,
      "Value of principal stress in direction of local axis e1 (pressure is positive). "
      "Often this is assumed to be the maximum principal stress.");
  params.addParam<Real>(
      "principal_stress_2",
      0.0,
      "Value of principal stress in direction of local axis e2 (pressure is positive).");
  params.addParam<Real>(
      "principal_stress_3",
      0.0,
      "Value of principal stress in direction of local axis e3 (pressure is positive). "
      "Often this is assumed to be the minimum principal stress.");

  params.addParam<Real>("stress_1_increment_x",
                        0.0,
                        "Increment of stress 1 in global x direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");
  params.addParam<Real>("stress_1_increment_y",
                        0.0,
                        "Increment of stress 1 in global y direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");
  params.addParam<Real>("stress_1_increment_z",
                        0.0,
                        "Increment of stress 1 in global z direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");

  params.addParam<Real>("stress_2_increment_x",
                        0.0,
                        "Increment of stress 2 in global x direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");
  params.addParam<Real>("stress_2_increment_y",
                        0.0,
                        "Increment of stress 2 in global y direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");
  params.addParam<Real>("stress_2_increment_z",
                        0.0,
                        "Increment of stress 2 in global z direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");

  params.addParam<Real>("stress_3_increment_x",
                        0.0,
                        "Increment of stress 3 in global x direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");
  params.addParam<Real>("stress_3_increment_y",
                        0.0,
                        "Increment of stress 3 in global y direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");
  params.addParam<Real>("stress_3_increment_z",
                        0.0,
                        "Increment of stress 3 in global z direction. The increment uses the "
                        "origin of the coordinate system provided via 'local_coordinate_system'.");

  return params;
}

ComputeEigenstrainFromGeostaticInitialStress::ComputeEigenstrainFromGeostaticInitialStress(
    const InputParameters & parameters)
  : ComputeEigenstrainBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _eigenstrain_old(getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_name)),
    _localCoordinateSystem(
        getUserObject<CartesianLocalCoordinateSystem>("local_coordinate_system")),
    _s1(parameters.get<Real>("principal_stress_1")),
    _s2(parameters.get<Real>("principal_stress_2")),
    _s3(parameters.get<Real>("principal_stress_3")),
    _s1_inc_x(parameters.get<Real>("stress_1_increment_x")),
    _s1_inc_y(parameters.get<Real>("stress_1_increment_y")),
    _s1_inc_z(parameters.get<Real>("stress_1_increment_z")),
    _s2_inc_x(parameters.get<Real>("stress_2_increment_x")),
    _s2_inc_y(parameters.get<Real>("stress_2_increment_y")),
    _s2_inc_z(parameters.get<Real>("stress_2_increment_z")),
    _s3_inc_x(parameters.get<Real>("stress_3_increment_x")),
    _s3_inc_y(parameters.get<Real>("stress_3_increment_y")),
    _s3_inc_z(parameters.get<Real>("stress_3_increment_z"))
{
}

void
ComputeEigenstrainFromGeostaticInitialStress::computeQpEigenstrain()
{
  if (_t_step == 1)
  {
    const auto x = _q_point[_qp](0);
    const auto y = _q_point[_qp](1);
    const auto z = _q_point[_qp](2);

    const auto bs1_inc_x = _s1_inc_x != 0;
    const auto bs1_inc_y = _s1_inc_y != 0;
    const auto bs1_inc_z = _s1_inc_z != 0;

    const auto bs2_inc_x = _s2_inc_x != 0;
    const auto bs2_inc_y = _s2_inc_y != 0;
    const auto bs2_inc_z = _s2_inc_z != 0;

    const auto bs3_inc_x = _s3_inc_x != 0;
    const auto bs3_inc_y = _s3_inc_y != 0;
    const auto bs3_inc_z = _s3_inc_z != 0;

    const auto x0 = (bs1_inc_x || bs2_inc_x || bs3_inc_x) ? _localCoordinateSystem.origin(0) : 0.0;
    const auto y0 = (bs1_inc_y || bs2_inc_y || bs3_inc_y) ? _localCoordinateSystem.origin(1) : 0.0;
    const auto z0 = (bs1_inc_z || bs2_inc_z || bs3_inc_z) ? _localCoordinateSystem.origin(2) : 0.0;

    auto s1 = _s1;
    if (bs1_inc_x)
      s1 += (x - x0) * _s1_inc_x;
    if (bs1_inc_y)
      s1 += (y - y0) * _s1_inc_y;
    if (bs1_inc_z)
      s1 += (z - z0) * _s1_inc_z;

    auto s2 = _s2;
    if (bs2_inc_x)
      s2 += (x - x0) * _s2_inc_x;
    if (bs2_inc_y)
      s2 += (y - y0) * _s2_inc_y;
    if (bs2_inc_z)
      s2 += (z - z0) * _s2_inc_z;

    auto s3 = _s3;
    if (bs3_inc_x)
      s3 += (x - x0) * _s3_inc_x;
    if (bs3_inc_y)
      s3 += (y - y0) * _s3_inc_y;
    if (bs3_inc_z)
      s3 += (z - z0) * _s3_inc_z;

    auto initial_stress = RankTwoTensor(s1, 0, 0,
                                        0, s2, 0,
                                        0, 0, s3);

    _localCoordinateSystem.rotateLocalToGlobal(&initial_stress);

    _eigenstrain[_qp] = _elasticity_tensor[_qp].invSymm() * initial_stress;
  }
  else
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
}
