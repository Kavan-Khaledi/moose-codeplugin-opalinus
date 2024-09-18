//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "OpalinusEigenstrainFromInitialStress.h"
#include "RankTwoTensor.h"
#include "Function.h"
#include "Conversion.h" // for stringify

registerMooseObject(MOOSEAPPNAME, OpalinusEigenstrainFromInitialStress);

InputParameters
OpalinusEigenstrainFromInitialStress::validParams()
{
  InputParameters params = ComputeEigenstrainBase::validParams();
  params.addClassDescription("Computes an eigenstrain from an initial stress");
  params.addParam<std::string>("base_name",
                               "The base_name for the elasticity tensor that will be "
                               "used to compute strain from stress.  Do not provide "
                               "any base_name if your elasticity tensor does not use "
                               "one.");
   params.addParam<Real>("principal_stress_1", 0.0, "magnitude of maximum principal stress");
   params.addParam<Real>("plunge_s1", 0.0, "Angle with respect to horizonal plane. The y-axis is assumed to point NORTH");
   params.addParam<Real>("trend_s1", 0.0, "Angle with respect to North. The y-axis is assumed to point NORTH");
   params.addParam<Real>("principal_stress_2", 0.0, "magnitude of second principal stress");
   params.addParam<Real>("plunge_s2", 0.0, "Angle with respect to horizonal plane. The y-axis is assumed to point NORTH");
   params.addParam<Real>("trend_s2", 0.0, "Angle with respect to North. The y-axis is assumed to point NORTH");
   params.addParam<Real>("principal_stress_3", 0.0, "magnitude of minimum principal stress");
   params.addParam<Real>("plunge_s3", 0.0, "Angle with respect to horizonal plane. The y-axis is assumed to point NORTH");
   params.addParam<Real>("trend_s3", 0.0, "Angle with respect to North. The y-axis is assumed to point NORTH");
  return params;
}

OpalinusEigenstrainFromInitialStress::OpalinusEigenstrainFromInitialStress(
    const InputParameters & parameters)
  : ComputeEigenstrainBase(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _eigenstrain_old(getMaterialPropertyOld<RankTwoTensor>(_eigenstrain_name)),
    _s1(parameters.get<Real>("principal_stress_1")),
    _s2(parameters.get<Real>("principal_stress_2")),
    _s3(parameters.get<Real>("principal_stress_3")),
    _tr1(parameters.get<Real>("trend_s1")),
    _tr2(parameters.get<Real>("trend_s2")),
    _tr3(parameters.get<Real>("trend_s3")),
    _pl1(parameters.get<Real>("plunge_s1")),
    _pl2(parameters.get<Real>("plunge_s2")),
    _pl3(parameters.get<Real>("plunge_s3"))
{

    Real phi_1 = _tr1 * (libMesh::pi / 180.0);
    Real Phi = _pl1 * (libMesh::pi / 180.0);
    Real c1 = std::cos(phi_1);
    Real c2 = std::cos(Phi);
    Real s1 = std::sin(phi_1);
    Real s2 = std::sin(Phi);
    _rotation_s1=RankTwoTensor(c2*c2*s1*s1,c2*c2*s1*c1,c2*s1*s2,c2*c2*s1*c1,c1*c1*c2*c2,c1*c2*s2,c2*s1*s2,c1*c2*s2,s2*s2);

    phi_1 = _tr2 * (libMesh::pi / 180.0);
    Phi = _pl2 * (libMesh::pi / 180.0);
    c1 = std::cos(phi_1);
    c2 = std::cos(Phi);
    s1 = std::sin(phi_1);
    s2 = std::sin(Phi);
    _rotation_s2=RankTwoTensor(c2*c2*s1*s1,c2*c2*s1*c1,c2*s1*s2,c2*c2*s1*c1,c1*c1*c2*c2,c1*c2*s2,c2*s1*s2,c1*c2*s2,s2*s2);
    
    phi_1 = _tr3 * (libMesh::pi / 180.0);
    Phi = _pl3 * (libMesh::pi / 180.0);
    c1 = std::cos(phi_1);
    c2 = std::cos(Phi);
    s1 = std::sin(phi_1);
    s2 = std::sin(Phi);
    _rotation_s3=RankTwoTensor(c2*c2*s1*s1,c2*c2*s1*c1,c2*s1*s2,c2*c2*s1*c1,c1*c1*c2*c2,c1*c2*s2,c2*s1*s2,c1*c2*s2,s2*s2);


}

void
OpalinusEigenstrainFromInitialStress::computeQpEigenstrain()
{
  if (_t_step == 1)
  {
    RankTwoTensor initial_stress;
    initial_stress=-_s1*_rotation_s1-_s2*_rotation_s2-_s3*_rotation_s3;

    _eigenstrain[_qp] = -_elasticity_tensor[_qp].invSymm() * initial_stress;
  }
  else
    _eigenstrain[_qp] = _eigenstrain_old[_qp];
}
