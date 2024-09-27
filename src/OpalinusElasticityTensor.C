//* This file calculates the elastic tensor for transversely isotropic rocks (Opalinus Clay)
//* Two geological rotation angles are needed to define the orientation of bedding planes.
//* dip_direction = rotation around global Z
//* dip= rotation around new x

#include "OpalinusElasticityTensor.h"
#include "RotationTensor.h"
#include "Function.h"
#include "CartesianLocalCoordinateSystem.h"

registerMooseObject(MOOSEAPPNAME, OpalinusElasticityTensor);

InputParameters
OpalinusElasticityTensor::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription(
      "Compute an elasticity tensor transversely isotropic material, such as shale.");

  params.addParam<FunctionName>(
      "elasticity_tensor_prefactor",
      "Optional function to use as a scalar prefactor on the elasticity tensor.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");

  params.addRequiredParam<UserObjectName>(
      "local_coordinate_system",
      "The UserObject that defines the local coordinate system. "
      "The local axis e1 and e2 of this coordinate system are considered in plane while "
      "the local axis e3 is assumed to be 'normal'.");

  params.addRequiredParam<Real>("youngs_modulus_in_plane",
                                "Module of elasticity P-sample (loading parallel to bedding)");
  params.addRequiredParam<Real>("youngs_modulus_normal",
                                "Module of elasticity S-sample (loading normal to bedding)");
  params.addRequiredParam<Real>("poisson_ratio_in_plane",
                                "Poisson's ratio P-sample (loading parallel to bedding)");
  params.addRequiredParam<Real>("poisson_ratio_normal",
                                "Poisson's ratio S-sample (loading normal to bedding)");
  params.addRequiredParam<Real>("shear_module_normal", "Shear module (out of plane)");

  return params;
}

OpalinusElasticityTensor::OpalinusElasticityTensor(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _elasticity_tensor_name(_base_name + "elasticity_tensor"),
    _elasticity_tensor(declareProperty<RankFourTensor>(_elasticity_tensor_name)),
    _prefactor_function(isParamValid("elasticity_tensor_prefactor")
                            ? &getFunction("elasticity_tensor_prefactor")
                            : nullptr),
    _Ep(parameters.get<Real>("youngs_modulus_in_plane")),
    _Es(parameters.get<Real>("youngs_modulus_normal")),
    _nu_p(parameters.get<Real>("poisson_ratio_in_plane")),
    _nu_s(parameters.get<Real>("poisson_ratio_normal")),
    _G_s(parameters.get<Real>("shear_module_normal")),
    _localCoordinateSystem(
        getUserObject<CartesianLocalCoordinateSystem>("local_coordinate_system")),
    _first_local_axis(declareProperty<std::vector<Real>>("first_local_vector")),
    _second_local_axis(declareProperty<std::vector<Real>>("second_local_vector")),
    _normal_local_axis(declareProperty<std::vector<Real>>("normal_local_vector"))

{
  // Check validity of input parameters
  if ((_nu_p > 1.0 - 2.0 * std::pow(_nu_s, 2) * _Ep / _Es))
    mooseError("The components of elastic tensor violate the condition "
               "nu_p <= 1 - 2 * (nu_s)^2 * (E_p / E_s)");

  // build Cijkl in local coordinate system
  std::vector<Real> input(9);
  input[0] = _Ep * (1.0 - (_Ep / _Es) * std::pow(_nu_s, 2.0)) / (1.0 + _nu_p) /
             (1.0 - _nu_p - 2.0 * (_Ep / _Es) * std::pow(_nu_s, 2.0));
  input[8] = _Ep / (2.0 * (1.0 + _nu_p));
  input[1] = input[0] - 2.0 * input[8];
  input[2] = _Ep * _nu_s / (1.0 - _nu_p - 2.0 * (_Ep / _Es) * std::pow(_nu_s, 2.0));
  input[3] = input[0];
  input[4] = input[2];
  input[5] = _Es * (1.0 - _nu_p) / (1.0 - _nu_p - 2.0 * (_Ep / _Es) * std::pow(_nu_s, 2.0));
  input[6] = input[7] = _G_s;

  _Cijkl.fillSymmetric9FromInputVector(input);

  // rotate Cijkl to global coordinates
  _localCoordinateSystem.rotateLocalToGlobal(&_Cijkl);

  // local copy of the base vectors to speed up computeQpProperties
  const auto e1 = _localCoordinateSystem.e1();
  const auto e2 = _localCoordinateSystem.e2();
  const auto e3 = _localCoordinateSystem.e3();
  _e1 = std::vector<Real>({e1(0), e1(1), e1(2)});
  _e2 = std::vector<Real>({e2(0), e2(1), e2(2)});
  _e3 = std::vector<Real>({e3(0), e3(1), e3(2)});
}

void
OpalinusElasticityTensor::initQpStatefulProperties()
{
  computeQpProperties();
}

void
OpalinusElasticityTensor::computeQpProperties()
{

  // Assign elasticity tensor at a given quad point
  // with and without prefactor
  if (_prefactor_function)
    _elasticity_tensor[_qp] = _Cijkl * _prefactor_function->value(_t, _q_point[_qp]);
  else
    _elasticity_tensor[_qp] = _Cijkl;

  // set local axis (for output)
  _first_local_axis[_qp] = _e1;
  _second_local_axis[_qp] = _e2;
  _normal_local_axis[_qp] = _e3;
}
