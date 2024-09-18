//* This file calculates the elastic tensor for transversely isotropic rocks (Opalinus Clay)
//* Two geological rotation angles are needed to define the orientation of bedding planes.
//* dip_direction = rotation around global Z
//* dip= rotation around new x

#include "OpalinusElasticityTensor.h"
#include "RotationTensor.h"
#include "Function.h"


registerMooseObject(MOOSEAPPNAME, OpalinusElasticityTensor);


InputParameters
OpalinusElasticityTensor::validParams()
{
  InputParameters params = Material::validParams();
  params.addParam<FunctionName>(
      "elasticity_tensor_prefactor",
      "Optional function to use as a scalar prefactor on the elasticity tensor.");
  params.addParam<std::string>("base_name",
                               "Optional parameter that allows the user to define "
                               "multiple mechanics material systems on the same "
                               "block, i.e. for multiple phases");
  params.addClassDescription("Compute an elasticity tensor transversely isotropic material, such as shale.");
  params.addParam<Real>("dip_direction", 0.0, "clock_wise rotation of around z-axis. For dip_direction =0, the y-axis is assumed to point NORTH");
  params.addParam<Real>("dip", 0.0, "counter clock-wise rotation around x-axis. The y-axis is assumed to point NORTH");
  params.addRequiredParam<Real>("youngs_modulus_in_plane", "Module of elasticity P-sample (loading parallel to bedding)");
  params.addRequiredParam<Real>("youngs_modulus_normal", "Module of elasticity S-sample (loading normal to bedding)");
  params.addRequiredParam<Real>("poisson_ratio_in_plane", "Poisson's ratio P-sample (loading parallel to bedding)");
  params.addRequiredParam<Real>("poisson_ratio_normal", "Poisson's ratio S-sample (loading normal to bedding)");
  params.addRequiredParam<Real>("shear_module_normal", "Shear module (out of plane)");
 
  return params;
}


OpalinusElasticityTensor::OpalinusElasticityTensor(
    const InputParameters & parameters)
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
  _geological_angles(this->template getParam<Real>("dip_direction"),
                       this->template getParam<Real>("dip")),
  _first_local_axis(declareProperty<std::vector<Real>>("first_local_vector")),
  _second_local_axis(declareProperty<std::vector<Real>>("second_local_vector")),
  _normal_local_axis(declareProperty<std::vector<Real>>("normal_local_vector"))


{


    std::vector<Real> input(9);
    input[0]=_Ep*(1.0-(_Ep/_Es)*std::pow(_nu_s,2.0))/(1.0+_nu_p)/
              (1.0-_nu_p-2.0*(_Ep/_Es)*std::pow(_nu_s,2.0));
    input[8]=_Ep/(2.0*(1.0+_nu_p));
    input[1]=input[0]-2.0*input[8];
    input[2]=_Ep*_nu_s/(1.0-_nu_p-2.0*(_Ep/_Es)*std::pow(_nu_s,2.0));
    input[3]=input[0];
    input[4]=input[2];
    input[5]=_Es*(1.0-_nu_p)/(1.0-_nu_p-2.0*(_Ep/_Es)*std::pow(_nu_s,2.0));
    input[6]=input[7]=_G_s;
    
    if ((_nu_p>1.0-2.0*std::pow(_nu_s,2)*_Ep/_Es))
    mooseError(
        "The components of elastic tensor violate the condition nu_p<=1-2*(nu_s)^2*(E_p/E_s)");

    _Cijkl.fillSymmetric9FromInputVector(input);
 
        const Real phi_1 = -1.0*(_geological_angles(0)) * (libMesh::pi / 180.0);
        const Real Phi = _geological_angles(1) * (libMesh::pi / 180.0);
        const Real c1 = std::cos(phi_1);
        const Real c2 = std::cos(Phi);
        const Real s1 = std::sin(phi_1);
        const Real s2 = std::sin(Phi);
     
      _geological_rotation(0,0)=c1;
      _geological_rotation(1,0)=s1;
      _geological_rotation(2,0)=0.0;
      _geological_rotation(0,1)=-c2*s1;
      _geological_rotation(1,1)=c1*c2;
      _geological_rotation(2,1)=s2;
      _geological_rotation(0,2)=s1*s2;
      _geological_rotation(1,2)=-c1*s2;
      _geological_rotation(2,2)=c2;

      _Cijkl.rotate(_geological_rotation);
  
  
}

void
OpalinusElasticityTensor::computeQpProperties()
{
  _first_local_axis[_qp].resize(3);
  _second_local_axis[_qp].resize(3);
  _normal_local_axis[_qp].resize(3);
  // Assign elasticity tensor at a given quad point
  _elasticity_tensor[_qp] = _Cijkl;

    // Multiply by prefactor
  if (_prefactor_function)
  {
    _elasticity_tensor[_qp] *= _prefactor_function->value(_t, _q_point[_qp]);
  }

  (_first_local_axis)[_qp][0]=_geological_rotation(0,0);
  (_first_local_axis)[_qp][1]=_geological_rotation(1,0);
  (_first_local_axis)[_qp][2]=_geological_rotation(2,0);
  (_second_local_axis)[_qp][0]=_geological_rotation(0,1);
  (_second_local_axis)[_qp][1]=_geological_rotation(1,1);
  (_second_local_axis)[_qp][2]=_geological_rotation(2,1);
  (_normal_local_axis)[_qp][0]=_geological_rotation(0,2);
  (_normal_local_axis)[_qp][1]=_geological_rotation(1,2);
  (_normal_local_axis)[_qp][2]=_geological_rotation(2,2);
}


