#include "ImplicitNonlocal.h"


registerMooseObject(MOOSEAPPNAME, ImplicitNonlocal);

//defineLegacyParams(ImplicitNonlocal);

InputParameters
ImplicitNonlocal::validParams()
{
  InputParameters params = Kernel::validParams();
  params.addClassDescription("Calculates the Non_local equation for "
                             "implicit gradient plasticity");
  params.addRequiredParam<Real>("length_scale", "length_scale");
  return params;
}

ImplicitNonlocal::ImplicitNonlocal(const InputParameters & parameters) : Kernel(parameters), 

_length_scale(parameters.get<Real>("length_scale")),
_intnl(getMaterialProperty<std::vector<Real>>
         ("plastic_internal_parameter")),
_plastic_strain(getMaterialProperty<RankTwoTensor>
               ("plastic_strain")),
_total_strain(getMaterialProperty<RankTwoTensor>
               ("total_strain"))

 {
}

Real
ImplicitNonlocal::computeQpResidual()
{

Real kesi=(_plastic_strain[_qp]).trace();
  return _u[_qp]*_test[_i][_qp]+_length_scale*_grad_u[_qp] * _grad_test[_i][_qp]
         -kesi*_test[_i][_qp];
}

Real
ImplicitNonlocal::computeQpJacobian()
{
  return _phi[_j][_qp]*_test[_i][_qp]+_length_scale*_grad_phi[_j][_qp] * _grad_test[_i][_qp];
}

