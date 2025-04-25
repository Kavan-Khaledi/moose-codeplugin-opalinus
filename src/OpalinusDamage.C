////* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
///** https://www.gnu.org/licenses/lgpl-2.1.html

#include "OpalinusDamage.h"

registerMooseObject(MOOSEAPPNAME, OpalinusDamage);

// defineLegacyParams(OpalinusDamage);

InputParameters
OpalinusDamage::validParams()
{
  InputParameters params = ScalarDamageBase::validParams();
  params.addClassDescription(
      "Scalar damage model for which the damage is prescribed by another material");
  params.addParam<Real>("parameter_damageI",
                        0.0,
                        "Onset of damage. This value is compared to the value of the "
                        "'nonlocal_variable' and if equal or exceeded, damage is considered.");
  params.addParam<Real>("parameter_damageF", 10000, "parameter_damageF");
  params.addParam<Real>("parameter_damageA", 1.0, "parameter_damageA");
  params.addParam<Real>("parameter_damageN", 1.0, "parameter_damageN");
  params.addParam<Real>("omega", 1.0, "omega");
  params.addRequiredCoupledVar("nonlocal_variable", "This is the non-local variable for damage.");

  return params;
}

OpalinusDamage::OpalinusDamage(const InputParameters & parameters)
  : ScalarDamageBase(parameters),
    _total_strain(getMaterialPropertyByName<RankTwoTensor>(_base_name + "total_strain")),
    _intnl(getMaterialPropertyByName<std::vector<Real>>(_base_name + "plastic_internal_parameter")),
    _dam_I(getParam<Real>("parameter_damageI")),
    _dam_F(getParam<Real>("parameter_damageF")),
    _dam_A(getParam<Real>("parameter_damageA")),
    _dam_N(getParam<Real>("parameter_damageN")),
    _omega(getParam<Real>("omega")),
    _nonlocal_var(coupledValue("nonlocal_variable"))

{
}

void
OpalinusDamage::updateQpDamageIndex()
{
  Real kesi = _nonlocal_var[_qp];
  // Real kesi=-_total_strain[_qp](2,2);
  Real dam_I = _dam_I;

  if (kesi <= dam_I)
  {
    _damage_index[_qp] = 0.0;
    _damage_index[_qp] = std::max(_damage_index[_qp], _damage_index_old[_qp]);
  }
  else
  {
    _damage_index[_qp] = std::min(
        _omega,
        _omega * std::max(0.0, 1.0 - std::exp(-_dam_A * pow((kesi - dam_I) / _dam_F, _dam_N))));
    _damage_index[_qp] = std::max(_damage_index[_qp], _damage_index_old[_qp]);
    ///_damage_index[_qp]=0.0;
  }

  if (MooseUtils::absoluteFuzzyLessThan(_damage_index[_qp], 0.0) ||
      MooseUtils::absoluteFuzzyGreaterThan(_damage_index[_qp], 1.0))
    mooseError(_base_name + "damage_index ",
               "must be between 0 and 1. Current value is: ",
               _damage_index[_qp]);
}
