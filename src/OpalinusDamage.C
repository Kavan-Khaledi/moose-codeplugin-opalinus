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

//defineLegacyParams(OpalinusDamage);

InputParameters
OpalinusDamage::validParams()
{
  InputParameters params = ScalarDamageBase::validParams();
  params.addClassDescription(
      "Scalar damage model for which the damage is prescribed by another material");
  params.addParam<Real>("pd1",0.0, "pd1");
  params.addParam<Real>("pd2",10000, "pd2");
  params.addParam<Real>("pd3",1.0, "pd3");
  params.addParam<Real>("pd4",1.0, "pd4");
  params.addParam<Real>("omega",1.0, "omega");
  params.addRequiredCoupledVar("nonlocal_variable",
                               "This is non local variable for damage");

  return params;
}

OpalinusDamage::OpalinusDamage(const InputParameters & parameters)
  : ScalarDamageBase(parameters),
  _total_strain(getMaterialPropertyByName<RankTwoTensor>
               (_base_name+"total_strain")),
  _intnl(getMaterialPropertyByName<std::vector<Real>>
         (_base_name+"plastic_internal_parameter")),
  _pd1(getParam<Real>("pd1")),
  _pd2(getParam<Real>("pd2")),
  _pd3(getParam<Real>("pd3")),
  _pd4(getParam<Real>("pd4")),
  _omega(getParam<Real>("omega")),
 _nonlocal_var(coupledValue("nonlocal_variable"))

{
}

void
OpalinusDamage::updateQpDamageIndex()
{
   Real kesi=_nonlocal_var[_qp];
  // Real kesi=-_total_strain[_qp](2,2);
   Real pd1=_pd1;

  if (kesi<=pd1)
  {
    _damage_index[_qp]=0.0;
     _damage_index[_qp] =std::max(_damage_index[_qp],_damage_index_old[_qp]);
  }else{
   _damage_index[_qp] = std::min(_omega, _omega*std::max(0.0, 1.0-std::exp(-_pd3*
                      pow((kesi -pd1)/_pd2,_pd4))));
   _damage_index[_qp] =std::max(_damage_index[_qp],_damage_index_old[_qp]);
    ///_damage_index[_qp]=0.0;
  }


  if (MooseUtils::absoluteFuzzyLessThan(_damage_index[_qp], 0.0) ||
      MooseUtils::absoluteFuzzyGreaterThan(_damage_index[_qp], 1.0))
    mooseError(_base_name + "damage_index ",
               "must be between 0 and 1. Current value is: ",
               _damage_index[_qp]);
}
