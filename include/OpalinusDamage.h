//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ScalarDamageBase.h"



/**
 * Scalar damage model for which the damage is prescribed by another material
 */
class OpalinusDamage : public ScalarDamageBase
{
public:
  static InputParameters validParams();

  OpalinusDamage(const InputParameters & parameters);

protected:
  virtual void updateQpDamageIndex() override;
  const MaterialProperty<RankTwoTensor> & _total_strain;
  ///@{ Material property that provides the damage index
  const MaterialProperty<std::vector<Real>> & _intnl;
  //const MaterialProperty<Real> & _dam_cri;
  /// Rotation increment material property
  ///@}
   const Real _pd1;
  const Real _pd2;
  const Real _pd3;
  const Real _pd4;
  const Real _omega;
  const VariableValue & _nonlocal_var;

};
