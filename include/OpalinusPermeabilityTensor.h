//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "PorousFlowPermeabilityBase.h"
#include "DerivativeMaterialInterface.h"
#include "Material.h"
#include "RankFourTensor.h"
#include "SymmetricRankFourTensor.h"
#include "GuaranteeProvider.h"
/**
 * Material designed to provide a constant permeability tensor
 */

class OpalinusPermeabilityTensor : public Material
{
public:
  static InputParameters validParams();

  OpalinusPermeabilityTensor(const InputParameters & parameters);

protected:
  void computeQpProperties() override;
   MaterialProperty<RealTensorValue> & _permeability_qp;
   MaterialProperty<std::vector<RealTensorValue>> & _dpermeability_qp_dvar;
   MaterialProperty<std::vector<std::vector<RealTensorValue>>> & _dpermeability_qp_dgradvar;
  /// Constant value of permeability tensor
  RealTensorValue _input_permeability;
  RealVectorValue _geological_angles;
  RankTwoTensor _geological_rotation;
  RankTwoTensor _unrotated_permeability;
  const Real _k_p;
  const Real _k_s;
  /// prefactor function to multiply the permeability tensor with
  const Function & _prefactor_f;

 
};


