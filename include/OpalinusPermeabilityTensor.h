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
#include "CartesianLocalCoordinateSystem.h"

/**
 * Material designed to provide a constant permeability tensor
 * rotated according to a local coordinate system
 */

class OpalinusPermeabilityTensor : public Material
{
public:
  static InputParameters validParams();

  OpalinusPermeabilityTensor(const InputParameters & parameters);

protected:
  void computeQpProperties() override;

  const CartesianLocalCoordinateSystem & _localCoordinateSystem;

  /// permeabilities
  const Real _permeability1;
  const Real _permeability2;
  const Real _permeability3;

  /// prefactor functor and material property to multiply the permeability tensor with
  const Moose::Functor<Real> * const _prefactor_functor;
  const MaterialProperty<Real> * const _prefactor_matprop;

  /// quadpoint permeability
  MaterialProperty<RealTensorValue> & _permeability_qp;

  /// d(quadpoint permeability)/d(PorousFlow variable)
  MaterialProperty<std::vector<RealTensorValue>> & _dpermeability_qp_dvar;

  /// d(quadpoint permeability)/d(grad(PorousFlow variable))
  MaterialProperty<std::vector<std::vector<RealTensorValue>>> & _dpermeability_qp_dgradvar;

private:
  /// Constant value of permeability tensor
  RankTwoTensor _input_permeability;
};
