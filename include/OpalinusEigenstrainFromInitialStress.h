//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ComputeEigenstrainBase.h"
#include "RankFourTensor.h"

/**
 * ComputeEigenstrain computes an Eigenstrain that results from an initial stress
 * The initial stress is defined in terms of Functions, which may be
 * multiplied by optional AuxVariables
 */
class OpalinusEigenstrainFromInitialStress : public ComputeEigenstrainBase
{
public:
  static InputParameters validParams();

  OpalinusEigenstrainFromInitialStress(const InputParameters & parameters);

protected:
  virtual void computeQpEigenstrain() override;

  /// base_name for elasticity tensor to use to convert stress to strain
  const std::string _base_name;

  /// elasticity tensor used to convert stress to strain
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  ///Stores the total eigenstrain in the previous step
  const MaterialProperty<RankTwoTensor> & _eigenstrain_old;

  const Real _s1;
  const Real _s2;
  const Real _s3;
  const Real _tr1;
  const Real _tr2;
  const Real _tr3;
  const Real _pl1;
  const Real _pl2;
  const Real _pl3;
  RankTwoTensor _rotation_s1;
  RankTwoTensor _rotation_s2;
  RankTwoTensor _rotation_s3;
};
