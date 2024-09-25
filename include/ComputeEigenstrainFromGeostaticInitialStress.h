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
#include "CartesianLocalCoordinateSystem.h"

/**
 * ComputeEigenstrain computes an Eigenstrain that results from an initial stress.
 * The initial stress is defined in terms of a local cartesian coordinate system
 * and local stress components.
 */
class ComputeEigenstrainFromGeostaticInitialStress : public ComputeEigenstrainBase
{
public:
  static InputParameters validParams();

  ComputeEigenstrainFromGeostaticInitialStress(const InputParameters & parameters);

protected:
  virtual void computeQpEigenstrain() override;

  /// base_name for elasticity tensor to use to convert stress to strain
  const std::string _base_name;

  /// elasticity tensor used to convert stress to strain
  const MaterialProperty<RankFourTensor> & _elasticity_tensor;

  /// Stores the total eigenstrain in the previous step
  const MaterialProperty<RankTwoTensor> & _eigenstrain_old;

private:

  const CartesianLocalCoordinateSystem & _localCoordinateSystem;

  const Real _s1;
  const Real _s2;
  const Real _s3;

  const Real _s1_inc_x;
  const Real _s1_inc_y;
  const Real _s1_inc_z;

  const Real _s2_inc_x;
  const Real _s2_inc_y;
  const Real _s2_inc_z;

  const Real _s3_inc_x;
  const Real _s3_inc_y;
  const Real _s3_inc_z;

};
