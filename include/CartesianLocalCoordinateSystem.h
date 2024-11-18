//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "UserObject.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "libmesh/point.h"

/**
 * Represents a cartesian local coordinate system.
 * By definition, the local axis e1, e2 and e3 are perpendicular to each other.
 * The local coordinate system uses the right hand rule.
 */
class CartesianLocalCoordinateSystem : public UserObject
{
public:
  static InputParameters validParams();

  CartesianLocalCoordinateSystem(const InputParameters & parameters);

  virtual void initialize() override {}
  virtual void finalize() override {}
  virtual void execute() override {}

  // This object is not threaded
  void threadJoin(const UserObject &) final {}

  //// Origin of the local coordinate system
  const Point origin;

  /// First base vector of the local coordinate system
  const RealVectorValue e1() const { return _e1; };

  /// Second base vector of the local coordinate system
  const RealVectorValue e2() const { return _e2; };

  //// Third base vector of the local coordinate system
  const RealVectorValue e3() const { return _e3; };

  void rotateGlobalToLocal(RankTwoTensor * t) const;
  void rotateGlobalToLocal(RankFourTensor * t) const;

  void rotateLocalToGlobal(RankTwoTensor * t) const;
  void rotateLocalToGlobal(RankFourTensor * t) const;

private:
  //// First base vector of the local coordinate system
  RealVectorValue _e1;

  //// Second base vector of the local coordinate system
  RealVectorValue _e2;

  //// Third base vector of the local coordinate system
  RealVectorValue _e3;

  //// rotation tensor rotating 'global to local'
  RankTwoTensor _rotate_global_to_local;

  //// rotation tensor rotating 'local to global'
  RankTwoTensor _rotate_local_to_global;

  // const RankTwoTensor _rotation_e1;
  // const RankTwoTensor _rotation_e2;
  // const RankTwoTensor _rotation_e3;

  void buildFromDip(const InputParameters & parameters);
  void buildFromTrendAndPlunge(const InputParameters & parameters);
  void buildFromVectorAndPoint(const InputParameters & parameters);
  void buildFromVectors(const InputParameters & parameters,
                        const RealVectorValue e1,
                        const RealVectorValue e2);
};
