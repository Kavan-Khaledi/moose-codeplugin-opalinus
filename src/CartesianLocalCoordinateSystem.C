//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "CartesianLocalCoordinateSystem.h"
#include "UserObject.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "libmesh/point.h"
#include "Conversion.h" // for stringify

registerMooseObject(MOOSEAPPNAME, CartesianLocalCoordinateSystem);

InputParameters
CartesianLocalCoordinateSystem::validParams()
{
  InputParameters params = UserObject::validParams();
  params.addClassDescription(
      "Represents an cartesian local coordinate system to be used by other objects.");

  params.addParam<Point>("origin",
                         Point(0.0, 0.0, 0.0),
                         "Origin of the local coordinate system in global coordinates.");

  params.addParam<Real>(
      "dip_direction_degree",
      0.0,
      "Dip direction of the plane defined by 'dip_option' "
      "given as clock-wise rotation in degree around the (vertical) global z-axis and "
      "the global y-axis (assumed to point NORTH).");
  params.addParam<Real>(
      "dip_angle_degree",
      0.0,
      "Dip Angle of the plane defined by 'dip_option' measured in degree to the horizontal plane.");
  MooseEnum dip_option("e1_e2_plane_e1_horizontal e1_e2_plane_e2_horizontal",
                       "e1_e2_plane_e1_horizontal");
  params.addParam<MooseEnum>("dip_option",
                             dip_option,
                             "Defines the interpretation of the parameters 'dip_direction_degree' "
                             "and 'dip_angle_degree'.");

  params.addParam<Real>(
      "e1_trend_degree",
      0.0,
      "Trend of vector of the local e1-axis. "
      "Clock-wise angle in degree with respect to North. The y-axis is assumed to point NORTH.");
  params.addParam<Real>("e1_plunge_degree",
                        0.0,
                        "Plunge of vector of the local e1-axis. "
                        "Angle in degree with respect to the horizontal plane.");

  params.addParam<Real>(
      "e2_trend_degree",
      0.0,
      "Trend of vector of the local e2-axis. "
      "Clock-wise angle in degree with respect to North. The y-axis is assumed to point NORTH. "
      "If the resulting e2 vector is not perpendicular to e1, it will be corrected.");
  params.addParam<Real>(
      "e2_plunge_degree",
      0.0,
      "Plunge of vector of the local e2-axis. "
      "Angle in degree with respect to the horizontal plane. "
      "If the resulting e2 vector is not perpendicular to e1, it will be corrected.");

  params.addParam<RealVectorValue>("e1", "Vector of the local e1-axis in global coordinates.");

  params.addParam<RealVectorValue>(
      "e2",
      "Vector of the local e2-axis in global coordinates. "
      "If this vector is not perpendicular to e1, it will be corrected. "
      "This parameter is an alternative to 'point_on_e1_e2'. "
      "This parameter and parameter 'point_on_e1_e2' are mutual exclusive.");

  params.addParam<Point>("point_on_e1_e2",
                         "Point on the e1-e2-plane (quadrant where e1 and e2 are positive). "
                         "This parameter is an alternative to 'e2'. "
                         "This point must not be located on the e1 vector. "
                         "This parameter and parameter 'e2' are mutual exclusive.");

  params.addParam<bool>("normalize",
                        true,
                        "If set to true, the vectors e1, e2 and e3 will be normalized "
                        "(e.g. set to a length of unity in global coordinates)");

  return params;
}

CartesianLocalCoordinateSystem::CartesianLocalCoordinateSystem(const InputParameters & parameters)
  : UserObject(parameters), origin(parameters.get<Point>("origin"))
{
  // std::cout << "CartesianLocalCoordinateSystem: '" << name() << "'\n";

  // how does the user plan to provide the local coordinate system?
  const auto from_dip = parameters.isParamSetByUser("dip_direction_degree") +
                        parameters.isParamSetByUser("dip_angle_degree") +
                        parameters.isParamSetByUser("dip_option");
  const auto from_trend_and_plunge = parameters.isParamSetByUser("e1_trend_degree") +
                                     parameters.isParamSetByUser("e1_plunge_degree") +
                                     parameters.isParamSetByUser("e2_trend_degree") +
                                     parameters.isParamSetByUser("e2_plunge_degree");
  const auto from_vectors = parameters.isParamSetByUser("e1") + parameters.isParamSetByUser("e2") +
                            parameters.isParamSetByUser("point_on_e1_e2");

  // check parameter set for inconsistencies
  const auto from_source_count =
      ((from_dip > 0) + (from_trend_and_plunge > 0) + (from_vectors > 0));
  if (from_source_count != 1)
    mooseError(
        "The local coordinate system must be provided by one of the following options: "
        "(1) 'dip_direction_degree', 'dip_angle_degree' and 'dip_option'), "
        "or (2) via input of trend and plunge of e1 and e2 ('e1_trend_degree', 'e1_plunge_degree', "
        "'e2_trend_degree' and 'e2_plunge_degree'). "
        "or (3) via input of vector and point ('e1' and 'e2' or 'point_on_e1_e2').");

  if (from_dip > 0)
  {
    if (from_dip != 3)
      mooseError("All parameters 'dip_direction_degree', 'dip_angle_degree' and 'dip_option' must "
                 "be specified.");
    buildFromDip(parameters);
  }
  else if (from_trend_and_plunge > 0)
  {
    if (from_trend_and_plunge != 4)
      mooseError("Parameters 'e1_trend_degree', 'e1_plunge_degree', 'e2_trend_degree' and "
                 "'e2_plunge_degree' must be specified.");
    buildFromTrendAndPlunge(parameters);
  }
  else if (from_vectors > 0)
  {
    if (from_vectors != 2 || !parameters.isParamSetByUser("e1"))
      mooseError(
          "Parameter 'e1' must be given and one of the parameters 'e2' or 'point_on_e1_e2'.");
    buildFromVectorAndPoint(parameters);
  }
  else
    mooseError("Unknown data source. Are you missing a parameter? Did you misspell one?");

  // build rotation matrix global to local
  _rotate_global_to_local(0, 0) = _e1(0);
  _rotate_global_to_local(1, 0) = _e1(1);
  _rotate_global_to_local(2, 0) = _e1(2);
  _rotate_global_to_local(0, 1) = _e2(0);
  _rotate_global_to_local(1, 1) = _e2(1);
  _rotate_global_to_local(2, 1) = _e2(2);
  _rotate_global_to_local(0, 2) = _e3(0);
  _rotate_global_to_local(1, 2) = _e3(1);
  _rotate_global_to_local(2, 2) = _e3(2);

  _rotate_local_to_global = _rotate_global_to_local.inverse();
}

void
CartesianLocalCoordinateSystem::buildFromDip(const InputParameters & parameters)
{
  const auto dip_direction = parameters.get<Real>("dip_direction_degree") * (libMesh::pi / 180.0);
  const auto dip_angle = parameters.get<Real>("dip_angle_degree") * (libMesh::pi / 180.0);
  const std::string dip_option = parameters.get<MooseEnum>("dip_option");

  const Real s1 = std::sin(dip_direction);
  const Real c1 = std::cos(dip_direction);
  const Real s2 = std::sin(dip_angle);
  const Real c2 = std::cos(dip_angle);

  if (dip_option == "e1_e2_plane_e1_horizontal")
  {
    // e1 is horizontal by definition; e2 is parallel to the dip
    _e1 = RealVectorValue(c1, -s1, 0.0);
    _e2 = RealVectorValue(s1 * c2, c1 * c2, -s2);
    _e3 = _e1.cross(_e2);
  }
  else if (dip_option == "e1_e2_plane_e2_horizontal")
  {
    // e1 is parallel to the dip; e2 is horizontal by definition
    _e1 = RealVectorValue(s1 * c2, c1 * c2, -s2);
    _e2 = RealVectorValue(-c1, s1, 0.0);
    _e3 = _e1.cross(_e2);
  }
  else
    mooseError("Unknown value for 'dip_option': '" + dip_option + "'");
}

void
CartesianLocalCoordinateSystem::buildFromTrendAndPlunge(const InputParameters & parameters)
{
  const auto e1_trend = parameters.get<Real>("e1_trend_degree") * (libMesh::pi / 180.0);
  const auto e1_plunge = parameters.get<Real>("e1_plunge_degree") * (libMesh::pi / 180.0);
  const auto e2_trend = parameters.get<Real>("e2_trend_degree") * (libMesh::pi / 180.0);
  const auto e2_plunge = parameters.get<Real>("e2_plunge_degree") * (libMesh::pi / 180.0);

  const Real e1_s1 = std::sin(e1_trend);
  const Real e1_c1 = std::cos(e1_trend);
  const Real e1_s2 = std::sin(e1_plunge);
  const Real e1_c2 = std::cos(e1_plunge);
  const auto e1 = RealVectorValue(e1_s1 * e1_c2, e1_c1 * e1_c2, -e1_s2);

  const Real e2_s1 = std::sin(e2_trend);
  const Real e2_c1 = std::cos(e2_trend);
  const Real e2_s2 = std::sin(e2_plunge);
  const Real e2_c2 = std::cos(e2_plunge);
  const auto e2 = RealVectorValue(e2_s1 * e2_c2, e2_c1 * e2_c2, -e2_s2);

  buildFromVectors(parameters, e1, e2);
}

void
CartesianLocalCoordinateSystem::buildFromVectorAndPoint(const InputParameters & parameters)
{
  // collect data given by the user on e1
  const auto user_e1 = parameters.get<RealVectorValue>("e1");

  // collect data given by the user on e2
  const auto user_e2 = (parameters.isParamSetByUser("e2"))
                           ? parameters.get<RealVectorValue>("e2")
                           : RealVectorValue(parameters.get<Point>("point_on_e1_e2") - origin);

  buildFromVectors(parameters, user_e1, user_e2);
}

void
CartesianLocalCoordinateSystem::buildFromVectors(const InputParameters & parameters,
                                                 const RealVectorValue e1,
                                                 const RealVectorValue e2)
{
  // some sanity checking on e1 and e2
  const auto e1_length = e1.norm();
  if (e1_length < 1E-6)
    mooseError("The length of vector 'e1' is too close to Zero.");
  const auto e2_length = e2.norm();
  if (e2_length < 1E-6)
    mooseError("The length of vector 'e2' is too close to Zero.");

  // shall we normalize?
  const auto normalize = parameters.get<bool>("normalize");

  // set e1
  _e1 = (normalize && e1_length != 1) ? RealVectorValue(e1 / e1_length) : e1;

  // set e3 (including some sanity checks)
  _e3 = e1.cross(e2);
  const auto e3_length = _e3.norm();
  if (e3_length < 1E-6)
    mooseError("The vectors 'e1' and 'e2' seem to be parallel or anti-parallel.");
  if (normalize && e3_length != 1)
    _e3 /= e3_length;

  // set e2
  _e2 = _e3.cross(_e1);
}

void
CartesianLocalCoordinateSystem::rotateGlobalToLocal(RankTwoTensor * t) const
{
  t->rotate(_rotate_global_to_local);
}

void
CartesianLocalCoordinateSystem::rotateGlobalToLocal(RankFourTensor * t) const
{
  t->rotate(_rotate_global_to_local);
}

void
CartesianLocalCoordinateSystem::rotateLocalToGlobal(RankTwoTensor * t) const
{
  t->rotate(_rotate_local_to_global);
}

void
CartesianLocalCoordinateSystem::rotateLocalToGlobal(RankFourTensor * t) const
{
  t->rotate(_rotate_local_to_global);
}
