// PerfectlyPlastic Desai Model beta1=0.0 beta=0.0 alfa=0.0 The model is basically equivalent to
// Drucker-Prager model

#include "OpalinusPerfectPlasticStressUpdate.h"
#include "libmesh/utility.h"
#include "ElasticityTensorTools.h"
#include "CartesianLocalCoordinateSystem.h"

registerMooseObject(MOOSEAPPNAME, OpalinusPerfectPlasticStressUpdate);

InputParameters
OpalinusPerfectPlasticStressUpdate::validParams()
{
  InputParameters params = MultiParameterPlasticityStressUpdate::validParams();

  params.addClassDescription("Compute stress for Non-associated Plasticity for Geomaterials");

  params.addParam<bool>("perfect_guess",
                        true,
                        "Provide a guess to the Newton-Raphson procedure "
                        "that is the result from perfect plasticity.  With "
                        "severe hardening/softening this may be "
                        "suboptimal.");
  params.addRequiredParam<Real>("gama_mean", "Slope of the failure line in p-q-diagram.");
  params.addRequiredParam<Real>("p_tensile",
                                "Intersection of the failure with the p-axis in the p-q-diagram.");
  params.addParam<Real>(
      "parameter_omega_1",
      0.0,
      "Anisotropic plasticity parameter omega_1 according to S. Pietruszczak and Z. Mraz.");
  params.addParam<Real>(
      "parameter_b_1",
      0.0,
      "Anisotropic plasticity parameter b_1 according to S. Pietruszczak and Z. Mraz.");
  params.addRequiredParam<UserObjectName>(
      "local_coordinate_system",
      "The UserObject that defines the local coordinate system. "
      "The local axis e1 and e2 of this coordinate system are considered in plane while "
      "the local axis e3 is assumed to be 'normal'.");
  params.addParam<Real>("psi_to_phi", 1.0, "psi_to_phi");
  params.addParam<Real>("tip_smoother", 0.0, "tip_smoother");
  params.addParam<Real>("curvature_yield", 0.0, "Curvature of plastic yield in p-q space");
  params.addParam<Real>(
      "lode_angle_coefficient", 0.0, "lode angle dependency coefficient, between 0 and 0.7");
  params.addParam<Real>("Fs_function_power", -0.25, "power of Fs function");
  params.addParam<Real>("parameter_n0", 3.0, "parameter_n0");
  params.addParam<Real>("hardening_a0", 0.0, "hardening_a0");
  params.addParam<Real>("hardening_eta", 1.0, "hardening_eta");

  return params;
}

OpalinusPerfectPlasticStressUpdate::OpalinusPerfectPlasticStressUpdate(
    const InputParameters & parameters)
  : MultiParameterPlasticityStressUpdate(parameters, 6, 1, 1),
    _psi_to_phi(getParam<Real>("psi_to_phi")),
    _St(parameters.get<Real>("p_tensile")),
    _perfect_guess(getParam<bool>("perfect_guess")),
    _small_smoother2(Utility::pow<2>(getParam<Real>("tip_smoother"))),
    _beta(getParam<Real>("lode_angle_coefficient")),
    _betta1(getParam<Real>("curvature_yield")),
    _mv(getParam<Real>("Fs_function_power")),
    _mean_gama(parameters.get<Real>("gama_mean") / std::sqrt(27.0)),
    _omega1(getParam<Real>("parameter_omega_1")),
    _b1(getParam<Real>("parameter_b_1")),
    _localCoordinateSystem(
        getUserObject<CartesianLocalCoordinateSystem>("local_coordinate_system")),
    _a0(getParam<Real>("hardening_a0")),
    _eta(getParam<Real>("hardening_eta")),
    _n0(getParam<Real>("parameter_n0")),
    _gama(declareProperty<Real>(_base_name + "gama"))

{
}

void
OpalinusPerfectPlasticStressUpdate::computeStressParams(const RankTwoTensor & stress,
                                                        std::vector<Real> & stress_params) const
{
  stress_params[0] = stress(0, 0);
  stress_params[1] = stress(1, 1);
  stress_params[2] = stress(2, 2);
  stress_params[5] = stress(1, 2);
  stress_params[4] = stress(0, 2);
  stress_params[3] = stress(0, 1);
  setGamaValue(stress_params, _gama[_qp]);
}

std::vector<RankTwoTensor>
OpalinusPerfectPlasticStressUpdate::dstress_param_dstress(const RankTwoTensor & stress) const
{
  (void)stress;

  std::vector<RankTwoTensor> dsp(_num_sp);
  dsp[0] = RankTwoTensor(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  dsp[1] = RankTwoTensor(0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0);
  dsp[2] = RankTwoTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  dsp[5] = RankTwoTensor(0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0);
  dsp[4] = RankTwoTensor(0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0);
  dsp[3] = RankTwoTensor(0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0);

  return dsp;
}

std::vector<RankFourTensor>
OpalinusPerfectPlasticStressUpdate::d2stress_param_dstress(const RankTwoTensor & stress) const
{
  (void)stress;

  std::vector<RankFourTensor> d2(_num_sp);
  for (unsigned i = 0; i < _num_sp; ++i)
    d2[i] = RankFourTensor();
  return d2;
}

void
OpalinusPerfectPlasticStressUpdate::preReturnMapV(const std::vector<Real> & /*trial_stress_params*/,
                                                  const RankTwoTensor & stress_trial,
                                                  const std::vector<Real> & /*intnl_old*/,
                                                  const std::vector<Real> & /*yf*/,
                                                  const RankFourTensor & Eijkl)
{
  (void)stress_trial; // silence compiler warning about unused variable
  (void)Eijkl;        // silence compiler warning about unused variable
}

void
OpalinusPerfectPlasticStressUpdate::setStressAfterReturnV(const RankTwoTensor & /*stress_trial*/,
                                                          const std::vector<Real> & stress_params,
                                                          Real /*gaE*/,
                                                          const std::vector<Real> & /*intnl*/,
                                                          const yieldAndFlow & /*smoothed_q*/,
                                                          const RankFourTensor & /*Eijkl*/,
                                                          RankTwoTensor & stress) const
{
  stress = RankTwoTensor(stress_params[0],
                         stress_params[3],
                         stress_params[4],
                         stress_params[3],
                         stress_params[1],
                         stress_params[5],
                         stress_params[4],
                         stress_params[5],
                         stress_params[2]);
  return;
}

void
OpalinusPerfectPlasticStressUpdate::yieldFunctionValuesV(const std::vector<Real> & stress_params,
                                                         const std::vector<Real> & intnl,
                                                         std::vector<Real> & yf) const
{
  (void)intnl; // silence compiler warning about unused variable

  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = stress_now.trace() - 3.0 * _St;
  const Real j2 = std::max(1e-2, stress_now.secondInvariant());
  const Real j3 = stress_now.thirdInvariant();

  const Real q = std::sqrt(j2);
  const Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);
  const Real fs = std::pow(std::exp(_betta1 * i1) + _beta * sr, _mv);

  const Real fb = _gama[_qp] * i1;

  yf[0] = std::sqrt(std::pow(q, 2.0) + _small_smoother2) + fb * fs;

  return;
}

void
OpalinusPerfectPlasticStressUpdate::computeAllQV(const std::vector<Real> & stress_params,
                                                 const std::vector<Real> & intnl,
                                                 std::vector<yieldAndFlow> & all_q) const
{
  (void)intnl; // silence compiler warning about unused variable

  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = stress_now.trace() - 3.0 * _St;
  const Real j2 = std::max(1e-2, stress_now.secondInvariant());
  const Real j3 = stress_now.thirdInvariant();

  const Real q = std::sqrt(j2);
  const Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);
  const Real fs = std::pow(std::exp(_betta1 * i1) + _beta * sr, _mv);

  const Real _gamaq = _psi_to_phi * _gama[_qp];
  const Real fb = _gama[_qp] * i1;
  const Real fbg = _gamaq * i1;

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj2ds_test = stress_now.dsecondInvariant();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();
  const RankFourTensor d2j3ds_t = stress_now.d2thirdInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = 1.0;
  di1ds[1] = 1.0;
  di1ds[2] = 1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 < 1e-2) // @Kavan-Khaledi: this condition is never met (we make sure j2 is at least 1e-2)
  {
    dj2ds[0] = 0.0;
    dj2ds[1] = 0.0;
    dj2ds[2] = 0.0;
    dj2ds[5] = 0.0;
    dj2ds[4] = 0.0;
    dj2ds[3] = 0.0;

    dj3ds[0] = 0.0;
    dj3ds[1] = 0.0;
    dj3ds[2] = 0.0;
    dj3ds[5] = 0.0;
    dj3ds[4] = 0.0;
    dj3ds[3] = 0.0;
  }
  else
  {
    dj2ds[0] = dj2ds_t(0, 0);
    dj2ds[1] = dj2ds_t(1, 1);
    dj2ds[2] = dj2ds_t(2, 2);
    dj2ds[5] = dj2ds_t(1, 2);
    dj2ds[4] = dj2ds_t(0, 2);
    dj2ds[3] = dj2ds_t(0, 1);

    dj3ds[0] = dj3ds_t(0, 0);
    dj3ds[1] = dj3ds_t(1, 1);
    dj3ds[2] = dj3ds_t(2, 2);
    dj3ds[5] = dj3ds_t(1, 2);
    dj3ds[4] = dj3ds_t(0, 2);
    dj3ds[3] = dj3ds_t(0, 1);
  }

  const Real dfbdi1 = _gama[_qp];
  const Real dfbgdi1 = _gamaq;

  // const Real dfsdi1 = 0.0;   // if this variable is always zero, let's comment it out
  const Real dfdi1 = dfbdi1 * fs;  // + fb * dfsdi1;
  const Real dgdi1 = dfbgdi1 * fs; // + fbg * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dfdj2 = (1.0 / 2.0) / (std::sqrt(pow(q, 2.0) + _small_smoother2)) + fb * dfsdj2;
  const Real dgdj2 = (1.0 / 2.0) / (std::sqrt(pow(q, 2.0) + _small_smoother2)) + fbg * dfsdj2;
  const Real dfdj3 = fb * dfsdj3;
  const Real dgdj3 = fbg * dfsdj3;

  // yield functions.  See comment in yieldFunctionValuesV
  all_q[0].f = std::sqrt(std::pow(q, 2.0) + _small_smoother2) + fb * fs;

  // d(yield function)/d(stress_params)
  for (unsigned yf = 0; yf < _num_yf; ++yf)
    for (unsigned a = 0; a < 3; ++a)
      all_q[yf].df[a] = dfdi1 * di1ds[a] + dfdj2 * dj2ds[a] + dfdj3 * dj3ds[a];

  all_q[0].df[3] = 2.0 * (dfdi1 * di1ds[3] + dfdj2 * dj2ds[3] + dfdj3 * dj3ds[3]);
  all_q[0].df[4] = 2.0 * (dfdi1 * di1ds[4] + dfdj2 * dj2ds[4] + dfdj3 * dj3ds[4]);
  all_q[0].df[5] = 2.0 * (dfdi1 * di1ds[5] + dfdj2 * dj2ds[5] + dfdj3 * dj3ds[5]);

  // the flow potential is just the yield function with phi->psi
  // d(flow potential)/d(stress_params)
  for (unsigned yf = 0; yf < _num_yf; ++yf)
    for (unsigned a = 0; a < 3; ++a)
      all_q[yf].dg[a] = dgdi1 * di1ds[a] + dgdj2 * dj2ds[a] + dgdj3 * dj3ds[a];

  all_q[0].dg[3] = 2.0 * (dgdi1 * di1ds[3] + dgdj2 * dj2ds[3] + dgdj3 * dj3ds[3]);
  all_q[0].dg[4] = 2.0 * (dgdi1 * di1ds[4] + dgdj2 * dj2ds[4] + dgdj3 * dj3ds[4]);
  all_q[0].dg[5] = 2.0 * (dgdi1 * di1ds[5] + dgdj2 * dj2ds[5] + dgdj3 * dj3ds[5]);

  // d(yield function)/d(intnl)

  all_q[0].df_di[0] = 0.0;

  // d(flow potential)/d(stress_params)/d(intnl)
  for (unsigned yf = 0; yf < _num_yf; ++yf)
    for (unsigned a = 0; a < _num_sp; ++a)
      for (unsigned i = 0; i < _num_intnl; ++i)
        all_q[yf].d2g_di[a][i] = 0.0;

  all_q[0].d2g_di[3][0] = 2.0 * all_q[0].d2g_di[3][0];
  all_q[0].d2g_di[4][0] = 2.0 * all_q[0].d2g_di[4][0];
  all_q[0].d2g_di[5][0] = 2.0 * all_q[0].d2g_di[5][0];

  std::vector<std::vector<Real>> d2j2ds(_num_sp, std::vector<Real>(_num_sp));
  d2j2ds[0][0] = 2.0 / 3.0;
  d2j2ds[1][1] = 2.0 / 3.0;
  d2j2ds[2][2] = 2.0 / 3.0;
  d2j2ds[5][5] = 0.5;
  d2j2ds[4][4] = 0.5;
  d2j2ds[3][3] = 0.5;
  d2j2ds[0][1] = -1.0 / 3.0;
  d2j2ds[0][2] = -1.0 / 3.0;
  d2j2ds[1][2] = -1.0 / 3.0;
  d2j2ds[1][0] = -1.0 / 3.0;
  d2j2ds[2][0] = -1.0 / 3.0;
  d2j2ds[2][1] = -1.0 / 3.0;
  d2j2ds[0][5] = d2j2ds[0][4] = d2j2ds[0][3] = 0.0;
  d2j2ds[1][5] = d2j2ds[1][4] = d2j2ds[1][3] = 0.0;
  d2j2ds[2][5] = d2j2ds[2][4] = d2j2ds[2][3] = 0.0;
  d2j2ds[5][0] = d2j2ds[4][0] = d2j2ds[3][0] = 0.0;
  d2j2ds[5][1] = d2j2ds[4][1] = d2j2ds[3][1] = 0.0;
  d2j2ds[5][2] = d2j2ds[4][2] = d2j2ds[3][2] = 0.0;
  d2j2ds[3][4] = d2j2ds[3][5] = d2j2ds[4][3] = 0.0;
  d2j2ds[4][5] = d2j2ds[5][3] = d2j2ds[5][4] = 0.0;

  std::vector<std::vector<Real>> d2j3ds(_num_sp, std::vector<Real>(_num_sp));

  d2j3ds[0][0] = d2j3ds_t(0, 0, 0, 0);
  d2j3ds[0][1] = d2j3ds_t(0, 0, 1, 1);
  d2j3ds[0][2] = d2j3ds_t(0, 0, 2, 2);
  d2j3ds[0][3] = d2j3ds_t(0, 0, 0, 1);
  d2j3ds[0][4] = d2j3ds_t(0, 0, 0, 2);
  d2j3ds[0][5] = d2j3ds_t(0, 0, 1, 2);

  d2j3ds[1][0] = d2j3ds_t(1, 1, 0, 0);
  d2j3ds[1][1] = d2j3ds_t(1, 1, 1, 1);
  d2j3ds[1][2] = d2j3ds_t(1, 1, 2, 2);
  d2j3ds[1][3] = d2j3ds_t(1, 1, 0, 1);
  d2j3ds[1][4] = d2j3ds_t(1, 1, 0, 2);
  d2j3ds[1][5] = d2j3ds_t(1, 1, 1, 2);

  d2j3ds[2][0] = d2j3ds_t(2, 2, 0, 0);
  d2j3ds[2][1] = d2j3ds_t(2, 2, 1, 1);
  d2j3ds[2][2] = d2j3ds_t(2, 2, 2, 2);
  d2j3ds[2][3] = d2j3ds_t(2, 2, 0, 1);
  d2j3ds[2][4] = d2j3ds_t(2, 2, 0, 2);
  d2j3ds[2][5] = d2j3ds_t(2, 2, 1, 2);

  d2j3ds[3][0] = d2j3ds[0][3] = d2j3ds_t(0, 1, 0, 0);
  d2j3ds[3][1] = d2j3ds[1][3] = d2j3ds_t(0, 1, 1, 1);
  d2j3ds[3][2] = d2j3ds[2][3] = d2j3ds_t(0, 1, 2, 2);
  d2j3ds[3][3] = d2j3ds_t(0, 1, 0, 1);
  d2j3ds[3][4] = d2j3ds[4][3] = d2j3ds_t(0, 1, 0, 2);
  d2j3ds[3][5] = d2j3ds[5][3] = d2j3ds_t(0, 1, 1, 2);

  d2j3ds[4][0] = d2j3ds[0][4] = d2j3ds_t(0, 2, 0, 0);
  d2j3ds[4][1] = d2j3ds[1][4] = d2j3ds_t(0, 2, 1, 1);
  d2j3ds[4][2] = d2j3ds[2][4] = d2j3ds_t(0, 2, 2, 2);
  d2j3ds[4][3] = d2j3ds[3][4] = d2j3ds_t(0, 2, 0, 1);
  d2j3ds[4][4] = d2j3ds_t(0, 2, 0, 2);
  d2j3ds[4][5] = d2j3ds[5][4] = d2j3ds_t(0, 2, 1, 2);

  d2j3ds[5][0] = d2j3ds[0][5] = d2j3ds_t(1, 2, 0, 0);
  d2j3ds[5][1] = d2j3ds[1][5] = d2j3ds_t(1, 2, 1, 1);
  d2j3ds[5][2] = d2j3ds[2][5] = d2j3ds_t(1, 2, 2, 2);
  d2j3ds[5][3] = d2j3ds[3][5] = d2j3ds_t(1, 2, 0, 1);
  d2j3ds[5][4] = d2j3ds[4][5] = d2j3ds_t(1, 2, 0, 2);
  d2j3ds[5][5] = d2j3ds_t(1, 2, 1, 2);

  // const Real d2fbgdi1i1 = 0.0;
  // const Real d2fsdi1i1 = 0.0;
  // const Real d2fsdi1j2 = 0.0;
  // const Real d2fsdi1j3 = 0.0;
  const Real d2fsdj2j2 = (15.0 * _beta * _mv * cof * j3 * std::pow(fs, (_mv - 1.0) / _mv)) /
                             (4.0 * std::pow(j2, 3.5)) +
                         (9.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(j3, 2.0) * std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                             (4.0 * std::pow(j2, 5.0));

  const Real d2fsdj2j3 =
      -(3.0 * _beta * _mv * cof * std::pow(fs, (_mv - 1.0) / _mv)) / (2.0 * std::pow(j2, 2.5)) -
      (3.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) * j3 *
       std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
          (2.0 * std::pow(j2, 4.0));

  const Real d2fsdj3j3 = (std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 3.0);

  const Real d2gdi1i1 =
      0.0; // d2fbgdi1i1 * fs + fbg * d2fsdi1i1 + dfbgdi1 * dfsdi1 + dfbgdi1 * dfsdi1
  const Real d2gdi1j2 = 0.0;              // dfbgdi1 * dfsdj2 + fbg * d2fsdi1j2;
  const Real d2gdi1j3 = dfbgdi1 * dfsdj3; // + fbg * d2fsdi1j3;
  const Real d2gdj2j2 =
      -(1.0 / 4.0) * std::pow(std::pow(q, 2.0) + _small_smoother2, -1.5) + fbg * d2fsdj2j2;
  const Real d2gdj2j3 = fbg * d2fsdj2j3;
  const Real d2gdj3j3 = fbg * d2fsdj3j3;

  std::vector<Real> d2gi1ds(_num_sp);
  std::vector<Real> d2gj2ds(_num_sp);
  std::vector<Real> d2gj3ds(_num_sp);

  for (unsigned i = 0; i < _num_sp; ++i)
  {
    d2gi1ds[i] = d2gdi1i1 * di1ds[i] + d2gdi1j2 * dj2ds[i] + d2gdi1j3 * dj3ds[i];
    d2gj2ds[i] = d2gdi1j2 * di1ds[i] + d2gdj2j2 * dj2ds[i] + d2gdj2j3 * dj3ds[i];
    d2gj3ds[i] = d2gdi1j3 * di1ds[i] + d2gdj2j3 * dj2ds[i] + d2gdj3j3 * dj3ds[i];
    for (unsigned j = 0; j < _num_sp; ++j)
      all_q[0].d2g[i][j] = d2gi1ds[i] * di1ds[j] + d2gj2ds[i] * dj2ds[j] + d2gj3ds[i] * dj3ds[j] +
                           dgdj2 * d2j2ds[i][j] + dgdj3 * d2j3ds[i][j];
  }

  for (unsigned yf = 0; yf < _num_yf; ++yf)
    for (unsigned a = 0; a < _num_sp; ++a)
      for (unsigned b = 0; b < _num_sp; ++b)
        if (a < 3 && b < 3)
          all_q[yf].d2g[a][b] = all_q[yf].d2g[a][b];
        else if (a >= 3 && b >= 3)
          all_q[yf].d2g[a][b] = 4.0 * all_q[yf].d2g[a][b];
        else
          all_q[yf].d2g[a][b] = 2.0 * all_q[yf].d2g[a][b];

  return;
}

void
OpalinusPerfectPlasticStressUpdate::setEffectiveElasticity(const RankFourTensor & Eijkl)
{

  _En = Eijkl(0, 0, 0, 0) + Eijkl(1, 1, 1, 1) + Eijkl(2, 2, 2, 2);
  _Eij[0][0] = 1.0 * Eijkl(0, 0, 0, 0);
  _Eij[0][1] = 1.0 * Eijkl(0, 0, 1, 1);
  _Eij[0][2] = 1.0 * Eijkl(0, 0, 2, 2);
  _Eij[0][5] = 1.0 * Eijkl(0, 0, 1, 2);
  _Eij[0][4] = 1.0 * Eijkl(0, 0, 0, 2);
  _Eij[0][3] = 1.0 * Eijkl(0, 0, 0, 1);

  _Eij[1][0] = 1.0 * Eijkl(1, 1, 0, 0);
  _Eij[1][1] = 1.0 * Eijkl(1, 1, 1, 1);
  _Eij[1][2] = 1.0 * Eijkl(1, 1, 2, 2);
  _Eij[1][5] = 1.0 * Eijkl(1, 1, 1, 2);
  _Eij[1][4] = 1.0 * Eijkl(1, 1, 0, 2);
  _Eij[1][3] = 1.0 * Eijkl(1, 1, 0, 1);

  _Eij[2][0] = 1.0 * Eijkl(2, 2, 0, 0);
  _Eij[2][1] = 1.0 * Eijkl(2, 2, 1, 1);
  _Eij[2][2] = 1.0 * Eijkl(2, 2, 2, 2);
  _Eij[2][5] = 1.0 * Eijkl(2, 2, 1, 2);
  _Eij[2][4] = 1.0 * Eijkl(2, 2, 0, 2);
  _Eij[2][3] = 1.0 * Eijkl(2, 2, 0, 1);

  _Eij[5][0] = 1.0 * Eijkl(1, 2, 0, 0);
  _Eij[5][1] = 1.0 * Eijkl(1, 2, 1, 1);
  _Eij[5][2] = 1.0 * Eijkl(1, 2, 2, 2);
  _Eij[5][5] = 1.0 * Eijkl(1, 2, 1, 2);
  _Eij[5][4] = 1.0 * Eijkl(1, 2, 0, 2);
  _Eij[5][3] = 1.0 * Eijkl(1, 2, 0, 1);

  _Eij[4][0] = 1.0 * Eijkl(0, 2, 0, 0);
  _Eij[4][1] = 1.0 * Eijkl(0, 2, 1, 1);
  _Eij[4][2] = 1.0 * Eijkl(0, 2, 2, 2);
  _Eij[4][5] = 1.0 * Eijkl(0, 2, 1, 2);
  _Eij[4][4] = 1.0 * Eijkl(0, 2, 0, 2);
  _Eij[4][3] = 1.0 * Eijkl(0, 2, 0, 1);

  _Eij[3][0] = 1.0 * Eijkl(0, 1, 0, 0);
  _Eij[3][1] = 1.0 * Eijkl(0, 1, 1, 1);
  _Eij[3][2] = 1.0 * Eijkl(0, 1, 2, 2);
  _Eij[3][5] = 1.0 * Eijkl(0, 1, 1, 2);
  _Eij[3][4] = 1.0 * Eijkl(0, 1, 0, 2);
  _Eij[3][3] = 1.0 * Eijkl(0, 1, 0, 1);

  return;
}

void
OpalinusPerfectPlasticStressUpdate::initializeVarsV(const std::vector<Real> & trial_stress_params,
                                                    const std::vector<Real> & intnl_old,
                                                    std::vector<Real> & stress_params,
                                                    Real & gaE,
                                                    std::vector<Real> & intnl) const
{

  if (!_perfect_guess)
  {
    for (unsigned i = 0; i < _num_sp; ++i)
      stress_params[i] = trial_stress_params[i];

    gaE = 0.0;
  }
  else
  {

    const RankTwoTensor stress_trial = RankTwoTensor(trial_stress_params[0],
                                                     trial_stress_params[3],
                                                     trial_stress_params[4],
                                                     trial_stress_params[3],
                                                     trial_stress_params[1],
                                                     trial_stress_params[5],
                                                     trial_stress_params[4],
                                                     trial_stress_params[5],
                                                     trial_stress_params[2]);

    const Real i1_trial = std::min(-1e-5, stress_trial.trace() - 3.0 * _St);

    Real j2_trial = stress_trial.secondInvariant();
    Real j3_trial = stress_trial.thirdInvariant();
    if (j2_trial < 1e-2)
    {
      j2_trial = 1e-2;
      j3_trial = 0.0;
    }

    const Real cof = std::sqrt(27.0) / 2.0;
    const Real q_trial = std::sqrt(j2_trial);

    const Real sr = std::clamp(cof * j3_trial / std::pow(j2_trial, 1.5), -1.0, +1.0);

    const Real fs_trial = std::pow(std::exp(_betta1 * i1_trial) + _beta * sr, _mv);

    const Real fb_trial = _gama[_qp] * i1_trial;

    const Real trial_desai_yf =
        std::sqrt(std::pow(q_trial, 2.0) + _small_smoother2) + fb_trial * fs_trial;

    const Real df_di = 0.0;

    bool found_solution = false;

    if (trial_desai_yf <= _f_tol)
    {
      for (unsigned i = 0; i < _num_sp; ++i)
        stress_params[i] = trial_stress_params[i];

      gaE = 0.0;
      found_solution = true;
    }

    if (!found_solution && trial_desai_yf > _f_tol)
    {

      std::vector<Real> dg(_num_sp);
      std::vector<Real> df(_num_sp);
      compute_dg(trial_stress_params, intnl_old, dg);
      compute_df(trial_stress_params, intnl_old, df);
      std::vector<Real> cdg(_num_sp);
      Real fcg = 0.0;
      for (unsigned i = 0; i < _num_sp; ++i)
      {
        cdg[i] = 0.0;

        for (unsigned j = 0; j < _num_sp; ++j)
          cdg[i] += _Eij[i][j] * dg[j];

        fcg += cdg[i] * df[i];
      }

      const Real ga = trial_desai_yf / (fcg - df_di);

      for (unsigned i = 0; i < _num_sp; ++i)
        stress_params[i] = trial_stress_params[i] - ga * cdg[i];

      gaE = ga * _En;
      const RankTwoTensor stress_test = RankTwoTensor(stress_params[0],
                                                      stress_params[3],
                                                      stress_params[4],
                                                      stress_params[3],
                                                      stress_params[1],
                                                      stress_params[5],
                                                      stress_params[4],
                                                      stress_params[5],
                                                      stress_params[2]);

      const Real i1_test = stress_test.trace() - 3.0 * _St;

      // Real j2_test = (stress_test.secondInvariant());
      // Real j3_test = (stress_test.thirdInvariant());
      // if (j2_test < 1e-2)
      // {
      //   j3_test = 0.0;
      //   j2_test = 1e-2;
      // }

      // const Real cof = std::sqrt(27.0) / 2.0;
      // const Real q_test = std::sqrt(j2_test);

      // const Real sr = std::clamp(cof * j3_test / std::pow(j2_test, 1.5), -1.0, +1.0);

      // const Real fs_test = std::pow(std::exp(_betta1 * i1_test) + _beta * sr, _mv);

      // const Real fb_test = _gama[_qp] * i1_test;

      // const Real test_desai_yf =
      //     std::sqrt(std::pow(q_test, 2.0) + _small_smoother2) + fb_test * fs_test;

      if (i1_test >= 3.0 * _St - std::sqrt(_small_smoother2) / (_mean_gama))
      {
        stress_params[0] = _St - std::sqrt(_small_smoother2) / (3.0 * _mean_gama);
        stress_params[1] = _St - std::sqrt(_small_smoother2) / (3.0 * _mean_gama);
        stress_params[2] = _St - std::sqrt(_small_smoother2) / (3.0 * _mean_gama);
        stress_params[3] = 0.0;
        stress_params[4] = 0.0;
        stress_params[5] = 0.0;
        // const Real ptrial =
        //     trial_stress_params[0] + trial_stress_params[1] + trial_stress_params[2];
        // const Real p = stress_params[0] + stress_params[1] + stress_params[2];

        gaE = _En * (trial_stress_params[0] - stress_params[0]) / (_Eij[0][0] * _mean_gama);
      }

      found_solution = true;
    }
  }

  setIntnlValuesV(trial_stress_params, stress_params, intnl_old, intnl);

  return;
}

void
OpalinusPerfectPlasticStressUpdate::setIntnlValuesV(const std::vector<Real> & trial_stress_params,
                                                    const std::vector<Real> & current_stress_params,
                                                    const std::vector<Real> & intnl_old,
                                                    std::vector<Real> & intnl) const
{

  std::vector<Real> cp3(_num_sp);
  std::vector<Real> dg(_num_sp);
  compute_dg(current_stress_params, intnl_old, dg);

  for (unsigned i = 0; i < _num_sp; ++i)
  {
    cp3[i] = 0.0;
    for (unsigned j = 0; j < _num_sp; ++j)
      cp3[i] += _Eij[i][j] * dg[j];
  }

  const Real ptrial = trial_stress_params[0] + trial_stress_params[1] + trial_stress_params[2];
  const Real p = current_stress_params[0] + current_stress_params[1] + current_stress_params[2];

  intnl[0] = intnl_old[0] + (ptrial - p) / (cp3[0] + cp3[1] + cp3[2]);

  return;
}

void
OpalinusPerfectPlasticStressUpdate::setIntnlDerivativesV(
    const std::vector<Real> & trial_stress_params,
    const std::vector<Real> & current_stress_params,
    const std::vector<Real> & intnl,
    std::vector<std::vector<Real>> & dintnl) const
{
  (void)trial_stress_params;

  // @Kavan-Khaledi: we do not use dg below. can we comment the next lines out?
  std::vector<Real> dg(_num_sp);
  compute_dg(current_stress_params, intnl, dg);

  // @Kavan-Khaledi: we do not use cp3 below. can we comment the next lines out?
  std::vector<Real> cp3(_num_sp);
  for (unsigned i = 0; i < _num_sp; ++i)
  {
    cp3[i] = 0.0;
    for (unsigned j = 0; j < _num_sp; ++j)
      cp3[i] += _Eij[i][j] * dg[j];
  }

  // @Kavan-Khaledi: we do not use d2g and cd2g below. can we comment the next lines out?
  std::vector<std::vector<Real>> d2g(_num_sp, std::vector<Real>(_num_sp));
  compute_d2g(current_stress_params, intnl, d2g);
  std::vector<Real> cd2g(_num_sp);
  for (unsigned i = 0; i < _num_sp; ++i)
  {
    cd2g[i] = 0.0;
    for (unsigned j = 0; j < _num_sp; ++j)
      cd2g[i] += _Eij[i][j] * d2g[j][i];
  }

  for (std::size_t i = 0; i < _num_sp; ++i)
    dintnl[0][i] = 0.0;
  //   ((-1.0 / cp3[i] + (trial_stress_params[i] - current_stress_params[i]) * cd2g[i]) /
  //    std::pow(cp3[i], 2.0));

  return;
}

void
OpalinusPerfectPlasticStressUpdate::consistentTangentOperatorV(
    const RankTwoTensor & stress_trial,
    const std::vector<Real> & trial_stress_params,
    const RankTwoTensor & /*stress*/,
    const std::vector<Real> & stress_params,
    Real gaE,
    const yieldAndFlow & /*smoothed_q*/,
    const RankFourTensor & elasticity_tensor,
    bool compute_full_tangent_operator,
    const std::vector<std::vector<Real>> & dvar_dtrial,
    RankFourTensor & cto)
{
  (void)stress_trial;
  (void)trial_stress_params;
  (void)dvar_dtrial;

  cto = elasticity_tensor;
  if (!compute_full_tangent_operator)
    return;

  Real ga = gaE / _En;

  RankTwoTensor dqdsig = RankTwoTensor();
  RankTwoTensor dFdsig = RankTwoTensor();
  RankFourTensor d2gds = RankFourTensor();
  RankTwoTensor gmatrix = RankTwoTensor();
  RankTwoTensor temp = RankTwoTensor();
  RankTwoTensor cdgdsig = RankTwoTensor();
  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);
  setGamaValue(stress_params, _gama[_qp]);
  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = stress_now.trace() - 3.0 * _St;
  const Real j2 = std::max(1e-2, stress_now.secondInvariant());
  const Real j3 = stress_now.thirdInvariant();

  const Real q = std::sqrt(j2);
  const Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  const Real fs = std::pow(std::exp(_betta1 * i1) + _beta * sr, _mv);

  ////
  const Real _gamaq = _psi_to_phi * _gama[_qp];
  const Real fb = _gama[_qp] * i1;
  const Real fbg = _gamaq * i1;

  const RankTwoTensor di1ds_t = stress_now.dtrace();
  RankTwoTensor dj2ds_t = stress_now.deviatoric();
  RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();
  RankFourTensor d2j3ds_t = stress_now.d2thirdInvariant();
  RankFourTensor d2j2ds_t = stress_now.d2secondInvariant();

  if (j2 < 1e-2) // @Kavan-Khaledi: this condition is never met (we make sure j2 is at least 1e-2)
  {
    dj2ds_t = RankTwoTensor();
    dj3ds_t = RankTwoTensor();
    d2j2ds_t = RankFourTensor();
    d2j3ds_t = RankFourTensor();
  }

  const Real dfbdi1 = _gama[_qp];

  const Real dfbgdi1 = _gamaq;
  // const Real dfsdi1 = 0.0;
  const Real dfdi1 = dfbdi1 * fs;  // + fb * dfsdi1;
  const Real dgdi1 = dfbgdi1 * fs; // + fbg * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dfdj2 = (1.0 / 2.0) / (std::sqrt(pow(q, 2.0) + _small_smoother2)) + fb * dfsdj2;
  const Real dgdj2 = (1.0 / 2.0) / (std::sqrt(pow(q, 2.0) + _small_smoother2)) + fbg * dfsdj2;
  const Real dfdj3 = fb * dfsdj3;
  const Real dgdj3 = fbg * dfsdj3;

  // const Real dfb_dalfa = 0.0;
  // const Real dalfa_di = 0.0;
  const Real df_di = 0.0; // (dfb_dalfa * fs) * dalfa_di;

  dFdsig = dfdi1 * di1ds_t + dfdj2 * dj2ds_t + dfdj3 * dj3ds_t;
  dqdsig = dgdi1 * di1ds_t + dgdj2 * dj2ds_t + dgdj3 * dj3ds_t;

  // const Real d2fbgdi1i1 = 0.0;

  // const Real d2fsdi1i1 = 0.0;
  // const Real d2fsdi1j2 = 0.0;
  // const Real d2fsdi1j3 = 0.0;
  const Real d2fsdj2j2 = (15.0 * _beta * _mv * cof * j3 * std::pow(fs, (_mv - 1.0) / _mv)) /
                             (4.0 * std::pow(j2, 3.5)) +
                         (9.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(j3, 2.0) * std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                             (4.0 * std::pow(j2, 5.0));

  const Real d2fsdj2j3 =
      -(3.0 * _beta * _mv * cof * std::pow(fs, (_mv - 1.0) / _mv)) / (2.0 * std::pow(j2, 2.5)) -
      (3.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) * j3 *
       std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
          (2.0 * std::pow(j2, 4.0));
  const Real d2fsdj3j3 = (std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 3.0);

  // const Real d2gdi1i1 = d2fbgdi1i1 * fs + dfbgdi1 * dfsdi1 + dfbgdi1 * dfsdi1 + fbg * d2fsdi1i1;
  const Real d2gdi1j2 = dfbgdi1 * dfsdj2; // + fbg * d2fsdi1j2;
  const Real d2gdi1j3 = dfbgdi1 * dfsdj3; // + fbg * d2fsdi1j3;
  const Real d2gdj2j2 =
      -(1.0 / 4.0) * std::pow(std::pow(q, 2.0) + _small_smoother2, -1.5) + fbg * d2fsdj2j2;
  const Real d2gdj2j3 = fbg * d2fsdj2j3;
  const Real d2gdj3j3 = fbg * d2fsdj3j3;

  d2gds = d2gdi1j2 * di1ds_t.outerProduct(dj2ds_t) + // d2gdi1i1 * di1ds_t.outerProduct(di1ds_t)
          d2gdi1j3 * di1ds_t.outerProduct(dj3ds_t) + d2gdi1j2 * dj2ds_t.outerProduct(di1ds_t) +
          d2gdj2j2 * dj2ds_t.outerProduct(dj2ds_t) + d2gdj2j3 * dj2ds_t.outerProduct(dj3ds_t) +
          d2gdi1j3 * dj3ds_t.outerProduct(di1ds_t) + d2gdj2j3 * dj3ds_t.outerProduct(dj2ds_t) +
          d2gdj3j3 * dj3ds_t.outerProduct(dj3ds_t) + dgdj2 * d2j2ds_t + dgdj3 * d2j3ds_t;

  RankFourTensor inv =
      RankFourTensor(RankFourTensor::initIdentityFour) + (ga)*elasticity_tensor * d2gds;

  inv = inv.transposeMajor().invSymm();
  inv = (elasticity_tensor.transposeMajor() * inv); //.transposeMajor();

  for (unsigned i = 0; i < _tensor_dimensionality; ++i)
    for (unsigned j = 0; j < _tensor_dimensionality; ++j)
      for (unsigned k = 0; k < _tensor_dimensionality; ++k)
        for (unsigned l = 0; l < _tensor_dimensionality; ++l)
        {
          cdgdsig(i, j) += inv(i, j, k, l) * dFdsig(k, l);
          temp(i, j) += inv(k, l, i, j) * dqdsig(k, l);
        }

  gmatrix = cdgdsig / (dFdsig.doubleContraction(temp) - df_di);

  cto = inv - (temp.outerProduct(gmatrix));

  return;
}

void
OpalinusPerfectPlasticStressUpdate::initializeReturnProcess()
{
}

void
OpalinusPerfectPlasticStressUpdate::compute_dg(const std::vector<Real> & stress_params,
                                               const std::vector<Real> & intnl,
                                               std::vector<Real> & dgg) const
{
  (void)intnl; // silence compiler warning about unused variable

  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = stress_now.trace() - 3.0 * _St;
  const Real j2 = std::max(1e-2, stress_now.secondInvariant());
  const Real j3 = stress_now.thirdInvariant();

  const Real q = std::sqrt(j2);
  const Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  const Real fs = std::pow(std::exp(_betta1 * i1) + _beta * sr, _mv);

  const Real _gamaq = _psi_to_phi * _gama[_qp];
  const Real fbg = _gamaq * i1;

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = 1.0;
  di1ds[1] = 1.0;
  di1ds[2] = 1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 <= 1e-8) // @Kavan-Khaledi: this condition is never met (we make sure j2 is at least 1e-2)
  {
    dj2ds[0] = 0.0;
    dj2ds[1] = 0.0;
    dj2ds[2] = 0.0;
    dj2ds[5] = 0.0;
    dj2ds[4] = 0.0;
    dj2ds[3] = 0.0;

    dj3ds[0] = 0.0;
    dj3ds[1] = 0.0;
    dj3ds[2] = 0.0;
    dj3ds[5] = 0.0;
    dj3ds[4] = 0.0;
    dj3ds[3] = 0.0;
  }
  else
  {
    dj2ds[0] = dj2ds_t(0, 0);
    dj2ds[1] = dj2ds_t(1, 1);
    dj2ds[2] = dj2ds_t(2, 2);
    dj2ds[5] = dj2ds_t(1, 2);
    dj2ds[4] = dj2ds_t(0, 2);
    dj2ds[3] = dj2ds_t(0, 1);

    dj3ds[0] = dj3ds_t(0, 0);
    dj3ds[1] = dj3ds_t(1, 1);
    dj3ds[2] = dj3ds_t(2, 2);
    dj3ds[5] = dj3ds_t(1, 2);
    dj3ds[4] = dj3ds_t(0, 2);
    dj3ds[3] = dj3ds_t(0, 1);
  }

  const Real dfbgdi1 = _gamaq;
  // const Real dfsdi1 = 0.0;
  const Real dgdi1 = dfbgdi1 * fs; // + fbg * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);
  const Real dgdj2 = (1.0 / 2.0) / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) + fbg * dfsdj2;
  const Real dgdj3 = fbg * dfsdj3;

  for (unsigned a = 0; a < 3; ++a)
    dgg[a] = dgdi1 * di1ds[a] + dgdj2 * dj2ds[a] + dgdj3 * dj3ds[a];

  dgg[3] = 2.0 * (dgdi1 * di1ds[3] + dgdj2 * dj2ds[3] + dgdj3 * dj3ds[3]);
  dgg[4] = 2.0 * (dgdi1 * di1ds[4] + dgdj2 * dj2ds[4] + dgdj3 * dj3ds[4]);
  dgg[5] = 2.0 * (dgdi1 * di1ds[5] + dgdj2 * dj2ds[5] + dgdj3 * dj3ds[5]);

  return;
}

void
OpalinusPerfectPlasticStressUpdate::compute_df(const std::vector<Real> & stress_params,
                                               const std::vector<Real> & intnl,
                                               std::vector<Real> & dff) const
{
  (void)intnl; // silence compiler warning about unused variable

  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = stress_now.trace() - 3.0 * _St;
  const Real j2 = std::max(1e-2, stress_now.secondInvariant());
  const Real j3 = stress_now.thirdInvariant();

  const Real q = std::sqrt(j2);
  const Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  const Real fs = std::pow(std::exp(_betta1 * i1) + _beta * sr, _mv);

  const Real fb = _gama[_qp] * i1;

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = 1.0;
  di1ds[1] = 1.0;
  di1ds[2] = 1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 <= 1e-8) // @Kavan-Khaledi: this condition is never met (we make sure j2 is at least 1e-2)
  {
    dj2ds[0] = 0.0;
    dj2ds[1] = 0.0;
    dj2ds[2] = 0.0;
    dj2ds[5] = 0.0;
    dj2ds[4] = 0.0;
    dj2ds[3] = 0.0;

    dj3ds[0] = 0.0;
    dj3ds[1] = 0.0;
    dj3ds[2] = 0.0;
    dj3ds[5] = 0.0;
    dj3ds[4] = 0.0;
    dj3ds[3] = 0.0;
  }
  else
  {
    dj2ds[0] = dj2ds_t(0, 0);
    dj2ds[1] = dj2ds_t(1, 1);
    dj2ds[2] = dj2ds_t(2, 2);
    dj2ds[5] = dj2ds_t(1, 2);
    dj2ds[4] = dj2ds_t(0, 2);
    dj2ds[3] = dj2ds_t(0, 1);

    dj3ds[0] = dj3ds_t(0, 0);
    dj3ds[1] = dj3ds_t(1, 1);
    dj3ds[2] = dj3ds_t(2, 2);
    dj3ds[5] = dj3ds_t(1, 2);
    dj3ds[4] = dj3ds_t(0, 2);
    dj3ds[3] = dj3ds_t(0, 1);
  }

  const Real dfbdi1 = _gama[_qp];
  // const Real dfsdi1 = 0.0;
  const Real dfdi1 = dfbdi1 * fs; // + fb * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dfdj2 = (1.0 / 2.0) / (std::sqrt(pow(q, 2.0) + _small_smoother2)) + fb * dfsdj2;
  const Real dfdj3 = fb * dfsdj3;

  for (unsigned a = 0; a < 3; ++a)
    dff[a] = dfdi1 * di1ds[a] + dfdj2 * dj2ds[a] + dfdj3 * dj3ds[a];

  dff[3] = 2.0 * (dfdi1 * di1ds[3] + dfdj2 * dj2ds[3] + dfdj3 * dj3ds[3]);
  dff[4] = 2.0 * (dfdi1 * di1ds[4] + dfdj2 * dj2ds[4] + dfdj3 * dj3ds[4]);
  dff[5] = 2.0 * (dfdi1 * di1ds[5] + dfdj2 * dj2ds[5] + dfdj3 * dj3ds[5]);

  return;
}

void
OpalinusPerfectPlasticStressUpdate::setGamaValue(const std::vector<Real> & stress_params,
                                                 Real & gama) const
{
  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);
  RankTwoTensor rotation = RankTwoTensor();
  std::vector<Real> lamda(3);
  RankTwoTensor eigv = RankTwoTensor();

  const auto e3 = _localCoordinateSystem.e3();

  stress_now.symmetricEigenvaluesEigenvectors(lamda, eigv);

  const Real cosb = eigv(0, 0) * e3(0) + eigv(1, 0) * e3(1) + eigv(2, 0) * e3(2);

  const Real lv2 = std::pow(cosb, 2.0);

  gama = (_mean_gama * (1.0 + _omega1 * (1.0 - 3.0 * lv2) +
                        _b1 * std::pow(_omega1, 2.0) * std::pow((1.0 - 3.0 * lv2), 2.0)));

  return;
}

void
OpalinusPerfectPlasticStressUpdate::compute_d2g(const std::vector<Real> & stress_params,
                                                const std::vector<Real> & intnl,
                                                std::vector<std::vector<Real>> & d2gg) const
{
  (void)intnl; // silence compiler warning about unused variable

  const RankTwoTensor stress_now = RankTwoTensor(stress_params[0],
                                                 stress_params[3],
                                                 stress_params[4],
                                                 stress_params[3],
                                                 stress_params[1],
                                                 stress_params[5],
                                                 stress_params[4],
                                                 stress_params[5],
                                                 stress_params[2]);

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = stress_now.trace() - 3.0 * _St;
  const Real j2 = std::max(1e-2, stress_now.secondInvariant());
  const Real j3 = stress_now.thirdInvariant();

  const Real q = std::sqrt(j2);

  const Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  const Real fs = std::pow(std::exp(_betta1 * i1) + _beta * sr, _mv);

  Real _gamaq = _psi_to_phi * _gama[_qp];

  const Real fbg = _gamaq * i1;

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();
  const RankFourTensor d2j3ds_t = stress_now.d2thirdInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = 1.0;
  di1ds[1] = 1.0;
  di1ds[2] = 1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 < 1e-2) // @Kavan-Khaledi: this condition is never met (we make sure j2 is at least 1e-2)
  {
    dj2ds[0] = 0.0;
    dj2ds[1] = 0.0;
    dj2ds[2] = 0.0;
    dj2ds[5] = 0.0;
    dj2ds[4] = 0.0;
    dj2ds[3] = 0.0;

    dj3ds[0] = 0.0;
    dj3ds[1] = 0.0;
    dj3ds[2] = 0.0;
    dj3ds[5] = 0.0;
    dj3ds[4] = 0.0;
    dj3ds[3] = 0.0;
  }
  else
  {
    dj2ds[0] = dj2ds_t(0, 0);
    dj2ds[1] = dj2ds_t(1, 1);
    dj2ds[2] = dj2ds_t(2, 2);
    dj2ds[5] = dj2ds_t(1, 2);
    dj2ds[4] = dj2ds_t(0, 2);
    dj2ds[3] = dj2ds_t(0, 1);

    dj3ds[0] = dj3ds_t(0, 0);
    dj3ds[1] = dj3ds_t(1, 1);
    dj3ds[2] = dj3ds_t(2, 2);
    dj3ds[5] = dj3ds_t(1, 2);
    dj3ds[4] = dj3ds_t(0, 2);
    dj3ds[3] = dj3ds_t(0, 1);
  }

  const Real dfbgdi1 = _gamaq;
  // const Real dfsdi1 = 0.0;

  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dgdj2 = (1.0 / 2.0) / (std::sqrt(pow(q, 2.0) + _small_smoother2)) + fbg * dfsdj2;

  const Real dgdj3 = fbg * dfsdj3;

  std::vector<std::vector<Real>> d2j2ds(_num_sp, std::vector<Real>(_num_sp));
  d2j2ds[0][0] = 2.0 / 3.0;
  d2j2ds[1][1] = 2.0 / 3.0;
  d2j2ds[2][2] = 2.0 / 3.0;
  d2j2ds[5][5] = 0.5;
  d2j2ds[4][4] = 0.5;
  d2j2ds[3][3] = 0.5;
  d2j2ds[0][1] = -1.0 / 3.0;
  d2j2ds[0][2] = -1.0 / 3.0;
  d2j2ds[1][2] = -1.0 / 3.0;
  d2j2ds[1][0] = -1.0 / 3.0;
  d2j2ds[2][0] = -1.0 / 3.0;
  d2j2ds[2][1] = -1.0 / 3.0;
  d2j2ds[0][5] = d2j2ds[0][4] = d2j2ds[0][3] = 0.0;
  d2j2ds[1][5] = d2j2ds[1][4] = d2j2ds[1][3] = 0.0;
  d2j2ds[2][5] = d2j2ds[2][4] = d2j2ds[2][3] = 0.0;
  d2j2ds[5][0] = d2j2ds[4][0] = d2j2ds[3][0] = 0.0;
  d2j2ds[5][1] = d2j2ds[4][1] = d2j2ds[3][1] = 0.0;
  d2j2ds[5][2] = d2j2ds[4][2] = d2j2ds[3][2] = 0.0;
  d2j2ds[3][4] = d2j2ds[3][5] = d2j2ds[4][3] = 0.0;
  d2j2ds[4][5] = d2j2ds[5][3] = d2j2ds[5][4] = 0.0;

  std::vector<std::vector<Real>> d2j3ds(_num_sp, std::vector<Real>(_num_sp));

  d2j3ds[0][0] = d2j3ds_t(0, 0, 0, 0);
  d2j3ds[0][1] = d2j3ds_t(0, 0, 1, 1);
  d2j3ds[0][2] = d2j3ds_t(0, 0, 2, 2);
  d2j3ds[0][3] = d2j3ds_t(0, 0, 0, 1);
  d2j3ds[0][4] = d2j3ds_t(0, 0, 0, 2);
  d2j3ds[0][5] = d2j3ds_t(0, 0, 1, 2);

  d2j3ds[1][0] = d2j3ds_t(1, 1, 0, 0);
  d2j3ds[1][1] = d2j3ds_t(1, 1, 1, 1);
  d2j3ds[1][2] = d2j3ds_t(1, 1, 2, 2);
  d2j3ds[1][3] = d2j3ds_t(1, 1, 0, 1);
  d2j3ds[1][4] = d2j3ds_t(1, 1, 0, 2);
  d2j3ds[1][5] = d2j3ds_t(1, 1, 1, 2);

  d2j3ds[2][0] = d2j3ds_t(2, 2, 0, 0);
  d2j3ds[2][1] = d2j3ds_t(2, 2, 1, 1);
  d2j3ds[2][2] = d2j3ds_t(2, 2, 2, 2);
  d2j3ds[2][3] = d2j3ds_t(2, 2, 0, 1);
  d2j3ds[2][4] = d2j3ds_t(2, 2, 0, 2);
  d2j3ds[2][5] = d2j3ds_t(2, 2, 1, 2);

  d2j3ds[3][0] = d2j3ds[0][3] = d2j3ds_t(0, 1, 0, 0);
  d2j3ds[3][1] = d2j3ds[1][3] = d2j3ds_t(0, 1, 1, 1);
  d2j3ds[3][2] = d2j3ds[2][3] = d2j3ds_t(0, 1, 2, 2);
  d2j3ds[3][3] = d2j3ds_t(0, 1, 0, 1);
  d2j3ds[3][4] = d2j3ds[4][3] = d2j3ds_t(0, 1, 0, 2);
  d2j3ds[3][5] = d2j3ds[5][3] = d2j3ds_t(0, 1, 1, 2);

  d2j3ds[4][0] = d2j3ds[0][4] = d2j3ds_t(0, 2, 0, 0);
  d2j3ds[4][1] = d2j3ds[1][4] = d2j3ds_t(0, 2, 1, 1);
  d2j3ds[4][2] = d2j3ds[2][4] = d2j3ds_t(0, 2, 2, 2);
  d2j3ds[4][3] = d2j3ds[3][4] = d2j3ds_t(0, 2, 0, 1);
  d2j3ds[4][4] = d2j3ds_t(0, 2, 0, 2);
  d2j3ds[4][5] = d2j3ds[5][4] = d2j3ds_t(0, 2, 1, 2);

  d2j3ds[5][0] = d2j3ds[0][5] = d2j3ds_t(1, 2, 0, 0);
  d2j3ds[5][1] = d2j3ds[1][5] = d2j3ds_t(1, 2, 1, 1);
  d2j3ds[5][2] = d2j3ds[2][5] = d2j3ds_t(1, 2, 2, 2);
  d2j3ds[5][3] = d2j3ds[3][5] = d2j3ds_t(1, 2, 0, 1);
  d2j3ds[5][4] = d2j3ds[4][5] = d2j3ds_t(1, 2, 0, 2);
  d2j3ds[5][5] = d2j3ds_t(1, 2, 1, 2);

  // const Real d2fbgdi1i1 = 0.0;

  // const Real d2fsdi1i1 = 0.0;
  // const Real d2fsdi1j2 = 0.0;
  // const Real d2fsdi1j3 = 0.0;
  const Real d2fsdj2j2 = (15.0 * _beta * _mv * cof * j3 * std::pow(fs, (_mv - 1.0) / _mv)) /
                             (4.0 * std::pow(j2, 3.5)) +
                         (9.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(j3, 2.0) * std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                             (4.0 * std::pow(j2, 5.0));

  const Real d2fsdj2j3 =
      -(3.0 * _beta * _mv * cof * std::pow(fs, (_mv - 1.0) / _mv)) / (2.0 * std::pow(j2, 2.5)) -
      (3.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) * j3 *
       std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
          (2.0 * std::pow(j2, 4.0));
  const Real d2fsdj3j3 = (std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 3.0);

  const Real d2gdi1i1 =
      0.0; // d2fbgdi1i1 * fs + dfbgdi1 * dfsdi1 + dfbgdi1 * dfsdi1 + fbg * d2fsdi1i1;
  const Real d2gdi1j2 = dfbgdi1 * dfsdj2; // + fbg * d2fsdi1j2;
  const Real d2gdi1j3 = dfbgdi1 * dfsdj3; // + fbg * d2fsdi1j3;
  const Real d2gdj2j2 =
      -(1.0 / 4.0) * std::pow(std::pow(q, 2.0) + _small_smoother2, -1.5) + fbg * d2fsdj2j2;
  const Real d2gdj2j3 = fbg * d2fsdj2j3;
  const Real d2gdj3j3 = fbg * d2fsdj3j3;

  std::vector<Real> d2gi1ds(_num_sp);
  std::vector<Real> d2gj2ds(_num_sp);
  std::vector<Real> d2gj3ds(_num_sp);

  for (unsigned i = 0; i < _num_sp; ++i)
  {
    d2gi1ds[i] = d2gdi1i1 * di1ds[i] + d2gdi1j2 * dj2ds[i] + d2gdi1j3 * dj3ds[i];
    d2gj2ds[i] = d2gdi1j2 * di1ds[i] + d2gdj2j2 * dj2ds[i] + d2gdj2j3 * dj3ds[i];
    d2gj3ds[i] = d2gdi1j3 * di1ds[i] + d2gdj2j3 * dj2ds[i] + d2gdj3j3 * dj3ds[i];
    for (unsigned j = 0; j < _num_sp; ++j)
    {
      d2gg[i][j] = d2gi1ds[i] * di1ds[j] + d2gj2ds[i] * dj2ds[j] + d2gj3ds[i] * dj3ds[j] +
                   dgdj2 * d2j2ds[i][j] + dgdj3 * d2j3ds[i][j];
    }
  }

  for (unsigned yf = 0; yf < _num_yf; ++yf)
    for (unsigned a = 0; a < _num_sp; ++a)
      for (unsigned b = 0; b < _num_sp; ++b)
        if (a < 3 && b < 3)
          d2gg[a][b] = d2gg[a][b];
        else if (a >= 3 && b >= 3)
          d2gg[a][b] = 4.0 * d2gg[a][b];
        else
          d2gg[a][b] = 2.0 * d2gg[a][b];
  return;
}
