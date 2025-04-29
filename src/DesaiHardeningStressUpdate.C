// PerfectlyPlastic Desai Model beta1=0.0 beta=0.0 alfa=0.0 The model is basically equivalent to
// Drucker-Prager model

#include "DesaiHardeningStressUpdate.h"
#include "libmesh/utility.h"
#include "ElasticityTensorTools.h"
#include "CartesianLocalCoordinateSystem.h"

registerMooseObject(MOOSEAPPNAME, DesaiHardeningStressUpdate);

InputParameters
DesaiHardeningStressUpdate::validParams()
{
  InputParameters params = MultiParameterPlasticityStressUpdate::validParams();

  params.addClassDescription("Compute stress for Non-associated Plasticity for Geomaterials");

  // === material coordinate system ===
  params.addRequiredParam<UserObjectName>(
      "local_coordinate_system",
      "The UserObject that defines the local coordinate system. "
      "The local axis e1 and e2 of this coordinate system are considered in plane while "
      "the local axis e3 is assumed to be 'normal'.");

  params.addRequiredParam<Real>("p_tensile",
                                "Intersection of the failure with the p-axis in the p-q-diagram "
                                "(physical unit: pressure).");

  // === Anisotropic Strength ===
  // parameters from S. Pietruszczak, Z. Mroz: "Formulation of anisotropic failure criteria
  // incorporating a microstructure tensor", Computers and Geotechnics 26 (2) (2000) 105-112.
  params.addRequiredParam<Real>(
      "gamma_mean",
      "Slope of the failure line in p-q-diagram. In the paper of Khaledi "
      "et al. (2023) this parameter is named gamma_m (physical unit: dimensionless).");
  params.addParam<Real>(
      "parameter_omega_1",
      0.0,
      "Anisotropic plasticity parameter omega_1 (physical unit: dimensionless) according to S. "
      "Pietruszczak and Z. Mróz (2000). "
      "In the paper of Khaledi et al. (2023) this parameter is named c_{gamma1}.");
  params.addParam<Real>("parameter_b_1",
                        0.0,
                        "Anisotropic plasticity parameter b_1 (physical unit: dimensionless) "
                        "according to S. Pietruszczak and Z. "
                        "Mróz (2000). In the "
                        "paper of Khaledi et al. (2023) this parameter is named c_{gamma2}.");

  // === shape of the plastic surface in the stress space ===
  params.addRangeCheckedParam<Real>(
      "tip_smoother",
      0.0,
      "tip_smoother>=0",
      "The plastic yield surface cone vertex at J2 = 0 will be smoothed by the given "
      "amount (physical unit: pressure). Typical value is 0.1*cohesion");
  params.addParam<Real>("curvature_yield",
                        0.0,
                        "Curvature of plastic yield in p-q space (physical unit: dimensionless). "
                        "In the paper of Khaledi et al. (2023) this parameter is named beta_1.");
  params.addParam<Real>("lode_angle_coefficient",
                        0.0,
                        "Lode angle dependency coefficient, between 0 and 0.7 (physical unit: "
                        "dimensionless). In the paper of "
                        "Khaledi et al. (2023) this parameter is named beta.");
  params.addParam<Real>(
      "Fs_function_power", -0.25, "Power of Fs function (physical unit: dimensionless).");
  params.addParam<Real>(
      "parameter_n0",
      3.0,
      "Anisotropy factor, i.e. the ratio between the undrained elastic moduli of P- and "
      "S-specimens. For Opalinus clay this parameter ranges between 2 and 3."
      "(physical unit: dimensionless). In the paper of Khaledi et al. (2023) this parameter is "
      "named n.");

  // === Hardening ===
  params.addParam<Real>("hardening_a0",
                        0.0,
                        "Strain hardening parameter (physical unit: dimensionless). "
                        "This parameter factorises the hardening. "
                        "In the paper of Khaledi et al. (2023) this parameter is named a_1.");

  params.addParam<Real>("hardening_eta",
                        1.0,
                        "Strain hardening parameter (physical unit: dimensionless). "
                        "This parameter is the exponent applied to the accumulated plastic strain. "
                        "In the paper of Khaledi et al. (2023) this parameter is named eta.");

  // === Damage ===
  params.addRequiredCoupledVar("nonlocal_variable", "This is the non-local variable for damage");
  params.addParam<Real>(
      "parameter_damageI",
      1.0,
      "Onset of damage. This value is compared to the value of the "
      "'nonlocal_variable' and if equal or exceeded, damage is considered."); // damage parameter
  params.addParam<Real>("parameter_damageA", 1.0, "parameter_damageA");       // damage parameter
  params.addParam<Real>("parameter_damageF", 1.0, "parameter_damageF");       // damage parameter
  params.addParam<Real>("parameter_damageN", 1.0, "parameter_damageN");       // damage parameter
  params.addParam<Real>("parameter_gammar",
                        0.3,
                        "Residual gamma (physical unit: dimensionless). "
                        "For Opalinus this value equals usually 0.045.");

  // === ??? ===
  params.addParam<Real>("psi_to_phi", 1.0, "psi_to_phi");

  // === ??? ===
  params.addParam<bool>("perfect_guess",
                        true,
                        "Provide a guess to the Newton-Raphson procedure "
                        "that is the result from perfect plasticity. With "
                        "severe hardening/softening this may be suboptimal.");

  // === parameter groups ===
  params.addParamNamesToGroup("gamma_mean parameter_omega_1 parameter_b_1", "Anisotropic Strength");
  params.addParamNamesToGroup(
      "tip_smoother curvature_yield lode_angle_coefficient Fs_function_power parameter_n0",
      "Shape of Plastic Surface in Stress Space");
  params.addParamNamesToGroup("hardening_a0 hardening_eta", "Strain Hardening");
  params.addParamNamesToGroup("nonlocal_variable parameter_damageI parameter_damageA "
                              "parameter_damageF parameter_damageN parameter_gammar",
                              "Post-Peak / Damage");

  return params;
}

DesaiHardeningStressUpdate::DesaiHardeningStressUpdate(const InputParameters & parameters)
  : MultiParameterPlasticityStressUpdate(parameters, 6, 1, 1),
    _psi_to_phi(getParam<Real>("psi_to_phi")),
    _St(parameters.get<Real>("p_tensile")),
    _perfect_guess(getParam<bool>("perfect_guess")),
    _small_smoother2(Utility::pow<2>(getParam<Real>("tip_smoother"))),
    _beta(getParam<Real>("lode_angle_coefficient")),
    _beta1(getParam<Real>("curvature_yield")),
    _mv(getParam<Real>("Fs_function_power")),
    _mean_gamma(parameters.get<Real>("gamma_mean")),
    _omega1(getParam<Real>("parameter_omega_1")),
    _b1(getParam<Real>("parameter_b_1")),
    _localCoordinateSystem(
        getUserObject<CartesianLocalCoordinateSystem>("local_coordinate_system")),
    _a0(getParam<Real>("hardening_a0")),
    _eta(getParam<Real>("hardening_eta")),
    _n0(getParam<Real>("parameter_n0")),
    _nonlocal_var(coupledValue("nonlocal_variable")),
    _gammar(getParam<Real>("parameter_gammar")),
    _dam_I(getParam<Real>("parameter_damageI")),
    _dam_A(getParam<Real>("parameter_damageA")),
    _dam_F(getParam<Real>("parameter_damageF")),
    _dam_N(getParam<Real>("parameter_damageN")),
    _alfa(declareProperty<Real>(_base_name + "alfa")),
    _gamma(declareProperty<Real>(_base_name + "gamma")),
    _gamma0(declareProperty<Real>(_base_name + "gamma0"))
//    _alfa_old(getMaterialPropertyOld<Real>(_base_name + "alfa")),
//    _gamma_old(getMaterialPropertyOld<Real>(_base_name + "gamma")),
//    _gamma0_old(getMaterialPropertyOld<Real>(_base_name + "gamma0"))
//  _gamma(parameters.get<Real>("gamma_mean")/std::sqrt(27.0)),
//  _gammaq(parameters.get<Real>("gamma_mean")*parameters.get<Real>("psi_to_phi")/std::sqrt(27.0)),
{
}

void
DesaiHardeningStressUpdate::computeStressParams(const RankTwoTensor & stress,
                                                std::vector<Real> & stress_params) const
{
  stress_params[0] = stress(0, 0);
  stress_params[1] = stress(1, 1);
  stress_params[2] = stress(2, 2);
  stress_params[5] = stress(1, 2);
  stress_params[4] = stress(0, 2);
  stress_params[3] = stress(0, 1);
}

std::vector<RankTwoTensor>
DesaiHardeningStressUpdate::dstress_param_dstress(const RankTwoTensor & stress) const
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
DesaiHardeningStressUpdate::d2stress_param_dstress(const RankTwoTensor & stress) const
{
  (void)stress;

  std::vector<RankFourTensor> d2(_num_sp);
  for (unsigned i = 0; i < _num_sp; ++i)
    d2[i] = RankFourTensor();
  return d2;
}

void
DesaiHardeningStressUpdate::preReturnMapV(const std::vector<Real> & /*trial_stress_params*/,
                                          const RankTwoTensor & stress_trial,
                                          const std::vector<Real> & /*intnl_old*/,
                                          const std::vector<Real> & /*yf*/,
                                          const RankFourTensor & Eijkl)
{
  (void)Eijkl; // silence compiler warning about unused variable

  std::vector<Real> stress_params(_num_sp);
  stress_params[0] = stress_trial(0, 0);
  stress_params[1] = stress_trial(1, 1);
  stress_params[2] = stress_trial(2, 2);
  stress_params[5] = stress_trial(1, 2);
  stress_params[4] = stress_trial(0, 2);
  stress_params[3] = stress_trial(0, 1);

  if (_t_step < 2)
  {
    const auto denom = std::sqrt(27.0);
    _gamma0[_qp] = _mean_gamma / denom;
    _gamma[_qp] = _mean_gamma / denom;
  }
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
  }
}

void
DesaiHardeningStressUpdate::setStressAfterReturnV(const RankTwoTensor & /*stress_trial*/,
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
DesaiHardeningStressUpdate::yieldFunctionValuesV(const std::vector<Real> & stress_params,
                                                 const std::vector<Real> & intnl,
                                                 std::vector<Real> & yf) const
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

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = std::max(1e-7, -(stress_now.trace()) + 3.0 * _St);
  Real j2 = stress_now.secondInvariant();
  Real q = std::pow(stress_now.secondInvariant(), 0.5);

  Real j3 = stress_now.thirdInvariant();
  if (j2 < 1e-8)
  {
    j3 = 0.0;
    j2 = 1e-8;
    // @Kavan-Khaledi do we have to update 'q' as well?
  }
  Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  Real fs = std::pow(std::exp(_beta1 * i1) + _beta * sr, _mv);

  _alfa[_qp] = _a0;
  if (_t_step < 2)
    _alfa[_qp] = _a0;
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
    _alfa[_qp] = _a0 * std::exp(-1.0 * intnl[0] / _eta);
  }
  // Real _gammaq = _psi_to_phi * _gamma[_qp];
  const Real fb = std::pow(
      std::pow(_gamma[_qp], 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);

  yf[0] = std::pow(std::pow(q, 2.0) + _small_smoother2, 0.5) - fb * fs;

  return;
}

void
DesaiHardeningStressUpdate::computeAllQV(const std::vector<Real> & stress_params,
                                         const std::vector<Real> & intnl,
                                         std::vector<yieldAndFlow> & all_q) const
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

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = std::max(1e-7, -(stress_now.trace()) + 3.0 * _St);
  Real j2 = stress_now.secondInvariant();
  Real q = std::pow(j2, 0.5);

  Real j3 = stress_now.thirdInvariant();
  if (j2 < 1e-8)
  {
    j3 = 0.0;
    j2 = 1e-8;
    // @Kavan-Khaledi do we have to update 'q' as well?
  }

  Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  Real fs = std::pow(std::exp(_beta1 * i1) + _beta * sr, _mv);

  _alfa[_qp] = _a0;
  if (_t_step < 2)
  {
    _alfa[_qp] = _a0;
    //   _alfa[_qp]xx=(j2+_small_smoother2)/pow(fs,2.0)/pow(i1,_n0)-pow(_gamma,2.0)*pow(i1,2.0-_n0);
  }
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
    _alfa[_qp] = _a0 * std::exp(-1.0 * intnl[0] / _eta);
  }

  Real _gammaq = _psi_to_phi * _gamma[_qp];
  const Real fb = std::pow(
      std::pow(_gamma[_qp], 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);
  const Real fbg =
      std::pow(std::pow(_gammaq, 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();
  const RankFourTensor d2j3ds_t = stress_now.d2thirdInvariant();
  // RankFourTensor d2j2ds_t = stress_now.d2secondInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = -1.0;
  di1ds[1] = -1.0;
  di1ds[2] = -1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 <= 1e-8)
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

  const Real dfbdi1 =
      (2.0 * std::pow(_gamma[_qp], 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fb);

  const Real dfbgdi1 =
      (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fbg);
  const Real dfsdi1 = _beta1 * _mv * std::exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv);
  const Real dfdi1 = -dfbdi1 * fs - fb * dfsdi1;
  const Real dgdi1 = -dfbgdi1 * fs - fbg * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dfdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fb * dfsdj2;
  const Real dgdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fbg * dfsdj2;
  const Real dfdj3 = -fb * dfsdj3;
  const Real dgdj3 = -fbg * dfsdj3;

  // yield functions.  See comment in yieldFunctionValuesV
  all_q[0].f = std::pow(std::pow(q, 2.0) + _small_smoother2, 0.5) - fb * fs;

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

  // @Kavan-Khaledi so dalfa_di is always Zero?
  Real dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * intnl[0] / _eta);
  if (_t_step < 2)
    dalfa_di = 0.0; //-(_a0/_eta)*std::exp(-1.0*intnl[0]/_eta);
  else
    dalfa_di = 0.0; //-(_a0/_eta)*std::exp(-1.0*intnl[0]/_eta);

  const Real dfb_dalfa = -std::pow(i1, _n0) / (2.0 * fb);
  const Real dfbg_dalfa = -std::pow(i1, _n0) / (2.0 * fbg);
  const Real d2fbg_di1_di =
      (-2.0 * _n0 * std::pow(i1, _n0 - 1.0) * fbg -
       2.0 * (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) *
           dfbg_dalfa) /
      (4.0 * std::pow(fbg, 2.0));
  const Real d2gdi1_alfa = -d2fbg_di1_di * fs - dfbg_dalfa * dfsdi1;
  const Real d2gdj2_alfa = -dfbg_dalfa * dfsdj2;
  const Real d2gdj3_alfa = -dfbg_dalfa * dfsdj3;

  // d(yield function)/d(intnl)

  all_q[0].df_di[0] = (-dfb_dalfa * fs) * dalfa_di;

  // d(flow potential)/d(stress_params)/d(intnl)
  for (unsigned yf = 0; yf < _num_yf; ++yf)
    for (unsigned a = 0; a < _num_sp; ++a)
      for (unsigned i = 0; i < _num_intnl; ++i)
        all_q[yf].d2g_di[a][i] =
            (d2gdi1_alfa * di1ds[a] + d2gdj2_alfa * dj2ds[a] + d2gdj3_alfa * dj3ds[a]) * dalfa_di;

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

  const Real d2fbgdi1i1 =
      ((2.0 * std::pow(_gammaq, 2.0) - _n0 * (_n0 - 1.0) * _alfa[_qp] * std::pow(i1, _n0 - 2.0)) *
           fbg -
       (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) * dfbgdi1) /
      (4.0 * std::pow(fbg, 2.0));
  const Real d2fsdi1i1 =
      std::pow(_beta1, 2.0) * _mv * exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv) +
      std::pow(_beta1, 2.0) * _mv * std::exp(2.0 * _beta1 * i1) * std::pow(fs, (_mv - 2.0) / _mv) *
          (_mv - 1.0);
  const Real d2fsdi1j2 = -(3.0 * _beta * _beta1 * _mv * cof * j3 * std::exp(_beta1 * i1) *
                           std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         (2.0 * std::pow(j2, 2.5));
  const Real d2fsdi1j3 = (_beta * _beta1 * _mv * cof * std::exp(_beta1 * i1) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 1.5);
  const Real d2fsdj2j2 = (15.0 * _beta * _mv * cof * j3 * std::pow(fs, (_mv - 1.0) / _mv)) /
                             (4.0 * std::pow(j2, 3.5)) +
                         (9.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(j3, 2.0) * std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                             (4.0 * std::pow(j2, 5.0));

  const Real d2fsdj2j3 =
      -(3.0 * _beta * _mv * cof * std::pow(fs, (_mv - 1.0) / _mv)) / (2.0 * pow(j2, 2.5)) -
      (3.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) * j3 *
       std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
          (2.0 * std::pow(j2, 4.0));
  const Real d2fsdj3j3 = (std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 3.0);

  const Real d2gdi1i1 = -(d2fbgdi1i1 * fs + dfbgdi1 * dfsdi1 + dfbgdi1 * dfsdi1 + fbg * d2fsdi1i1);
  const Real d2gdi1j2 = -(dfbgdi1 * dfsdj2 + fbg * d2fsdi1j2);
  const Real d2gdi1j3 = -(dfbgdi1 * dfsdj3 + fbg * d2fsdi1j3);
  const Real d2gdj2j2 =
      -(1.0 / (4.0 * std::pow(q, 2.0))) * std::pow(std::pow(q, 2.0) + _small_smoother2, -0.5) +
      (1.0 / (4.0 * std::pow(q, 2.0))) *
          (_small_smoother2 / std::pow(pow(q, 2.0) + _small_smoother2, 1.5)) -
      fbg * d2fsdj2j2;
  const Real d2gdj2j3 = -fbg * d2fsdj2j3;
  const Real d2gdj3j3 = -fbg * d2fsdj3j3;

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

  //   if(j2<=1e-8)
  //   {
  //    for (unsigned a = 0; a < _num_sp; ++a)
  //       for (unsigned b = 0; b < _num_sp; ++b)
  //          all_q[0].d2g[a][b] = 0.0;
  //   }

  return;
}

void
DesaiHardeningStressUpdate::setEffectiveElasticity(const RankFourTensor & Eijkl)
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
DesaiHardeningStressUpdate::initializeVarsV(const std::vector<Real> & trial_stress_params,
                                            const std::vector<Real> & intnl_old,
                                            std::vector<Real> & stress_params,
                                            Real & gaE,
                                            std::vector<Real> & intnl) const
{

  // Real _gammaq = _psi_to_phi * _gamma[_qp];

  if (_t_step < 2)
    _alfa[_qp] = _a0;
  else
    _alfa[_qp] = _a0 * std::exp(-1.0 * intnl_old[0] / _eta);

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
    const Real i1_trial =
        std::max(1e-7, -(stress_trial.trace()) + 3.0 * _St); //-std::sqrt(_small_smoother2)/_gamma;
    const Real q_trial = std::pow(stress_trial.secondInvariant(), 0.5);
    Real j2_trial = (stress_trial.secondInvariant());
    Real j3_trial = (stress_trial.thirdInvariant());

    if (j2_trial < 1e-8)
    {
      j3_trial = 0.0;
      j2_trial = 1e-8;
    }
    const Real cof = std::sqrt(27.0) / 2.0;

    Real sr = std::clamp(cof * j3_trial / std::pow(j2_trial, 1.5), -1.0, +1.0);

    Real fs_trial = std::pow(std::exp(_beta1 * i1_trial) + _beta * sr, _mv);

    const Real fb_trial = std::pow(std::pow(_gamma[_qp], 2.0) * std::pow(i1_trial, 2.0) -
                                       _alfa[_qp] * std::pow(i1_trial, _n0),
                                   0.5);

    const Real trial_desai_yf =
        std::pow(std::pow(q_trial, 2.0) + _small_smoother2, 0.5) - fb_trial * fs_trial;
    const Real dfb_dalfa = -std::pow(i1_trial, _n0) / (2.0 * fb_trial);

    Real dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * intnl[0] / _eta);
    if (_t_step < 2)
      dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * intnl[0] / _eta);

    else
      dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * intnl[0] / _eta);

    const Real df_di = (-dfb_dalfa * fs_trial) * dalfa_di;

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
      compute_dg(stress_params, intnl_old, dg);
      compute_df(stress_params, intnl_old, df);
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

      gaE = ga * (_Eij[0][0] + _Eij[1][1] + _Eij[2][2]);
      const RankTwoTensor stress_test = RankTwoTensor(stress_params[0],
                                                      stress_params[3],
                                                      stress_params[4],
                                                      stress_params[3],
                                                      stress_params[1],
                                                      stress_params[5],
                                                      stress_params[4],
                                                      stress_params[5],
                                                      stress_params[2]);
      const Real i1_test =
          -(stress_test.trace()) + 3.0 * _St; //-std::sqrt(_small_smoother2)/_gamma;
      const Real q_test = std::pow(stress_test.secondInvariant(), 0.5);

      if (i1_test <= 1e-7 && q_test < 1.0e-5)
      {

        for (unsigned i = 0; i < _num_sp; ++i)
          stress_params[i] = trial_stress_params[i];

        gaE = 0.0;

        // stress_params[0] = _St-std::sqrt(_small_smoother2) / (3.0 * _gamma);
        // stress_params[1] = _St-std::sqrt(_small_smoother2) / (3.0 * _gamma);
        // stress_params[2] = _St-std::sqrt(_small_smoother2) / (3.0 * _gamma);
        // stress_params[3] = 0.0;
        // stress_params[4] = 0.0;
        // stress_params[5] = 0.0;
        // gaE=(_Eij[0][0]+_Eij[1][1]+_Eij[2][2])*(trial_stress_params[0]-stress_params[0])/(_gammaq*(_Eij[0][0]+_Eij[1][1]+_Eij[2][2]));
      }
      found_solution = true;
    }
  }

  setIntnlValuesV(trial_stress_params, stress_params, intnl_old, intnl);

  return;
}

void
DesaiHardeningStressUpdate::setIntnlValuesV(const std::vector<Real> & trial_stress_params,
                                            const std::vector<Real> & current_stress_params,
                                            const std::vector<Real> & intnl_old,
                                            std::vector<Real> & intnl) const
{

  std::vector<Real> cp3(_num_sp);
  std::vector<Real> dg(_num_sp);
  compute_dg(current_stress_params, intnl_old, dg);

  const Real lamdaF = std::pow(dg[0], 2.0) + std::pow(dg[1], 2.0) + std::pow(dg[2], 2.0) +
                      0.5 * std::pow(dg[3], 2.0) + 0.5 * std::pow(dg[4], 2.0) +
                      0.5 * std::pow(dg[5], 2.0);

  for (unsigned i = 0; i < _num_sp; ++i)
  {
    cp3[i] = 0.0;
    for (unsigned j = 0; j < _num_sp; ++j)
      cp3[i] += _Eij[i][j] * dg[j];
  }

  intnl[0] = intnl_old[0] +
             ((trial_stress_params[0] - current_stress_params[0]) / cp3[0]) * std::sqrt(lamdaF);

  return;
}

void
DesaiHardeningStressUpdate::setIntnlDerivativesV(const std::vector<Real> & trial_stress_params,
                                                 const std::vector<Real> & current_stress_params,
                                                 const std::vector<Real> & intnl,
                                                 std::vector<std::vector<Real>> & dintnl) const
{
  (void)trial_stress_params;

  std::vector<Real> cp3(_num_sp);
  std::vector<Real> dg(_num_sp);
  compute_dg(current_stress_params, intnl, dg);

  // Real lamdaF = std::pow(dg[0], 2.0) + std::pow(dg[1], 2.0) + std::pow(dg[2], 2.0) +
  //               0.5 * std::pow(dg[3], 2.0) + 0.5 * std::pow(dg[4], 2.0) +
  //               0.5 * std::pow(dg[5], 2.0);

  for (unsigned i = 0; i < _num_sp; ++i)
  {
    cp3[i] = 0.0;
    for (unsigned j = 0; j < _num_sp; ++j)
      cp3[i] += _Eij[i][j] * dg[j];
  }

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
    dintnl[0][i] =
        0.0; //((-1.0*cp3[i]-(trial_stress_params[i]-current_stress_params[i])*cd2g[i])/std::pow(cp3[i],2.0));

  return;
}

void
DesaiHardeningStressUpdate::consistentTangentOperatorV(
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

  RankTwoTensor dqdsig;
  RankTwoTensor dFdsig;
  RankFourTensor d2gds;
  RankTwoTensor gmatrix;
  RankTwoTensor temp;
  RankTwoTensor cdgdsig;
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
  const Real i1 = std::max(1e-7, -(stress_now.trace()) + 3.0 * _St);
  Real j2 = stress_now.secondInvariant();
  Real q = std::pow(stress_now.secondInvariant(), 0.5);
  Real j3 = stress_now.thirdInvariant();

  if (j2 < 1e-8)
  {
    j3 = 0.0;
    j2 = 1e-8;
    // @Kavan-Khaledi do we have to update 'q' as well?
  }

  Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  Real fs = std::pow(std::exp(_beta1 * i1) + _beta * sr, _mv);

  if (_t_step < 2)
  {
  }
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
  }
  Real _gammaq = _psi_to_phi * _gamma[_qp];
  const Real fb = std::pow(
      std::pow(_gamma[_qp], 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);
  const Real fbg =
      std::pow(std::pow(_gammaq, 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);

  const RankTwoTensor di1ds_t = -1.0 * stress_now.dtrace();
  RankTwoTensor dj2ds_t = stress_now.deviatoric();
  RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();
  RankFourTensor d2j3ds_t = stress_now.d2thirdInvariant();
  RankFourTensor d2j2ds_t = stress_now.d2secondInvariant();
  if (j2 <= 1e-8)
  {
    dj2ds_t = RankTwoTensor();
    dj3ds_t = RankTwoTensor();
    d2j2ds_t = RankFourTensor();
    d2j3ds_t = RankFourTensor();
  }

  const Real dfbdi1 =
      (2.0 * std::pow(_gamma[_qp], 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fb);
  const Real dfbgdi1 =
      (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fbg);
  const Real dfsdi1 = _beta1 * _mv * std::exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv);
  const Real dfdi1 = -(dfbdi1 * fs + fb * dfsdi1);
  const Real dgdi1 = -(dfbgdi1 * fs + fbg * dfsdi1);
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dfdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fb * dfsdj2;
  const Real dgdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fbg * dfsdj2;
  const Real dfdj3 = -fb * dfsdj3;
  const Real dgdj3 = -fbg * dfsdj3;

  const Real dfb_dalfa = -std::pow(i1, _n0) / (2.0 * fb);

  Real dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * _intnl[_qp][0] / _eta);
  if (_t_step < 2)
    dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * _intnl[_qp][0] / _eta);
  else
    dalfa_di = -(_a0 / _eta) * std::exp(-1.0 * _intnl[_qp][0] / _eta);

  const Real df_di = (-dfb_dalfa * fs) * dalfa_di;

  dFdsig = dfdi1 * di1ds_t + dfdj2 * dj2ds_t + dfdj3 * dj3ds_t;
  dqdsig = dgdi1 * di1ds_t + dgdj2 * dj2ds_t + dgdj3 * dj3ds_t;

  const Real d2fbgdi1i1 =
      ((2.0 * std::pow(_gammaq, 2.0) - _n0 * (_n0 - 1.0) * _alfa[_qp] * std::pow(i1, _n0 - 2.0)) *
           fbg -
       (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) * dfbgdi1) /
      (4.0 * std::pow(fbg, 2.0));
  const Real d2fsdi1i1 =
      std::pow(_beta1, 2.0) * _mv * exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv) +
      std::pow(_beta1, 2.0) * _mv * std::exp(2.0 * _beta1 * i1) * std::pow(fs, (_mv - 2.0) / _mv) *
          (_mv - 1.0);
  const Real d2fsdi1j2 = -(3.0 * _beta * _beta1 * _mv * cof * j3 * std::exp(_beta1 * i1) *
                           std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         (2.0 * std::pow(j2, 2.5));
  const Real d2fsdi1j3 = (_beta * _beta1 * _mv * cof * std::exp(_beta1 * i1) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 1.5);
  const Real d2fsdj2j2 = (15.0 * _beta * _mv * cof * j3 * std::pow(fs, (_mv - 1.0) / _mv)) /
                             (4.0 * std::pow(j2, 3.5)) +
                         (9.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(j3, 2.0) * std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                             (4.0 * std::pow(j2, 5.0));

  const Real d2fsdj2j3 =
      -(3.0 * _beta * _mv * cof * std::pow(fs, (_mv - 1.0) / _mv)) / (2.0 * pow(j2, 2.5)) -
      (3.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) * j3 *
       std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
          (2.0 * std::pow(j2, 4.0));
  const Real d2fsdj3j3 = (std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 3.0);

  const Real d2gdi1i1 = -(d2fbgdi1i1 * fs + dfbgdi1 * dfsdi1 + dfbgdi1 * dfsdi1 + fbg * d2fsdi1i1);
  const Real d2gdi1j2 = -(dfbgdi1 * dfsdj2 + fbg * d2fsdi1j2);
  const Real d2gdi1j3 = -(dfbgdi1 * dfsdj3 + fbg * d2fsdi1j3);
  const Real d2gdj2j2 =
      -(1.0 / (4.0 * std::pow(q, 2.0))) * std::pow(std::pow(q, 2.0) + _small_smoother2, -0.5) +
      (1.0 / (4.0 * std::pow(q, 2.0))) *
          (_small_smoother2 / std::pow(pow(q, 2.0) + _small_smoother2, 1.5)) -
      fbg * d2fsdj2j2;
  const Real d2gdj2j3 = -fbg * d2fsdj2j3;
  const Real d2gdj3j3 = -fbg * d2fsdj3j3;

  d2gds = d2gdi1i1 * di1ds_t.outerProduct(di1ds_t) + d2gdi1j2 * di1ds_t.outerProduct(dj2ds_t) +
          d2gdi1j3 * di1ds_t.outerProduct(dj3ds_t) + d2gdi1j2 * dj2ds_t.outerProduct(di1ds_t) +
          d2gdj2j2 * dj2ds_t.outerProduct(dj2ds_t) + d2gdj2j3 * dj2ds_t.outerProduct(dj3ds_t) +
          d2gdi1j3 * dj3ds_t.outerProduct(di1ds_t) + d2gdj2j3 * dj3ds_t.outerProduct(dj2ds_t) +
          d2gdj3j3 * dj3ds_t.outerProduct(dj3ds_t) + dgdj2 * d2j2ds_t + dgdj3 * d2j3ds_t;

  RankFourTensor inv =
      RankFourTensor(RankFourTensor::initIdentityFour) + (ga)*elasticity_tensor * d2gds;

  inv = inv.transposeMajor().invSymm();
  inv = (elasticity_tensor.transposeMajor() * inv).transposeMajor();
  //  temp=inv*bximatrix;

  for (unsigned i = 0; i < _tensor_dimensionality; ++i)
    for (unsigned j = 0; j < _tensor_dimensionality; ++j)
      for (unsigned k = 0; k < _tensor_dimensionality; ++k)
        for (unsigned l = 0; l < _tensor_dimensionality; ++l)
          cdgdsig(i, j) += inv(i, j, k, l) * dFdsig(k, l);

  for (unsigned i = 0; i < _tensor_dimensionality; ++i)
    for (unsigned j = 0; j < _tensor_dimensionality; ++j)
      for (unsigned k = 0; k < _tensor_dimensionality; ++k)
        for (unsigned l = 0; l < _tensor_dimensionality; ++l)
          temp(i, j) += inv(k, l, i, j) * dqdsig(k, l);

  gmatrix = cdgdsig / (dFdsig.doubleContraction(temp) - df_di);

  cto = inv - (temp.outerProduct(gmatrix));

  return;
}

void
DesaiHardeningStressUpdate::compute_dg(const std::vector<Real> & stress_params,
                                       const std::vector<Real> & intnl,
                                       std::vector<Real> & dgg) const
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

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = std::max(1e-7, -(stress_now.trace()) + 3.0 * _St);
  Real j2 = stress_now.secondInvariant();
  Real q = std::pow(j2, 0.5);
  Real j3 = stress_now.thirdInvariant();
  if (j2 < 1e-8)
  {
    j3 = 0.0;
    j2 = 1e-8;
    // @Kavan-Khaledi do we have to update 'q' as well?
  }

  Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  Real fs = std::pow(std::exp(_beta1 * i1) + _beta * sr, _mv);

  if (_t_step < 2)
    _alfa[_qp] = _a0;
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
    _alfa[_qp] = _a0 * std::exp(-1.0 * intnl[0] / _eta);
  }
  Real _gammaq = _psi_to_phi * _gamma[_qp];
  const Real fbg =
      std::pow(std::pow(_gammaq, 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = -1.0;
  di1ds[1] = -1.0;
  di1ds[2] = -1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 <= 1e-8)
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

  const Real dfbgdi1 =
      (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fbg);
  const Real dfsdi1 = _beta1 * _mv * std::exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv);
  const Real dgdi1 = -dfbgdi1 * fs - fbg * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);
  const Real dgdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fbg * dfsdj2;
  const Real dgdj3 = -fbg * dfsdj3;

  for (unsigned a = 0; a < 3; ++a)
    dgg[a] = dgdi1 * di1ds[a] + dgdj2 * dj2ds[a] + dgdj3 * dj3ds[a];

  dgg[3] = 2.0 * (dgdi1 * di1ds[3] + dgdj2 * dj2ds[3] + dgdj3 * dj3ds[3]);
  dgg[4] = 2.0 * (dgdi1 * di1ds[4] + dgdj2 * dj2ds[4] + dgdj3 * dj3ds[4]);
  dgg[5] = 2.0 * (dgdi1 * di1ds[5] + dgdj2 * dj2ds[5] + dgdj3 * dj3ds[5]);

  return;
}

void
DesaiHardeningStressUpdate::compute_df(const std::vector<Real> & stress_params,
                                       const std::vector<Real> & intnl,
                                       std::vector<Real> & dff) const
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

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = std::max(1e-7, -(stress_now.trace()) + 3.0 * _St);
  Real j2 = stress_now.secondInvariant();
  Real q = std::pow(j2, 0.5);
  Real j3 = stress_now.thirdInvariant();
  if (j2 < 1e-8)
  {
    j3 = 0.0;
    j2 = 1e-8;
    // @Kavan-Khaledi do we have to update 'q' as well?
  }

  Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  Real fs = std::pow(std::exp(_beta1 * i1) + _beta * sr, _mv);

  if (_t_step < 2)
    _alfa[_qp] = _a0;
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
    _alfa[_qp] = _a0 * std::exp(-1.0 * intnl[0] / _eta);
  }
  // Real _gammaq = _psi_to_phi * _gamma[_qp];
  const Real fb = std::pow(
      std::pow(_gamma[_qp], 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = -1.0;
  di1ds[1] = -1.0;
  di1ds[2] = -1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 <= 1e-8)
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

  const Real dfbdi1 =
      (2.0 * std::pow(_gamma[_qp], 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fb);
  const Real dfsdi1 = _beta1 * _mv * std::exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv);
  const Real dfdi1 = -dfbdi1 * fs - fb * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  const Real dfdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fb * dfsdj2;
  const Real dfdj3 = -fb * dfsdj3;

  for (unsigned a = 0; a < 3; ++a)
    dff[a] = dfdi1 * di1ds[a] + dfdj2 * dj2ds[a] + dfdj3 * dj3ds[a];

  dff[3] = 2.0 * (dfdi1 * di1ds[3] + dfdj2 * dj2ds[3] + dfdj3 * dj3ds[3]);
  dff[4] = 2.0 * (dfdi1 * di1ds[4] + dfdj2 * dj2ds[4] + dfdj3 * dj3ds[4]);
  dff[5] = 2.0 * (dfdi1 * di1ds[5] + dfdj2 * dj2ds[5] + dfdj3 * dj3ds[5]);

  return;
}

void
DesaiHardeningStressUpdate::setGammaValue(const std::vector<Real> & stress_params,
                                          Real & gamma) const
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

  // const auto e1 = _localCoordinateSystem.e1();
  // const auto e2 = _localCoordinateSystem.e2();
  const auto e3 = _localCoordinateSystem.e3();

  (stress_now).symmetricEigenvaluesEigenvectors(lamda, eigv);

  const Real cosb = eigv(0, 0) * e3(2) + eigv(1, 0) * e3(1) + eigv(2, 0) * e3(0);
  // const Real sinb = std::sqrt(1.0-cosb*cosb);
  // const Real cosa = eigv(0,1)*rotation(0,2)+eigv(1,1)*rotation(1,2)+eigv(2,1)*rotation(2,2);
  // const Real sina = std::sqrt(1.0-cosa*cosa);

  // const Real lv2 =
  //   (pow(lamda[0], 2.0) * pow(cosb, 2.0) + pow(lamda[1], 2.0) * pow(cosa, 2.0) * pow(sinb, 2.0) +
  //    pow(lamda[2], 2.0) * pow(sina, 2.0) * pow(sinb, 2.0)) /
  //   (pow(lamda[0], 2.0) + pow(lamda[1], 2.0) + pow(lamda[2], 2.0));
  const Real lv2 = pow(cosb, 2.0);

  gamma = (_mean_gamma * (1.0 + _omega1 * (1.0 - 3.0 * lv2) +
                          _b1 * pow(_omega1, 2.0) * pow((1.0 - 3.0 * lv2), 2.0))) /
          std::sqrt(27.0);

  return;
}

void
DesaiHardeningStressUpdate::compute_d2g(const std::vector<Real> & stress_params,
                                        const std::vector<Real> & intnl,
                                        std::vector<std::vector<Real>> & d2gg) const
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

  const Real cof = std::sqrt(27.0) / 2.0;
  const Real i1 = std::max(1e-7, -(stress_now.trace()) + 3.0 * _St);
  Real j2 = stress_now.secondInvariant();
  Real q = std::pow(stress_now.secondInvariant(), 0.5);

  Real j3 = stress_now.thirdInvariant();
  if (j2 < 1e-8)
  {
    j3 = 0.0;
    j2 = 1e-8;
    // @Kavan-Khaledi do we have to update 'q' as well?
  }

  Real sr = std::clamp(cof * j3 / std::pow(j2, 1.5), -1.0, +1.0);

  Real fs = std::pow(std::exp(_beta1 * i1) + _beta * sr, _mv);

  _alfa[_qp] = _a0;
  if (_t_step < 2)
  {
    _alfa[_qp] = _a0;
    //   _alfa[_qp]xx=(j2+_small_smoother2)/pow(fs,2.0)/pow(i1,_n0)-pow(_gamma,2.0)*pow(i1,2.0-_n0);
  }
  else
  {
    setGammaValue(stress_params, _gamma0[_qp]);
    setGammaDamage(_gamma0[_qp], _gamma[_qp]);
    _alfa[_qp] = _a0 * std::exp(-1.0 * intnl[0] / _eta);
  }
  Real _gammaq = _psi_to_phi * _gamma[_qp];
  // const Real fb =
  //     std::pow(std::pow(_gamma[_qp], 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0),
  //     0.5);
  const Real fbg =
      std::pow(std::pow(_gammaq, 2.0) * std::pow(i1, 2.0) - _alfa[_qp] * std::pow(i1, _n0), 0.5);

  const RankTwoTensor dj2ds_t = stress_now.deviatoric();
  const RankTwoTensor dj3ds_t = stress_now.dthirdInvariant();
  const RankFourTensor d2j3ds_t = stress_now.d2thirdInvariant();
  // RankFourTensor d2j2ds_t = stress_now.d2secondInvariant();

  std::vector<Real> di1ds(_num_sp);
  std::vector<Real> dj2ds(_num_sp);
  std::vector<Real> dj3ds(_num_sp);
  di1ds[0] = -1.0;
  di1ds[1] = -1.0;
  di1ds[2] = -1.0;
  di1ds[5] = 0.0;
  di1ds[4] = 0.0;
  di1ds[3] = 0.0;

  if (j2 <= 1e-8)
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

  // const Real dfbdi1 =
  //     (2.0 * std::pow(_gamma[_qp], 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
  //     (2.0 * fb);

  const Real dfbgdi1 =
      (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) /
      (2.0 * fbg);
  const Real dfsdi1 = _beta1 * _mv * std::exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv);
  // const Real dfdi1 = -dfbdi1 * fs - fb * dfsdi1;
  // const Real dgdi1 = -dfbgdi1 * fs - fbg * dfsdi1;
  const Real dfsdj2 =
      (-1.5 * _beta * cof * j3 * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 2.5);
  const Real dfsdj3 = (_beta * cof * _mv * std::pow(fs, (_mv - 1.0) / _mv)) / std::pow(j2, 1.5);

  // const Real dfdj2 =
  //     (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fb * dfsdj2;
  const Real dgdj2 =
      (1.0 / (2.0 * q)) * q / (std::pow(pow(q, 2.0) + _small_smoother2, 0.5)) - fbg * dfsdj2;
  // const Real dfdj3 = -fb * dfsdj3;
  const Real dgdj3 = -fbg * dfsdj3;

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

  const Real d2fbgdi1i1 =
      ((2.0 * std::pow(_gammaq, 2.0) - _n0 * (_n0 - 1.0) * _alfa[_qp] * std::pow(i1, _n0 - 2.0)) *
           fbg -
       (2.0 * std::pow(_gammaq, 2.0) * i1 - _n0 * _alfa[_qp] * std::pow(i1, _n0 - 1.0)) * dfbgdi1) /
      (4.0 * std::pow(fbg, 2.0));
  const Real d2fsdi1i1 =
      std::pow(_beta1, 2.0) * _mv * exp(_beta1 * i1) * std::pow(fs, (_mv - 1.0) / _mv) +
      std::pow(_beta1, 2.0) * _mv * std::exp(2.0 * _beta1 * i1) * std::pow(fs, (_mv - 2.0) / _mv) *
          (_mv - 1.0);
  const Real d2fsdi1j2 = -(3.0 * _beta * _beta1 * _mv * cof * j3 * std::exp(_beta1 * i1) *
                           std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         (2.0 * std::pow(j2, 2.5));
  const Real d2fsdi1j3 = (_beta * _beta1 * _mv * cof * std::exp(_beta1 * i1) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 1.5);
  const Real d2fsdj2j2 = (15.0 * _beta * _mv * cof * j3 * std::pow(fs, (_mv - 1.0) / _mv)) /
                             (4.0 * std::pow(j2, 3.5)) +
                         (9.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(j3, 2.0) * std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                             (4.0 * std::pow(j2, 5.0));

  const Real d2fsdj2j3 =
      -(3.0 * _beta * _mv * cof * std::pow(fs, (_mv - 1.0) / _mv)) / (2.0 * pow(j2, 2.5)) -
      (3.0 * std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) * j3 *
       std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
          (2.0 * std::pow(j2, 4.0));
  const Real d2fsdj3j3 = (std::pow(_beta, 2.0) * _mv * std::pow(cof, 2.0) *
                          std::pow(fs, (_mv - 2.0) / _mv) * (_mv - 1.0)) /
                         std::pow(j2, 3.0);

  const Real d2gdi1i1 = -(d2fbgdi1i1 * fs + dfbgdi1 * dfsdi1 + dfbgdi1 * dfsdi1 + fbg * d2fsdi1i1);
  const Real d2gdi1j2 = -(dfbgdi1 * dfsdj2 + fbg * d2fsdi1j2);
  const Real d2gdi1j3 = -(dfbgdi1 * dfsdj3 + fbg * d2fsdi1j3);
  const Real d2gdj2j2 =
      -(1.0 / (4.0 * std::pow(q, 2.0))) * std::pow(std::pow(q, 2.0) + _small_smoother2, -0.5) +
      (1.0 / (4.0 * std::pow(q, 2.0))) *
          (_small_smoother2 / std::pow(pow(q, 2.0) + _small_smoother2, 1.5)) -
      fbg * d2fsdj2j2;
  const Real d2gdj2j3 = -fbg * d2fsdj2j3;
  const Real d2gdj3j3 = -fbg * d2fsdj3j3;

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

void
DesaiHardeningStressUpdate::setGammaDamage(const Real & gamma0, Real & gamma) const
{
  // here we define damage-dependent material parameters. If damage is zero, then they are equal to
  // their initial value
  Real softVar = _nonlocal_var[_qp];

  if (softVar < _dam_I)
    gamma = gamma0;
  else
    gamma = _gammar +
            (gamma0 - _gammar) * (exp(-1.0 * _dam_A * pow((softVar - _dam_I) / _dam_F, _dam_N)));

  return;
}
