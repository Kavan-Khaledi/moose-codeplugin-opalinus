//PerfectlyPlastic Desai Model beta1=0.0 beta=0.0 alfa=0.0 The model is basically equivalent to Drucker-Prager model

#pragma once

#include "MultiParameterPlasticityStressUpdate.h"
#include "TensorMechanicsHardeningModel.h"
#include <cmath>
#include "libmesh/libmesh.h"

class DesaiHardeningStressUpdate : public MultiParameterPlasticityStressUpdate
{
public:
  static InputParameters validParams();

  DesaiHardeningStressUpdate(const InputParameters & parameters);

  /**
   * Does the model require the elasticity tensor to be isotropic?
   */
  bool requiresIsotropicTensor() override { return false; }

  bool isIsotropic() override { return false; };

protected:
  const Real _psi_to_phi;
  const Real _St;
  const bool _perfect_guess;
  const Real _small_smoother2;
  const Real _beta;
  const Real _betta1;
  const Real _mv;
  const Real _mean_gama;
  const Real _omega1;
  const Real _b1;
  const Real _dip;
  const Real _strike;
  const Real _a0;
  const Real _eta;
  const Real _n0;
  const Real _kesi0;
  const VariableValue & _nonlocal_var;
  const Real  _gamar;
  const Real _dam_I;
  const Real _dam_A;
  const Real _dam_F;
  const Real _dam_N;
  /// plastic strain
  MaterialProperty<Real> & _alfa;

  /// Old value of plastic strain
  const MaterialProperty<Real> & _alfa_old;
  
    /// plastic strain
  MaterialProperty<Real> & _gama;

  /// Old value of plastic strain
  const MaterialProperty<Real> & _gama_old;
      /// plastic strain
  MaterialProperty<Real> & _gama0;

  /// Old value of plastic strain
  const MaterialProperty<Real> & _gama0_old;

  void computeStressParams(const RankTwoTensor & stress,
                           std::vector<Real> & stress_params) const override;

  std::vector<RankTwoTensor> dstress_param_dstress(const RankTwoTensor & stress) const override;

  std::vector<RankFourTensor> d2stress_param_dstress(const RankTwoTensor & stress) const override;

  virtual void setStressAfterReturnV(const RankTwoTensor & stress_trial,
                                     const std::vector<Real> & stress_params,
                                     Real gaE,
                                     const std::vector<Real> & intnl,
                                     const yieldAndFlow & smoothed_q,
                                     const RankFourTensor & Eijkl,
                                     RankTwoTensor & stress) const override;

  virtual void preReturnMapV(const std::vector<Real> & trial_stress_params,
                             const RankTwoTensor & stress_trial,
                             const std::vector<Real> & intnl_old,
                             const std::vector<Real> & yf,
                             const RankFourTensor & Eijkl) override;

  void setEffectiveElasticity(const RankFourTensor & Eijkl) override;

  void yieldFunctionValuesV(const std::vector<Real> & stress_params,
                            const std::vector<Real> & intnl,
                            std::vector<Real> & yf) const override;

  void computeAllQV(const std::vector<Real> & stress_params,
                    const std::vector<Real> & intnl,
                    std::vector<yieldAndFlow> & all_q) const override;

  void initializeVarsV(const std::vector<Real> & trial_stress_params,
                       const std::vector<Real> & intnl_old,
                       std::vector<Real> & stress_params,
                       Real & gaE,
                       std::vector<Real> & intnl) const override;

  void setIntnlValuesV(const std::vector<Real> & trial_stress_params,
                       const std::vector<Real> & current_stress_params,
                       const std::vector<Real> & intnl_old,
                       std::vector<Real> & intnl) const override;

  void setIntnlDerivativesV(const std::vector<Real> & trial_stress_params,
                            const std::vector<Real> & current_stress_params,
                            const std::vector<Real> & intnl,
                            std::vector<std::vector<Real>> & dintnl) const override;

  virtual void consistentTangentOperatorV(const RankTwoTensor & stress_trial,
                                          const std::vector<Real> & trial_stress_params,
                                          const RankTwoTensor & stress,
                                          const std::vector<Real> & stress_params,
                                          Real gaE,
                                          const yieldAndFlow & smoothed_q,
                                          const RankFourTensor & Eijkl,
                                          bool compute_full_tangent_operator,
                                          const std::vector<std::vector<Real>> & dvar_dtrial,
                                          RankFourTensor & cto) override;
                                          
virtual void compute_dg(const std::vector<Real> & stress_params,
                                             const std::vector<Real> & intnl,
                                             std::vector<Real> & dgg) const;
virtual void compute_df(const std::vector<Real> & stress_params,
                                             const std::vector<Real> & intnl,
                                             std::vector<Real> & dff) const;
virtual void setGamaValue(const std::vector<Real> & stress_params,
                                                  Real & gama) const;

virtual void compute_d2g(const std::vector<Real> & stress_params,
                                             const std::vector<Real> & intnl,
                                              std::vector<std::vector<Real>> & d2gg) const;
virtual void  setGamaDamage(const Real & gama0,Real & gama) const;

private:


};