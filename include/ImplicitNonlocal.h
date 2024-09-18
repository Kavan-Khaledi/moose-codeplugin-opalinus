#pragma once

#include "Kernel.h"
//


/**
 * This kernel implements the Laplacian operator:
 * $\nabla u \cdot \nabla \phi_i$
 */
class ImplicitNonlocal : public Kernel
{
public:
  static InputParameters validParams();

  ImplicitNonlocal(const InputParameters & parameters);

private:
 Real _length_scale;



protected:

 const MaterialProperty<std::vector<Real>> & _intnl;
 
 const MaterialProperty<RankTwoTensor> & _plastic_strain;
  const MaterialProperty<RankTwoTensor> & _total_strain;

  virtual Real computeQpResidual() override;

  virtual Real computeQpJacobian() override;
};