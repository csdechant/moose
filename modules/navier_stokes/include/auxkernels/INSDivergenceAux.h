<<<<<<< HEAD
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

// Forward Declarations
=======
#ifndef INSDIVERGENCEAUX_H
#define INSDIVERGENCEAUX_H

#include "AuxKernel.h"

//Forward Declarations
class INSDivergenceAux;

template<>
InputParameters validParams<INSDivergenceAux>();
>>>>>>> Merging Modules into MOOSE #2460

/**
 * Computes h_min / |u|
 */
class INSDivergenceAux : public AuxKernel
{
public:
<<<<<<< HEAD
  static InputParameters validParams();

  INSDivergenceAux(const InputParameters & parameters);
=======
  INSDivergenceAux(const std::string & name, InputParameters parameters);
>>>>>>> Merging Modules into MOOSE #2460

  virtual ~INSDivergenceAux() {}

protected:
  virtual Real computeValue();

  // Velocity gradients
<<<<<<< HEAD
  const VariableGradient & _grad_u_vel;
  const VariableGradient & _grad_v_vel;
  const VariableGradient & _grad_w_vel;
};
=======
  VariableGradient& _grad_u_vel;
  VariableGradient& _grad_v_vel;
  VariableGradient& _grad_w_vel;
};

#endif //VELOCITYAUX_H
>>>>>>> Merging Modules into MOOSE #2460
