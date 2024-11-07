//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernel.h"

/**
 *  A base class that supplies a field vector that is defined
 *  depending on the suppled variable:
 *  
 *  1. If user supplies a scalar variable type:
 *     field = scaler * grad_scalar_variable, where scaler is a scaling factor defined
 *     in the inheriting object
 * 
 *  2. If user supplies a vector variable type:
 *     field = vector_variable
 */
class ADVectorSupplier : public ADKernel
{
public:
  static InputParameters validParams();

  ADVectorSupplier(const InputParameters & parameters);

protected:
  /*
   *  NOTE: Not overriding computeQpResidual() causes an error,
   *  so computeQpResidual() = 0.0 in body file
   */
  virtual ADReal computeQpResidual() override;

  /*
   *  Function the defines the field depending on supplied variable type.
   *  scaler is the scaling factor to multiplied by the gradient of the variable
   *  if a scalar variable is supplied.
   */
  virtual ADRealVectorValue computeQpFieldValue(const Real & scaler);

  /*
   *  The following of made protected instead of private the in cases the the inheriting
   *  object needs to know if the supplied variable is a vector or scalar.
   */
  /// The variable data of the supplied variable for the field
  const MooseVariableFieldBase & _field_var;
  /// True is the supplied variable is a vector
  const bool _is_vector;

private:
  /// The field defined from a vector variable
  const ADVectorVariableValue & _vector_field;
  /// The field defined from the gradient of a scalar variable
  const ADVariableGradient & _grad_scalar_field;
};
