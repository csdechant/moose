//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADVectorSupplier.h"

registerMooseObject("MooseApp", ADVectorSupplier);

InputParameters
ADVectorSupplier::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription(
      "Base object to supply a field vector that is defined eiter by a suppled vector variable or "
      "the gradient of a suppled scalar variable.");
  params.addCoupledVar("field", "The coupled vector or scalar variable to defined the field");
  return params;
}

ADVectorSupplier::ADVectorSupplier(const InputParameters & parameters)
  : ADKernel(parameters),
    _field_var(*getFieldVar("field", 0)),
    _is_vector(_field_var.isVector()),
    _vector_field(_is_vector ? adCoupledVectorValue("field") : _ad_grad_zero),
    _grad_scalar_field(_is_vector ? _ad_grad_zero : adCoupledGradient("field"))
{
}

/*
 *  NOTE: Not overriding computeQpResidual() causes an error,
 *  so computeQpResidual() = 0.0
 */
ADReal
ADVectorSupplier::computeQpResidual()
{
  return 0.0;
}

ADRealVectorValue
ADVectorSupplier::computeQpFieldValue(const Real & scaler)
{
  if (_is_vector)
    return _vector_field[_qp];
  else
    return scaler * _grad_scalar_field[_qp];
}
