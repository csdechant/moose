//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADJouleHeatingSource_WithMaterialField.h"

registerMooseObject("ElectromagneticsApp", ADJouleHeatingSource_WithMaterialField);

InputParameters
ADJouleHeatingSource_WithMaterialField::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addCoupledVar("elec", "Electrostatic potential for joule heating.");
  params.addParam<MaterialPropertyName>(
      "heating_term", "field_heating", "Material property providing the Joule Heating.");
  params.addParam<MaterialPropertyName>(
      "electrical_conductivity",
      "electrical_conductivity",
      "Material property providing electrical conductivity of the material.");
  params.addClassDescription("Calculates the heat source term corresponding to electrostatic Joule "
                             "heating, with Jacobian contributions calculated using the automatic "
                             "differentiation system.");
  return params;
}

ADJouleHeatingSource_WithMaterialField::ADJouleHeatingSource_WithMaterialField(const InputParameters & parameters)
  : ADKernelValue(parameters),

    _supplied_potential(isParamValid("elec")),

    _grad_potential(_supplied_potential ? adCoupledGradient("elec") : _ad_grad_zero),
    _elec_cond(_supplied_potential ? getADMaterialProperty<Real>("electrical_conductivity")
                                   : getGenericZeroMaterialProperty<Real, true>()),

    _heating_residual(_supplied_potential ? getGenericZeroMaterialProperty<Real, true>()
                                 : getADMaterialProperty<Real>("heating_term"))
{
}

ADReal
ADJouleHeatingSource_WithMaterialField::precomputeQpResidual()
{
  if (_supplied_potential)
    return -_elec_cond[_qp] * _grad_potential[_qp] * _grad_potential[_qp];
  else
    return -_heating_residual[_qp];
}
