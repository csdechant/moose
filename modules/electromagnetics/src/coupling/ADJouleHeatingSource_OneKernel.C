//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ADJouleHeatingSource_OneKernel.h"
#include "FormulationEnums.h"

registerMooseObject("ElectromagneticsApp", ADJouleHeatingSource_OneKernel);

InputParameters
ADJouleHeatingSource_OneKernel::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addCoupledVar("elec", "Electrostatic potential for joule heating.");
  params.deprecateParam("elec", "field", "12/31/2024");
  params.addCoupledVar(
      "field", "The electic field vector or electrostatic potential scalar to produce the field.");
  params.addCoupledVar(
      "complex_field", "The complex component of the electic field vector for the harmonic formulation.");
  params.addParam<MaterialPropertyName>(
      "electrical_conductivity",
      "electrical_conductivity",
      "Material property providing electrical conductivity of the material.");
  params.addParam<Real>("value", 1.0, "Coefficient to multiply by heating term.");
  MooseEnum formulation("standard harmonic", "standard");
  params.addParam<MooseEnum>("formulation", formulation, "The formulation of the Joule heating (standard or harmonic).");
  params.addClassDescription("Calculates the heat source term corresponding to Joule "
                             "heating, with Jacobian contributions calculated using the automatic "
                             "differentiation system.");
  return params;
}

ADJouleHeatingSource_OneKernel::ADJouleHeatingSource_OneKernel(const InputParameters & parameters)
  : ADKernelValue(parameters),
    _field_var(*getFieldVar("field", 0)),
    _is_vector(_field_var.isVector()),
    _efield(_is_vector ? adCoupledVectorValue("field") : _ad_grad_zero),
    _efield_complex(_is_vector ? adCoupledVectorValue("complex_field") : _ad_grad_zero),
    _grad_potential(_is_vector ?  _ad_grad_zero : adCoupledGradient("field")),
    _elec_cond(getADMaterialProperty<Real>("electrical_conductivity")),
    _scale(getParam<Real>("value")),
    _formulation(getParam<MooseEnum>("formulation"))
{
  if ((_formulation == FM::HARMONIC) && !_is_vector)
  {
    mooseError("The harmonic formulation is selected, but only a scalar potential is provided! Please check input file.");
  }
}

ADReal
ADJouleHeatingSource_OneKernel::precomputeQpResidual()
{
  ADReal value;

  if (_is_vector)
    value = _efield[_qp] * _efield[_qp];
  else
    value = _grad_potential[_qp] * _grad_potential[_qp];
  
  if (_formulation == FM::HARMONIC)
    value = 0.5 * (value + _efield_complex[_qp] * _efield_complex[_qp]);

  return -_elec_cond[_qp] * value;
}
