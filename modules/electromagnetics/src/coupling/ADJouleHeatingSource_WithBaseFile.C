//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html


#include "ADJouleHeatingSource_WithBaseFile.h"
#include "FormulationEnums.h"

registerMooseObject("ElectromagneticsApp", ADJouleHeatingSource_WithBaseFile);

InputParameters
ADJouleHeatingSource_WithBaseFile::validParams()
{
  InputParameters params = ADVectorSupplier::validParams();
  params.addCoupledVar("elec", "Electrostatic potential for joule heating.");
  params.deprecateParam("elec", "field", "12/31/2024");
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

ADJouleHeatingSource_WithBaseFile::ADJouleHeatingSource_WithBaseFile(const InputParameters & parameters)
  : ADVectorSupplier(parameters),
    _efield_complex(_is_vector ? adCoupledVectorValue("complex_field") : _ad_grad_zero),
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
ADJouleHeatingSource_WithBaseFile::computeQpResidual()
{
  ADRealVectorValue field = computeQpFieldValue(-1.0);

  ADReal value = field * field;

  if (_formulation == FM::HARMONIC)
    value = 0.5 * (value + _efield_complex[_qp] * _efield_complex[_qp]);

  return -_test[_i][_qp] * _elec_cond[_qp] * value;
}

