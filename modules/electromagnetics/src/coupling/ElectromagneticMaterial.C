//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ElectromagneticMaterial.h"
#include "FormulationEnums.h"

registerMooseObject("ElectromagneticsApp", ElectromagneticMaterial);

InputParameters
ElectromagneticMaterial::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription(
      "Material class used to provide the electric field as a material property and compute some "
      "residuals for the electromagnetic equations that differ between standard and harmonic "
      "formulation.");
  params.addCoupledVar(
      "electric_field_or_potential",
      "The electic field vector or electrostatic potential scalar to produce the field.");
  params.addCoupledVar(
      "complex_field",
      "The complex component of the electic field vector for the harmonic formulation.");
  params.addParam<std::string>(
      "field_name", "field", "User-specified material property name for the field.");
  params.addParam<std::string>("field_heating_name",
                               "field_heating",
                               "User-specified material property name for the Joule Heating.");
  params.addParam<Real>("heating_scaling", 1.0, "Coefficient to multiply by heating term.");
  params.addParam<MaterialPropertyName>(
      "conductivity",
      "conductivity",
      "Material property providing electrical conductivity of the material.");
  MooseEnum formulation("standard harmonic", "standard");
  MooseEnum solver("electrostatic electromagnetic", "electrostatic");
  params.addParam<MooseEnum>(
      "formulation", formulation, "The formulation of the Joule heating (standard or harmonic).");
  params.addParam<MooseEnum>(
      "solver", solver, "Electrostatic or electromagnetic field solver (default = electrostatic).");
  return params;
}

ElectromagneticMaterial::ElectromagneticMaterial(const InputParameters & parameters)
  : ADMaterial(parameters),
    _field_var(*getFieldVar("electric_field_or_potential", 0)),
    _is_vector(_field_var.isVector()),
    _efield(_is_vector ? adCoupledVectorValue("electric_field_or_potential") : _ad_grad_zero),
    _efield_complex(_is_vector ? adCoupledVectorValue("complex_field") : _ad_grad_zero),
    _grad_potential(_is_vector ? _ad_grad_zero : adCoupledGradient("electric_field_or_potential")),
    _field(declareADProperty<RealVectorValue>(getParam<std::string>("field_name"))),
    _field_complex(
        declareADProperty<RealVectorValue>(getParam<std::string>("field_name") + "_complex")),
    _field_heating(declareADProperty<Real>(getParam<std::string>("field_heating_name"))),
    _heating_scaling(getParam<Real>("heating_scaling")),
    _elec_cond(getADMaterialProperty<Real>("conductivity")),
    _formulation(getParam<MooseEnum>("formulation")),
    _solver(getParam<MooseEnum>("solver"))
{
  if ((_formulation == FM::HARMONIC) && (_solver == FM::ELECTROSTATIC))
  {
    mooseError("The harmonic formulation is selected, but the solver type is electrostatic! Please "
               "check input file.");
  }

  if ((_solver == FM::ELECTROMAGNETIC) && !_is_vector)
  {
    mooseError("The solver type is electromagnetic, but only a scalar potential is provided! "
               "Please check input file.");
  }

  if ((_formulation == FM::HARMONIC) && !_is_vector)
  {
    mooseError("The harmonic formulation is selected, but only a scalar potential is provided! "
               "Please check input file.");
  }
}

void
ElectromagneticMaterial::computeQpProperties()
{
  computeFieldValue();
  computeJouleHeating();
}

void
ElectromagneticMaterial::computeFieldValue()
{
  if (_solver == FM::ELECTROSTATIC)
    _field[_qp] = -_grad_potential[_qp];
  else
    _field[_qp] = _efield[_qp];

  if (_formulation == FM::HARMONIC)
    _field_complex[_qp] = _efield_complex[_qp];
}

void
ElectromagneticMaterial::computeJouleHeating()
{
  if (_formulation == FM::HARMONIC)
    _field_heating[_qp] = _heating_scaling * 0.5 * _elec_cond[_qp] *
                          (_field[_qp] * _field[_qp] + _field_complex[_qp] * _field_complex[_qp]);
  else
    _field_heating[_qp] = _heating_scaling * _elec_cond[_qp] * _field[_qp] * _field[_qp];
    
}