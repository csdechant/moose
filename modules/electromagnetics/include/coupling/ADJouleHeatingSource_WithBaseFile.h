//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADVectorSupplier.h"

/**
 * This kernel calculates the heat source term corresponding to joule heating,
 * in the three following definitions, depend on EM formulations:
 *
 * 1. If EM formulation is set to STANDARD:
 *    Q = conductivity * field * field, where E is the electric field
 *
 * 2. If EM formulation is set to HARMONIC:
 *    Q = 0.5 * conductivity * field * E^*, where E^* is the complex conjugate of the electric field
 *
 * The field variable is suppled by ADVectorSupplier and is depend by the supplied
 * variable type, as:
 * 
 * 1. If user supplies a scalar variable type:
 *    field = -1.0 * grad_potential, where potential is the suppled scalar variable
 * 
 * 2. If user supplies a vector variable:
 *    field = electric_field, where electric_field is the suppled vector variable
 */
class ADJouleHeatingSource_WithBaseFile : public ADVectorSupplier
{
public:
  static InputParameters validParams();

  ADJouleHeatingSource_WithBaseFile(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

private:
  /// The complex component of the electric field, needed for harmonic formulation
  const ADVectorVariableValue & _efield_complex;
  /// Real component of the material conductivity (in S/m)
  const ADMaterialProperty<Real> & _elec_cond;
  /// Coefficient to multiply by heating term
  const Real & _scale;
  /// The formulation the the EM Joule heating (either STANDARD or HARMONIC)
  MooseEnum _formulation;
};
