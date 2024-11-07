//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernelValue.h"

/**
 * This kernel calculates the heat source term corresponding to joule heating,
 * in the three following definitions, depend on variable type and EM formulations:
 *
 * 1. If user supplies a scalar variable type:
 *    Q = conductivity * grad_potential * grad_potential
 *
 * 2. If user supplies a vector variable type and EM formulation is set to STANDARD:
 *    Q = conductivity * E * E, where E is the electric field
 *
 * 3. If user supplies a vector variable type and EM formulation is set to HARMONIC:
 *    Q = 0.5 * conductivity * E * E^*, where E^* is the complex conjugate of the electric field
 */
class ADJouleHeatingSource_OneKernel : public ADKernelValue
{
public:
  static InputParameters validParams();

  ADJouleHeatingSource_OneKernel(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

private:
  /// The variable data of the supplied variable for the field
  const MooseVariableFieldBase & _field_var;
  /// True is the supplied variable is a vector
  const bool _is_vector;
  /// The field defined from a vector variable
  const ADVectorVariableValue & _efield;
  /// The complex component of the electric field, needed for harmonic formulation
  const ADVectorVariableValue & _efield_complex;
  /// The field defined from the gradient of a scalar variable
  const ADVariableGradient & _grad_potential;
  /// Real component of the material conductivity (in S/m)
  const ADMaterialProperty<Real> & _elec_cond;
  /// Coefficient to multiply by heating term
  const Real & _scale;
  /// The formulation the the EM Joule heating (either STANDARD or HARMONIC)
  MooseEnum _formulation;
};
