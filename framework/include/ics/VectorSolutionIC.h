//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "VectorInitialCondition.h"

class SolutionUserObjectBase;

/**
 * Class for reading an initial condition for vector field variable
 * from a solution user object
 */
class VectorSolutionIC : public VectorInitialCondition
{
public:
  VectorSolutionIC(const InputParameters & parameters);

  virtual void initialSetup() override;

  virtual RealVectorValue value(const Point & p) override;

protected:
  /// SolutionUserObject containing the solution of interest
  const SolutionUserObjectBase & _solution_object;

  /// The variable name extracted from the SolutionUserObject
  const VariableName & _solution_object_var_name;

public:
  static InputParameters validParams();
};
