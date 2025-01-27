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

class SolutionUserObject;

/**
 * Class for reading an initial condition from a solution user object
 */
class VectorSolutionIC : public VectorInitialCondition
{
public:
  VectorSolutionIC(const InputParameters & parameters);

  virtual void initialSetup() override;

  virtual RealVectorValue value(const Point & p) override;

protected:
  /// SolutionUserObject containing the solution of interest
  const SolutionUserObject & _solution_object;

  /// The variable name extracted from the SolutionUserObject
  const VariableName & _solution_object_var_name;
  //const VariableName & _solution_object_var_x_name;
  //const VariableName & _solution_object_var_y_name;
  //std::string _solution_object_var_x_name;
  //std::string _solution_object_var_y_name;
  //std::string _solution_object_var_z_name;

  /// Remapped IDs from the current mesh to the ExodusII mesh
  std::set<SubdomainID> _exo_block_ids;

public:
  static InputParameters validParams();
};
