//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VectorSolutionIC.h"
#include "SolutionUserObjectBase.h"
#include "MooseMesh.h"
#include "SystemBase.h"

registerMooseObject("MooseApp", VectorSolutionIC);

InputParameters
VectorSolutionIC::validParams()
{
  InputParameters params = VectorInitialCondition::validParams();
  params.addRequiredParam<UserObjectName>("solution_uo",
                                          "The SolutionUserObject to extract data from.");
  params.addRequiredParam<VariableName>(
      "from_variable", "The name of the variable in the file that is to be extracted");
  params.addClassDescription("Sets the initial condition from a vector field variable "
                             "retrieved by a SolutionUserObject");
  return params;
}

VectorSolutionIC::VectorSolutionIC(const InputParameters & parameters)
  : VectorInitialCondition(parameters),
    _solution_object(getUserObject<SolutionUserObjectBase>("solution_uo")),
    _solution_object_var_name(getParam<VariableName>("from_variable"))
{
}

void
VectorSolutionIC::initialSetup()
{
  if (_solution_object.getSolutionFileType() == "exodusII")
  {
    // NOTE: Maybe this should just be an error, since combining Lagrange scalar field components
    //       for any vector field data set is the incorrect way to go...
    mooseWarning(
        "You are reading from an Exodus file for vector field data! This is ill-advised because "
        "Exodus files supply vector field data as a set of Lagrange scalar field components, "
        "regardless of the family type of the vector field. Please use XDA output data instead.");
  }
}

RealVectorValue
VectorSolutionIC::value(const Point & p)
{
  if (_solution_object.getSolutionFileType() == "exodusII")
  {
    RealVectorValue output(0., 0., 0.);

    output(0) = _solution_object.pointValue(0., p, _solution_object_var_name + "_x");
    if (_dim > 1)
      output(1) = _solution_object.pointValue(0., p, _solution_object_var_name + "_y");
    if (_dim > 2)
      output(2) = _solution_object.pointValue(0., p, _solution_object_var_name + "_z");

    return output;
  }

  return _solution_object.pointVectorValue(0., p, _solution_object_var_name);
}
