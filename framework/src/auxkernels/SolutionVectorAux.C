//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "MooseError.h"
#include "SolutionVectorAux.h"
#include "SolutionUserObjectBase.h"

registerMooseObject("MooseApp", SolutionVectorAux);

InputParameters
SolutionVectorAux::validParams()
{
  InputParameters params = VectorAuxKernel::validParams();
  params.addClassDescription(
      "Creates vector fields by using information from a SolutionUserObject.");
  params.addRequiredParam<UserObjectName>("solution", "The name of the SolutionUserObject");
  params.addParam<std::string>("from_variable",
                               "The name of the variable to extract from the file");

  params.addParam<bool>(
      "direct",
      false,
      "If true the meshes must be the same and then the values are simply copied over.");
  params.addParam<Real>("scale_factor",
                        1.0,
                        "Scale factor (s)  to be applied to the solution (x,y,z): "
                        "s(x,y,z)+(a,b,c), where b is the 'add_factor'");
  params.addParam<RealVectorValue>("add_factor",
                                   RealVectorValue(0, 0, 0),
                                   "Add this vector value (a,b,c) to the solution (x,y,z): "
                                   "s(x,y,z)+(a,b,c), where a is the 'scale_factor'");
  return params;
}

SolutionVectorAux::SolutionVectorAux(const InputParameters & parameters)
  : VectorAuxKernel(parameters),
    _solution_object(getUserObject<SolutionUserObjectBase>("solution")),
    _direct(getParam<bool>("direct")),
    _scale_factor(getParam<Real>("scale_factor")),
    _add_factor(getParam<RealVectorValue>("add_factor"))
{
}

void
SolutionVectorAux::initialSetup()
{
  // If 'from_variable' is supplied, use the value
  if (isParamValid("from_variable"))
    _var_name = getParam<std::string>("from_variable");

  // If not, get the value from the SolutionUserObject
  else
  {
    // Get all the variables from the SolutionUserObject
    const std::vector<std::string> & vars = _solution_object.variableNames();

    // If there are more than one, throw an error
    if (vars.size() > 1)
      mooseError("The SolutionUserObject contains multiple variables, please specify the desired "
                 "variables in the input file with 'from_variable' parameter.");

    // Define the variable
    _var_name = vars[0];
  }
}

RealVectorValue
SolutionVectorAux::computeValue()
{
  // The value to output
  RealVectorValue output;

  // _direct=true, extract the values using the dof
  if (_direct)
  {
    if (isNodal())
      output = _solution_object.directVectorValue(_current_node, _var_name);

    else
      output = _solution_object.directVectorValue(_current_elem, _var_name);
  }

  // _direct=false, extract the values using time and point
  else
  {
    /*
    if (isNodal())
      output = _solution_object.pointValue(_t, *_current_node, _var_name);

    else
      output = _solution_object.pointValue(_t, _current_elem->vertex_average(), _var_name);
    */
    mooseError("Currently  'pointValue' in SolutionUserObject is broken for vector field variable "
               "due to libMesh's MeshFuction operator. This should be fixed soon.");
  }

  // Apply factors and return the value
  return _scale_factor * output + _add_factor;
}
