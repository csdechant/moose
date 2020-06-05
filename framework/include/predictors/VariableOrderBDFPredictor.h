//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "Predictor.h"

// Forward declarations
class VariableOrderBDFPredictor;
//class SystemBase;

namespace libMesh
{
template <typename T>
class NumericVector;
}

template <>
InputParameters validParams<VariableOrderBDFPredictor>();

/**
 * Implements an explicit Adams predictor based on two old solution
 * vectors.
 */
class VariableOrderBDFPredictor : public Predictor
{
public:
  static InputParameters validParams();

  VariableOrderBDFPredictor(const InputParameters & parameters);

  virtual int order() override { return _order; }
  virtual void timestepSetup() override;
  virtual bool shouldApply() override;
  virtual void apply(NumericVector<Number> & sln) override;
  virtual NumericVector<Number> & solutionPredictor() override { return _solution_predictor; }

protected:

  //SystemBase & _sys;

  Real & _time;
  Real & _time_old;

  //New older time variables. They are currently calculated in Transients and
  //used by both BDF Time Integrator, Stepper, and Predictor to calculate coefficients
  Real & _time_older;
  Real & _time_old3;
  Real & _time_old4;
  Real & _time_old5;
  Real & _time_old6;

  int _user_order;
  //New order variable. It is initiated in FEProblem and calculated in the
  // BDF Time Stepper.
  int & _k_order;

  int & _order;
  int & _setup_step;

  //An array for the time variables
  std::vector<Real> & _time_array;

  //Variable for the Fix Leading Coeff. Variable Order BDF
  std::vector<Real> & _psi;
  std::vector<Real> & _psi_old;
  std::vector<Real> & _beta;
  std::vector<Real> & _phi_coeff;

  //Array for the coeff. in the computeTimeDerivativeHelper
  std::vector<Real> & _weight;

  //Calling older solutions.
  const NumericVector<Number> & _solution_old3;
  const NumericVector<Number> & _solution_old4;
  const NumericVector<Number> & _solution_old5;
  const NumericVector<Number> & _solution_old6;
};
