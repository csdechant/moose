//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "TimeIntegrator.h"
#include "MathUtils.h"

class VariableOrderBDF;

template <>
InputParameters validParams<VariableOrderBDF>();

/**
 * VariableOrderBDF time integrator
 */
class VariableOrderBDF : public TimeIntegrator
{
public:
  static InputParameters validParams();

  VariableOrderBDF(const InputParameters & parameters);

  virtual int order() override { return _order; }
  virtual void preStep() override;
  virtual void computeTimeDerivatives() override;
  void computeADTimeDerivatives(DualReal & ad_u_dot, const dof_id_type & dof) const override;
  virtual void postResidual(NumericVector<Number> & residual) override;

protected:
  /**
   * Helper function that actually does the math for computing the time derivative
   */
  template <typename T, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
  void
  computeTimeDerivativeHelper(T & u_dot,
                              const T2 & u,
                              const T3 & u_old,
                              const T4 & u_older,
                              const T5 & u_old3,
                              const T6 & u_old4,
                              const T7 & u_old5,
                              const T8 & u_old6) const;

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
  std::vector<Real> & _alpha;
  std::vector<Real> & _beta;
  std::vector<Real> & _phi_coeff;
  std::vector<Real> & _gamma;
  Real & _alpha_s;

  //Array for the coeff. in the computeTimeDerivativeHelper
  std::vector<Real> & _weight;

  //Calling older solutions.
  const NumericVector<Number> & _solution_old3;
  const NumericVector<Number> & _solution_old4;
  const NumericVector<Number> & _solution_old5;
  const NumericVector<Number> & _solution_old6;

};

namespace VariableOrderBDFHelper
{
}

template <typename T, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
void
VariableOrderBDF::computeTimeDerivativeHelper(T & u_dot,
                                              const T2 & u,
                                              const T3 & u_old,
                                              const T4 & u_older,
                                              const T5 & u_old3,
                                              const T6 & u_old4,
                                              const T7 & u_old5,
                                              const T8 & u_old6) const
{
  if (_t_step <= _setup_step)
  {
    u_dot -= u_old;
    u_dot *= 1 / _dt;
  }
  else
  {
    MathUtils::addScaled(_weight[0], u, u_dot);
    MathUtils::addScaled(_weight[1], u_old, u_dot);
    MathUtils::addScaled(_weight[2], u_older, u_dot);
    MathUtils::addScaled(_weight[3], u_old3, u_dot);
    MathUtils::addScaled(_weight[4], u_old4, u_dot);
    MathUtils::addScaled(_weight[5], u_old5, u_dot);
    MathUtils::addScaled(_weight[6], u_old6, u_dot);
  }
}
