//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableOrderBDF.h"
#include "NonlinearSystem.h"
#include "FEProblem.h"

#include "libmesh/nonlinear_solver.h"

registerMooseObject("MooseApp", VariableOrderBDF);

defineLegacyParams(VariableOrderBDF);

InputParameters
VariableOrderBDF::validParams()
{
  InputParameters params = TimeIntegrator::validParams();

  params.addParam<int>("order", 0, "If used, sets the time integrator to a constant user "
                       "specified order.");
  params.addClassDescription("Time Integrator based on the Fixed-Leading Coefficient "
                             "Variable Order Backwards Difference Formulation.");

  return params;
}

VariableOrderBDF::VariableOrderBDF(const InputParameters & parameters)
  : TimeIntegrator(parameters),

    _time(_fe_problem.time()),
    _time_old(_fe_problem.timeOld()),
    _time_older(_fe_problem.timeOlder()),
    _time_old3(_fe_problem.timeOld3()),
    _time_old4(_fe_problem.timeOld4()),
    _time_old5(_fe_problem.timeOld5()),
    _time_old6(_fe_problem.timeOld6()),

    _user_order(getParam<int>("order")),
    _k_order(_fe_problem.kOrder()),
    _order(declareRestartableData<int>("true_order")),
    _setup_step(declareRestartableData<int>("setup_step")),

    _time_array(declareRestartableData<std::vector<Real>>("time_array")),
    _psi(declareRestartableData<std::vector<Real>>("psi")),
    _psi_old(declareRestartableData<std::vector<Real>>("psi_old")),
    _alpha(declareRestartableData<std::vector<Real>>("alpha")),
    _beta(declareRestartableData<std::vector<Real>>("beta")),
    _phi_coeff(declareRestartableData<std::vector<Real>>("theta_coeff")),
    _gamma(declareRestartableData<std::vector<Real>>("gamma")),
    _alpha_s(declareRestartableData<Real>("alpha_s")),
    _weight(declareRestartableData<std::vector<Real>>("weight")),

    _solution_old3(_sys.solutionState(3)),
    _solution_old4(_sys.solutionState(4)),
    _solution_old5(_sys.solutionState(5)),
    _solution_old6(_sys.solutionState(6))
{
  _time_array.resize(7);
  _psi.resize(6);
  _psi_old.resize(6);
  _alpha.resize(6);
  _beta.resize(6);
  _phi_coeff.resize(6);
  _gamma.resize(6);
  _weight.resize(7);
}

void
VariableOrderBDF::preStep()
{
  /**
   * During the preStep, the coefficients are calculated based on formulation
   * found in Ch. 5 of "Numerical Solution of Initial-Vaglue Problems in
   * Differential-Algebral Equations" by K.E. Brenan, S. L. Campbell,
   * and L.R. Petzold. Instead on holding vectors of phi and the predictors,
   * a single coefficient for each older solution is determinded instead
   * (e.g y' = y'^(0)-a[y-y^(0)] = sum[C_i*y_i])
   */

  //Setting time array
  _time_array[6] = _time_old6;
  _time_array[5] = _time_old5;
  _time_array[4] = _time_old4;
  _time_array[3] = _time_old3;
  _time_array[2] = _time_older;
  _time_array[1] = _time_old;
  _time_array[0] = _time;

  /**
   * If the order is not specified, then the order is determined in the
   * Time Stepper.
   */
  if (_user_order > 0)
  {
    _order = _user_order;
    _setup_step = _user_order;
  }
  else
  {
    _order = _k_order;
    _setup_step = 1;
  }

  if (_t_step > _setup_step)
  {
    if (_order == 0)
    {
      mooseError("Order K is still zero after first time step."
                 "Please make sure VariableOrderBDFStepper is being used.");
    }

    std::fill(_psi.begin(), _psi.end(), 0.0);
    std::fill(_psi_old.begin(), _psi_old.end(), 0.0);
    std::fill(_alpha.begin(), _alpha.end(), 0.0);
    std::fill(_beta.begin(), _beta.end(), 0.0);
    std::fill(_phi_coeff.begin(), _phi_coeff.end(), 0.0);
    std::fill(_gamma.begin(), _gamma.end(), 0.0);
    std::fill(_weight.begin(), _weight.end(), 0.0);

    _alpha_s = 0.0;

    for (unsigned int i = 0; i < _order; ++i)
    {
      _alpha_s += (-1.0 / (1.0 + i));
    }

    Real fix_coeff = _alpha_s / _dt;
    _weight[0] = -fix_coeff;

    for (unsigned int i = 0; i < (_order + 1); ++i)
    {
      _psi[i] = _time_array[0] - _time_array[i+1];
      if (i < 5)
      {
        _psi_old[i] = _time_array[1] - _time_array[i+2];
      }
      _alpha[i] = _dt / _psi[i];
      if (i==0)
      {
        _beta[i] = 1.0;
        _phi_coeff[i] = 1.0;
        _gamma[i] = 0.0;
      }
      else
      {
        _beta[i] = _psi[i-1] / _psi_old[i-1] * _beta[i-1];
        _phi_coeff[i] = _psi_old[i-1] * _phi_coeff[i-1];
        _gamma[i] = _gamma[i-1] + (_alpha[i-1] / _dt);
      }

      Real var_coeff = _phi_coeff[i]*_beta[i]*(_gamma[i] + fix_coeff);
      for (unsigned int j = 0; j <= i; ++j)
      {
        Real DD_coeff = 1.0;
        for (unsigned int k = 0; k <= i ; ++k)
        {
          if (k == j)
          {
            DD_coeff = DD_coeff * 1.0;
          }
          else
          {
            DD_coeff = DD_coeff * (_time_array[j+1] - _time_array[k+1]);
          }
        }

        if (j == 0)
        {
          _weight[1] += (var_coeff / DD_coeff);
        }
        if (j == 1)
        {
          _weight[2] += (var_coeff / DD_coeff);
        }
        if (j == 2)
        {
          _weight[3] += (var_coeff / DD_coeff);
        }
        if (j == 3)
        {
          _weight[4] += (var_coeff / DD_coeff);
        }
        if (j == 4)
        {
          _weight[5] += (var_coeff / DD_coeff);
        }
        if (j == 5)
        {
          _weight[6] += (var_coeff / DD_coeff);
        }
      }
    }
  }
}

void
VariableOrderBDF::computeTimeDerivatives()
{
  if (!_sys.solutionUDot())
    mooseError("VariableOrderBDF: Time derivative of solution (`u_dot`) is not stored. Please set "
               "uDotRequested() to true in FEProblemBase befor requesting `u_dot`.");

  NumericVector<Number> & u_dot = *_sys.solutionUDot();

  if (_t_step <= _setup_step)
  {
    u_dot = *_solution;
    _du_dot_du = 1. / _dt;
  }
  else
  {
    u_dot.zero();
    _du_dot_du = -_alpha_s / _dt;
  }
  computeTimeDerivativeHelper(u_dot, *_solution, _solution_old, _solution_older, _solution_old3, _solution_old4, _solution_old5, _solution_old6);
  u_dot.close();
}

void
VariableOrderBDF::computeADTimeDerivatives(DualReal & ad_u_dot, const dof_id_type & dof) const
{
  auto ad_sln = ad_u_dot;
  if (_t_step > _setup_step)
    ad_u_dot != 0;
  computeTimeDerivativeHelper(ad_u_dot, ad_sln, _solution_old(dof), _solution_older(dof), _solution_old3(dof), _solution_old4(dof), _solution_old5(dof), _solution_old6(dof));
}

void
VariableOrderBDF::postResidual(NumericVector<Number> & residual)
{
  residual += _Re_time;
  residual += _Re_non_time;
  residual.close();
}
