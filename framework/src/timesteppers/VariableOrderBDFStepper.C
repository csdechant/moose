//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableOrderBDFStepper.h"
#include "Problem.h"
#include "FEProblem.h"
#include "MooseApp.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"
#include "TimeIntegrator.h"
#include "Conversion.h"
#include "Transient.h"
#include "MathUtils.h"

#include "libmesh/nonlinear_solver.h"
#include "libmesh/numeric_vector.h"

// C++ Includes
#include <iomanip>
#include <iostream>
#include <fstream>

registerMooseObject("MooseApp", VariableOrderBDFStepper);

defineLegacyParams(VariableOrderBDFStepper);

InputParameters
VariableOrderBDFStepper::validParams()
{
  InputParameters params = TimeStepper::validParams();
  params.addRequiredParam<Real>("dt", "Initial time step size");
  params.addParam<int>("starting_step", 1, "The time step to start the increasing the "
                                            "order of the BDF time integrator.");
  params.addParam<int>("max_order", 5, "The max order the BDF time integrator is "
                                        "allow to reach.");
  //Currently the Absolute and Relative are user defined but they should be the
  //'Nonlinear Absolute Tolerance' and 'Nonlinear Relative Tolerance'
  params.addParam<Real>("nl_abs_tol", 1.0e-50, "Absolute Tolerance");
  params.addParam<Real>("nl_rel_tol", 1.0e-8, "Relative Tolerance");
  return params;
}

VariableOrderBDFStepper::VariableOrderBDFStepper(const InputParameters & parameters)
  : TimeStepper(parameters),

    _tmp_vector(_fe_problem.getNonlinearSystemBase().addVector("tmp_vector", true, GHOSTED)),
    _norm_weight(_fe_problem.getNonlinearSystemBase().addVector("norm_weight", true, GHOSTED)),

    _diff(_fe_problem.getNonlinearSystemBase().addVector("diff", true, GHOSTED)),
    _diff_old(_fe_problem.getNonlinearSystemBase().addVector("diff_old", true, GHOSTED)),
    _TERKM1_vector(_fe_problem.getNonlinearSystemBase().addVector("TERKM1_vector", true, GHOSTED)),
    _TERKM2_vector(_fe_problem.getNonlinearSystemBase().addVector("TERKM2_vector", true, GHOSTED)),
    _TERKP1_vector(_fe_problem.getNonlinearSystemBase().addVector("TERKP1_vector", true, GHOSTED)),

    _k_order(_fe_problem.kOrder()),
    _k_old(declareRestartableData<int>("k_old", 0)),
    _k_new(declareRestartableData<int>("k_new", 1)),

    _dt_old(_fe_problem.dtOld()),

    _time_older(_fe_problem.timeOlder()),
    _time_old3(_fe_problem.timeOld3()),
    _time_old4(_fe_problem.timeOld4()),
    _time_old5(_fe_problem.timeOld5()),
    _time_old6(_fe_problem.timeOld6()),

    _time_array(declareRestartableData<std::vector<Real>>("time_array")),
    _psi(declareRestartableData<std::vector<Real>>("psi")),
    _psi_old(declareRestartableData<std::vector<Real>>("psi_old")),
    _alpha(declareRestartableData<std::vector<Real>>("alpha")),
    _beta(declareRestartableData<std::vector<Real>>("beta")),
    _phi_coeff(declareRestartableData<std::vector<Real>>("theta_coeff")),
    _sigma(declareRestartableData<std::vector<Real>>("sigma")),
    _alpha_s(declareRestartableData<Real>("alpha_s")),
    _alpha_0(declareRestartableData<Real>("alpha_0")),

    _errk(declareRestartableData<Real>("errk", 0.0)),
    _terk(declareRestartableData<Real>("terk", 0.0)),
    _errkm1(declareRestartableData<Real>("errkm1", 0.0)),
    _terkm1(declareRestartableData<Real>("terkm1", 0.0)),
    _errkm2(declareRestartableData<Real>("errkm2", 0.0)),
    _terkm2(declareRestartableData<Real>("terkm2", 0.0)),
    _errkp1(declareRestartableData<Real>("errkp1", 0.0)),
    _terkp1(declareRestartableData<Real>("terkp1", 0.0)),
    _error(declareRestartableData<Real>("error", 0.0)),

    _err_counter(declareRestartableData<int>("err_counter", 0)),

    _is_beginning_settup(declareRestartableData<bool>("is_beginning_settup", true)),

    _abs_tol(getParam<Real>("nl_abs_tol")),
    _rel_tol(getParam<Real>("nl_rel_tol")),
    _stepup_step(getParam<int>("starting_step")),
    _k_max(getParam<int>("max_order"))

{
  _time_array.resize(7);
  _psi.resize(6);
  _psi_old.resize(6);
  _alpha.resize(6);
  _beta.resize(6);
  _phi_coeff.resize(6);
  _sigma.resize(6);
}

void
VariableOrderBDFStepper::preSolve()
{
  /**
   * During the preSolve, the coefficients are calculated. Some of these
   * terms are already calculated in the Time Integrator. I should found a
   * way for the Time Integrator and Stepper to share values.
   */

  _time_array[6] = _time_old6;
  _time_array[5] = _time_old5;
  _time_array[4] = _time_old4;
  _time_array[3] = _time_old3;
  _time_array[2] = _time_older;
  _time_array[1] = _time_old;
  _time_array[0] = _time_old + _dt;

  std::fill(_psi.begin(), _psi.end(), 0.0);
  std::fill(_psi_old.begin(), _psi_old.end(), 0.0);
  std::fill(_alpha.begin(), _alpha.end(), 0.0);
  std::fill(_beta.begin(), _beta.end(), 0.0);
  std::fill(_phi_coeff.begin(), _phi_coeff.end(), 0.0);
  std::fill(_sigma.begin(), _sigma.end(), 0.0);

  _alpha_s = 0.0;
  _alpha_0 = 0.0;

  for (unsigned int i = 0; i < _k_order; ++i)
  {
    _alpha_s += (-1.0 / (1.0 + i));
  }

  for (unsigned int i = 0; i < (_k_order + 1); ++i)
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
      _sigma[i] = 1.0;
    }
    else
    {
      _beta[i] = _psi[i-1] / _psi_old[i-1] * _beta[i-1];
      _phi_coeff[i] = _psi_old[i-1] * _phi_coeff[i-1];
      _sigma[i] = _sigma[i-1] * _dt * 1.0 * i / _psi[i];
    }
  }

  for (unsigned int i = 0; i < _k_order; ++i)
  {
    _alpha_0 += (-1.0 * _alpha[i]);
  }
}

void
VariableOrderBDFStepper::step()
{
  if(!_at_sync_point)
  {
    /**
     * Determining for the the startup phase of the simulation is over.
     * During the startup phase, the order increase by 1 and the dt is
     * doubled after each step. This is true until the order is maxed,
     * the order is lowered, or until converged == false.
     */
    if((_is_beginning_settup) && (_k_order == _k_max))
    {
      _is_beginning_settup = false;
    }

    NonlinearSystemBase & nl = _fe_problem.getNonlinearSystemBase();

    _fe_problem.solve();
    _converged = _fe_problem.converged();

    if ((_converged) && (_t_step <= _stepup_step))
    {
      _k_new = 1;

      _k_old = _k_order;
      _k_order = _k_new;
    }

    if ((_converged) && (_t_step > _stepup_step))
    {
      /**
       * Calling solution and older solutions
       */
      const NumericVector<Number> & _u0 = *nl.currentSolution();
      NumericVector<Number> & _u1 = nl.solutionOld();
      _u1.close();
      NumericVector<Number> & _u2 = nl.solutionOlder();
      _u2.close();
      NumericVector<Number> & _u3 = nl.solutionState(3);
      _u3.close();
      NumericVector<Number> & _u4 = nl.solutionState(4);
      _u4.close();
      NumericVector<Number> & _u5 = nl.solutionState(5);
      _u5.close();
      NumericVector<Number> & _u6 = nl.solutionState(6);
      _u6.close();

      /**
       * Forming the weight for the weighted norm calculations.
       */
      _norm_weight.zero();
      _norm_weight.add(_rel_tol, _u1);
      _norm_weight.abs();
      _norm_weight.add(_abs_tol);
      _norm_weight.close();

      _diff.zero();
      _diff = DifferenceCalculation(_u0,_u1,_u2,_u3,_u4,_u5,_u6);
      Real _diff_norm = WeightedNorm(_diff);
      //Real _diff_norm = _diff.l2_norm();

      _errk = _sigma[_k_order] * _diff_norm;
      _terk = (1.*_k_order + 1.) * _errk;

      /**
       * Determining for the order should be lowered before error calculation.
       */
      if (_k_order == 2)
      {
        _TERKM1_vector.zero();
        _TERKM1_vector = TERKM1VectorCalculation(_diff,_u1,_u2,_u3,_u4,_u5,_u6);
        Real _TERKM1_norm = WeightedNorm(_TERKM1_vector);

        _errkm1 = _sigma[_k_order-1] * _TERKM1_norm;
        _terkm1 = (1.*_k_order) * _errkm1;

        if (_terkm1 <= _terk/2.0)
        {
          _k_new = _k_order - 1;
          _is_beginning_settup = false;
        }
      }
      else if (_k_order > 2)
      {
        _TERKM1_vector.zero();
        _TERKM1_vector = TERKM1VectorCalculation(_diff,_u1,_u2,_u3,_u4,_u5,_u6);
        Real _TERKM1_norm = WeightedNorm(_TERKM1_vector);

        _TERKM2_vector.zero();
        _TERKM2_vector = TERKM2VectorCalculation(_TERKM1_vector,_u1,_u2,_u3,_u4,_u5,_u6);
        Real _TERKM2_norm = WeightedNorm(_TERKM2_vector);

        _errkm1 = _sigma[_k_order-1] * _TERKM1_norm;
        _terkm1 = (1.*_k_order) * _errkm1;

        _errkm2 = _sigma[_k_order-2] * _TERKM2_norm;
        _terkm2 = (1.*_k_order - 1.) * _errkm2;

        if (std::max(_terkm1,_terkm2) <= _terk)
        {
          _k_new = _k_order - 1;
          _is_beginning_settup = false;
        }
      }

      /**
       * Calculating error that determines convergence (converged for error < 1)
       */
      Real _m = std::max(_alpha[_k_order],std::abs(_alpha[_k_order]+_alpha_s-_alpha_0));
      _error = _m * _diff_norm;

      if (_error <= 1.0)
      {
        /**
         * During setup, the order is increase by 1 after each step after the first
         * time step.
         */
        if (_is_beginning_settup)
        {
          _k_new = _k_order + 1;
        }
        /**
         * Determing to raise the order if converged, last order was constant,
         * order is less than 5, last dt was constant, and the order was not decided
         * to be lowered the the steps above.
         */
        else if ((_k_new == _k_order) && (_k_order < _k_max) &&
                 (_dt == _dt_old) && (_k_order == _k_old))
        {
          _TERKP1_vector.zero();
          _TERKP1_vector.add(1.0, _diff);
          _TERKP1_vector.add(-1.0, _diff_old);
          _TERKP1_vector.close();

          Real _TERKP1_norm = WeightedNorm(_TERKP1_vector);

          _errkp1 = _TERKP1_norm / (1.*_k_order + 1.);
          _terkp1 = (1.*_k_order + 1.) * _errkp1;

          if ((_k_order == 1) && (_terkp1 <= _terk/2.0) )
          {
            _k_new = _k_order + 1;
          }
          if (_k_order > 1)
          {
            if (_terkm1 <= std::min(_terk,_terkp1))
            {
              _k_new = _k_order - 1;
            }
            else if ((_k_new == _k_order) && (_terkp1 < _terk))
            {
              _k_new = _k_order + 1;
            }
          }
        }

        /**
         * Saving the difference between the conveged solution and predictor, in
         * case the order needs to be rised next time step.
         */
        _diff_old = _diff;
      }
      else
      {
        _err_counter = _err_counter + 1;

        if (_err_counter > 2)
        {
          _k_new = 1;
        }
      }

      /**
       * Defining next order and saving current order.
       */
      _k_old = _k_order;
      _k_order = _k_new;
    }
  }
  else
  {
    _converged = _executioner.picardSolve().solve();

    _diff.zero();
    _TERKM1_vector.zero();
    _TERKM2_vector.zero();
  }
}

bool
VariableOrderBDFStepper::converged() const
{
  if (!_converged)
    return false;
  if (_error <= 1.0)
    return true;
  else
    return false;
}

void
VariableOrderBDFStepper::postSolve()
{
  if(!_at_sync_point)
    {
      if (converged())
      {
        _console << "Marking last solve converged with BDF order of " << _k_order << std::endl;
        _console << "TERM of " << _errk << std::endl;
        _console << "TERM-1 of " << _errkm1 << std::endl;
        _console << "TERM-2 of " << _errkm2 << std::endl;
        _console << "TERM+1 of " << _errkp1 << std::endl;
      }
      else if (_error > 1.0 && !converged())
      {
        _console << "Marking last solve not converged with local truncation error of " << _error << std::endl;
        _console << "BDF order of " << _k_order << std::endl;
        _console << "TERM of " << _errk << std::endl;
        _console << "TERM-1 of " << _errkm1 << std::endl;
        _console << "TERM-2 of " << _errkm2 << std::endl;
        _console << "TERM+1 of " << _errkp1 << std::endl;
        //for (unsigned int i = 0; i < _diff.size(); ++i)
        //{
        //  _console << "Diff " << i << " of " << _diff(i) << std::endl;
        //}
      }
    }
}

Real
VariableOrderBDFStepper::computeDT()
{
  if(!_at_sync_point)
  {
    /**
     * Resetting error counter.
     */
    _err_counter = 0;

    if (_t_step <= _stepup_step)
    {
      Real _int_dt = computeInitialDT();
      if (_dt < _int_dt)
      {
        return std::min(_int_dt, 2.0 * _dt);
      }
      else
      {
        return _int_dt;
      }
    }
    else
    {
      if ((_is_beginning_settup))
      {
        /**
         * During setup, the time step is doubled after each step.
         */
         return std::min(_dt_max, 2.0 * _dt);
      }
      else
      {
        /**
         * The new dt is determined by the local error of the new order.
         * New dt can be either double, reduced between 0.5 - 0.9 of the current dt
         * or remain the same.
         */
        Real _est;
        if(_k_order == (_k_old - 1)) {
          _est = _errkm1; }
        if(_k_order == (_k_old + 1)) {
          _est = _errkp1; }
        if(_k_order == _k_old) {
          _est = _errk; }

        Real _tau = std::pow(2.0 * _est + 1e-4, (-1./(_k_order+1.)));
        if(_tau >= 2.)
        {
          return std::min(_dt_max, 2.0 * _dt);
        }
        if(_tau <= 1.)
        {
          return std::max(0.5,std::min(0.9,_tau)) * _dt;
        }
        else
        {
          return _dt;
        }
      }
    }
  }
  else
  {
    if (_k_order == 0)
    {
      _k_order = 1;
    }

    return _dt;
  }
}

Real
VariableOrderBDFStepper::computeFailedDT()
{
  if (_dt <= _dt_min)
  {
    mooseError("Solve failed and timestep already at or below dtmin, cannot continue!");
  }

  if (_t_step <= _stepup_step)
  {
    if (1./4. * _dt >= _dt_min)
    {
      return 1./4. * _dt;
    }
    else
    {
      return _dt_min;
    }
  }
  else
  {
    /**
     * The setup is over due to not convergence.
     */
    _is_beginning_settup = false;

    /**
     * Failure due to the error calculation.
     * For the first failure, the dt is lower between 0.25 - 0.9 of
     * the current dt based of the local error.
     * For the second failure, the dt is lower by a fourth of the current dt.
     * For later failures, the order is lowered to 1 and the dt is lower
     * by a fourth of the current dt.
     */
    if(_error > 1.0)
    {
      if(_err_counter == 1)
      {
        Real _est;
        if(_k_order < _k_old)
        {
          _est = _errkm1;
        }
        else
        {
          _est = _errk;
        }

        Real _tau = std::pow(2.0 * _est + 1e-4, (-1./(_k_order+1.)));
        return std::max(0.25,std::min(0.9,0.9*_tau)) * _dt;
      }
      else
      {
        return 1./4. * _dt;
      }
    }

    /**
     * Failure due to the nonlinear solver. The dt is lowered by a fourth.
     */
    else
    {
      if (1./4. * _dt >= _dt_min)
      {
        return 1./4. * _dt;
      }
      else
      {
        return _dt_min;
      }
    }
  }
}

Real
VariableOrderBDFStepper::computeInitialDT()
{
  return getParam<Real>("dt");
}

/**
 * Computes the l2 weighted norm
 */
Real
VariableOrderBDFStepper::WeightedNorm(NumericVector<Number> & vector)
{
  Real norm = 0.0;
  Real dof = _norm_weight.size();
  for (unsigned int i = 0; i < dof; ++i)
  {
    norm += std::pow(vector(i) / _norm_weight(i), 2.0);
  }
  norm = std::sqrt(norm / dof);


  return norm;
}

/**
 * Computes the difference between the solution and the predictor
 */
NumericVector<Number> &
VariableOrderBDFStepper::DifferenceCalculation(const NumericVector<Number> & _u0,
                                                     NumericVector<Number> & _u1,
                                                     NumericVector<Number> & _u2,
                                                     NumericVector<Number> & _u3,
                                                     NumericVector<Number> & _u4,
                                                     NumericVector<Number> & _u5,
                                                     NumericVector<Number> & _u6)
{
  std::vector<Real> weight;
  weight.resize(7);

  weight[0] = 1.0;
  for (unsigned int i = 0; i < (_k_order + 1); ++i)
  {
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

      Real coeff = _beta[i] * _phi_coeff[i];
      if (j == 0)
      {
        weight[1] -= (coeff / DD_coeff);
      }
      if (j == 1)
      {
        weight[2] -= (coeff / DD_coeff);
      }
      if (j == 2)
      {
        weight[3] -= (coeff / DD_coeff);
      }
      if (j == 3)
      {
        weight[4] -= (coeff / DD_coeff);
      }
      if (j == 4)
      {
        weight[5] -= (coeff / DD_coeff);
      }
      if (j == 5)
      {
        weight[6] -= (coeff / DD_coeff);
      }
    }
  }

  NumericVector<Number> & diff = _tmp_vector;
  diff.zero();
  MathUtils::addScaled(weight[0], _u0, diff);
  MathUtils::addScaled(weight[1], _u1, diff);
  MathUtils::addScaled(weight[2], _u2, diff);
  MathUtils::addScaled(weight[3], _u3, diff);
  MathUtils::addScaled(weight[4], _u4, diff);
  MathUtils::addScaled(weight[5], _u5, diff);
  MathUtils::addScaled(weight[6], _u6, diff);
  diff.close();

  return diff;
}

/**
 * Computes the vector for the order-1 local error
 */
NumericVector<Number> &
VariableOrderBDFStepper::TERKM1VectorCalculation(NumericVector<Number> & diff,
                                                 NumericVector<Number> & _u1,
                                                 NumericVector<Number> & _u2,
                                                 NumericVector<Number> & _u3,
                                                 NumericVector<Number> & _u4,
                                                 NumericVector<Number> & _u5,
                                                 NumericVector<Number> & _u6)
{
  std::vector<Real> weight;
  weight.resize(7);

  weight[0] = 1.0;
  for (unsigned int j = 0; j <= (_k_order); ++j)
  {
    Real DD_coeff = 1.0;
    for (unsigned int k = 0; k <= (_k_order) ; ++k)
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

    Real coeff = _beta[_k_order] * _phi_coeff[_k_order];
    if (j == 0)
    {
      weight[1] += (coeff / DD_coeff);
    }
    if (j == 1)
    {
      weight[2] += (coeff / DD_coeff);
    }
    if (j == 2)
    {
      weight[3] += (coeff / DD_coeff);
    }
    if (j == 3)
    {
      weight[4] += (coeff / DD_coeff);
    }
    if (j == 4)
    {
      weight[5] += (coeff / DD_coeff);
    }
    if (j == 5)
    {
      weight[6] += (coeff / DD_coeff);
    }
  }

  NumericVector<Number> & TERKM1_vector = _tmp_vector;
  TERKM1_vector.zero();
  MathUtils::addScaled(weight[0], diff, TERKM1_vector);
  MathUtils::addScaled(weight[1], _u1, TERKM1_vector);
  MathUtils::addScaled(weight[2], _u2, TERKM1_vector);
  MathUtils::addScaled(weight[3], _u3, TERKM1_vector);
  MathUtils::addScaled(weight[4], _u4, TERKM1_vector);
  MathUtils::addScaled(weight[5], _u5, TERKM1_vector);
  MathUtils::addScaled(weight[6], _u6, TERKM1_vector);
  TERKM1_vector.close();

  return TERKM1_vector;
}

/**
 * Computes the vector for the order-2 local error
 */
NumericVector<Number> &
VariableOrderBDFStepper::TERKM2VectorCalculation(NumericVector<Number> & diff,
                                                 NumericVector<Number> & _u1,
                                                 NumericVector<Number> & _u2,
                                                 NumericVector<Number> & _u3,
                                                 NumericVector<Number> & _u4,
                                                 NumericVector<Number> & _u5,
                                                 NumericVector<Number> & _u6)
{
  std::vector<Real> weight;
  weight.resize(7);

  weight[0] = 1.0;
  for (unsigned int j = 0; j <= (_k_order-1); ++j)
  {
    Real DD_coeff = 1.0;
    for (unsigned int k = 0; k <= (_k_order-1) ; ++k)
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

    Real coeff = _beta[_k_order-1] * _phi_coeff[_k_order-1];
    if (j == 0)
    {
      weight[1] += (coeff / DD_coeff);
    }
    if (j == 1)
    {
      weight[2] += (coeff / DD_coeff);
    }
    if (j == 2)
    {
      weight[3] += (coeff / DD_coeff);
    }
    if (j == 3)
    {
      weight[4] += (coeff / DD_coeff);
    }
    if (j == 4)
    {
      weight[5] += (coeff / DD_coeff);
    }
    if (j == 5)
    {
      weight[6] += (coeff / DD_coeff);
    }
  }

  NumericVector<Number> & TERKM2_vector = _tmp_vector;
  TERKM2_vector.zero();
  MathUtils::addScaled(weight[0], diff, TERKM2_vector);
  MathUtils::addScaled(weight[1], _u1, TERKM2_vector);
  MathUtils::addScaled(weight[2], _u2, TERKM2_vector);
  MathUtils::addScaled(weight[3], _u3, TERKM2_vector);
  MathUtils::addScaled(weight[4], _u4, TERKM2_vector);
  MathUtils::addScaled(weight[5], _u5, TERKM2_vector);
  MathUtils::addScaled(weight[6], _u6, TERKM2_vector);
  TERKM2_vector.close();

  return TERKM2_vector;
}
