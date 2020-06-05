//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "VariableOrderBDFPredictor.h"
#include "NonlinearSystem.h"
//#include "SystemBase.h"
#include "FEProblem.h"

#include "libmesh/numeric_vector.h"

registerMooseObject("MooseApp", VariableOrderBDFPredictor);

defineLegacyParams(VariableOrderBDFPredictor);

InputParameters
VariableOrderBDFPredictor::validParams()
{
  InputParameters params = Predictor::validParams();
  params.addParam<int>("order", 0, "If used, sets the time integrator to a constant user "
                       "specified order.");
  return params;
}

VariableOrderBDFPredictor::VariableOrderBDFPredictor(const InputParameters & parameters)
  : Predictor(parameters),

  //_sys(*getCheckedPointerParam<SystemBase *>("_sys")),

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
  _beta(declareRestartableData<std::vector<Real>>("beta")),
  _phi_coeff(declareRestartableData<std::vector<Real>>("theta_coeff")),
  _weight(declareRestartableData<std::vector<Real>>("weight")),

  _solution_old3(_nl.solutionState(3)),
  _solution_old4(_nl.solutionState(4)),
  _solution_old5(_nl.solutionState(5)),
  _solution_old6(_nl.solutionState(6))
{
  _time_array.resize(7);
  _psi.resize(6);
  _psi_old.resize(6);
  _beta.resize(6);
  _phi_coeff.resize(6);
  _weight.resize(6);
}

void
VariableOrderBDFPredictor::timestepSetup()
{
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
    std::fill(_beta.begin(), _beta.end(), 0.0);
    std::fill(_phi_coeff.begin(), _phi_coeff.end(), 0.0);
    std::fill(_weight.begin(), _weight.end(), 0.0);

    for (unsigned int i = 0; i < (_order + 1); ++i)
    {
      _psi[i] = _time_array[0] - _time_array[i+1];
      if (i < 5)
      {
        _psi_old[i] = _time_array[1] - _time_array[i+2];
      }

      if (i==0)
      {
        _beta[i] = 1.0;
        _phi_coeff[i] = 1.0;
      }
      else
      {
        _beta[i] = _psi[i-1] / _psi_old[i-1] * _beta[i-1];
        _phi_coeff[i] = _psi_old[i-1] * _phi_coeff[i-1];
      }

      Real var_coeff = _phi_coeff[i]*_beta[i];
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
          _weight[0] += (var_coeff / DD_coeff);
        }
        if (j == 1)
        {
          _weight[1] += (var_coeff / DD_coeff);
        }
        if (j == 2)
        {
          _weight[2] += (var_coeff / DD_coeff);
        }
        if (j == 3)
        {
          _weight[3] += (var_coeff / DD_coeff);
        }
        if (j == 4)
        {
          _weight[4] += (var_coeff / DD_coeff);
        }
        if (j == 5)
        {
          _weight[5] += (var_coeff / DD_coeff);
        }
      }
    }
  }
}

bool
VariableOrderBDFPredictor::shouldApply()
{
  if (!Predictor::shouldApply())
    return false;

  if (_t_step <= _setup_step)
    return false;
  else
    return true;
}

void
VariableOrderBDFPredictor::apply(NumericVector<Number> & sln)
{
  // localize current solution to working vec
  sln.localize(_solution_predictor);

  _solution_predictor.zero();
  _solution_predictor.add(_weight[0], _solution_old);
  _solution_predictor.add(_weight[1], _solution_older);
  _solution_predictor.add(_weight[2], _solution_old3);
  _solution_predictor.add(_weight[3], _solution_old4);
  _solution_predictor.add(_weight[4], _solution_old5);
  _solution_predictor.add(_weight[5], _solution_old6);

  _solution_predictor.localize(sln);
}
