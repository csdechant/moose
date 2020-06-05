//* This file is part of Zapdos, an open-source
//* application for the simulation of plasmas
//* https://github.com/shannon-lab/zapdos
//*
//* Zapdos is powered by the MOOSE Framework
//* https://www.mooseframework.org
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

// MOOSE includes
#include "TimeStepper.h"

// C++ includes
#include <fstream>

// Forward Declarations
class VariableOrderBDFStepper;

namespace libMesh
{
template <typename T>
class NumericVector;
}

template <>
InputParameters validParams<VariableOrderBDFStepper>();

/**
 * A TimeStepper based on the the Fix Leading Coeff. Variable Order BDF method.
 * Forumlation can be found in Ch. 5 of "Numerical Solution of Initial-Vaglue
 * Problems in Differential-Algebral Equations" by K.E. Brenan, S. L. Campbell,
 * and L.R. Petzold. IDA of the Lawrence Livermore National Laboratory
 * was referenced as an example of implementation.
 */
class VariableOrderBDFStepper : public TimeStepper
{
public:
  static InputParameters validParams();

  VariableOrderBDFStepper(const InputParameters & parameters);

  virtual void step() override;
  virtual void preSolve() override;
  virtual bool converged() const override;
  virtual void postSolve() override;

  virtual int & KOrder() const { return _k_order; }

protected:
  virtual Real computeDT() override;
  virtual Real computeInitialDT() override;
  virtual Real computeFailedDT() override;

  virtual Real WeightedNorm(NumericVector<Number> & vector);
  virtual NumericVector<Number> & DifferenceCalculation(const NumericVector<Number> & _u0,
                                                        NumericVector<Number> & _u1,
                                                        NumericVector<Number> & _u2,
                                                        NumericVector<Number> & _u3,
                                                        NumericVector<Number> & _u4,
                                                        NumericVector<Number> & _u5,
                                                        NumericVector<Number> & _u6);
  virtual NumericVector<Number> & TERKM1VectorCalculation(NumericVector<Number> & diff,
                                                          NumericVector<Number> & _u1,
                                                          NumericVector<Number> & _u2,
                                                          NumericVector<Number> & _u3,
                                                          NumericVector<Number> & _u4,
                                                          NumericVector<Number> & _u5,
                                                          NumericVector<Number> & _u6);
  virtual NumericVector<Number> & TERKM2VectorCalculation(NumericVector<Number> & diff,
                                                          NumericVector<Number> & _u1,
                                                          NumericVector<Number> & _u2,
                                                          NumericVector<Number> & _u3,
                                                          NumericVector<Number> & _u4,
                                                          NumericVector<Number> & _u5,
                                                          NumericVector<Number> & _u6);

  NumericVector<Number> & _tmp_vector;
  NumericVector<Number> & _norm_weight;

  NumericVector<Number> & _diff;
  NumericVector<Number> & _diff_old;
  NumericVector<Number> & _TERKM1_vector;
  NumericVector<Number> & _TERKM2_vector;
  NumericVector<Number> & _TERKP1_vector;

  int & _k_order;
  int & _k_old;
  int & _k_new;

  Real & _dt_old;

  Real & _time_older;
  Real & _time_old3;
  Real & _time_old4;
  Real & _time_old5;
  Real & _time_old6;

  std::vector<Real> & _time_array;
  std::vector<Real> & _psi;
  std::vector<Real> & _psi_old;
  std::vector<Real> & _alpha;
  std::vector<Real> & _beta;
  std::vector<Real> & _phi_coeff;
  std::vector<Real> & _sigma;

  Real & _alpha_s;
  Real & _alpha_0;

  Real & _errk;
  Real & _terk;
  Real & _errkm1;
  Real & _terkm1;
  Real & _errkm2;
  Real & _terkm2;
  Real & _errkp1;
  Real & _terkp1;
  Real & _error;

  int & _err_counter;

  bool & _is_beginning_settup;

  Real _abs_tol;
  Real _rel_tol;
  int _stepup_step;
  int _k_max;
};
