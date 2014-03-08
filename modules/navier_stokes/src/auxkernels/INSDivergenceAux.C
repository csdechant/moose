<<<<<<< HEAD
//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "INSDivergenceAux.h"
#include "MooseMesh.h"

registerMooseObject("NavierStokesApp", INSDivergenceAux);

InputParameters
INSDivergenceAux::validParams()
{
  InputParameters params = AuxKernel::validParams();

  params.addClassDescription("Computes h_min / |u|.");
=======
#include "INSDivergenceAux.h"

template<>
InputParameters validParams<INSDivergenceAux>()
{
  InputParameters params = validParams<AuxKernel>();

>>>>>>> Merging Modules into MOOSE #2460
  // Coupled variables
  params.addRequiredCoupledVar("u", "x-velocity");
  params.addCoupledVar("v", "y-velocity"); // only required in 2D and 3D
  params.addCoupledVar("w", "z-velocity"); // only required in 3D

  return params;
}

<<<<<<< HEAD
INSDivergenceAux::INSDivergenceAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _grad_u_vel(coupledGradient("u")),
    _grad_v_vel(_mesh.dimension() >= 2 ? coupledGradient("v") : _grad_zero),
    _grad_w_vel(_mesh.dimension() == 3 ? coupledGradient("w") : _grad_zero)
{
}
=======
INSDivergenceAux::INSDivergenceAux(const std::string & name, InputParameters parameters)
  :AuxKernel(name, parameters),
  _grad_u_vel(coupledGradient("u")),
  _grad_v_vel(_mesh.dimension() >= 2 ? coupledGradient("v") : _grad_zero),
  _grad_w_vel(_mesh.dimension() == 3 ? coupledGradient("w") : _grad_zero)
{}
>>>>>>> Merging Modules into MOOSE #2460

Real
INSDivergenceAux::computeValue()
{
  // div U = du/dx + dv/dy + dw/dz
  return _grad_u_vel[_qp](0) + _grad_v_vel[_qp](1) + _grad_w_vel[_qp](2);
}
