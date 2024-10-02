# AuxComplexHeating

!syntax description /AuxKernels/AuxComplexHeating

## Overview

!style halign=left
The AuxComplexHeating object calculates the heating term imparted to the medium based on the conduction current. The term is define as:

\begin{equation}
  0.5 * Re \left(\sigma \; \vec{E} \cdot \vec{E}^{*} \right)
\end{equation}

where

- $\sigma$ is the conductivity of the medium, 
- $\vec{E}$ is the electric field, and
- $\vec{E}^{*}$ is the complex conjugate of the electric field.

## Example Input File Syntax

!listing aux_microwave_heating.i block=AuxKernels/aux_microwave_heating

!syntax parameters /AuxKernels/AuxComplexHeating

!syntax inputs /AuxKernels/AuxComplexHeating

!syntax children /AuxKernels/AuxComplexHeating
