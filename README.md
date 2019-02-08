# FluxQon
This repository contains the core classes and functions I wrote to produce the simulations presented in my PhD thesis. Some basic quantum-modelling tools, such as [operator matrix construction](/Operator), [solving unitary time evolution](/+Solve/UTE.m), [solving Lindblad master equation](/+Solve/LME.m), and [calculating radiation spectrum](/Readout/RadiationSpectrum.m), are provided in this repository. Although the codes provided here can be tailored for more general uses, they were written specifically to do quantum simulations of flux qubits, YSO:Er<sup>3+</sup> ions, and combinations of the two. Coupled systems a flux qubit and paramagnetic ions are referred to in the thesis as *flux qons*. A flux qon with YSO:Er<sup>3+</sup> could be used to perform quantum frequency conversions between the microwave and the fiber-optic domain. For more information on the theory and dynamics of flux qons, please refer to [my PhD thesis](http://hdl.handle.net/11343/220532). Some [examples](https://github.com/heriantolim/FluxQon#examples) are given below to show how the codes in this repository can be implemented. These examples are simplified forms of the directives I used to produce the figures contained in my PhD thesis.

## Academic Use
If you use any of the materials in this repository, please consider citing my PhD thesis:

> Lim, H. (2017). *Erbium for optical modulation and quantum computation*. PhD thesis, School of Physics, The University of Melbourne. Retrieved from <http://hdl.handle.net/11343/220532>.

## Licensing
This software is licensed under the GNU General Public License (version 3).

## Tested On
- MATLAB R2017a - R2018a

## Requirements
- [MatCommon](https://github.com/heriantolim/MatCommon)
- [PhysConst](https://github.com/heriantolim/PhysConst)
- [MatGraphics](https://github.com/heriantolim/MatGraphics) - required for doing the examples

## Setting Up
1. Download or git-clone this repository and other repositories listed in the [Requirements](https://github.com/heriantolim/PeakFit#requirements).
2. Add the repositories to the MATLAB's search path via `addpath(genpath( ... ))` OR this [version control system](https://github.com/heriantolim/MatlabVerCon).

## Examples
### YSO:Er<sup>3+</sup> Energies vs Rotation Angle
[Examples/YSOEr_energy_vs_angle.m](/Examples/YSOEr_energy_vs_angle.m)

![YSO:Er3+ Energies vs Rotation Angle](/Examples/YSOEr_energy_vs_angle.png)

### Wavefunctions of a 3JJ Flux Qubit
[Examples/a3JJ_wavefunction.m](/Examples/a3JJ_wavefunction.m)

![Wavefunctions of a 3JJ Flux Qubit](/Examples/a3JJ_wavefunction.png)

### Solving Schrödinger Equation
[Examples/qon_mw_dynamics.m](/Examples/qon_mw_dynamics.m)

![Microwave Dynamics of a Flux Qon](/Examples/qon_mw_dynamics.png)

### Solving Lindblad Master Equation
[Examples/unitary_downconversion.m](/Examples/unitary_downconversion.m)

![Unitary Frequency Downconversion](/Examples/unitary_downconversion.png)

### Transmission Spectrum
[Examples/transmission_spectrum.m](/Examples/transmission_spectrum.m)

![Transmission Spectrum](/Examples/transmission_spectrum.png)
