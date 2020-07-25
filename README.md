# FluxQon
This repository contains the core classes and functions I wrote to produce the simulations presented in my PhD thesis. Some basic quantum-modelling tools, such as [operator matrix construction](/Operator), [solving unitary time evolution](/+Solve/UTE.m), [solving Lindblad master equation](/+Solve/LME.m), and [calculating radiation spectrum](/Readout/RadiationSpectrum.m), are provided in this repository. Although the codes provided here can be tailored for more general uses, they were written specifically to do quantum simulations of flux qubits, YSO:Er<sup>3+</sup> ions, and combinations of the two. Coupled systems a flux qubit and paramagnetic ions are referred to in the thesis as *flux qons*. A flux qon with YSO:Er<sup>3+</sup> could be used to perform quantum frequency conversions between the microwave and the fiber-optic domain. For more information on the theory and dynamics of flux qons, please refer to [my PhD thesis](http://hdl.handle.net/11343/220532). Some [examples](https://github.com/heriantolim/FluxQon#examples) are given below to show how the codes in this repository can be implemented. These examples are simplified forms of the directives I used in my PhD work. The codes used to generate the figures in my PhD thesis are given in [PhDWork folder](/PhDWork) and explained in the [last section](https://github.com/heriantolim/FluxQon#phd-work) of this page.

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
- [MatVerCon](https://github.com/heriantolim/MatVerCon) - required for doing the PhD work
- [ResFileSys](https://github.com/heriantolim/ResFileSys) - required for doing the PhD work

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

## PhD Work
### Setting Up
In order to run the scripts in the [PhDWork folder](/PhDWork), the following needs to be set up first in your system.
1. Download or git-clone this repository and all the repositories listed in the [Requirements](https://github.com/heriantolim/PeakFit#requirements).
2. Place the `.m` files from **MatVerCon** in the default working directory of your MATLAB.
3. For each downloaded repository except **MatVerCon**, place the contents in a directory structure as follows:
```
{MATLAB_default_working_dir}/Packages/FluxQon/v1.0.0/
{MATLAB_default_working_dir}/Packages/MatCommon/v1.0.0/
...
```
4. Edit `librarypath.m` in **MatVerCon** so that the function returns a string that points to `{MATLAB_default_working_dir}/Packages`.
5. Move the contents in `FluxQon/v1.0.0/PhDWork` to a folder, say `{Your_home_folder}/Projects/`.
6. Edit `basepath.m` in **ResFileSys** so that the function returns a string that points to the project folder: `{Your_home_folder}/Projects`.

### Contents
The scripts in the [PhDWork folder](/PhDWork) perform numerical simulations that produce the figures in my PhD thesis Chapter 5 and 6 as outlined below.

| Date ID | Figure number |
|---------|---------------|
| 20160625 | 5.2.4 |
| 20160715 | 5.2.6, 5.2.7 |
| 20160925 | 5.1.1 |
| 20161215 | 5.3.1, 5.3.2 |
| 20170621 | 6.2.1 -- 6.2.3 |
| 20170626 | 6.3.13 |
| 20170627 | 6.3.14, 6.3.15 |
| 20170628 | 6.3.16 |
| 20180501 | 6.3.1 -- 6.3.8, 6.3.12 |
| 20180502 | 6.3.9 -- 6.3.11 |
| 20180521 | 6.4.4, 6.4.7 |
| 20180522 | 6.4.2 |
| 20180525 | 6.4.5, 6.4.6, 6.4.8, 6.4.9 |
| 20180603 | 6.4.3 |

### Running the Simulations
Some of the simulations require about 100GB of RAM. You will also need hundreds of CPU cores, otherwise the simulations will take years to complete. To run the simulations, run the following on your MATLAB console.
```
addpackage('ResFileSys');
openproject('hlim_FluxQon');
runall;
```
