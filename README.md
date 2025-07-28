# Thermo-diffusive
This is code to solve thermo-diffusive combustion equations i.e Cobvectiv-Diffusive-Reactive partial differential equations in two-dimensional written in mordern Fortran 90.

This code solves the unsteady molecular transport and heat trasport equation along with source/sink term from combustion reaction to study the                            
thermo-diffusive instability in laminar flames with wide range of premixiedness. The chemical kinetics is onse-step overall reaction mechanism is used. The code is nondimensionalized.
So the parametric study with varying parameters Llike Lewis number, Damkohler number and premixedness is possible.
The speciality of this code is it uses the higher order compact schemes (6th order in this case) devised by Lele. 
This was written by Mr. Ganesh Vijaykumar and subsequently modified by Mr. David Bhatt at IIT Madras. (All glory to God)

  !More details can be found at following publication:                                                                   
  !David S. Bhatt, S.R. Chakravarthy (2012) Nonlinear dynamical behaviour of intrinsic thermal-diffusive oscillations of
laminar flames with varying premixedness, Combustion and Flame, 159(6), 2115-2125. https://doi.org/10.1016/j.combustflame.2012.01.025

In order to run, A list of input paramters are required and is described in inputs2d_info.txt
Also a sample input file inputs2d.txt is also included.
to compile use 
gfortran -o 2d tdcode_2d.f90
./compile.sh
to run use
./2d<inputs2d.txt
inlet.dat contains the mass fraction of fuel and oxidiser at the inlet and is generated using a mathematica code using the formula for premxideness
a sample output field is required tostart the code. 
