# Thermo-diffusive
This is code to solve thermo-diffusive combustion equations.

This code solves the unsteady molecular transport and heat trasport equation  to study the                            
thermo-diffusive instability in laminar flames with wide range of premixiedness. The speciality of this code is it use\
s the higher order (6th order) compact schemes devised by Lele. This was written by Mr. Ganesh Vijaykumar and subsequently modified by Mr. David Bhatt at IIT Madras. (All glory to God)

  !More details can be found at following publication:                                                                   
  !David S. Bhatt, S.R. Chakravarthy (2012) Nonlinear dynamical behaviour of intrinsic thermal-diffusive oscillations of
laminar flames with varying premixedness, Combustion and Flame, 159(6), 2115-2125. https://doi.org/10.1016/j.combustflame
.2012.01.025

In order to run, A list of input paramters are required and is described in inputs2d_info.txt
Also a sample input file inputs2d.txt is also included
inlet.dat contains the mass fraction of fuel and oxidiser at the inlet and is generated using a mathematica code using the formula for premxideness
a sample output field is required tostart the code. 