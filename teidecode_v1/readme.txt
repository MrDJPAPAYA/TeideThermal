# Thermal Subsystem - TEIDESAT-I

Leyre Hernández Palacios && Javier González Vilar

## How to use the code
In order to run the full simulation
only running the plotter.m script is enough.
When run, the script will prompt the user to enter
the desired orbit case to run. Valid inputs are:
* 1 -> theta_SC = 0 and phi_SC = 0
* 2 -> theta_SC = 45deg and phi_SC = 0
If other value is input the program will run case 1.
When finished the script will ask the user whether
to run another case or not, repeating the process.
If the answer is no, the program will ask whether to
close the plots and clear results and cmd. Answer
these questions with true/1 or false/0.

## data.m
This script contains:
* the data of the particular problem
* constants
* functions for the orbiting problem

## solver.m
This scripts calls the data.m script
It implements the solver for 10 orbital periods
and stores the results for the last period.

## plotter.m
This scripts calls the solver.m script
It generates the plots that answer the questions
for the laboratory report.