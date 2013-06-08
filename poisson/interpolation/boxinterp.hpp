// Jan Becker, AVT.PT, RWTH Aachen University
// 2012-06-19

#pragma once

//Enumerate the squares (needed for the bilinear interpolation) in x- resp. in y-direction
int boxnumber(double x, double delta, double offset);

//Calculate a phase-value within [0,L], that is equivalent to the original one.
double phaseModuloWavelenghtLiangData(double phi, double L);

//Calculate space- ant time-depentent phase (c designates the phase-velocity)
double phase(double x, double t, double c);

// definition interpolation-function
double BilinearInterpol(double UQ11, double UQ21, double UQ12, double UQ22, double x1, double y1, double x2, double y2, double X, double Y);
