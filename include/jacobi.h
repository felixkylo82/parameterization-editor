#ifndef JACOBI_H
#define JACOBI_H

void   jacobi(double a[][3], int n, double d[], double v[][3], int &nrot);
double findNormalEstimate(double covarianceMatrix[][3], double normal[]);

#endif