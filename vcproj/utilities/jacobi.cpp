#include "jacobi.h"

#include <cfloat>
#include <cmath>
#include <iostream>

double findNormalEstimate(double covarianceMatrix[][3], double normal[])
{
	double d[3];
	double Q[3][3];
	int nrot;
	short unsigned int lambda1, lambda3;

	jacobi(covarianceMatrix, 3, d, Q, nrot);

	if (d[0] < d[1])
	{
		lambda1 = 0;
		lambda3 = 1;
	}
	else
	{
		lambda1 = 1;
		lambda3 = 0;
	}

	if (d[2] < d[lambda1])
	{
		lambda1 = 2;
	}
	else if (d[2] > d[lambda3])
	{
		lambda3 = 2;
	}

	normal[0] = Q[0][lambda1];
	normal[1] = Q[1][lambda1];
	normal[2] = Q[2][lambda1];

	return d[lambda1] / d[lambda3];
}

/* jocobi source from Numerical Recipes in C
* Copyright(C) 1988 - 1992 by Cambridge University Press
* Programs Copyright(C) 1988 - 1992 by Numerical Recipes Software
* memory leakage fixed by Felix Lo
*/

#define ROTATE(a, i, j, k, l) \
    g = a[i][j]; \
    h = a[k][l]; \
    a[i][j] = g - s * (h + g * tau); \
    a[k][l] = h + s * (g - h * tau);

void jacobi(double a[][3], int n, double d[], double v[][3], int &nrot)
{
	const double epsilon = 10.0f * FLT_EPSILON;
	const double zeroVector3[3] = { 0.0f, 0.0f, 0.0f, };
	const double Ident3[3][3] =
	{
		{ 1.0f, 0.0f, 0.0f, },
		{ 0.0f, 1.0f, 0.0f, },
		{ 0.0f, 0.0f, 1.0f, },
	};

	int   j, iq, ip, i;
	double tresh, theta, tau, t, sm, s, h, g, c, *b, *z;

	b = new double[n];
	z = new double[n];

	if (n == 3)
	{
		memcpy(v, Ident3, 9 * sizeof(double));
		memcpy(z, zeroVector3, 3 * sizeof(double));

		for (ip = 0; ip < 3; ++ip) {
			b[ip] = d[ip] = a[ip][ip];
		}

	}
	else
	{
		for (ip = 0; ip < n; ++ip)
		{
			for (iq = 0; iq < n; ++iq)
				v[ip][iq] = 0.0f;
			v[ip][ip] = 1.0f;

			z[ip] = 0.0f;
			b[ip] = d[ip] = a[ip][ip];
		}
	}

	nrot = 0;
	for (i = 0; i < 50; ++i)
	{
		sm = 0.0f;
		for (ip = 0; ip < n - 1; ++ip)
		{
			for (iq = ip + 1; iq < n; ++iq)
				sm += fabs(a[ip][iq]);
		}

		if (sm <= epsilon)
			break;

		if (i < 3)
			tresh = 0.2f * sm / (n*n);
		else
			tresh = 0.0f;

		for (ip = 0; ip < n - 1; ++ip)
		{
			for (iq = ip + 1; iq < n; ++iq)
			{
				g = 100.0f * fabs(a[ip][iq]);
				if (i > 3 && (fabs(d[ip]) + g) == fabs(d[ip]) &&
					(fabs(d[iq]) + g) == fabs(d[iq]))
				{
					a[ip][iq] = 0.0;
				}
				else if (fabs(a[ip][iq]) > tresh)
				{
					h = d[iq] - d[ip];
					if ((fabs(h) + g) == fabs(h))
					{
						t = (a[ip][iq]) / h;
					}
					else
					{
						theta = 0.5f * h / (a[ip][iq]);
						t = 1.0f / (fabs(theta) + sqrt(1.0f + theta*theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0f / sqrt(1 + t*t);
					s = t * c;
					tau = s / (1.0f + c);
					h = t * a[ip][iq];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;

					a[ip][iq] = 0.0f;
					for (j = 0; j <= ip - 1; ++j)
					{
						ROTATE(a, j, ip, j, iq);
					}
					for (j = ip + 1; j <= iq - 1; ++j)
					{
						ROTATE(a, ip, j, j, iq);
					}
					for (j = iq + 1; j < n; ++j)
					{
						ROTATE(a, ip, j, iq, j);
					}
					for (j = 0; j < n; ++j)
					{
						ROTATE(v, j, ip, j, iq);
					}
					++nrot;
				}
			}
		}

		for (ip = 0; ip < n; ++ip)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0f;
		}
	}

	if (sm > epsilon)
		std::cerr << "Too many iterations in routine jacobi" << std::endl;

	if (n > 32)
	{
		delete[] z;
		delete[] b;
	}

	return;
}
