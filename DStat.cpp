#include <iostream>
#include <cmath>
#include <math.h>

#include "DStat.h"

using namespace std;

namespace DFinMath
{
	double DStat::NormDist(double x)
	{
		// The cumulative normal distribution function            
		double z;
		if (x == 0)
			z = 0.5;
		else
		{
			double L, k;
			const double a1 = 0.31938153;
			const double a2 = -0.356563782;
			const double a3 = 1.781477937;
			const double a4 = -1.821255978;
			const double a5 = 1.330274429;
			
			L = abs(x);
			k = 1 / (1 + 0.2316419 * L);
			z = 1 - 1 / sqrt(2 * DStat::pi) * exp(-L * L / 2) * (k * (a1 + k * (a2 + k * (a3 + k * (a4 + k * a5)))));

			if (x < 0)
				z = 1 - z;
		}

		return z;
	}

	double DStat::N_Inv(double x)
	{
		//const double SQRT_TWO_PI = 2.506628274631;
		const double e_1 = -39.6968302866538;
		const double e_2 = 220.946098424521;
		const double e_3 = -275.928510446969;
		const double e_4 = 138.357751867269;
		const double e_5 = -30.6647980661472;
		const double e_6 = 2.50662827745924;

		const double f_1 = -54.4760987982241;
		const double f_2 = 161.585836858041;
		const double f_3 = -155.698979859887;
		const double f_4 = 66.8013118877197;
		const double f_5 = -13.2806815528857;

		const double g_1 = -0.00778489400243029;
		const double g_2 = -0.322396458041136;
		const double g_3 = -2.40075827716184;
		const double g_4 = -2.54973253934373;
		const double g_5 = 4.37466414146497;
		const double g_6 = 2.93816398269878;
		
		const double h_1 = 0.00778469570904146;
		const double h_2 = 0.32246712907004;
		const double h_3 = 2.445134137143;
		const double h_4 = 3.75440866190742;

		const double x_l = 0.02425;
		const double x_u = 0.97575;
		
		double z, r;

		// Lower region: 0 < x < x_l
		if (x < x_l)
		{
			z = sqrt(-2.0 * log(x));
			z = (((((g_1 * z + g_2) * z + g_3) * z + g_4) * z + g_5) * z + g_6) / ((((h_1 * z + h_2) * z + h_3) * z + h_4) * z + 1.0);
		}
		// Central region: x_l <= x <= x_u
		else if (x <= x_u)
		{
			z = x - 0.5;
			r = z * z;
			z = (((((e_1 * r + e_2) * r + e_3) * r + e_4) * r + e_5) * r + e_6) * z / (((((f_1 * r + f_2) * r + f_3) * r + f_4) * r + f_5) * r + 1.0);
		}
		// Upper region. ( x_u < x < 1 )
		else
		{
			z = sqrt(-2.0 * log(1.0 - x));
			z = -(((((g_1 * z + g_2) * z + g_3) * z + g_4) * z + g_5) * z + g_6) / ((((h_1 * z + h_2) * z + h_3) * z + h_4) * z + 1.0);
		}

		// Now |relative error| < 1.15e-9.  One iteration of Halley's third
		// order zero finder gives full machine precision:
		//
		//r = (N(z) - x) * SQRT_TWO_PI * exp( 0.5 * z * z );  //  f(z)/df(z)
		//z -= r/(1+0.5*z*r);

		return z;
	}
}
