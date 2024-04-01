#ifndef DFinMath_DStat_h
#define DFinMath_DStat_h

#pragma once

namespace DFinMath
{
	class DStat
	{
	public:
		const double pi = 3.14159265358979323846;
		double NormDist(double x);
		double N_Inv(double x);
	};
}

#endif
