#pragma once
#define _HAS_STD_BYTE 0
#include<iostream>
#include<string>
#include<math.h>
#include "meshengine.h"
#include "Timer.hpp"
#include<complex>
#include "utilities.h"


using namespace std;
class EFEM{
protected:
	struct Portmodes {
		dmatrix <double> mode;
		double kz;
		double fc;
	};
	double epsr[30];
	double mur[30];
	double sigma[30];
	MeshEngine me; 

	
public:
	Timer tmr;
 	const double mu0=4*pi*1e-7;
    const double eps0=8.85e-12;
     double light_velocity = 299792458;
    const double pi=3.141592653589793238462643383279;
	Portmodes ipportmode;
	
	EFEM();

	void load(string fname);

	void start(double f) ;
	
	/* Extracts single dominant mode in xy plane */
	Portmodes portModexy(int fno,double f);

	
	
};