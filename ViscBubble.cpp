// ViscBubble.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <math.h>
#include <time.h>

//#define _CRT_SECURE_NO_WARNINGS

#include "Data.h"

double D_dihot(double eps)
{
	double D_min, D_max, D_int, D_tmp, dXi, Xi, Eq_tmp;
	int N, i;

	N = 10000000;
	D_min = 0.0;
	D_max = 1.0e10;
	dXi = 1.0e-4;

	D_tmp = (D_max+D_min) / 2.0;
	D_tmp = D_max;
	do{
		D_int = 0.0;
		for (i = 0; i < N; i++)
		{
			Xi = 1.0 + i*dXi;
			D_int = D_int + (exp(-(Xi*Xi / 2.0 + 1.0 / Xi - 3.0 / 2.0)*D_tmp / 2.0) / (Xi*Xi))*dXi;
		}
		Eq_tmp = D_tmp / 2.0*D_int - eps;
		if (Eq_tmp < 0.0) { D_min = D_tmp; }
		else { D_max = D_tmp; };
		D_tmp = (D_max + D_min) / 2.0;
	} while (D_max - D_min > 1.0e-10);

	return D_tmp;
}

int _tmain(int argc, _TCHAR* argv[])
{
	Nu_flag = 0;
	Cp_flag = 1;

	Tm = 10000.0;
	T_i = 0.0;
	tau = 1.0e-5;
	T_out = 10.0 * tau;

	P_i = 1500.0*p0;
	P_f = 150.0*p0;
	P = P_i;
	dP = P_i - P_f;
	Rb = 2.0*sigma / (dP);
	Pg = P_i;
	Mg = 4.0*Pi * Pg * Rb * Rb * Rb * MH2O / (3.0* Rg * T);
	Rho_g = P_f * MH2O / (Rg * T);

	Nu_0 = 1.0e4;

	Cp_i = Kh*sqrt(P_i);
	Cp_f = Kh*sqrt(P_f);
	Cp = Kh*sqrt(Pg);


	eps = rho0*(Cp_i-Cp_f) / Rho_g;
	if (eps > 10.0) { D_eff = 12.0*De*eps*eps / Pi;}
	else { D_eff = 2.0*De*eps;};
	
	D_eff = De*D_dihot(eps);


	Im = 250;
	dr = 1.0e-9;

	r = new double[Im];
	for (i = 0; i <= Im; i++){
		r[i] = Rb+i*dr;
	}
	
	CpR = new double[Im];	
	NuR = new double[Im];
	for (i = 0; i <= Im; i++){
		CpR[i] = Cp_i;
		E = E0*(1 - 12.0 * CpR[i]);
		NuR[i] = mu0*exp(E / (kB*T));
	}

	if (Nu_flag == 0) { Nu = Nu_0; }
	else { Nu = NuR[0]; };
	
	RRb = (Pg - P)*Rb / (4.0*Nu);

	out_file = fopen("ViscBubble.dat", "wt");
	fprintf(out_file, "# P_i = %10.8lf \n", P_i / p0);
	fprintf(out_file, "# P_f = %10.8lf \n", P_f / p0);
	fprintf(out_file, "# Cp_i = %10.8lf \n", Cp_i);
	fprintf(out_file, "# Cp_f = %10.8lf \n", Cp_f);
	fprintf(out_file, "# eps = %10.8lf \n", eps);
	fprintf(out_file, "# D_eff = %e \n", D_eff);
	fprintf(out_file, "# D_min = %e \n", 2.0*De*eps);
	fprintf(out_file, "# D_max = %e \n", 12.0*De*eps*eps / Pi);
	fprintf(out_file, "# Exp = %10.8lf \n", dP/(4.0*Nu));

	Cp_file = fopen("Cp.dat", "wt");
	out_num = 0;

	do {
		printf("Time %3.8lf s \n", T_i);

		if (T_i > out_num * T_out) {

			fprintf(out_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %e \t %10.8lf \n",
				T_i, P / p0, Pg / p0, Rb*1.0e6, Mg, Nu);

			if (out_num < 100000) { sprintf(out_name, "Data/%d.dat", out_num); };
			if (out_num < 10000) { sprintf(out_name, "Data/0%d.dat", out_num); };
			if (out_num < 1000) { sprintf(out_name, "Data/00%d.dat", out_num); };
			if (out_num < 100) { sprintf(out_name, "Data/000%d.dat", out_num); };
			if (out_num < 10) { sprintf(out_name, "Data/0000%d.dat", out_num); };
			cut_file = fopen(out_name, "wt");
			for (i = 0; i <= Im; i++){
				fprintf(cut_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n", r[i] * 1.0e6, r[i] / r[0], CpR[i], NuR[i]);
				fprintf(Cp_file, "%10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \t %10.8lf \n", T_i, r[i] * 1.0e6, r[i]/r[0], CpR[i], NuR[i]);
			}
			fclose(cut_file);
			out_num = out_num + 1;
		}
		dr = 49.0*Rb / Im;
		for (i = 0; i <= Im; i++){
			r[i] = Rb + i*dr;
		}

		Cp_int = 0.0;
		for (i = Im; i >= 0; i--)
		{
			Cp_int = Cp_int + dr*exp(-(r[i] * r[i] / 2.0 + 1.0 / r[i])*Rb*RRb) / (r[i] * r[i]);
			CpR[i] = Cp_int;
		}


		for (i = Im; i >= 0; i--)
		{
			CpR[i] = CpR[i] / CpR[0];
			CpR[i] = Cp_i - (Cp_i - Cp)* CpR[i];
			E = E0*(1 - 12.0 * CpR[i]);
			NuR[i] = mu0*exp(E / (kB*T));
		}
		
		Nu = 0.0;
		for (i = 0; i <Im; i++)
		{
			E = E0*(1 - 12.0 * CpR[i]);
			NuR[i] = mu0*exp(E / (kB*T));
		//	NuR[i] = 1.0;
		//	NuR[i+1] = 1.0;
			r_tmp = (r[i + 1] + r[i]) / 2.0;
			Nu = Nu + (NuR[i + 1] + NuR[i]) / 2.0 * dr / (r_tmp * r_tmp * r_tmp * r_tmp);
		}
		
		if (Nu_flag == 0) { Nu = Nu_0; }
		else { Nu = 3.0* Rb*Rb*Rb*Nu;}

		T_i = T_i + tau;


		P = P_f;

		if (Cp_flag == 0) { dCp = Cp_i - Cp_f;}
		else { dCp = (CpR[1] - CpR[0]) / dr; };
		if (dCp < 0.0) { dCp = 0.0; }

		RRb = (Pg - P)*Rb / (4.0*Nu);
		Rb = Rb + RRb *tau;
		Mg = Mg + 4.0*Pi* Rb * Rb * rho0 * De * dCp *tau;
		Pg = Mg * 3.0* Rg * T / (4.0*Pi *  Rb * Rb * Rb * MH2O);
		Cp = Kh*sqrt(Pg);

	

	} while (T_i < Tm+tau);


	fclose(out_file);	
	fclose(Cp_file);
	return 0;
}

