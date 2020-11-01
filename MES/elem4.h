#pragma once
#include <cmath>

struct elem4
{
	double ksi[4];
	double eta[4];
	double ip = 1.0 / sqrt(3);
	elem4()
	{
		ksi[0] = - ip;
		ksi[1] = ip;
		ksi[2] = ip;
		ksi[3] = -ip;

		eta[0] = - ip;
		eta[1] = -ip;
		eta[2] = ip;
		eta[3] = ip;
	}
};

void calculate_H()
{
	int x[4] = { 0,4,4,0 };
	int y[4] = { 0,0,4,4 };

	double N_ksi[4][4];
	double N_eta[4][4];

	double dx_dksi[4] = { 0.0,0.0,0.0,0.0 };
	double dy_dksi[4] = { 0.0,0.0,0.0,0.0 };
	double dx_deta[4] = { 0.0,0.0,0.0,0.0 };
	double dy_deta[4] = { 0.0,0.0,0.0,0.0 };

	double J[4][2][2];
	double det_J[4];


	elem4 element;

	//Assigning values to the arrays dN1/dksi , dN2/dksi ... dN4/dksi
	//N_ksi [integration point][N1,N2,N3,N4]
	for (auto i = 0; i < 4; i++)
	{
		N_ksi[i][0] = -1.0 / 4.0 * (1.0 - element.eta[i]);
		N_ksi[i][1] = 1.0 / 4.0 * (1 - element.eta[i]);
		N_ksi[i][2] = 1.0 / 4.0 * (1 + element.eta[i]);
		N_ksi[i][3] = -1.0 / 4.0 * (1 + element.eta[i]);

		N_eta[i][0] = -1.0 / 4.0 * (1 - element.ksi[i]);
		N_eta[i][1] = -1.0 / 4.0 * (1 + element.ksi[i]);
		N_eta[i][2] = 1.0 / 4.0* (1 + element.ksi[i]);
		N_eta[i][3] = 1.0 / 4.0 * (1 - element.ksi[i]);
	}
	
	for (auto ip = 0; ip < 4; ip++)
	{
		for (auto i=0;i<4;i++)
		{
			dx_dksi[ip] += N_ksi[ip][i] * x[i];
			dy_dksi[ip] += N_ksi[ip][i] * y[i];
			dx_deta[ip] += N_eta[ip][i] * x[i];
			dy_deta[ip] += N_eta[ip][i] * y[i];
		}
	}

	//J[integration_point][row][column]
	for (auto ip=0;ip<4;ip++)
	{
		J[ip][0][0] = dy_deta[ip];
		J[ip][0][1] = -dy_dksi[ip];
		J[ip][1][0] = -dx_deta[ip];
		J[ip][1][1] = dx_dksi[ip];

		det_J[ip] = J[ip][0][0] * J[ip][1][1] - J[ip][0][1] * J[ip][1][0];
	}



	//Lab4:
	//[integration_point][dN1,dN2,dN3,dN4]
	double dN_dx[4][4];
	double dN_dy[4][4];


	for ( auto ip = 0; ip<4; ip++)
	{
		for ( auto j=0;j<4;j++)
		{
			dN_dx[ip][j] = 1.0 / det_J[ip] * (N_ksi[ip][j] * dy_deta[ip] + N_eta[ip][j] * (-dy_dksi[ip]));
			dN_dy[ip][j] = 1.0 / det_J[ip] * (N_ksi[ip][j] * (-dx_deta[ip]) + N_eta[ip][j] * dx_dksi[ip]);
		}
	}



	
	
	//[integration_point][column][row]
	double dN_dx_dN_dx_T[4][4][4];
	double dN_dy_dN_dy_T[4][4][4];
	double H_point[4][4][4];
	double H[4][4] = {0.0};

	for(auto ip=0;ip<4;ip++)
	{
		for (auto i=0;i<4;i++ )
		{
			for (auto j = 0; j < 4; j++)
			{
				dN_dx_dN_dx_T[ip][i][j] = dN_dx[ip][i] * dN_dx[ip][j];
				dN_dy_dN_dy_T[ip][i][j] = dN_dy[ip][i] * dN_dy[ip][j];
				H_point[ip][i][j] = 30 * (dN_dx_dN_dx_T[ip][i][j] + dN_dy_dN_dy_T[ip][i][j]) * det_J[ip];
				H[i][j] += H_point[ip][i][j];
			}
		}
	}
	
}