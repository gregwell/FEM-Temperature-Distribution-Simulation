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

void calculate_jacobian()
{
	int x[4] = { 0,4,4,0 };
	int y[4] = { 0,0,4,4 };

	double N_ksi[4][4];
	double N_eta[4][4];

	double dx_dksi = 0.0;
	double dy_dksi = 0.0;
	double dx_deta = 0.0;
	double dy_deta = 0.0;

	double J[2][2];
	double det_J;


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

		//edit to powinno byc dla 4 pkt calkowania
		dx_dksi += N_ksi[0][i] * x[i];
		dy_dksi += N_ksi[0][i] * y[i];
		dx_deta += N_eta[0][i] * x[i];
		dy_deta += N_eta[0][i] * y[i];
	}

	//[row][column]
	J[0][0] = dy_deta;
	J[0][1] = -dy_dksi;
	J[1][0] = -dx_deta;
	J[1][1] = dx_dksi;

	det_J = J[0][0] * J[1][1] - J[0][1] * J[1][0];

	//cout << "Det J: " << det_J << endl;

	//Lab4:
	//[integration_point][dN1,dN2,dN3,dN4]
	double dN_dx[4][4];
	double dN_dy[4][4];


	for ( auto ip = 0; ip<4; ip++)
	{
		for ( auto j=0;j<4;j++)
		{
			dN_dx[ip][j] = 1.0 / det_J * (N_ksi[ip][j] * dy_deta + N_eta[ip][j] * (-dy_dksi));
			dN_dy[ip][j] = 1.0 / det_J * (N_ksi[ip][j] * (-dx_deta) + N_eta[ip][j] * dx_dksi);
		}
	}

	//[inetgration_point][column][row]
	double dN_dx_dN_dx_T[4][4][4];
	double dN_dy_dN_dy_T[4][4][4];

	for(auto ip=0;ip<4;ip++)
	{
		for (auto i=0;i<4;i++ )
		{
			for (auto j = 0; j < 4; j++)
			{
				dN_dx_dN_dx_T[ip][i][j] = dN_dx[ip][i] * dN_dx[ip][j];
				dN_dy_dN_dy_T[ip][i][j] = dN_dy[ip][i] * dN_dy[ip][j];
			}
		}
	}

	cout << ""


	
}