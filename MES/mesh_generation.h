#pragma once
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;

struct node {
	double x; //coordinates
	double y;
	node()
	{
		x = 0.0;
		y = 0.0;
	}
};

//struct element_coordinates {
//	double x[4]; //coordinates
//	double y[4];
//};

struct element {
	int id[4]; //id of all four nodes from each side
	double H[4][4] = {0.0}; //H matrix of element
};

struct global_data {
	double w; //width
	double h; //height
	int n_w; //number of nodes in width length
	int n_h; //number of nodes in height length
	int n_n; // m_n = m_w * m_h // number of nodes at all
	int n_e; // m_e = (m_w-1)*(m_h-1) // number of elements?
	friend std::istream& operator>>(std::istream& is, global_data& global_data)
	{
		std::string line;
		std::getline(is, line);
		std::istringstream iss(line);

		std::string temp;
		iss >> global_data.w;
		iss >> global_data.h;
		iss >> global_data.n_w;
		iss >> global_data.n_h;

		return is;
	}
	double return_w()
	{
		return w;
	}
	double return_h()
	{
		return h;
	}
};

struct elem4
{
	double ksi[4];
	double eta[4];
	double ip = 1.0 / sqrt(3);
	elem4()
	{
		ksi[0] = -ip;
		ksi[1] = ip;
		ksi[2] = ip;
		ksi[3] = -ip;

		eta[0] = -ip;
		eta[1] = -ip;
		eta[2] = ip;
		eta[3] = ip;
	}
	double H[4][4] = { 0.0 };
};



void calculate_H(element input_element[], int n_El, node ND[])
{
	for (auto iterator = 1 ; iterator < n_El ; iterator++ )
	{
		double x[4];
		double y[4];
		int temp;

		for (auto i = 0; i < 4; i++)
		{
			temp = input_element[iterator].id[i];
			x[i] = ND[temp].x;
			y[i] = ND[temp].y;
		}

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
			for (auto i = 0; i < 4; i++)
			{
				dx_dksi[ip] += N_ksi[ip][i] * x[i];
				dy_dksi[ip] += N_ksi[ip][i] * y[i];
				dx_deta[ip] += N_eta[ip][i] * x[i];
				dy_deta[ip] += N_eta[ip][i] * y[i];
			}
		}

		//J[integration_point][row][column]
		for (auto ip = 0; ip < 4; ip++)
		{
			J[ip][0][0] = dy_deta[ip];
			J[ip][0][1] = -dy_dksi[ip];
			J[ip][1][0] = -dx_deta[ip];
			J[ip][1][1] = dx_dksi[ip];

			det_J[ip] = J[ip][0][0] * J[ip][1][1] - J[ip][0][1] * J[ip][1][0];
		}

		//[integration_point][dN1,dN2,dN3,dN4]
		double dN_dx[4][4];
		double dN_dy[4][4];

		for (auto ip = 0; ip < 4; ip++)
		{
			for (auto j = 0; j < 4; j++)
			{
				dN_dx[ip][j] = 1.0 / det_J[ip] * (N_ksi[ip][j] * dy_deta[ip] + N_eta[ip][j] * (-dy_dksi[ip]));
				dN_dy[ip][j] = 1.0 / det_J[ip] * (N_ksi[ip][j] * (-dx_deta[ip]) + N_eta[ip][j] * dx_dksi[ip]);
			}
		}


		//[integration_point][column][row]
		double dN_dx_dN_dx_T[4][4][4];
		double dN_dy_dN_dy_T[4][4][4];
		double H_point[4][4][4];

		for (auto ip = 0; ip < 4; ip++)
		{
			for (auto i = 0; i < 4; i++)
			{
				for (auto j = 0; j < 4; j++)
				{
					dN_dx_dN_dx_T[ip][i][j] = dN_dx[ip][i] * dN_dx[ip][j];
					dN_dy_dN_dy_T[ip][i][j] = dN_dy[ip][i] * dN_dy[ip][j];
					H_point[ip][i][j] = 30 * (dN_dx_dN_dx_T[ip][i][j] + dN_dy_dN_dy_T[ip][i][j]) * det_J[ip];
					element.H[i][j] += H_point[ip][i][j];
				}
			}
		}

		for (auto i = 0; i < 4; i++)
		{
			for (auto j = 0; j < 4; j++)
			{
				input_element[iterator].H[i][j] = element.H[i][j];
			}
		}
	}
}

void inline generate_mesh()
{

	ifstream input_file("data.txt");
	global_data gdata;
	input_file >> gdata;

	//TEST
	cout << "gdata w:" << gdata.return_w() << endl;
	cout << "gdata h:" << gdata.return_h() << endl;
	//TEST/End

	//func

	int n_ND = 1; //node array number.
	double x = 0;
	double y = 0;
	double delta_x = gdata.w / (gdata.n_w - 1); //x-axis element length
	double delta_y = gdata.h / (gdata.n_h - 1); //y-axis element length
	gdata.n_n = gdata.n_w*gdata.n_h; //number of nodes
	gdata.n_e = (gdata.n_w - 1)*(gdata.n_h - 1); //number of elements

	//Assigning coordinates (x,y) to every node.
	node *ND = new node[gdata.n_n + 1];

	for (auto i = 0; i < gdata.n_w; i++) //width
	{
		for (auto j = 0; j < gdata.n_h; j++) //height
		{
			ND[n_ND].x = i * delta_x;
			ND[n_ND].y = j * delta_y;
			cout << "I'm inside node n_Nd=" << n_ND << "/ ND[n_ND].x = " << i * delta_x << "/ ND[n_ND].y=" << j * delta_y << endl;
			n_ND++;

		}
	}

	//Assigning nodes to elements.
	int n_El = 1; //node array number
	int k = 0;
	element *Elem = new element[gdata.n_e + 1];

	for (auto i = 0; i < gdata.n_w - 1; i++) //width
	{
		for (auto j = 0; j < gdata.n_h - 1; j++) //height
		{
			Elem[n_El].id[0] = n_El + k;
			Elem[n_El].id[1] = n_El + gdata.n_h + k;
			Elem[n_El].id[2] = n_El + gdata.n_h + 1 + k;
			Elem[n_El].id[3] = n_El + 1 + k;
			cout << "I'm inside n_el=" << n_El << "/ id1=  " << n_El + k << "/// id2= " << n_El + gdata.n_h + k << endl;
			n_El++;
		}
		k++;
	}

	calculate_H(Elem, n_El, ND);
	
	cout << "to jest elem8 h00: " << Elem[8].H[1][0] << endl;
	cout << "to jest elem8 h01: " << Elem[8].H[1][1] << endl;
	cout << "to jest elem8 h00: " << Elem[8].H[1][2] << endl;
	cout << "to jest elem8 h01: " << Elem[8].H[1][3] << endl;
}