#pragma once
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
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
	int n_e; // n_e = (m_w-1)*(m_h-1) // number of elements
	int order_of_integration; //number of Gauss points
	double thermal_conductivity;
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
		iss >> global_data.order_of_integration;
		iss >> global_data.thermal_conductivity;

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

//ELEMENT IN LOCAL COORDINATE SYSTEM
struct elem4
{
	double ksi[16];
	double eta[16];
	double ip;
	double ip_greater;
	double weight[4];
	double multiplier[16];
	
	elem4(int order_of_integration)
	{
		switch(order_of_integration)
		{
		case 2:
			ip = 1.0 / sqrt(3);
			ksi[0] = -ip;
			ksi[1] = ip;
			ksi[2] = ip;
			ksi[3] = -ip;

			eta[0] = -ip;
			eta[1] = -ip;
			eta[2] = ip;
			eta[3] = ip;
			for(auto n =0; n<4; n++) multiplier[n] = 1.0;
			break;
		case 3:
			ip = 1.0*sqrt(3.0 / 5.0);
			weight[0] = weight[2] = 5.0 / 9.0;
			weight[1] = 8.0 / 9.0;
			for ( auto n=0; n<3;n++)
			{
				ksi[0+n*3] = eta[0 + n]  = -ip;
				ksi[1+n*3] = eta[3 + n] = 0;
				ksi[2+n*3] = eta[6 + n]= ip;
				multiplier[0+n*3] = weight[0] * weight[n];
				multiplier[1+n*3] = weight[1] * weight[n];
				multiplier[2+n*3] = weight[2] * weight[n];
			}
			break;
		case 4:
			ip = 1.0 * sqrt((3.0 / 7.0) - (2.0 / 7.0)*sqrt((6.0 / 5.0)));
			ip_greater = 1.0 * sqrt((3.0 / 7.0) + (2.0 / 7.0)*sqrt((6.0 / 5.0)));
			weight[0] = weight[3] = (18.0 - sqrt(30.0)) / 36;
			weight[1] = weight[2] = (18.0 + sqrt(30.0)) / 36;
			for (auto n = 0; n < 4; n++)
			{
				ksi[0 + n * 4] = eta[0 + n] = -ip_greater;
				ksi[1 + n * 4] = eta[4 + n] = -ip;
				ksi[2 + n * 4] = eta[8 + n] = ip;
				ksi[3 + n * 4] = eta[12 + n] = ip_greater;
				multiplier[0 + n * 4] = weight[0] * weight[n];
				multiplier[1 + n * 4] = weight[1] * weight[n];
				multiplier[2 + n * 4] = weight[2] * weight[n];
				multiplier[3 + n * 4] = weight[3] * weight[n];
			}
			break;
		default:
			cout << "";
		}

	}
};


void calculate_H(element input_element[], int n_El, node ND[], int order_of_integration, double thermal_conductivity)
{
	

	int points = order_of_integration * order_of_integration;
	const int max_points = 16;
	
	//ITERATOR == NO OF CURRENT ELEMENT 
	for (auto iterator = 1 ; iterator < n_El ; iterator++ )
	{
		//ASSIGNING X,Y COORDINATES OF CURRENT ELEMENT 
		double x[4];
		double y[4];
		
		for (auto i = 0; i < 4; i++)
		{
			int temp = input_element[iterator].id[i];
			x[i] = ND[temp].x;
			y[i] = ND[temp].y;
		}
		
		//ASSIGNING VALUES TO SHAPE FUNCTIONS IN LOCAL COORDINATE SYSTEM
		double dN_dksi[max_points][4];
		double dN_deta[max_points][4];
		elem4 element(order_of_integration);

		//dN_dksi[integration point][number of shape function]
		for (auto ip = 0; ip < points; ip++)
		{
			dN_dksi[ip][0] = -1.0 / 4.0 * (1.0 - element.eta[ip]);
			dN_dksi[ip][1] = 1.0 / 4.0 * (1 - element.eta[ip]);
			dN_dksi[ip][2] = 1.0 / 4.0 * (1 + element.eta[ip]);
			dN_dksi[ip][3] = -1.0 / 4.0 * (1 + element.eta[ip]);

			dN_deta[ip][0] = -1.0 / 4.0 * (1 - element.ksi[ip]);
			dN_deta[ip][1] = -1.0 / 4.0 * (1 + element.ksi[ip]);
			dN_deta[ip][2] = 1.0 / 4.0* (1 + element.ksi[ip]);
			dN_deta[ip][3] = 1.0 / 4.0 * (1 - element.ksi[ip]);
		}

		//CALCULATING VALUES OF JACOBI MATRIX COMPONENTS
		double dx_dksi[max_points] = { 0.0 };
		double dy_dksi[max_points] = { 0.0 };
		double dx_deta[max_points] = { 0.0 };
		double dy_deta[max_points] = { 0.0 };
		
		for (auto ip = 0; ip < points; ip++)
		{
			for (auto i = 0; i < 4; i++)
			{
				dx_dksi[ip] += dN_dksi[ip][i] * x[i];
				dy_dksi[ip] += dN_dksi[ip][i] * y[i];
				dx_deta[ip] += dN_deta[ip][i] * x[i];
				dy_deta[ip] += dN_deta[ip][i] * y[i];
			}
		}

		//ASSIGNING VALUES TO JACOBI MATRIX COMPONENTS
		double J[max_points][2][2];
		double det_J[max_points];

		//J[integration_point][row][column]
		for (auto ip = 0; ip < points; ip++)
		{
			J[ip][0][0] = dy_deta[ip];
			J[ip][0][1] = -dy_dksi[ip];
			J[ip][1][0] = -dx_deta[ip];
			J[ip][1][1] = dx_dksi[ip];

			det_J[ip] = J[ip][0][0] * J[ip][1][1] - J[ip][0][1] * J[ip][1][0];
		}

		//CALCULATING DERIVATIVES OF SHAPE FUNCTIONS WITH RESPECT TO X,Y
		
		//[integration_point][dN1,dN2,dN3,dN4]
		double dN_dx[max_points][4];
		double dN_dy[max_points][4];

		for (auto ip = 0; ip < points; ip++)
		{
			for (auto j = 0; j < 4; j++)
			{
				dN_dx[ip][j] = 1.0 / det_J[ip] * (dN_dksi[ip][j] * dy_deta[ip] + dN_deta[ip][j] * (-dy_dksi[ip]));
				dN_dy[ip][j] = 1.0 / det_J[ip] * (dN_dksi[ip][j] * (-dx_deta[ip]) + dN_deta[ip][j] * dx_dksi[ip]);
			}
		}

		//CALCULATING H MATRIX
		
		//[integration_point][column][row]

		double dN_dx_dN_dx_T[max_points][4][4];
		double dN_dy_dN_dy_T[max_points][4][4];

		double H_point[max_points][4][4];

		for (auto ip = 0; ip < points; ip++)
		{
			for (auto i = 0; i < 4; i++)
			{
				for (auto j = 0; j < 4; j++)
				{
					dN_dx_dN_dx_T[ip][i][j] = dN_dx[ip][i] * dN_dx[ip][j];
					dN_dy_dN_dy_T[ip][i][j] = dN_dy[ip][i] * dN_dy[ip][j];
					H_point[ip][i][j] = thermal_conductivity * (dN_dx_dN_dx_T[ip][i][j] + dN_dy_dN_dy_T[ip][i][j]) * det_J[ip];
					input_element[iterator].H[i][j] += H_point[ip][i][j]*element.multiplier[ip];
					
				}
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

	calculate_H(Elem, n_El, ND, gdata.order_of_integration, gdata.thermal_conductivity);
	
	cout << "to jest elem8 h00: " << Elem[2].H[0][0] << endl;
	cout << "to jest elem8 h01: " << Elem[2].H[0][1] << endl;
	cout << "to jest elem8 h00: " << Elem[2].H[0][2] << endl;
	cout << "to jest elem8 h01: " << Elem[2].H[0][3] << endl;

	//delete the last row/column TODO

	double** HG = new double*[gdata.n_n];
	
	for (int i = 0; i < gdata.n_n; i++) 
	{
		HG[i] = new double[gdata.n_n];
	}
	
	for (int i = 0; i < gdata.n_n; i++)
	{
		for (int j = 0; j < gdata.n_n; j++)
		{
			HG[i][j] = 0.0;
		}
	}

	for (int k = 1; k < gdata.n_e; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				HG[Elem[k].id[i]-1][Elem[k].id[j]-1] += Elem[k].H[i][j];
			}
		}
	}

	for (int i = 0; i < gdata.n_n; i++)
	{
		for (int j = 0; j < gdata.n_n; j++)
		{
			cout << setfill(' ') << setw(5) << setprecision(4) << HG[i][j] << "\t";
		}
		cout << endl;
	}
}