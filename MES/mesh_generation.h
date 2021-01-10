#pragma once
#include <string>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <iostream>


using namespace std;

struct node
{
	double x; //coordinates
	double y;
	int BC;
	double t0;

	node()
	{
		x = 0.0;
		y = 0.0;
		BC = 0;
		t0 = 0.0;
	}
};

struct element
{
	int id[4]; //id of all four nodes from each side
	double H[4][4] = {0.0}; //H matrix of element
	double C[4][4] = {0.0}; //C matrix of element
	double P[4] = {0.0}; // P vector
	double thermal_conductivity; //(k)
	double density; //(ro)
	double specific_heat; //(c)
};

struct global_data
{
	double w; //width
	double h; //height
	int n_w; //number of nodes in width length
	int n_h; //number of nodes in height length
	int n_n; // m_n = m_w * m_h // number of nodes at all
	int n_e; // n_e = (m_w-1)*(m_h-1) // number of elements
	int order_of_integration; //number of Gauss points
	double thermal_conductivity[2]; //(k)
	double density[2]; //(ro)
	double specific_heat[2]; //(c)
	double convective_heat_transfer_coefficient; //(alfa)
	double ambient_temperature;
	double initial_temperature;
	double simulation_time;
	double simulation_step_time;

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
		iss >> global_data.thermal_conductivity[0];
		iss >> global_data.density[0];
		iss >> global_data.specific_heat[0];
		iss >> global_data.thermal_conductivity[1];
		iss >> global_data.density[1];
		iss >> global_data.specific_heat[1];
		iss >> global_data.convective_heat_transfer_coefficient;
		iss >> global_data.ambient_temperature;
		iss >> global_data.initial_temperature;
		iss >> global_data.simulation_time;
		iss >> global_data.simulation_step_time;
		return is;
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
		switch (order_of_integration)
		{
		case 2: //2 (2D)
			ip = 1.0 / sqrt(3);
			for (auto n = 0; n < 2; n++)
			{
				ksi[0 + n * 3] = eta[0 + n] = -ip;
				ksi[1 + n] = eta[2 + n] = ip;
			}
			for (auto n = 0; n < 4; n++) multiplier[n] = 1.0;
			break;
		case 3: //3 (2D)
			ip = 1.0 * sqrt(3.0 / 5.0);
			weight[0] = weight[2] = 5.0 / 9.0;
			weight[1] = 8.0 / 9.0;
			for (auto n = 0; n < 3; n++)
			{
				ksi[0 + n * 3] = eta[0 + n] = -ip;
				ksi[1 + n * 3] = eta[3 + n] = 0;
				ksi[2 + n * 3] = eta[6 + n] = ip;
				for (auto i = 0; i < 3; i++) multiplier[i + n * 3] = weight[i] * weight[n];
			}
			break;
		case 4: //4 (2D)
			ip = 1.0 * sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt((6.0 / 5.0)));
			ip_greater = 1.0 * sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt((6.0 / 5.0)));
			weight[0] = weight[3] = (18.0 - sqrt(30.0)) / 36;
			weight[1] = weight[2] = (18.0 + sqrt(30.0)) / 36;
			for (auto n = 0; n < 4; n++)
			{
				ksi[0 + n * 4] = eta[0 + n] = -ip_greater;
				ksi[1 + n * 4] = eta[4 + n] = -ip;
				ksi[2 + n * 4] = eta[8 + n] = ip;
				ksi[3 + n * 4] = eta[12 + n] = ip_greater;
				for (auto i = 0; i < 4; i++) multiplier[i + n * 4] = weight[i] * weight[n];
			}
			break;
		case 20: // 2 (1D)  for surfaces..
			//0,1 - first surface // 2,3 - second ... and so on
			ip = 1.0 / sqrt(3);
			ksi[0] = eta[2] = ksi[5] = eta[7] = -ip;
			ksi[1] = eta[3] = ksi[4] = eta[6] = ip;
			ksi[2] = ksi[3] = eta[4] = eta[5] = 1;
			eta[0] = eta[1] = ksi[6] = ksi[7] = -1;
			for (auto i = 0; i < 2; i++) multiplier[i] = 1;
			break;
		case 30: // 3 (1D)
			//0,1,2 - first surface // 3,4,5 - second ... and so on
			ip = 1.0 * sqrt(3.0 / 5.0);
			weight[0] = weight[2] = 5.0 / 9.0;
			weight[1] = 8.0 / 9.0;
			for (auto i = 0; i < 3; i++)
			{
				eta[i] = ksi[9 + i] = -1;
				eta[6 + i] = ksi[3 + i] = 1;
			}
			ksi[0] = ksi[6] = eta[3] = eta[9] = -ip;
			ksi[1] = ksi[7] = eta[4] = eta[10] = 0;
			ksi[2] = ksi[8] = eta[5] = eta[11] = ip;
			for (auto i = 0; i < 3; i++) multiplier[i] = weight[i];
			break;

		case 40:
			//0,1,2,3 - first surface // 4,5,6,7 - second ... and so on
			ip = 1.0 * sqrt((3.0 / 7.0) - (2.0 / 7.0) * sqrt((6.0 / 5.0)));
			ip_greater = 1.0 * sqrt((3.0 / 7.0) + (2.0 / 7.0) * sqrt((6.0 / 5.0)));
			weight[0] = weight[3] = (18.0 - sqrt(30.0)) / 36;
			weight[1] = weight[2] = (18.0 + sqrt(30.0)) / 36;
			for (auto i = 0; i < 4; i++)
			{
				eta[i] = ksi[12 + i] = -1;
				eta[8 + i] = ksi[4 + i] = 1;
			}
			ksi[0] = eta[4] = ksi[11] = eta[15] = -ip_greater;
			ksi[1] = eta[5] = ksi[10] = eta[14] = -ip;
			ksi[2] = eta[6] = ksi[9] = eta[13] = ip;
			ksi[3] = eta[7] = ksi[8] = eta[12] = ip_greater;
			for (auto i = 0; i < 4; i++) multiplier[i] = weight[i];
			break;
		default:
			cout << "";
		}
	}
};

void calculate_H_C(element input_element[], int n_El, node ND[], int order_of_integration);
void calculate_HBC_P(element input_element[], int order_of_integration, int n_El, node ND[],
                  double convective_heat_transfer_coefficient, double delta_x, double delta_y,
                  double ambient_temperature);
void print_square_matrix(double** CG, int n_n, int style);
double max(double* arr, int n);
double min(double* arr, int n);
double* solve_sof(int n_n, double** HR, double* PR);
void generate_mesh();
