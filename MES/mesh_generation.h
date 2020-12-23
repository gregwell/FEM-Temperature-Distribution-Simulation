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

struct element {
	int id[4]; //id of all four nodes from each side
	double H[4][4] = { 0.0 }; //H matrix of element
	double C[4][4] = { 0.0 }; //C matrix of element
	double P[4] = { 0.0 }; // P vector
};

struct global_data {
	double w; //width
	double h; //height
	int n_w; //number of nodes in width length
	int n_h; //number of nodes in height length
	int n_n; // m_n = m_w * m_h // number of nodes at all
	int n_e; // n_e = (m_w-1)*(m_h-1) // number of elements
	int order_of_integration; //number of Gauss points
	double thermal_conductivity; //(k)
	double density; //(ro)
	double specific_heat; //(c)
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
		iss >> global_data.thermal_conductivity;
		iss >> global_data.density;
		iss >> global_data.specific_heat;
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
		switch(order_of_integration)
		{
		case 2: //2 (2D)
			ip = 1.0 / sqrt(3);
			for(auto n =0; n<2;n++)
			{
				ksi[0 + n * 3] = eta[0 + n] = -ip;
				ksi[1 + n] = eta[2 + n] = ip;
			}
			for(auto n =0; n<4; n++) multiplier[n] = 1.0;
			break;
		case 3: //3 (2D)
			ip = 1.0*sqrt(3.0 / 5.0);
			weight[0] = weight[2] = 5.0 / 9.0;
			weight[1] = 8.0 / 9.0;
			for ( auto n=0; n<3;n++)
			{
				ksi[0+n*3] = eta[0 + n]  = -ip;
				ksi[1+n*3] = eta[3 + n] = 0;
				ksi[2+n*3] = eta[6 + n]= ip;
				for (auto i = 0; i < 3; i++) multiplier[i + n * 3] = weight[i] * weight[n];
			}
			break;
		case 4: //4 (2D)
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
				for (auto i=0; i<4;i++) multiplier[i + n * 4] = weight[i] * weight[n];
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
			ip = 1.0*sqrt(3.0 / 5.0);
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
			ip = 1.0 * sqrt((3.0 / 7.0) - (2.0 / 7.0)*sqrt((6.0 / 5.0)));
			ip_greater = 1.0 * sqrt((3.0 / 7.0) + (2.0 / 7.0)*sqrt((6.0 / 5.0)));
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

void calculate_H(element input_element[], int n_El, node ND[], int order_of_integration, double thermal_conductivity, double density, double specific_heat)
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

		//24.11 added (1)
		double N[max_points][4];
		for (auto ip = 0; ip < points; ip++)
		{
			N[ip][0] = 1.0 / 4.0 * (1.0 - element.ksi[ip]) * (1.0 - element.eta[ip]);
			N[ip][1] = 1.0 / 4.0 * (1.0 + element.ksi[ip]) * (1.0 - element.eta[ip]);
			N[ip][2] = 1.0 / 4.0 * (1.0 + element.ksi[ip]) * (1.0 + element.eta[ip]);
			N[ip][3] = 1.0 / 4.0 * (1.0 - element.ksi[ip]) * (1.0 + element.eta[ip]);
		}

		//CALCULATING H MATRIX +//24.11 added C MATRIX
		//[integration_point][column][row]
		double dN_dx_dN_dx_T[max_points][4][4];
		double dN_dy_dN_dy_T[max_points][4][4];
		double H_point[max_points][4][4];

		double NNT[max_points][4][4];
		double C_point[max_points][4][4];

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

					NNT[ip][i][j] = N[ip][i] * N[ip][j];
					C_point[ip][i][j] = specific_heat * density * NNT[ip][i][j] * det_J[ip];
					input_element[iterator].C[i][j] += C_point[ip][i][j]*element.multiplier[ip];
				}
			}
		}
	}
}

void calculateHBC(element input_element[], int order_of_integration, int n_El, node ND[], double convective_heat_transfer_coefficient, double delta_x, double delta_y, double ambient_temperature)
{
	// TODO: 3,4 point Gauss integration method solutions
	const int max_points = 16;
	int points = order_of_integration;
	int points_bc; //if there is boundary condition points_bc = {points}{0} (switch in elem4 constructor)
	
	if (points == 2) points_bc = 20;
	else if (points == 3) points_bc = 30;
	else points_bc = 40;
	
	elem4 element(points_bc);

	for (auto i = 0 ; i <12 ; i++)
	{
		cout << "ksi[" << i << "]= " << element.ksi[i] << endl;
		cout << "eta[" << i << "]= " << element.eta[i] << endl;
	}

	for ( auto i = 0; i<16; i++) cout << "multiplier[" << i << "]= " << element.multiplier[i] << endl;
	
	// TODO: rethink array sizes
	double N[max_points][4];
	double HBC_point[max_points][4][4];
	double NNT[max_points][4][4];
	double HBC[16][4][4] = {0.0};
	double L;

	double P_point[max_points][4];
	double P[max_points][4];
	
	for (auto iterator = 1; iterator < n_El; iterator++)
	{
		//boundary = 0 if the first surface, 2 if the right surface, 4i if the top surface, 6 if the left surface
		for (auto id_el = 0; id_el < 4; id_el++)
		{
			int id_el2;
			if (id_el < 3) id_el2 = id_el+1;
			else id_el2 = 0;
			
			if (ND[input_element[iterator].id[id_el]].BC == 1 && ND[input_element[iterator].id[id_el2]].BC == 1)
			{
				if(iterator==1) cout << "input element true id= " << input_element[iterator].id[id_el] << endl;
				if(iterator==1) cout << "input element true id2= " << input_element[iterator].id[id_el2] << endl;
				int boundary = id_el * points;
				//0 do liczby punktow, 2 do liczby punktow +2

				int ip_1d = 0;
				for (auto ip = boundary; ip < boundary + points; ip++)
				{
					N[ip_1d][0] = 1.0 / 4.0 * (1.0 - element.ksi[ip]) * (1.0 - element.eta[ip]);
					N[ip_1d][1] = 1.0 / 4.0 * (1.0 + element.ksi[ip]) * (1.0 - element.eta[ip]);
					N[ip_1d][2] = 1.0 / 4.0 * (1.0 + element.ksi[ip]) * (1.0 + element.eta[ip]);
					N[ip_1d][3] = 1.0 / 4.0 * (1.0 - element.ksi[ip]) * (1.0 + element.eta[ip]);
					if(iterator==1) cout << " N[ip_1d][0,1,2,3]= " << N[ip_1d][0]<< ", " << N[ip_1d][1] << ", " << N[ip_1d][2] << ", " << N[ip_1d][3] << ", " << endl;
					if (iterator == 1) cout << "ip: " << ip << ", ip_1d=" << ip_1d << endl;
					ip_1d++;
				}

				for (auto ip = 0; ip < points; ip++) // integration point 1d
				{
					for (auto i = 0; i < 4; i++) // surface
					{
						if (boundary == 0 || boundary == 4) L = delta_x;
						else L = delta_y;
						//TODO: patrzec po wspolrzednych a nie tak
						
						for (auto j = 0; j < 4; j++)
						{
							NNT[ip][i][j] = N[ip][i] * N[ip][j];
							
							HBC_point[ip][i][j] = convective_heat_transfer_coefficient * NNT[ip][i][j] * L/2; 
							input_element[iterator].H[i][j] += HBC_point[ip][i][j] * element.multiplier[ip];
						}

						P_point[ip][i] = -convective_heat_transfer_coefficient * N[ip][i] *ambient_temperature * L / 2;
						input_element[iterator].P[i] += P_point[ip][i] * element.multiplier[ip];
						
					}
				}			
			}
		}
	}
}


void print_square_matrix(double** CG, int n_n) {
	for (int i = 0; i < n_n; i++)
	{
		for (int j = 0; j < n_n; j++)
		{
			cout << setfill(' ') << setw(5) << setprecision(4) << CG[i][j] << "\t";
		}
		cout << endl;
	}
}

void inline generate_mesh()
{
	//READING FROM FILE AND GENERATING FEM MESH
	ifstream input_file("data.txt");
	global_data gdata;
	input_file >> gdata;
	int n_ND = 1; //node array starting index.
	double x = 0;
	double y = 0;
	double delta_x = gdata.w / (gdata.n_w - 1); //x-axis element length
	double delta_y = gdata.h / (gdata.n_h - 1); //y-axis element length
	gdata.n_n = gdata.n_w*gdata.n_h; //number of nodes
	gdata.n_e = (gdata.n_w - 1)*(gdata.n_h - 1); //number of elements

	//ASSIGNING COORDINATES(X,Y) TO EACH NODE (+ INITIAL TEMPERATURE)
	node *ND = new node[gdata.n_n + 1];
	for (auto i = 0; i < gdata.n_w; i++) //width
	{
		for (auto j = 0; j < gdata.n_h; j++) //height
		{
			ND[n_ND].x = i * delta_x;
			ND[n_ND].y = j * delta_y;
			ND[n_ND].t0 = gdata.initial_temperature; 
			n_ND++;
		}
	}

	//SETTING BOUNDARY CONDITION (BC) FLAG DEPENDING WHETHER THE NODE LIES ON A SIDE EDGE
	n_ND = 1;
	for (auto i = 0 ; i<gdata.n_n ; i++)
	{
		if (ND[n_ND].x == 0 || ND[n_ND].y == 0 || ND[n_ND].y == gdata.w || ND[n_ND].x == gdata.h)
		{
			ND[n_ND].BC = 1;
		}
		n_ND++;
	}

	//ASSINGING NODES TO ELEMENTS
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
			n_El++;
		}
		k++;
	}


	//DEFINING GLOBAL MATRICES.
	double PG[17] = { 0.0 }; // TODO: dynamic sizing
	double PR[17] = { 0.0 }; // replacement P for SOE
	
	double** HG = new double*[gdata.n_n]; 	// TODO : delete the last row/column 
	double** CG = new double*[gdata.n_n];
	double** HR = new double*[gdata.n_n]; //replacement H for SOE
	
	//ASSIGNING 0.0 TO NEWLY CREATED VALUES
	for (auto i = 0; i < gdata.n_n; i++)
	{
		HG[i] = new double[gdata.n_n];
		CG[i] = new double[gdata.n_n];
		HR[i] = new double[gdata.n_n]; //replacement H
		for (auto j = 0; j < gdata.n_n; j++)
		{
			HG[i][j] = { 0.0 };
			CG[i][j] = { 0.0 };
			HR[i][j] = { 0.0 }; // replacement H
		}
	}

	
	//HEAT TRANSFER PRCOESS SIMULATION: CALCULATIONS FOR EACH ITERATION
	double nt = gdata.simulation_time / gdata.simulation_step_time; // no. of iterations
	for (int iteration_no= 0; iteration_no < 1; iteration_no++)
	{
		cout << "ITERATION NO: " << iteration_no << "... " << endl;
		// 1. CALCULATE H (WITHOUT BC) AND C MATRIX
		calculate_H(Elem, n_El, ND, gdata.order_of_integration, gdata.thermal_conductivity, gdata.density, gdata.specific_heat);
		// 2. CALCULATE H(BC party only) AND SUM IT UP TO H, CALCULATE P MATRIX
		calculateHBC(Elem, gdata.order_of_integration, n_El, ND, gdata.convective_heat_transfer_coefficient, delta_x, delta_y, gdata.ambient_temperature);

		// 3. HG, CG, PG MATRICES AGGREGATION
		for (auto k = 1; k < gdata.n_e+1; k++)
		{
			for (auto i = 0; i < 4; i++)
			{
				cout << "P: elem nr: " <<k <<" "<< Elem[k].P[i] << endl;
				for (auto j = 0; j < 4; j++)
				{
					HG[Elem[k].id[i] - 1][Elem[k].id[j] - 1] += Elem[k].H[i][j];
					CG[Elem[k].id[i] - 1][Elem[k].id[j] - 1] += Elem[k].C[i][j];
				}
				PG[Elem[k].id[i]] += Elem[k].P[i];
			}
		}

		cout << "PRINTING CG MATRIX: " << endl;
		print_square_matrix(CG, gdata.n_n);
		cout << "PRINTING HG MATRIX: " << endl;
		print_square_matrix(HG, gdata.n_n);
		cout << "PRINTING PG VECTOR: " << endl;
		for (int j = 1; j < gdata.n_n+1; j++) cout << fixed << setprecision(2) << PG[j] << endl;

		// CALCULATE REPLACEMENT H (HR):
		for (auto i = 0; i < gdata.n_n; i++)
		{
			for (auto j = 0; j < gdata.n_n; j++)
			{
				HR[i][j] = HG[i][j] + CG[i][j] / gdata.simulation_step_time;
			}
		}

		//CALCULATE REPLACEMENT P (PR)
		for (auto i = 1; i < gdata.n_n+1; i++)
		{
			PR[i] = PG[i];
			for (auto j = 1; j < gdata.n_n+1; j++)
			{
				PR[i] -= (CG[i-1][j-1] / gdata.simulation_step_time) * ND[j].t0;
			}
		}
		
		cout << "PRINTING HR MATRIX: " << endl;
		print_square_matrix(HR, gdata.n_n);
		cout << "PRINTING PR VECTOR: " << endl;
		for (int j = 1; j < gdata.n_n+1; j++) cout << fixed << setprecision(2) << PR[j] << endl;

		//SOLVING THE SYSTEM OF LINEAR EQUATIONS: HR * {t} = PR


		

	}
	//???? reset element.H values?

}

