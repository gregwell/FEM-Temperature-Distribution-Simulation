#include "pch.h"
#include <iostream>
#include <fstream>
#include "mesh_generation.h"

using namespace std;

void calculate_H(element input_element[], int n_El, node ND[], int order_of_integration)
{
	int points = order_of_integration * order_of_integration;
	const int max_points = 16;



	//ITERATOR == NO OF CURRENT ELEMENT 
	for (auto iterator = 1; iterator < n_El; iterator++)
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
					H_point[ip][i][j] = input_element[iterator].thermal_conductivity * (dN_dx_dN_dx_T[ip][i][j] + dN_dy_dN_dy_T[ip][i][j]) * det_J[ip];
					input_element[iterator].H[i][j] += H_point[ip][i][j] * element.multiplier[ip];

					NNT[ip][i][j] = N[ip][i] * N[ip][j];
					C_point[ip][i][j] = input_element[iterator].specific_heat * input_element[iterator].density * NNT[ip][i][j] * det_J[ip];
					input_element[iterator].C[i][j] += C_point[ip][i][j] * element.multiplier[ip];
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

	//for (auto i = 0 ; i <12 ; i++)
	//{
	//	cout << "ksi[" << i << "]= " << element.ksi[i] << endl;
	//	cout << "eta[" << i << "]= " << element.eta[i] << endl;
	//}

	//for ( auto i = 0; i<16; i++) cout << "multiplier[" << i << "]= " << element.multiplier[i] << endl;

	// TODO: rethink array sizes
	double N[max_points][4];
	double HBC_point[max_points][4][4];
	double NNT[max_points][4][4];
	double HBC[16][4][4] = { 0.0 };
	double L;

	double P_point[max_points][4];
	double P[max_points][4];

	for (auto iterator = 1; iterator < n_El; iterator++)
	{
		//boundary = 0 if the first surface, 2 if the right surface, 4i if the top surface, 6 if the left surface
		for (auto id_el = 0; id_el < 4; id_el++)
		{
			int id_el2;
			if (id_el < 3) id_el2 = id_el + 1;
			else id_el2 = 0;

			if (ND[input_element[iterator].id[id_el]].BC == 1 && ND[input_element[iterator].id[id_el2]].BC == 1)
			{
				//if(iterator==1) cout << "input element true id= " << input_element[iterator].id[id_el] << endl;
				//if(iterator==1) cout << "input element true id2= " << input_element[iterator].id[id_el2] << endl;
				int boundary = id_el * points;
				//0 do liczby punktow, 2 do liczby punktow +2

				int ip_1d = 0;
				for (auto ip = boundary; ip < boundary + points; ip++)
				{
					N[ip_1d][0] = 1.0 / 4.0 * (1.0 - element.ksi[ip]) * (1.0 - element.eta[ip]);
					N[ip_1d][1] = 1.0 / 4.0 * (1.0 + element.ksi[ip]) * (1.0 - element.eta[ip]);
					N[ip_1d][2] = 1.0 / 4.0 * (1.0 + element.ksi[ip]) * (1.0 + element.eta[ip]);
					N[ip_1d][3] = 1.0 / 4.0 * (1.0 - element.ksi[ip]) * (1.0 + element.eta[ip]);
					//if(iterator==1) cout << " N[ip_1d][0,1,2,3]= " << N[ip_1d][0]<< ", " << N[ip_1d][1] << ", " << N[ip_1d][2] << ", " << N[ip_1d][3] << ", " << endl;
					//if (iterator == 1) cout << "ip: " << ip << ", ip_1d=" << ip_1d << endl;
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

							HBC_point[ip][i][j] = convective_heat_transfer_coefficient * NNT[ip][i][j] * L / 2;
							input_element[iterator].H[i][j] += HBC_point[ip][i][j] * element.multiplier[ip];
						}

						P_point[ip][i] = -convective_heat_transfer_coefficient * N[ip][i] * ambient_temperature * L / 2;
						input_element[iterator].P[i] += P_point[ip][i] * element.multiplier[ip];

					}
				}
			}
		}
	}
}

void print_square_matrix(double** CG, int n_n, int style) {
	for (int i = 0; i < n_n; i++)
	{
		for (int j = 0; j < n_n; j++)
		{
			if (style == 1) cout << fixed << setfill(' ') << setw(5) << setprecision(3) << CG[i][j] << "\t";
			else if (style == 2) cout << setfill(' ') << setw(5) << setprecision(4) << CG[i][j] << "\t";
			else if (style == 3) cout << setfill(' ') << setw(5) << setprecision(2) << CG[i][j] << "\t";
		}
		cout << endl;
	}
}

double max(double* arr, int n)
{
	double max = arr[0];
	for (auto i = 1; i < n; i++)
		if (arr[i] > max)
			max = arr[i];
	return max;
}

double min(double* arr, int n)
{
	double min = arr[0];
	for (auto i = 1; i < n; i++)
		if (arr[i] < min)
			min = arr[i];
	return min;
}

double* solve_sof(int n_n, double** HR, double* PR)
{
	//CREATING TEMPORARY MATRIX TO MERGE HR MATRIX WITH PR VECTOR
	double **a;
	a = new double*[n_n];
	for (auto i = 0; i < n_n; ++i) a[i] = new double[n_n + 1];


	for (auto i = 0; i < n_n; i++)
	{
		for (auto j = 0; j < n_n + 1; j++)
		{
			if (j < n_n) a[i][j] = HR[i][j];
			else a[i][j] = PR[i];
		}
	}

	/*cout << "THE a MATRIX BEFORE PIVOTISATION" << endl;
	for (int i = 0; i < n_n; i++) {
		for (int j = 0; j < n_n + 1; j++) {
			if (j > n_n -1 ) cout << "| ";
			cout << a[i][j] << "   ";
		}
		cout << "\n";
	}*/

	//DEFINING VECTOR TO RETURN SOLUTION
	double *x;
	x = new double[n_n];

	//PIVOTISATION
	for (int i = 0; i < n_n; i++)
		for (int k = i; k < n_n; k++)
			if (a[i][i] < a[k][i])
				for (int j = 0; j <= n_n; j++) {
					double temp = a[i][j];
					a[i][j] = a[k][j];
					a[k][j] = temp;
				}

	//cout << "\nTHE a MATRIX AFTER PIVOTISATION:"<<endl;
	//for (int i = 0; i < n_n; i++) {
	//	for (int j = 0; j < n_n + 1; j++) {
	//		if (j > n_n - 1) cout << "| ";
	//		cout << a[i][j] << "   ";
	//	}
	//	cout << "\n";
	//}

	//PERFORMING GAUSS ELIMINATION
	for (int i = 0; i < n_n - 1; i++) // not going to the last row
		for (int k = i + 1; k < n_n; k++)
		{
			double mult = a[k][i] / a[i][i];
			for (int j = 0; j <= n_n; j++)
				a[k][j] = a[k][j] - mult * a[i][j]; //make the elements below the pivot equal to zero.
		}

	//cout << "\nTHE MATRIX AFTER GAUSS ELIMINATION:\n";
	//for (int i = 0; i < n_n; i++) {
	//	for (int j = 0; j < n_n + 1; j++) {
	//		if (j > n_n - 1) cout << "| ";
	//		cout << a[i][j] << "   ";
	//	}
	//	cout << "\n";
	//}

	for (int i = n_n - 1; i >= 0; i--)
	{
		x[i] = a[i][n_n];                //x[i] is now the right side of equation in line [i]
		for (int j = i + 1; j < n_n; j++)
			if (j != i)            //substract all left side numbers except the coefficient of the variable whose value is being calculated
				x[i] = x[i] - a[i][j] * x[j];
		x[i] = x[i] / a[i][i];            //divide the right side of equation by the coefficient of the variable to be calculated
	}

	//cout << "\nTHE SOLUTION VECTOR :\n";
	//for (int i = 0; i < n_n; i++) cout << x[i] << " ";
	//cout << endl;


	return x;
}

void generate_mesh()
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
	for (auto i = 0; i < gdata.n_n; i++)
	{
		if (ND[n_ND].x == 0 || ND[n_ND].y == 0 || ND[n_ND].y == gdata.h || ND[n_ND].x == gdata.w)
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

	//ASSIGNINGS PHYSICAL PROPERTIES TO ELEMENTS
	int m;
	for (int i = 1; i < gdata.n_e + 1; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (ND[Elem[i].id[j]].x < delta_x + 0.001 || ND[Elem[i].id[j]].x > gdata.w - delta_x - 0.001 || ND[Elem[i].id[j]].y < delta_y + 0.001 || ND[Elem[i].id[j]].y > gdata.h - delta_y - 0.001) m = 0;
			else m = 1;

			Elem[i].thermal_conductivity = gdata.thermal_conductivity[m];
			Elem[i].density = gdata.density[m];
			Elem[i].specific_heat = gdata.specific_heat[m];
		}
	}




	//DEFINING GLOBAL MATRICES.
	double PG[1000] = { 0.0 }; // TODO: dynamic sizing
	double PR[1000] = { 0.0 }; // replacement P for SOE

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

	double* t0;
	t0 = new double[gdata.n_n];
	for (auto i = 0; i < gdata.n_n; i++) t0[i] = ND[i + 1].t0;
	double * t1;



	int nt = gdata.simulation_time / gdata.simulation_step_time; // no. of iterations
	double* t_min;
	double* t_max;
	t_min = new double[nt];
	t_max = new double[nt];


	//HEAT TRANSFER PRCOESS SIMULATION: CALCULATIONS FOR EACH ITERATION

	for (int iteration_no = 0; iteration_no < nt; iteration_no++)
	{
		//cout << "ITERATION NO: " << iteration_no << "... _____________________---------------_____________--------__________---__--_-" << endl;
		// 1. CALCULATE H (WITHOUT BC) AND C MATRIX
		calculate_H(Elem, n_El, ND, gdata.order_of_integration);
		// 2. CALCULATE H(BC party only) AND SUM IT UP TO H, CALCULATE P MATRIX
		calculateHBC(Elem, gdata.order_of_integration, n_El, ND, gdata.convective_heat_transfer_coefficient, delta_x, delta_y, gdata.ambient_temperature);

		// 3. HG, CG, PG MATRICES AGGREGATION
		for (auto k = 1; k < gdata.n_e + 1; k++)
		{
			for (auto i = 0; i < 4; i++)
			{
				//cout << "P: elem nr: " <<k <<" "<< Elem[k].P[i] << endl;
				for (auto j = 0; j < 4; j++)
				{
					HG[Elem[k].id[i] - 1][Elem[k].id[j] - 1] += Elem[k].H[i][j];
					CG[Elem[k].id[i] - 1][Elem[k].id[j] - 1] += Elem[k].C[i][j];
				}
				PG[Elem[k].id[i]] += Elem[k].P[i];
			}
		}

		//cout << "PRINTING CG MATRIX: " << endl;
		//print_square_matrix(CG, gdata.n_n,2);
		//cout << "PRINTING HG MATRIX: " << endl;
		//print_square_matrix(HG, gdata.n_n,2);
		//cout << "PRINTING PG VECTOR: " << endl;
		//for (int j = 1; j < gdata.n_n+1; j++) cout << fixed << setprecision(2) << PG[j] << " ";
		//cout << endl;

		// 4. CALCULATE REPLACEMENT H (HR):
		for (auto i = 0; i < gdata.n_n; i++)
		{
			for (auto j = 0; j < gdata.n_n; j++)
			{
				HR[i][j] = HG[i][j] + CG[i][j] / gdata.simulation_step_time;
			}
		}

		// 5. CALCULATE REPLACEMENT P (PR)
		for (auto i = 0; i < gdata.n_n; i++)
		{
			PR[i] = -PG[i + 1];
			for (auto j = 0; j < gdata.n_n; j++)
			{
				PR[i] += (CG[i][j] / gdata.simulation_step_time) * t0[j];
			}
		}

		//cout << "PRINTING HR MATRIX: " << endl;
		//print_square_matrix(HR, gdata.n_n,1);
		//cout << "PRINTING PR VECTOR: " << endl;
		//for (int j = 0; j < gdata.n_n; j++) cout << fixed << setprecision(2) << PR[j] << " ";
		//cout << endl;

		//6. SOLVE THE SYSTEM OF LINEAR EQUATIONS: HR * {t} = PR
		t1 = solve_sof(gdata.n_n, HR, PR);

		double temp_t_max = max(t1, gdata.n_n);
		double temp_t_min = min(t1, gdata.n_n);
		//cout << "[" << (iteration_no+1)*gdata.simulation_step_time<<" sec] MAX: " << temp_t_max << "   MIN:" << temp_t_min << endl;

		cout << "[" << (iteration_no + 1)*gdata.simulation_step_time << " sec] " << t1[19] << " " << t1[64] << " " << t1[109] << " " << t1[154] << " " << t1[199] << " " << t1[244] << " " << t1[289] << " " << t1[334] << " " <<
			t1[379] << " " << t1[424] << " " << t1[469] << " " << t1[514] << " " << t1[559] << endl;


		t_min[iteration_no] = temp_t_min;
		t_max[iteration_no] = temp_t_max;


		t0 = t1;

		// 0.0 to the values:
		for (auto i = 0; i < gdata.n_n; i++)
		{
			for (auto j = 0; j < gdata.n_n; j++)
			{
				HG[i][j] = { 0.0 };
				CG[i][j] = { 0.0 };
			}
		}
		for (auto i = 0; i < gdata.n_n + 1; i++)
		{
			PG[i] = 0.0;
		}

		for (auto iterator = 1; iterator < n_El; iterator++)
		{
			for (auto i = 0; i < 4; i++)
			{
				for (auto j = 0; j < 4; j++)
				{
					Elem[iterator].C[i][j] = 0.0;
					Elem[iterator].H[i][j] = 0.0;

				}
				Elem[iterator].P[i] = 0.0;
			}
		}
	}
	int g = 0;
	for (int j = 0; j < nt; j++)
	{
		g += 50;
		cout << fixed << setprecision(2) << "   " << g << "  |";
	}
	cout << " [after seconds]" << endl;
	for (int j = 0; j < nt; j++) cout << fixed << setprecision(2) << t_min[j] << " | ";
	cout << " [t_min]" << endl;
	for (int j = 0; j < nt; j++) cout << fixed << setprecision(2) << t_max[j] << " | ";
	cout << " [t_max]" << endl;


}

