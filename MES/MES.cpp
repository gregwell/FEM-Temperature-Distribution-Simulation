#include "pch.h"
#include <iostream>
#include "Header.h"
#include <fstream>


using namespace std;

double delta_x;
double delta_y;

int main()
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
	delta_x = gdata.w / (gdata.n_w - 1); //x-axis element length
	delta_y = gdata.h / (gdata.n_h - 1); //y-axis element length
	gdata.n_n = gdata.n_w*gdata.n_h; //number of nodes
	gdata.n_e = (gdata.n_w - 1)*(gdata.n_h - 1); //number of elements

	//Assigning coordinates (x,y) to every node.
	node *ND = new node[gdata.n_n+1];
	
	for (auto i=0; i<gdata.n_w; i++) //width
	{
		for(auto j=0; j<gdata.n_h;j++) //height
		{
			ND[n_ND].x = i * delta_x;
			ND[n_ND].y = j * delta_y;
			//cout << "Jestem w etapie: n_Nd=" << n_ND << "/ ND[n_ND].x = " << i * delta_x << "/ ND[n_ND].y=" << j * delta_y << endl;
			n_ND++;

		}
	}

	//Assigning nodes to elements.
	int n_El = 1; //node array number
	int k = 0;
	element *Elem = new element[gdata.n_e+1];

	for (auto i = 0; i < gdata.n_w-1; i++) //width
	{
		for (auto j = 0; j < gdata.n_h-1; j++) //height
		{
			Elem[n_El].id[0] = n_El+k;
			Elem[n_El].id[1] = n_El + gdata.n_h+k;
			Elem[n_El].id[2] = n_El + gdata.n_h + 1 + k;
			Elem[n_El].id[3] = n_El + 1 + k;
			//cout << "Jestem w etapie: n_el=" << n_El << "/ id1=  " << n_El + k << "/// id2= " << n_El + gdata.n_h + k << endl;
			n_El++;
		}
		k++;
	}

	

	system("pause");
	return 0;
}
