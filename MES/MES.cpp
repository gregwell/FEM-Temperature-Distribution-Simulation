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

	double n_ND = 1;
	delta_x = gdata.w / (gdata.n_w - 1);
	delta_y = gdata.h / (gdata.n_h - 1);

	gdata.n_n = gdata.n_w*gdata.n_h;
	gdata.n_e = (gdata.n_w - 1)*(gdata.n_h - 1);
	
	
	for (auto i=1; i<gdata.n_h; i++)
	{
		for(auto j=1; j<gdata.n_w;j++)
		{
			node *ND = new node[gdata.n_n];
			ND[j].x = j * delta_x;
			ND[j].x = j * delta_x;
			n_ND++;
		}
	}
	
	

	system("pause");
	return 0;
}
