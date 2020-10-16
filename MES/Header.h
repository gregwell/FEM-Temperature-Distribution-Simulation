#pragma once
#include <string>
#include <sstream>

struct node {
	double x; //coordinates
	double y;
};

struct element {
	double id[4]; //id of all four nodes from each side
};

struct global_data {
public:
	double w; //width
	double h; //height
	double n_w; //number of nodes in width length
	double n_h; //number of nodes in height length
	double n_n; // m_n = m_w * m_h // number of nodes at all
	double n_e; // m_e = (m_w-1)*(m_h-1) // number of elements?
public:
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
