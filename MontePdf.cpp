// MontePdf.cpp 

#include <iostream>
#include <armadillo>
#include <cmath>

//Used for easier reading
using namespace arma;

// A simple function to calculate area under straight line y=x using monte carlo method
double monteStraightLine(double& lower, double& upper, int& n);

// Use transform density to uniform density to calculate area of y=x
double monteUniformArea(double& lower, double& upper, int& n);


int main()
{
	int n{ 100000 };
	double low{ 0 }, up{ 16 }, ans{0};
	
	// Make hypothesis test to check number of points is sufficient
	for (int i = 0; i < 10; i++)
	{
		ans = monteStraightLine(low, up, n);
		cout << ans << endl;
	}

	cout << endl << "Alternative test:" << endl;
	for (int i = 0; i < 10; i++)
	{
		ans = monteUniformArea(low, up, n);
		cout << ans << endl;
	}

	return 0;
}

double monteStraightLine(double& lower, double& upper, int& n)
// Add for lower =/= 0
// Add for additional constant y = x + b
// Add for different slope y = ax
{
	
	double area = (upper-lower)* (upper - lower); //first test for lower = 0
	mat A(n,2, fill::randu); 
	A = A * (upper - lower);
	double under{ 0 };

	for (int i = 0; i < n; i++)
	{
		if (A.at(i, 0) > A.at(i, 1)) // Since y=x
			under++;
	}

	area = area * under / n;

	return area;
}

double monteUniformArea(double& lower, double& upper, int& n)
{
	double uni_value = upper - lower; 
	uni_value *= uni_value;
	vec A(n, fill::randu);
	double area = accu(A) * uni_value/n;

	return area;
}