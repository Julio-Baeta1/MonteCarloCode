// MontePdf.cpp 

#include <iostream>
#include <armadillo>
#include <cmath>

//Used for easier reading
using namespace arma;

// A simple function to calculate area under straight line y=x using monte carlo method
double monteStraightLine(double& lower, double& upper, int& n, vec& coeffs);

// Use transform density to uniform density to calculate area of y=x
double monteUniformArea(double& lower, double& upper, int& n);

//double polyYvalue(double& x, vec& coeffs);
int interval_contained(vec& partition, double& val);


int main()
{
	int n{ 100 }; //10000
	double low{ 0 }, up{ 1 }, ans{0};
	vec straight_no_const("1 0");
	

	//// Make hypothesis test to check number of points is sufficient
	for (int i = 0; i < 2; i++)
	{
		ans = monteStraightLine(low, up, n, straight_no_const);
		cout << ans << endl;
	}

	//cout << endl << "Alternative test:" << endl;
	//for (int i = 0; i < 10; i++)
	//{
	//	ans = monteUniformArea(low, up, n);
	//	cout << ans << endl;
	//}

	return 0;
}

int interval_contained(vec& partition, double& val)
////test code
//vec x = linspace(low, up, 9);
//double val = 0.5;
//int indx = interval_contained(x, val);
//x.print("x:");
//cout << "val: " << val << endl << "index: " << indx << endl;
{
	for (int i = 1; i < partition.n_elem+1; i++)
	{
		if (val <= partition.at(i))
			return i;
	}
	
	return 0;
	// make code for when not in partition set
}

double monteStraightLine(double& lower, double& upper, int& n, vec& coeffs)
//// coeffs [x^n, x^(n-1), ..., x^2, x, 1]
// Add for lower =/= 0
// Add for additional constant y = x + b
// Add for different slope y = ax
{

	//For lower = 0
	double area = (upper-lower)* (upper - lower); 
	double under_count{ 0 };
	int indx{ 0 };

	vec x = linspace(lower, upper, n);
	vec y = polyval(coeffs, x);

	mat A(n,2, fill::randu); 
	A = A * (upper - lower);

	cout << "start" << endl;
	for (int i = 0; i < n; i++)
	{
		indx = interval_contained(x, A.at(i, 0));

		if (A.at(i, 1) <= y.at(indx))
			under_count++;
	}
	cout << "end" << endl;

	//for (int i = 0; i < n; i++)
	//{
		//if (A.at(i, 0) > A.at(i, 1)) // Since y=x
		//	under++;
	//}

	area = area * under_count / n;

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