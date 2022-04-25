// MontePdf.cpp 

#include <iostream>
#include <armadillo>
#include <cmath>

//Used for easier reading
using namespace arma;

// A simple function to calculate area under straight line y=x using monte carlo method
double monteStraightLine(double& lower, double& upper, int& n, vec& coeffs);

// A simple function to calculate area under straight line y=x using monte carlo method
double monteStraightLineIntervalMethod(double& lower, double& upper, int& n, vec& coeffs);

// Use transform density to uniform density to calculate area of y=x
double monteUniformArea(double& lower, double& upper, int& n);

//double polyYvalue(double& x, vec& coeffs);
int interval_contained(vec& partition, double& val);


int main()
{
	int n{ 10000 }; //10000
	double low{ 0 }, up{ 1 }, ans{ 0 }, tot{0};
	vec straight_no_const("1 1");
	

	// Make hypothesis test to check number of points is sufficient
	for (int i = 0; i < 10; i++)
	{
		ans = monteStraightLine(low, up, n, straight_no_const);
		cout << ans << endl;
		tot += ans;
	}
	cout << tot / 10 << endl;


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
// Add for area below 0 option 
// Check for higher order polynomial
{

	//For lower = 0
	vec x = linspace(lower, upper, n);
	vec y = polyval(coeffs, x);

	double area = y.max() * x.max(); //(upper-lower)* (upper - lower); 
	double under_count{ 0 }, y_val{0};

	mat A(n,2, fill::randu); 
	A.col(1) *= y.max();
	A.col(0) *= x.max();

	vec y_for_Ax = polyval(coeffs, A.col(0));

	for (int i = 0; i < n; i++)
	{
		if (A.at(i,1) <= y_for_Ax.at(i))
			under_count++;
	}

	//for (int i = 0; i < n; i++)
	//{
		//if (A.at(i, 0) > A.at(i, 1)) // Since y=x
		//	under++;
	//}

	area = area * under_count / n;

	return area;
}

double monteStraightLineIntervalMethod(double& lower, double& upper, int& n, vec& coeffs)
//// coeffs [x^n, x^(n-1), ..., x^2, x, 1]
// Add for lower =/= 0
// Add for area below 0 option 
// Check for higher order polynomial
{

	//For lower = 0
	vec x = linspace(lower, upper, n);
	vec y = polyval(coeffs, x);

	double area = y.max() * x.max(); //(upper-lower)* (upper - lower); 
	double under_count{ 0 }, a{ 0 }, b{ 0 };
	int indx{ 0 };

	mat A(n, 2, fill::randu);
	A.col(1) *= y.max();
	A.col(0) *= x.max();
	//A = A * (upper - lower);

	for (int i = 0; i < n; i++)
	{
		indx = interval_contained(x, A.at(i, 0));
		//if (A.at(i, 1) <= y.at(indx)) //Assume y value constant in interval
		a = (y.at(indx) - y.at(indx - 1)) / (x.at(indx) - x.at(indx - 1));
		b = y.at(indx) - a * x.at(indx);
		if (A.at(i, 1) <= a * A.at(i, 0) + b)
			under_count++;
	}

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