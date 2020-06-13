//////////////////////////////////////////////////
//  Aviraj Singh Bevli, 18MA20009			         	//
//	Solving differential equations of the form	//
// 	y' = f(x,y)								                 	//
//	using the Backward Euler Method			      	//
//////////////////////////////////////////////////

#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<cmath>
using namespace std;

long double e = 2.718281828459045235360287471352;

//This is the accuracy required to be acheived
//before terminating the Newton Raphson method
long double NewtonRaphson_accuracy = 0.000000000001;

//This is the step size
//Constant to be changed
long double h = 0.1;

//The question you intend to solve
//Constant to be changed
int question_to_solve = 3;

long double func(long double x, long double y){
	if(question_to_solve == 1)
		return x*x*y;
	else if(question_to_solve == 2)
		return x/y;
	else if(question_to_solve == 3)
		return x-y;
}

long double funcAnalytical(long double x){
	if(question_to_solve == 1)
		return pow( e, (pow(x,3)/3) );
	else if(question_to_solve == 2)
		return sqrt( pow(x,2) + 1 );
	else if(question_to_solve == 3)
		return ( 2*pow(e,-1*x) + x - 1 );
}

//The function defined for the Newton Raphson method
long double F(long double x1, long double y0, long double y1){
	long double temp = y1 - ( y0 + h*func(x1,y1) );
	return temp;
}

//Takes in the given value of y(0) and sets it as the initial value
//Takes in the value of x_n+1
long double NewtonRaphson(long double x1, long double y0){
	long double y_curr = y0;
	long double y_prev = y_curr - 1;
	long double dy_prev = 0.0000000000001;
	int i=0;
	while(y_curr - y_prev > NewtonRaphson_accuracy ){
		y_prev = y_curr;	
		long double Fd = ( F(x1,y0,y_prev+dy_prev) - F(x1,y0,y_prev) )/dy_prev ; 
		y_curr = y_prev - ( F(x1,y0,y_prev)/Fd );
	}
	return y_curr;
}

int main()
{
	int num_steps = 1+(1/h);
	long double y0 = 1;
	long double y1 = 1;
	long double y1_analytical = 1;
	long double x1=0;
	for(int i=0;i<num_steps;i++){
		x1+=h;
		y0 = y1;
		y1 = NewtonRaphson(x1,y0);
		y1_analytical = funcAnalytical(x1);
		printf("x1: %0.15Lf    y1: %0.15Lf   y1_analytical: %0.15Lf   Error: %0.15Lf \n", x1, y1, y1_analytical, y1 - y1_analytical);
	}
	return 0;
}
