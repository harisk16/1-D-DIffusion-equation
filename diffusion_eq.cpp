#include<iostream>
#include<math.h>
#include <fstream>
#include <string>
using namespace std;
int main()
{
	// Define Grid points and time step //
	int Nx= 128; // Number of grid points
	double dt= 0.00001; // Time step 
	double Xmax = 1;
	double Xmin = 0;
	int Tmax=1001;
	
	// Discritize mesh & Time
	double Mesh[Nx];
	double Time[Tmax];
	double dx;
	dx = (Xmax - Xmin)/(Nx-1) ;
	for(int i=0; i<Nx; i++){
		Mesh[i] = Xmin + i*dx; 	
	}
	for(int j=0; j<Tmax; j++){
		Time[j] = j*dt; 	
	}
	
	// Define variable array
	double U[Nx][Tmax]={ }; // Initialize all elements to zero 
	double Uexact[Nx][Tmax]={ };
	
	// Initial Conditions and exact solution
	for(int i=0; i<Nx; i++){
		double x;
		x = Mesh[i];
		U[i][0] = sin(2*M_PI*x);
		for(int j=0; j<Tmax; j++){
			double t;
			t = Time[j];
			Uexact[i][j] = sin(2*M_PI*x)*exp(-4*M_PI*M_PI*t);		
		}
		
	}
	
	// Solver 
	for(int j=0; j < Tmax; j++){
		double Error = 0;
		for(int i=1; i< Nx-1; i++){
			U[Nx-1][j] = U[0][j]; // Periodic Boundary Condition
			U[i][j+1] = U[i][j] + dt*( U[i+1][j] -2*U[i][j] + U[i-1][j])/(2*dx*dx);
			double A[Nx];
			if(j == 0){
				std::ofstream File_1;
				File_1.open("t_0s.csv");
				for(int p=0; p<Nx; p++){
					A[p] = U[p][j];
					File_1 << Mesh[p] << "," << A[p] << std::endl;
				}
			}
			if(j == 99){
				std::ofstream File_2;
				File_2.open("t_100s.csv");
				for(int p=0; p<Nx; p++){
					A[p] = U[p][j];
					File_2 << Mesh[p] << "," << A[p] << std::endl;
				}
			}
			if(j == 499){
				std::ofstream File_3;
				File_3.open("t_500s.csv");
				for(int p=0; p<Nx; p++){
					A[p] = U[p][j];
					File_3 << Mesh[p] << "," << A[p] << std::endl;
				}
			}
			if(j == 998){
				std::ofstream File_4;
				File_4.open("t_1000s.csv");
				for(int p=0; p<Nx; p++){
					A[p] = U[p][j+1];
					File_4 << Mesh[p] << "," << A[p] << std::endl;
				}
			}
			
			 	Error += Uexact[i][j] - U[i][j] ; // Algebraic average is taken in this case 
			 	Error = Error/Nx;
			 	
			}
			double Er[Tmax];
			Er[j] = Error;  
			cout << "Error at timestep " << Time[j] << " is " << Error << endl;
			std::ofstream File_5;
			File_5.open("Error.csv");
			for(int p=0; p<Tmax; p++){
			File_5 << Time[p] << "," << Er[p] << std::endl;
		    }
}
return 0;
}
