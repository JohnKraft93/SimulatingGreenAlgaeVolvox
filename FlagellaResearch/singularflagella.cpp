#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <math.h>  
#include <fstream>

#ifndef M_PI
	#define M_PI 3.14159265358979323846
#endif

using namespace std;

double const upperLambda = 0.1;
double const T = (1.0 / 33.0); //time scale
double const d = 10.0; //length scale


vector<vector<double>> A(double r0, double a, double z) {
	double pos = 1 + ((9 * a) / (16 * z));
	double pos3x3 = 1 + ((9 * a) / (8 * z));
	double val = pos * (2 * M_PI * r0);
	double val3x3 = pos3x3 * (2 * M_PI * r0);
	vector<vector<double>> retval = { {val, 0, 0},
								      {0, val, 0},
									  {0, 0, val3x3} };
	return retval;
}


vector<double> B(double r, double r0, vector<double> erhat, vector<double> enothat) {
	double retval = upperLambda * (r - r0);
	vector<double> tmp;
	
	for (double val : erhat) {
		tmp.push_back(val * retval);
	}
	vector<double> results;

	transform(tmp.begin(), tmp.end(), enothat.begin(), back_inserter(results), std::minus<double>());

	return results;
}


vector<vector<double>> inverse(vector<vector<double>> A) {
	double determinant = A[0][0] * (A[1][1] * A[2][2] - A[2][1] * A[1][2]) -
		A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
		A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
	vector<vector<double>> A_inverse;

	for (int i = 0; i < 3; i++) {
		vector<double> tmp;
		for (int j = 0; j < 3; j++) {
			if (i != j) {
				tmp.push_back(0.0);
			}
			else {
				tmp.push_back(A[(i + 1) % 3][(j + 1) % 3] * A[(i + 2) % 3][(j + 2) % 3]/determinant);
			}
		}
		A_inverse.push_back(tmp);
	}
	return A_inverse;
}


void printVectors(vector<double> erarrow, vector<double> erhat, vector<double> enotarrow, vector<double> enothat) {
	cout << "ERARROW: " << endl;
	for (double val : erarrow) {
		cout << val << " ";
	}
	cout << endl;

	cout << "ERHAT: " << endl;
	for (double val : erhat) {
		cout << val << " ";
	}
	cout << endl;

	cout << "ENOTARROW: " << endl;
	for (double val : enotarrow) {
		cout << val << " ";
	}
	cout << endl;
	cout << "ENOTHAT: " << endl;
	for (double val : enothat) {
		cout << val << " ";
	}
	cout << endl;
}


int main() {
	double r0s = 5.0; //preferred radius of the sphere's prbit
	double as = 1.0; //radius of sphere

	double r0 = r0s / d;
	double a = as / d;

	vector<double> center = { 0.0, 0.0, 1.0 };
	vector<double> position = { 0.0, 0.0, 1.6 };
	double z = position[2];

	vector<double> erArrow;
	vector<double> erHat;
	vector<double> enotArrow;
	vector<double> enotHat;
	vector<double> vArrow;
	vector<double> newVArrow;
	vector<double> initB;
	vector<double> newPos;
	vector<vector<double>> initA;
	vector<vector<double>> A_inverse;

	double radius;
	ofstream outfile;
	outfile.open("output.txt");
	

	for (double timeStep = (1.0 / 30.0); timeStep <= 1.0; timeStep += (1.0 / 30.0)) {
		//clear vectors
		initA.clear();
		A_inverse.clear();
		erArrow.clear();
		erHat.clear();
		enotArrow.clear();
		enotHat.clear();
		vArrow.clear();
		newVArrow.clear();
		initB.clear();
		newPos.clear();
		

		initA = A(r0, a, z);
		A_inverse = inverse(initA);

		transform(position.begin(), position.end(), center.begin(), back_inserter(erArrow), std::minus<double>());
		
		radius = 0;
		for (double val : erArrow) {
			radius += pow(val, 2);
		}
		radius = pow(radius, .5);

		for (double val : erArrow) {
			erHat.push_back(val / radius);
		}

		//calculating enotarrow/enothat
		enotArrow = { -erHat[2], 0, erHat[0] };

		double unit = 0;
		for (double val : enotArrow) {
			unit += pow(val, 2);
		}
		unit = pow(unit, .5);


		for (double val : enotArrow) {
			enotHat.push_back(val / unit);
		}

		initB = B(radius, r0, erHat, enotHat);

		//find V vector

		//A^-1 B = varrow

		vArrow = {A_inverse[0][0] * initB[0],
			A_inverse[1][1] * initB[1],
			A_inverse[2][2] * initB[2]};

		// multiply varrow by timestep
		
		for (double val : vArrow) {
			newVArrow.push_back(val * timeStep);
		}

		
		transform(position.begin(), position.end(), newVArrow.begin(), back_inserter(newPos), plus<double>());
		cout << "TIMESTEP: " << timeStep << endl;
		for (int i = 0; i < 3; i++) {
			cout << "POSITIONS : " << newPos[i] << " ";
		}
		cout << endl;

		for (int i = 0; i < 3; i++) {
			outfile << position[i];
			if (i < 2) { outfile << ","; }
		}
		outfile << endl;

		position = newPos;
		z = position[2];
		//printVectors(erarrow, erhat, enotarrow, enothat);
		/*
		if (timeStep > 1.0/30.0) {
			break;
		}
		*/
	}

	return 0;
}