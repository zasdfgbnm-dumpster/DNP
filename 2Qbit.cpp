#include "spin.hpp"
#include <fstream>
#include <iostream>
#include <sstream>

/* electron resonance with a H nuclear
 * H = w1*Sx1 + w2*Sz2 + Sz1*(A*Sx2 + B*Sy2)
 */

int main() {
	double delta_f = -0.5_MHz;
	double f2 = -14.5_MHz;
	double A = -5.8293_MHz;
	double B = 0.7154_MHz;
	double T = 2000_ns;
	int n = 20000;
	double delta_t = T/n;
	for(int i=0;i<60;i++) {
		double f1 = delta_f*i;		
		Operator H = 2*pi*(f1*Sx(0)+f2*Sz(1)+Sz(0)*(A*Sx(1)+B*Sy(1)));
		auto U = H.U();
		
		double alpha_i = -1.1579e-6;
		Operator rho0 = Op<2>(0,0.5,0.5,0.5,0.5)*Op<2>(1,0.5*(1-alpha_i),0,0,0.5*(1+alpha_i));
		//Operator rho0 = Op<2>(0,0.5,0.5,0.5,0.5)*Op<2>(1,0,0,0,1);
		
		stringstream fnstream;
		fnstream << "2Qbit_" << f1/1_MHz << "MHz.txt";
		string fn = fnstream.str();
		ofstream out(fn);
		
		for(int j=0;j<=n;j++) {
			double t = delta_t*j;
			Operator rho = U(t)*rho0*U(-t);
			Operator rhoe = rho.tr(1);
			Operator rhoi = rho.tr(0);
			out << t/1_ns << "\t" << real(tr(rhoe*Sx(0))) << "\t" << real(tr(rhoi*Sz(1))) << endl;
		}
		out.close();
		cout << fn << endl;
	}
}