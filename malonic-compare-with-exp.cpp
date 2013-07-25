#include "spin.hpp"
#include <iostream>
#include <fstream>
#include <vector>

double fi = -14.5_MHz;
double fe = -28.3775_MHz;
constexpr int N = 8;
double Azx[N] = { -5.8293_MHz,1.0625_MHz,-2.5745_MHz,0.1657_MHz,1.2042_MHz,0.0968_MHz,0.7024_MHz,-2.6575_MHz };
double Azy[N] = { 0.7154_MHz,0.4234_MHz,-0.5842_MHz,-1.8497_MHz,-1.6795_MHz,-1.9991_MHz,-0.6635_MHz,-1.9556_MHz };
double Azz[N] = { -56.3421_MHz,-3.9900_MHz,-1.3463_MHz,-0.5889_MHz,-0.5885_MHz,1.3987_MHz,-1.2610_MHz,7.4265_MHz };

int main() {
	int n = 2000;
	double T = 2000_ns;
	double delta_t = T/n;
	vector<double> f1_list = { 7.445_MHz,9.236_MHz,10.271_MHz,16.113_MHz,11.719_MHz,20.509_MHz };
	for(auto k=f1_list.begin();k!=f1_list.end();k++) {
		/* Hamiltonian */
		double f1 = *k;
		Operator H = 2*pi*(fe*Sz(0)+f1*Sx(0));
		for(int i=0;i<N;i++)
			H += 2*pi*(fi*Sz(i+1)+Sz(0)*(Azx[i]*Sx(i+1)+Azy[i]*Sy(i+1)+Azz[i]*Sz(i+1)));
		auto U = H.U();
		/* initial state */
		double alpha_e = 7.6216e-4;
		double alpha_i = -1.1579e-6;
		Operator rho0 = Op<2>(0,0.5*(1-alpha_e),0,0,0.5*(1+alpha_e));
		//Operator rho0 = Op<2>(0,0.5,0.5,0.5,0.5);
		for(int i=1;i<=N;i++)
			rho0 *= Op<2>(i,0.5*(1-alpha_i),0,0,0.5*(1+alpha_i));
		/* output stream */
		stringstream fnstream;
		fnstream << "malonic_" << f1/1_MHz << "MHz.txt";
		string fn = fnstream.str();
		ofstream out(fn);
		cout << fn << ":" << endl;
		/* calculate and output */
		for(int i=0;i<=n;i++) {
			double t = delta_t*i;
			cout << fn << ":\t" << t/T*100 << "%" << endl;
			out << t/1_ns << '\t';
			Operator rho = U(t)*rho0*U(-t);
			Operator rhoe = rho;
			for(int j=1;j<=N;j++)
				rhoe = rhoe.tr(j);
			out << real(tr(rhoe*Sz(0))) << endl;
		}
		out.close();
	}
}