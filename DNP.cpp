#include "spin.hpp"
#include <iostream>
#include <fstream>

double fi = -14.5_MHz;
double fe = -28.3775_MHz;
double f1 = 14.5_MHz;
double t_half_pi = (pi/2)/(2*pi*f1);
double t_cycle = 200_ns;
constexpr int N = 8;
int N_cycle = 100;
double Azx[N] = { -5.8293_MHz,1.0625_MHz,-2.5745_MHz,0.1657_MHz,1.2042_MHz,0.0968_MHz,0.7024_MHz,-2.6575_MHz };
double Azy[N] = { 0.7154_MHz,0.4234_MHz,-0.5842_MHz,-1.8497_MHz,-1.6795_MHz,-1.9991_MHz,-0.6635_MHz,-1.9556_MHz };
double Azz[N] = { -56.3421_MHz,-3.9900_MHz,-1.3463_MHz,-0.5889_MHz,-0.5885_MHz,1.3987_MHz,-1.2610_MHz,7.4265_MHz };

void output_data(double t,Operator &rho,ofstream &out) {
	cout << "time: " << t/1_ns << "ns" << endl;
	out << t/1_ns << '\t';
	/* electron */
	Operator rhoe = rho;
	for(int j=1;j<=N;j++)
		rhoe = rhoe.tr(j);
	out << real(tr(rhoe*Sx(0))) << '\t';
	out << real(tr(rhoe*Sy(0))) << '\t';
	out << real(tr(rhoe*Sz(0))) << '\t';
	/* nuclei */
	for(int j=1;j<=N;j++) {
		Operator rhoj = rho;
		for(int k=0;k<=N;k++) {
			if(k==j)
				continue;
			rhoj = rhoj.tr(k);
		}
		out << real(tr(rhoj*Sz(j))) << '\t';
	}
	out << endl;
}

int main() {
	//output
	ofstream out("DNP.txt");
	cout << "time of pi/2 pulse: " << t_half_pi/1_ns << "ns" << endl;
	
	//Static Hamiltonian
	Operator Hstatic = 2*pi*fe*Sz(0);
	for(int i=0;i<N;i++)
		Hstatic += 2*pi*(fi*Sz(i+1)+Sz(0)*(Azx[i]*Sx(i+1)+Azy[i]*Sy(i+1)+Azz[i]*Sz(i+1)));
	auto Ustatic = Hstatic.U();
	
	//x nutation
	Operator Hx = 2*pi*f1*Sx(0) + Hstatic;
	auto Ux = Hx.U();
	
	//y nutation
	Operator Hy = 2*pi*f1*Sy(0) + Hstatic;
	auto Uy = Hy.U();
	
	//Initial state is thermal equilibrium
	double alpha_e = 7.6216e-4;
	double alpha_i = -1.1579e-6;
	//Operator rho0e = Op<2>(0,0.5*(1-alpha_e),0,0,0.5*(1+alpha_e));
	Operator rho0e = Op<2>(0,0.5,0.5,0.5,0.5);
	Operator rho0 = rho0e;
	for(int i=1;i<=N;i++)
		rho0 *= Op<2>(i,0.5*(1-alpha_i),0,0,0.5*(1+alpha_i));
	
	//do cycle
	for(int cycle=0;cycle<100;cycle++) {
		cout << "cycle: " << cycle << endl;
		double t = 0;
		//apply an x direction pi/2 pulse
		while(t<t_half_pi){
			Operator rho = Ux(t)*rho0*Ux(-t);
			output_data(t+t_cycle*cycle,rho,out);
			t += 1_ns;
		}
		t += 1_ns;
		Operator rho1 = Ux(t_half_pi)*rho0*Ux(-t_half_pi);
	
		//apply an y direction spin lock
		while(t<t_cycle) {
			Operator rho = Uy(t-t_half_pi)*rho1*Uy(-(t-t_half_pi));
			output_data(t+t_cycle*cycle,rho,out);
			t += 1_ns;
		}
		
		//reset the electron's state to thermal equilibrium
		Operator rho_cycle_end = Uy(t_cycle-t_half_pi)*rho1*Uy(-(t_cycle-t_half_pi));
		rho0 = rho0e*rho_cycle_end.tr(0);
	}
	
	//clean up
	out.close();
	return 0;
}