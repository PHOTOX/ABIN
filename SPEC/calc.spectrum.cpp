#include <stdio.h>
#include <math.h>
#include <iostream>
using namespace std;

const double  pi=4*atan(1.0);
int MyRound( double value )  {
	return floor( value + 0.5 );
	}

int Prevod(int nsample,double aint,double E);

 int main(int argc, char* argv[])
 {	
	int nsample,nexc,l;
	double t1,t2,t3,exc,trans;
	double prop[1000];
	for (int i=0;i<=801;i++) {
		prop[i]=0;
	}
	cin >> nsample;
	cin >> nexc;
	for (int j=1;j<=nsample;j++) {
		for (int k=1;k<=nexc;k++) {
                        // excitation energies in electronvolts
			cin >> exc;
			if ( exc > 10.0 ) {
				cout << "Excitation greater than 10 eV! Modify source code accordingly. Exiting..." << endl;
				return 1;
			}
                        // transition electric dipole moment in atomic units
			cin >> t1 >> t2 >> t3;
			trans=t1*t1+t2*t2+t3*t3;
			l=MyRound((10.-exc)/0.02);
			if (l >= 0 && l <=500) {
				prop[l]=prop[l]+trans;
			}
			t1=t2=t3=0;
		}
	}
	for (int i=0;i<=499;i++) {
	Prevod(nsample,prop[i],10.00-0.02*i);
	}
	return 0;
}
 
 int Prevod(int nsample,double aint,double E) { 
 const double  deltaE=0.02;
 const double  evtoj=1.602e-19;
 const double  eps=8.854e-12;
 const double  hprime=6.626e-34/(2*pi);
 const double  c=299792e3;
 const double  deb=2.5*3.34e-30;
       double sigma,cross_sec;

       aint=aint/nsample;
       sigma=(pi*E*aint)/(3*hprime*eps*c*deltaE);
       cross_sec=sigma*deb*deb*1e4;
       double nm=1239.8/E;
       if (nm < 1000) {
	       cout << nm << "\t" << cross_sec << endl;
       }
       return 0;
 }
