/* IMPLEMENTATION OF THE SVJ MODEL FOR A EUROCALL OPTION 6/6/2011 */

#include<iostream>
#include<complex>
#include<iomanip>
#include<boost/function.hpp>
#include<boost/bind.hpp>
#include<boost/mem_fn.hpp>

typedef std::complex<double> dcomp;
double const PI = 3.141592653589793;
dcomp const i(0,1.0); //Complex i
double SimpsonsRule(boost::function<double (double u)> f, double a, double b, int n);
class Euro_Call_Option;

using namespace std;

class Euro_Call_Option
{
/*	This class represents a Euro option with stochastic volatility and jumps
	in the underlying asset only.  The implementation is based on the 
	closed form solution provided in "The Volatility Surface", by Jim
	Gatheral.
*/
	private:
/*		Option paramteters: K is the strike price, r is the interest rate,
		v_bar is the long-term vol, rho is the correlation of the Brownian
		Motions, eta is the vol of the vol and lambda is the speed of 
		reversion of the vol.  The remainding parameters are for the 
		lognormally distributed jump process: alpha_J is the mean log-jump,
		delta_J is the standard deviation and lambda_J is the jump intensity.
*/
		double K, r, v_bar, rho, eta, lambda, alpha_J, delta_J, lambda_J;

/*		Integrand returns the integrand required to compute option prices from the
		characteristic function
*/
		double Integrand(double u, double k, double v, double T);
public:
//		Option constructor. Initiate all the option parameters.
		Euro_Call_Option::Euro_Call_Option(double xK, double xv_bar=0, double xeta=0, 
				double xrho=0, double xl=0,	double xalpha_J=0, double xdelta_J=0, double xl_J=0):
					K(xK), v_bar(xv_bar), eta(xeta), rho(xrho), lambda(xl), 
						alpha_J(xalpha_J), delta_J(xdelta_J), lambda_J(xl_J){};

/*		phi returns the characteristic function of the SVJ process at the value
		u - I*complex_shift.	Implemented such a way so that the characteristic 
		function is independant of the integration routine. This function is declare
		as a public function for future use.	*/
		dcomp phi(double u, double complex_shift, double v, double T);

/*		Value returns the value of a Euro call option with x=ln(S/E),
		v = current vol and T = time to expiry.	*/
		double Value(double S, double v, double T);
};

dcomp Euro_Call_Option::phi(double u, double complex_shift, double v, double T)
{
/*	This function computes the characteristic function at u-i*complex_shift required 
	for the valuation formula.  Most of the definitions come directly from the 
	"The Volatility Surface", pages 16-21 and page 66
	(except evaluated at u-i*complex_shift rather than u). Formula (4) in report.
*/
	dcomp u_shift(u,-complex_shift);

	dcomp alpha= -u_shift*u_shift*0.5-u_shift*i*0.5;
	dcomp beta = lambda-rho*eta*u_shift*i;
	double gamma = (eta*eta)/2;
	dcomp d = sqrt((beta*beta)-4.0*alpha*gamma);
	dcomp rp =0.5*(beta+d)/gamma;
	dcomp rm =0.5*(beta-d)/gamma;

	dcomp g = (rm)/(rp);
	dcomp D = rm*(1.0-exp(-d*T))/(1.0-g*exp(-d*T));
	dcomp C = lambda*(rm*T-(1/gamma)*log((1.0-g*exp(-d*T))/(1.0-g)));

/*	psi(u) is the jump contribution as defined on pg 66 of "The Volatility Surface". 
*///
	dcomp complex_exponent = i*u_shift*alpha_J-u_shift*u_shift*delta_J*delta_J*0.5;
	dcomp psi = exp(complex_exponent)-1.0 - u_shift * i * (exp(alpha_J+delta_J*delta_J*0.5)-1.0);
	//cout << "CHAR: " << psi <<endl;
	return exp(C*v_bar+D*v+psi*lambda_J*T);
}

double Euro_Call_Option::Value(double S, double v, double T)
{
//	Valuation function. Formula (3) in report.
	boost::function<double (double u)> f;
	f = boost::bind(boost::mem_fn(&Euro_Call_Option::Integrand),this,_1,log(K/S),v,T);
	//cout<< "f(0): " << f(0.876) << endl;
	//cout<< "Integrand(0,log(S/K),v,T): " << Integrand(0.876,log(S/K),v,T)<< endl;
	return S-sqrt(S*K)*(1/PI)*(SimpsonsRule(f, 0.0, 1.0,10000)+SimpsonsRule(f,1.0,100000.0,5000000));
}

double Euro_Call_Option::Integrand(double u, double k, double v, double T)
{
	dcomp complexIntegrand =  exp(-i*u*k)*phi(u,0.5,v,T)/(u*u+0.25);
	return complexIntegrand.real();
}

double SimpsonsRule(boost::function<double (double u)> f, double a, double b, int n)
{
//	Standard Simpons rule for numerical integration. Formula (5) in report.
	double h = (b-a)/n;
	double runningSum = 0;
	int j = 1;
	
	runningSum = f(a);
	if (n%2==1) n+=1;

	for(j; j<n; j++){
		(j % 2 ==1) ? runningSum+=4*f(a+j*h) : runningSum+=2*f(a+j*h);
	}
	runningSum+= f(b);

	return (h/3)*runningSum;
}

int main(int argc, char *argv[])
{
	//	Main function simply constructs an option, with parameters as below,
	//	prompts the user the strike and outputs the option value.

	double K,vbar = 0.04, lambda = 1.15, eta = 0.39, rho = -0.64,u=1;
	double S0 = 1, v0=0.04, T=1, alpha_J=-0.1151, delta_J=0.0967, lambda_J=0.1308;

	printf("Enter the strike: \n");
	scanf("%lf",&K);
	
	cout << setprecision(6);
	Euro_Call_Option myOption(K,vbar,eta,rho,lambda,alpha_J,delta_J,lambda_J);
	
	cout<<"The option value is: " << myOption.Value(S0,v0,T)<<endl;

	return EXIT_SUCCESS;
}