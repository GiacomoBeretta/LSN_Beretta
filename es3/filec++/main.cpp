#include <iostream>
#include <fstream>
#include "random.hpp"
#include <cmath> //exp, erf, log
#include <algorithm> //max
#include <numeric>
#include <vector>

using namespace std;

/*per calcolare esattamente i valori di call e put option
double N_x(double x);
double d1(double T_f, double t, double r, double sigma, double S_t, double k);//T_f è l'expiry date, t è il tempo, r è l'interest rate, sigma è la volatilità, S_t è il prezzo al tempo t, k è il prezzo concordato
double d2(double d1, double sigma, double T_f, double t);
double EuroCallOpt(double S_t, double N_d1, double k, double r, double T_f, double t, double N_d2);
double EuroPutOpt(double S_t, double N_d1, double k, double r, double T_f, double t, double N_d2);*/
double S_t_MC(double S_t1, double mu, double sigma, double t2, double t1, double z);
double S_t_MC_discrete(double S_t1, double mu, double sigma, double t_f, double t_0, Random rnd, int N_int ); //N_int è il numero di intervalli
double EuroCallOpt_MC(double S_t, double k);
double EuroPutOpt_MC(double S_t, double k);

int main (){
	
	Random rnd;
	rnd.SetRandomPrimes();

	int N = 100000; //number of throws
	int M = 1000; //number of blocks
	int L = N/M;
	

	double t_0 = 0;
	double S_0 = 100; //asset price at t=0
	double t_f = 1; //delivery time
	double k = 100; //strike price
	double r = 0.1; //risk-free interest rate
	double sigma = 0.25; //volatility

	double discount = exp( -r * ( t_f - t_0 ) );


	//variables for call option (direct method)
	double ave1_c_dir;
	double ave1_squared_c_dir;
	double ave2_c_dir;
	double ave2_squared_c_dir;
	double StDevMean_c_dir;
	double sum1_c_dir;
	double sum2_c_dir = 0;
	double sum3_c_dir = 0;

	//variables for put option (direct method)
	double ave1_p_dir;
	double ave1_squared_p_dir;
	double ave2_p_dir;
	double ave2_squared_p_dir;
	double StDevMean_p_dir;
	double sum1_p_dir;
	double sum2_p_dir = 0;
	double sum3_p_dir = 0;



//variables for call option (discrete method)
	int N_int = 100;

	double ave1_c_dis;
	double ave1_squared_c_dis;
	double ave2_c_dis;
	double ave2_squared_c_dis;
	double StDevMean_c_dis;
	double sum1_c_dis;
	double sum2_c_dis = 0;
	double sum3_c_dis = 0;

	//variables for put option (discrete method)
	double ave1_p_dis;
	double ave1_squared_p_dis;
	double ave2_p_dis;
	double ave2_squared_p_dis;
	double StDevMean_p_dis;
	double sum1_p_dis;
	double sum2_p_dis = 0;
	double sum3_p_dis = 0;	

	ofstream output[5];
	output[0].open("/mnt/c/Users/giaco/esercizi_python/3/N_and_M_values.txt");
	output[1].open("/mnt/c/Users/giaco/esercizi_python/3/call_option_direct.txt");
	output[2].open("/mnt/c/Users/giaco/esercizi_python/3/put_option_direct.txt");
	output[3].open("/mnt/c/Users/giaco/esercizi_python/3/call_option_discrete.txt");
	output[4].open("/mnt/c/Users/giaco/esercizi_python/3/put_option_discrete.txt");

	output[0] << N << " " << M <<endl;

	double z;
	double S_t_dir;
	double c_dir;
	double p_dir;

	double S_t_dis;
	double c_dis;
	double p_dis;

	for(int i=0; i < M; i++){

		sum1_c_dir = 0;
		sum1_p_dir = 0;

		sum1_c_dis = 0;
		sum1_p_dis = 0;

		for(int j=0; j<L; j++){

			//direct
			z = rnd.Gauss(0, 1);
			S_t_dir = S_t_MC(S_0, r, sigma, t_f, t_0, z);//(double S_t1, double mu, double sigma, double t2, double t1, double z);
			
			c_dir = EuroCallOpt_MC(S_t_dir, k);
			p_dir = EuroPutOpt_MC(S_t_dir, k);

			sum1_c_dir += c_dir;
			sum1_p_dir += p_dir;

			//discrete

			S_t_dis = S_t_MC_discrete(S_0, r, sigma, t_f, 0, rnd, N_int);
			
			c_dis = EuroCallOpt_MC(S_t_dis, k);
			p_dis = EuroPutOpt_MC(S_t_dis, k);
						
			sum1_c_dis += c_dis;
			sum1_p_dis += p_dis;
		}

		//calcolo le medie dei blocchi	
		
		//direct
		ave1_c_dir = sum1_c_dir / double(L);
		ave1_squared_c_dir = ave1_c_dir * ave1_c_dir;

		ave1_p_dir = sum1_p_dir / double(L);
		ave1_squared_p_dir = ave1_p_dir * ave1_p_dir;

		//discrete
		ave1_c_dis = sum1_c_dis / double(L);
		ave1_squared_c_dis = ave1_c_dis * ave1_c_dis;
		
		ave1_p_dis = sum1_p_dis / double(L);
		ave1_squared_p_dis = ave1_p_dis * ave1_p_dis;


		//sommo i valori delle medie ottenute

		//direct
		sum2_c_dir += ave1_c_dir;
		sum3_c_dir += ave1_squared_c_dir;

		sum2_p_dir += ave1_p_dir;
		sum3_p_dir += ave1_squared_p_dir;

		//discrete
		sum2_c_dis += ave1_c_dis;
		sum3_c_dis += ave1_squared_c_dis;
		
		sum2_p_dis += ave1_p_dis;
		sum3_p_dis += ave1_squared_p_dis;


		//calcolo la media sui valori delle medie dei blocchi

		//direct
		ave2_c_dir = sum2_c_dir / double(i+1);//i+1 altrimenti sarebbe N-1 al denominatore
		ave2_squared_c_dir = sum3_c_dir / double(i+1);

		ave2_p_dir = sum2_p_dir / double(i+1);
		ave2_squared_p_dir = sum3_p_dir / double(i+1);

		StDevMean_c_dir = sqrt( (ave2_squared_c_dir - ave2_c_dir*ave2_c_dir) / double(i+1) );
		StDevMean_p_dir = sqrt( (ave2_squared_p_dir - ave2_p_dir*ave2_p_dir) / double(i+1) );

		//discrete
		ave2_c_dis = sum2_c_dis / double(i+1);
		ave2_squared_c_dis = sum3_c_dis / double(i+1);
		
		ave2_p_dis = sum2_p_dis / double(i+1);
		ave2_squared_p_dis = sum3_p_dis / double(i+1);
		
		StDevMean_c_dis = sqrt( (ave2_squared_c_dis - ave2_c_dis*ave2_c_dis) / double(i+1) );
		StDevMean_p_dis = sqrt( (ave2_squared_p_dis - ave2_p_dis*ave2_p_dis) / double(i+1) );
		
		//scrivo su un file i valori ottenuti moltiplicando per e^-r(t_f - t_0) che lo sconto dovuto al fatto che se il writer investisse in banca il denaro dell'opzione così ricaverebbe al delivery time un ammonto pari a quello che guadagna l'holder
		
		output[1] << ave2_c_dir * discount << " " << StDevMean_c_dir * discount << endl;
		output[2] << ave2_p_dir * discount << " " << StDevMean_p_dir * discount << endl;
		output[3] << ave2_c_dis * discount << " " << StDevMean_c_dis * discount << endl;
		output[4] << ave2_p_dis * discount << " " << StDevMean_p_dis * discount << endl;
			
	}

	for(int i=0; i<5; i++){
			output[i].close();
	}
	
	return 0;
}


/*
double N_x(double x){
	return 0.5*( 1 + erf( x/sqrt(2.0) ) );
}

double d1(double T_f, double t, double r, double sigma, double S_t, double k){
	return ( log(S_t/k) + r*sigma*sigma*(T_f - t) ) / ( sigma * sqrt(T_f - t) );
}

double d2(double d1, double sigma, double T_f, double t){
	return d1 - sigma*sqrt(T_f - t);
}

double EuroCallOpt(double S_t, double N_d1, double k, double r, double T_f, double t, double N_d2){
	return S_t * N_d1 - k * N_d2 * exp( -r * ( T_f - t ) );
}

double EuroPutOpt(double S_t, double N_d1, double k, double r, double T_f, double t, double N_d2){
	return S_t * ( N_d1 - 1 ) - k * ( N_d2 - 1 ) * exp( -r * ( T_f - t ) );
}

*/

double S_t_MC(double S_t1, double mu, double sigma, double t2, double t1, double z){
	return S_t1 * exp( ( mu - sigma * sigma / 2.0 ) * ( t2 - t1 ) + sigma * z * sqrt( t2 - t1 ) );
}

double S_t_MC_discrete(double S_t1, double mu, double sigma, double t_f, double t_0, Random rnd, int N_int ){

	double t_2;
	double t_1 = t_0;
	double z;
	
	for(int i=0; i<N_int; i++){

		t_2 = t_1 + (t_f - t_0) / double(N_int);
		z = rnd.Gauss(0, 1);
		S_t1 = S_t_MC(S_t1, mu, sigma, t_2, t_1, z);
		t_1 = t_2;

	}

	return S_t1;
}


double EuroCallOpt_MC(double S_t, double k){
	return max(0.0, S_t - k);
}

double EuroPutOpt_MC(double S_t, double k){
	return max(0.0, k - S_t);	
}
