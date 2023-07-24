#include "random.hpp"
#include <vector>

using namespace std;

double ProbGauss(double a, double b);//calcola la differenza F(b)-F(a) dove F è la funzione di ripartizione della gaussiana
double Chi(vector<int> O, vector<double> E);
 
int main (){
	
	Random rnd;
	rnd.SetRandomPrimes();
	
	int M = 10000;//10000 iterazioni o realizzazioni
	
	int N[4] = {1,2,10,100};

	int sum_dice = 0;
	double sum_exp = 0;
	double sum_cauchy = 0;
	double S_N_dice;
	double S_N_exp;
	double S_N_cauchy;

	int bins = 101;	
	
	vector<int> O_dice(bins);//vettore che contiene i conteggi osservati
	vector<double> E_dice(bins);//contiene i conteggi aspettati usando una gaussiana
	
	double exp_min = 0.6;
	double exp_max = 1.4;
	int N_exp = 0;
	vector<int> O_exp(bins);
	vector<double> E_exp(bins);
	
	//vector<int> O_cauchy(bins);
	//vector<double> E_cauchy(bins);
	
	double mean_dice = 3.5; //   =(1+2+3+4+5+6)/6
	double devStMean_dice = sqrt( ( (1+4+9+16+25+36)/double(6) - mean_dice*mean_dice ) / double(N[3]) ) ;
	
	double mean_exp = 1;//    =1/lambda
	double devStMean_exp = 1/sqrt(N[3]); //  =1/ ( lambda * sqrt(N) )
	
	ofstream output[4];
	output[0].open("/mnt/c/Users/giaco/esercizi_python/1/Standard_Dice.txt");
	output[1].open("/mnt/c/Users/giaco/esercizi_python/1/Exponential.txt");
	output[2].open("/mnt/c/Users/giaco/esercizi_python/1/Cauchy_Lorentz.txt");
	output[3].open("/mnt/c/Users/giaco/esercizi_python/1/chi_quadro_dice.txt");
	
	for(int i=0; i<4; i++){
		
		for(int j=0; j<M; j++){
			
			sum_dice   = 0;
			sum_exp    = 0;
			sum_cauchy = 0;
			
			for(int k=0; k<N[i]; k++){//somma sugli N=1,2,10,100
				
				sum_dice   += rnd.Integer(1,6);
				sum_exp    += rnd.Exp(1);
				sum_cauchy += rnd.CauchyLorentz(0,1);		
				
			}
			
			S_N_dice   = sum_dice   / double( N[i] );
			S_N_exp    = sum_exp    / double( N[i] );
			S_N_cauchy = sum_cauchy / double( N[i] );
			
			output[0] << S_N_dice   << endl;
			output[1] << S_N_exp    << endl;
			output[2] << S_N_cauchy << endl;

			//calculate observed values for dice and exp
			if(i==3){//N=100
			
				int x_dice = (S_N_dice - 1) * N[i] / double(5); // ( (S_N - 1)/5 ) * 100
				O_dice[x_dice]++;	
	
				if(S_N_exp > exp_min && S_N_exp < exp_max){
					
					N_exp++;
					int x_exp = (S_N_exp - exp_min) * N[i] / double(exp_max - exp_min); //  ( (S_N - 0.6) / (1.4 - 0.6) ) * 100
					O_exp[x_exp]++;
				}
			}
		
		}
	}
	double chi_dice = 0;
	double chi_exp = 0;
	//double chi_cauchy = 0;
	
	for(int i=0; i<bins; i++){
		
		double a_dice = 1 + i*(6-1)/double(bins-1) - (6-1)/double( (bins-1)*2 );
		double b_dice = 1 + i*(6-1)/double(bins-1) + (6-1)/double( (bins-1)*2 );
		
		double z_a_dice = (a_dice - mean_dice) / double(devStMean_dice);
		double z_b_dice = (b_dice - mean_dice) / double(devStMean_dice);
		
		double p = ProbGauss(z_a_dice, z_b_dice);
		E_dice[i] = M * p; //valore atteso = N*p		
		
		double a_exp = exp_min + i*(exp_max - exp_min)/double(bins-1) - (exp_max - exp_min)/double( (bins-1)*2 );
		double b_exp = exp_min + i*(exp_max - exp_min)/double(bins-1) + (exp_max - exp_min)/double( (bins-1)*2 );
		
		double z_a_exp = (a_exp - mean_exp) / double(devStMean_exp);
		double z_b_exp = (b_exp - mean_exp) / double(devStMean_exp);
		
		p = ProbGauss(z_a_exp, z_b_exp);
		E_exp[i] = N_exp * p;
		//cout<<"i="<<i<<" -> "<<E_exp[i]<<endl;
	}
	
	chi_dice    = Chi(O_dice, E_dice);
	chi_exp     = Chi(O_exp, E_exp);
	//chi_cauchy  = chi<double>(O_cauchy, E_cauchy);
	
	output[3] << chi_dice << endl;
	output[3] << chi_exp  << endl;
	//output[3] << chi_cauchy <<endl;
	
	for(int i=0; i<4; i++){
		output[i].close();
	}

   return 0;
   
}


double ProbGauss(double a, double b){//da vedere su wikipedia per calcolare la cumulativa a partire dalla erf
	
	if(a > b){
		cout<<"a > b, dev'essere a < b, riprovare"<<endl;
		exit(-1);
	}
	
	double t_a = a / sqrt(2);
	double t_b = b / sqrt(2);
	
	return 0.5 * ( erf(t_b) - erf(t_a) );
}

double Chi(vector<int> O, vector<double> E){//ho scritto una funzione apposita per evitare la divisione per 0
	
	if( O.size() != E.size() ){
		cout<<"si è cercato di calcolare il chi quadro ma i due vettori hanno dimensione diversa"<<endl;
		exit(-1);
	}
	
	double chi = 0;
	
	for(int i=0; i<O.size(); i++){
		
		if( E[i] == 0){
			
			if( O[i] - E[i] != 0){
				chi += (O[i] - E[i]) * (O[i] - E[i]) / double(E[i] + 1.0e-13);
			}//else chi += 0
			
		}else{
			chi += (O[i] - E[i]) * (O[i] - E[i]) / double(E[i]);
		}
			
	}
	
	return chi;
}
