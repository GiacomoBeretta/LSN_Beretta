#include "random.hpp"
#include "StatsCalculator.hpp"
#include <vector>

using namespace std;

class TestRandomGen : public StatsCalculator, public Random {

	public:

		//TestRandom();
		TestRandomGen(int N, int M) : StatsCalculator(N, M), Random(){}
		//~TestRandom();

		double GetValue(){
			cout<<"this function is not implemented you maybe want to use GetVectorValues"<<endl;
			exit(-1);
		}
		
		vector<double> GetVectorValues(){
		
			vector<double> values(2);
			
			double r = Rannyu();
			
			values[0] = r;
			values[1] = ( r - 0.5 ) * ( r - 0.5 );
			return values;
		}
		
};


int main (){

	int N = 100000;//numero di random totali
	int M = 100; //numero di blocchi

	TestRandomGen r(N, M);
	r.SetRandomPrimes();

	vector<string> Filenames;
	Filenames.push_back("/mnt/c/Users/giaco/esercizi_python/1/random_mean_dev.txt");
	Filenames.push_back("/mnt/c/Users/giaco/esercizi_python/1/random_sigma2.txt");

	r.Print_Parameters("/mnt/c/Users/giaco/esercizi_python/1/test_parameters.txt");
	r.Print_Mean_StDevMean( Filenames );

	//calcolo del chi quadro
	int bins = 100; //numero di intervalli in cui divido (0,1)
	int N_throws = 10000; 
	int N_chi = 100000;//numero di volte con cui ripeto l'esperimento

	double E_chi = N_throws/bins; //valore di aspettazione E=n*p (numero di tentavi per probabilità p -> E = 10^4/M)
	double chi;
	ofstream output_chi;
	output_chi.open("/mnt/c/Users/giaco/esercizi_python/1/random_chi_quadro.txt");

	output_chi << bins << endl;
	output_chi << N_chi << endl;
	
	for(int i=0; i<N_chi; i++){

		chi = 0;
		vector<int> n(bins);//vettore che contiene i conteggi, n[0] � il numero di numeri casuali usciti nel primo intervallo (0, 0.01), n[1] quelli usciti nel secondo intervallo (0.01, 0.02)
				
		for(int j=0; j < N_throws; j++){ //ciclo sui numeri casuali ogni 10^4

			int k = r.Rannyu() * bins;			
			n[k]++;			
				
		}
		
		for(int j=0; j<bins; j++){
			
			chi += (n[j] - E_chi) * (n[j] - E_chi) / double(E_chi);

		}
		output_chi<<chi<<endl;

	}
	output_chi.close();
	return 0;
}
