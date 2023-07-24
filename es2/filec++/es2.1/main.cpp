#include "random.hpp"
#include "StatsCalculator.hpp"

using namespace std;

class CosineGen : public StatsCalculator, public Random{
	public:
		CosineGen(int N, int M) : StatsCalculator(N, M), Random(){};

		double GetValue(){
			cout<<"this function is not implemented you maybe want to use GetVectorValues"<<endl;
			exit(-1);
		}
		
		vector<double> GetVectorValues(){

			vector<double> values;
			
			double r = Rannyu();

			//uniform sampling
			double x_uni = M_PI*cos(M_PI*r/2.0) / 2.0;
			values.push_back(x_uni);
						
			//importance sampling with inverted cumulative function
			double r2 = 1 - sqrt(1 - r);
			double x_imp = M_PI*cos(M_PI*r2/2.0) / (4 - 4*r2);
			values.push_back(x_imp);
			
			//importance sampling with Accept Reject
			double p_r = 3*(1-r*r)/2.0;
			double r3 = Rannyu();
			
			while(r3 > p_r*2/3.0){
				r = Rannyu();
				p_r = 3*(1-r*r)/2.0;
				r3 = Rannyu();
			} 

			double x_imp_AR = M_PI*cos( M_PI*r/2.0 ) / ( 3*( 1 - r*r ) );
			values.push_back(x_imp_AR);

			return values;
		}
	private:
};
 
int main (){

	int N = 100000;
	int M = 100;//number of blocks
	//int L = N/M;//length of a block
	CosineGen c(N, M);
	c.SetRandomPrimes();

	c.Print_Parameters("/mnt/c/Users/giaco/esercizi_python/2/integral_parameters.txt");

	vector<string> filenames;
	filenames.push_back("/mnt/c/Users/giaco/esercizi_python/2/integral_uniform.txt");
	filenames.push_back("/mnt/c/Users/giaco/esercizi_python/2/integral_importance_sampling.txt");
	filenames.push_back("/mnt/c/Users/giaco/esercizi_python/2/integral_importance_sampling_AR.txt");
	c.Print_Mean_StDevMean(filenames);
	
	return 0;
}
