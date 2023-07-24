#include "random.hpp"
#include "StatsCalculator.hpp"

using namespace std;

class piCalculator : public StatsCalculator, public Random {

	public:

		piCalculator(int N, int M, double L_needle, double d, int N_throws) : StatsCalculator(N, M), Random(){
			L_needle_ = L_needle;
			d_ = d;
			N_throws_ = N_throws;
		}

		void SetL_needle(double L_needle){L_needle_ = L_needle;}
		void SetD(double d){d_ = d;}
		void SetN_throws(double N_throws){N_throws_ = N_throws;}

		double GetValue(){
			cout<<"this function is not implemented you maybe want to use GetVectorValues"<<endl;
			exit(-1);
		}
		
		vector<double> GetVectorValues(){

			vector<double> values;
		
			int N_hit = 0;

			double y0;//coordinata y di un estremo dell'ago		
			double x1;//coordinate di un punto all'interno di un quadrato di lato 2
			double y1;
		
			for(int i=0; i<N_throws_; i++){
				
				//metodo accept-reject per determinare un punto dentro il cerchio di raggio L_needle
				
				y0 = Rannyu(0, d_);
				
				x1 = Rannyu(-1, 1);
				y1 = Rannyu(-1, 1);
				
				while( y1 < -sqrt( 1 - x1*x1 ) || y1 > sqrt( 1 - x1*x1 ) ){
					x1 = Rannyu(-1, 1);
					y1 = Rannyu(-1, 1);
				}
				//ora (x1, y1) sono le coordinate di un punto a caso dentro il cerchio di raggio 1
				//normalizzo (x1, y1) e moltiplico per L_needle
				//cout<<"ok2"<<i<<endl;
				y1 = y1 * L_needle_ / sqrt( x1*x1 + y1*y1 );
				
				y1 = y1 + y0;//ora y1 è la coordinata y dell'altro estremo dell'ago
				
				if(y1 > d_ || y1 < 0){
					N_hit++;
				}
	
			}
			values.push_back( 2 * L_needle_ * N_throws_ / double( N_hit * d_ ) );
			
			return values;
			
		}

		void Print_Parameters(const string Filename){
			ofstream output;
			output.open(Filename);
			output << N_ << " " << M_ << " " << L_needle_ << " " << d_ << " " << N_throws_ << endl;
		}

	private:
		double L_needle_;
		double d_;
		int N_throws_;
};
 
int main (){
	
	int N = 10000;//number on which mean is calculated
	int M = 100;//number of blocks
	double L_needle = 1;
	double d = 2;
	int N_throws = 10000;
	
	piCalculator p(N, M, L_needle, d, N_throws);
	p.SetRandomPrimes();
	
	p.Print_Parameters("/mnt/c/Users/giaco/esercizi_python/1/pi_greco_parameters.txt");

	vector<string> Filenames;
	Filenames.push_back("/mnt/c/Users/giaco/esercizi_python/1/pi_greco.txt");
	p.Print_Mean_StDevMean(Filenames);
		
	return 0;
}




/*
				double x1 = rnd.Rannyu();//coordinate di un secondo punto all'interno di un quadrato
				double y1 = rnd.Rannyu();
				
				double t = rnd.sign() * L_needle / sqrt( (x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0) );
		
				double y2 = (y1 - y0)*t + y0;//coordinata y del secondo estremo dell'ago
				double x2 = (x1 - x0)*t + x0;
	
				if(y2 > d || y2 < 0){
					N_hit++;
				}
				*/
