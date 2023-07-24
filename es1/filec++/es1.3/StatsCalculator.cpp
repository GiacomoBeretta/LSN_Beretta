#include "StatsCalculator.hpp"

using namespace std;

StatsCalculator :: StatsCalculator(){}

StatsCalculator :: StatsCalculator(int N, int M){

	if( N<0 || M<0 || N%M != 0 ){
		cout<<"N and M must be positive and N/M must be an integer"<<endl;
		exit(-1);
	}else{
		N_ = N;
		M_ = M;
		L_ = N/M;
	}

	mean_ = 0;	
	StDevMean_ = 0;
}

StatsCalculator :: ~StatsCalculator(){}

void StatsCalculator :: Print_Mean_StDevMean(const string Filename){

	double value;
	double ave1;
	double ave1_squared;
	double ave2;
	double ave2_squared;
	double StDevMean;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	
	ofstream out;
	out.open(Filename);

	for(int i=0; i<M_; i++){
	
		sum1 = 0;
				
		for(int j=0; j<L_; j++){
		
			value = GetValue();
			sum1 += value;	
		}
				
		//calcolo le medie dei blocchi
		ave1 = sum1 / double(L_);
		ave1_squared = ave1 * ave1;
	
		//calcolo la media sui valori delle medie dei blocchi
		sum2 += ave1;
		sum3 += ave1_squared;
	
		ave2 = sum2 / double(i+1); //i+1 altrimenti sarebbe N-1 al denominatore
		
		ave2_squared = sum3 / double(i+1);
			
		StDevMean = sqrt( (ave2_squared - ave2*ave2) / double(i+1) );
	
		//scrivo su un file i valori ottenuti
		out << ave2 << " " << StDevMean << endl;
		
	}				
	
	out.close();
	
}

void StatsCalculator :: Print_Mean_StDevMean(const vector<string> Filenames){

	int dim = Filenames.size();

	vector<double> value(dim);
	vector<double> ave1(dim);
	vector<double> ave1_squared(dim);
	vector<double> ave2(dim);
	vector<double> ave2_squared(dim);
	vector<double> StDevMean(dim);
	vector<double> sum1(dim);
	vector<double> sum2(dim);
	vector<double> sum3(dim);

	vector<double> null(dim);

	ofstream *out = new ofstream[ dim ];
	
	for(int i=0; i<dim; i++){
		
		out[i].open(Filenames[i]);

	}
	
	for(int i=0; i<M_; i++){

		sum1 = null;
			
		for(int j=0; j<L_; j++){
	
			value = GetVectorValues();
			transform(sum1.begin(), sum1.end(), value.begin(), sum1.begin(), plus<double>() ); //sum1 += value

		}
			
		//calcolo le medie dei blocchi
		//ave1 = sum1 / double(L_);
		//ave1_squared = ave1 * ave1;
		transform(sum1.begin(), sum1.end(), ave1.begin(), [ this ](const double &c ){ return c / double( this->L_ ); } );	
		transform(ave1.begin(), ave1.end(), ave1_squared.begin(), [](const double &c){ return c*c; } ); 

		//calcolo la media sui valori delle medie dei blocchi
		//sum2 += ave1;
		//sum3 += ave1_squared;
		transform(sum2.begin(), sum2.end(), ave1.begin(), sum2.begin(), plus<double>() );
		transform(sum3.begin(), sum3.end(), ave1_squared.begin(), sum3.begin(), plus<double>() ); 

		//ave2 = sum2 / double(i+1); //i+1 altrimenti sarebbe N-1 al denominatore
		transform(sum2.begin(), sum2.end(), ave2.begin(), [ i ]( const double &c ){ return c / double(i+1); } ); 
		
		//ave2_squared = sum3 / double(i+1);
		transform(sum3.begin(), sum3.end(), ave2_squared.begin(), [ i ]( const double &c ){ return c / double(i+1); } );
	
		//StDevMean = sqrt( (ave2_squared - ave2*ave2) / double(i+1) );
		transform(ave2_squared.begin(), ave2_squared.end(), ave2.begin(), StDevMean.begin(), [ i ]( const double &a, const double &b ){ return sqrt( (a - b*b) / double(i+1) ); } );
		
		//scrivo su un file i valori ottenuti
		for(int j=0; j<dim; j++){
			out[j] << ave2[j] << " " << StDevMean[j] << endl;
		}
		
			
	}

	for(int i=0; i<dim; i++){
			
		out[i].close();
		
	}	
		
	//mean_ = ave2;
	//StDevMean_ = StDevMean;
}

void StatsCalculator :: Print_Parameters(const std::string Filename){
	ofstream output;
	output.open(Filename);
	output << N_ << " " << M_ << endl;
}
