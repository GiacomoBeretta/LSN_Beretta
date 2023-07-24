#include "random.hpp"
#include <vector>

using namespace std;

template <typename T>
double squared_norm(T* walker);
 
int main (){
	
	Random rnd;
	rnd.SetRandomPrimes();
		
	int N = 10000; //realizzazioni
	int M = 100; // numero di blocchi
	int L = N/M;
	int N_step = 100; // numero di step del random walk
	int a = 1; //lattice constant

	vector<double> ave2_dis( N_step + 1 );
	vector<double> ave2_squared_dis( N_step + 1 );
	vector<double> StDevMean_dis( N_step + 1 );
		
	vector<double> ave2_con( N_step + 1 );
	vector<double> ave2_squared_con( N_step + 1 );
	vector<double> StDevMean_con( N_step + 1 );
	
	
	vector<double> sum2_dis(N_step+1);
	vector<double> sum3_dis(N_step+1);//sum on the squared values
	
	vector<double> sum2_con(N_step+1);
	vector<double> sum3_con(N_step+1);

	ofstream output[2];
	output[0].open("/mnt/c/Users/giaco/esercizi_python/2/discrete_random_walk.txt");
	output[1].open("/mnt/c/Users/giaco/esercizi_python/2/continuous_random_walk.txt");
	
	for(int i=0; i<M; i++){
		
		vector<double> sum1_dis( N_step + 1 ); //inizializzo a zero le componenti
		vector<double> ave1_dis( N_step + 1 ); // N_step + 1 perchè i random walk partono da 0
		vector<double> ave1_squared_dis( N_step + 1 );
		
		vector<double> sum1_con( N_step + 1 );
		vector<double> ave1_con( N_step + 1 );
		vector<double> ave1_squared_con( N_step + 1 );
			
		for(int j=0; j<L; j++){
		
			int walker_dis[3]    = {0, 0, 0};//reset walker
			double walker_con[3] = {0, 0, 0};
		
			for(int k=1; k<N_step+1; k++){//così fa N_step e il primo rimane 0
			
				int dir  = rnd.Rannyu(0, 3);	
				int step = a*rnd.sign();
			
				walker_dis[dir] += step;
			
				sum1_dis[k] += squared_norm<int>(walker_dis);//c'è il += così sommo L volte

				//continuous walker
				double theta = acos( 1 - 2 * rnd.Rannyu() );
				double phi   = rnd.Rannyu( 0, 2*M_PI );
				
				walker_con[0] += a*sin(theta)*cos(phi);
				walker_con[1] += a*sin(theta)*sin(phi);
				walker_con[2] += a*cos(theta);

				sum1_con[k] += squared_norm<double>(walker_con);
				
			}			
		}

		
		for(int j=0; j<N_step+1; j++){
			
			ave1_dis[j]         = sum1_dis[j] / double(L);//calcolo le medie nei blocchi per ogni passo j (ora i valori ave1 rappresentano i miei nuovi valori di partenza su cui calcolo la media e l'incertezza)
			ave1_squared_dis[j] = ave1_dis[j]*ave1_dis[j];

			sum2_dis[j] += ave1_dis[j];// c'è += così sommo M volte
			sum3_dis[j] += ave1_squared_dis[j];
						
			//continuous walker
			ave1_con[j]         = sum1_con[j] / double(L);//calcolo le medie nei blocchi per ogni passo j (ora i valori ave1 rappresentano i miei nuovi valori di partenza su cui calcolo la media e l'incertezza)
			ave1_squared_con[j] = ave1_con[j]*ave1_con[j];
			
			sum2_con[j] += ave1_con[j];// c'è += così sommo M volte
			sum3_con[j] += ave1_squared_con[j];

		}
	}

	for(int i=0; i<N_step+1; i++){

		ave2_dis[i] = sum2_dis[i] / double(M);
		ave2_squared_dis[i] = sum3_dis[i] / double(M);

		StDevMean_dis[i] = sqrt( ( ave2_squared_dis[i] - ave2_dis[i]*ave2_dis[i] ) / double(M - 1) );


		ave2_con[i] = sum2_con[i] / double(M);
		ave2_squared_con[i] = sum3_con[i] / double(M);
		
				StDevMean_con[i] = sqrt( ( ave2_squared_con[i] - ave2_con[i]*ave2_con[i] ) / double(M - 1) );

		output[0] << sqrt( ave2_dis[i] ) <<" "<< StDevMean_dis[i] << endl;
		output[1] << sqrt( ave2_con[i] ) <<" "<< StDevMean_con[i] << endl;
		
	}

	output[0].close();
	output[1].close();

	//chi quadro
	
	return 0;
}

template <typename T>
double squared_norm(T* walker){
	return walker[0]*walker[0] + walker[1]*walker[1] + walker[2]*walker[2];
}
	/*
	vector<double> ave2_dis( N_step + 1 );
	vector<double> ave2_squared_dis( N_step + 1 );
	vector<double> StDevMean_dis( N_step + 1 );
		
	vector<double> ave2_con( N_step + 1 );
	vector<double> ave2_squared_con( N_step + 1 );
	vector<double> StDevMean_con( N_step + 1 );
	
	
	vector<double> sum2_dis(N_step+1);
	vector<double> sum3_dis(N_step+1);//sum on the squared values
	
	vector<double> sum2_con(N_step+1);
	vector<double> sum3_con(N_step+1);

	ofstream output[2];
	output[0].open("/mnt/c/Users/giaco/esercizi_python/2/discrete_random_walk.txt");
	output[1].open("/mnt/c/Users/giaco/esercizi_python/2/continuous_random_walk.txt");
	
	for(int i=0; i<M; i++){
		
		vector<double> sum1_dis( N_step + 1 ); //inizializzo a zero le componenti
		vector<double> ave1_dis( N_step + 1 ); // N_step + 1 perchè i random walk partono da 0
		vector<double> ave1_squared_dis( N_step + 1 );
		
		vector<double> sum1_con( N_step + 1 );
		vector<double> ave1_con( N_step + 1 );
		vector<double> ave1_squared_con( N_step + 1 );
			
		for(int j=0; j<L; j++){
		
			int walker_dis[3]    = {0, 0, 0};//reset walker
			double walker_con[3] = {0, 0, 0};
		
			for(int k=1; k<N_step+1; k++){//così fa N_step e il primo rimane 0
			
				int dir  = rnd.Rannyu(0, 3);	
				int step = a*rnd.sign();
			
				walker_dis[dir] += step;
			
				sum1_dis[k] += squared_norm<int>(walker_dis);//c'è il += così sommo L volte

				//continuous walker
				double theta = acos( 1 - 2 * rnd.Rannyu() );
				double phi   = rnd.Rannyu( 0, 2*M_PI );
				
				walker_con[0] += a*sin(theta)*cos(phi);
				walker_con[1] += a*sin(theta)*sin(phi);
				walker_con[2] += a*cos(theta);

				sum1_con[k] += squared_norm<double>(walker_con);
				
			}			
		}

		
		for(int j=0; j<N_step+1; j++){
			
			ave1_dis[j]         = sum1_dis[j] / double(L);//calcolo le medie nei blocchi per ogni passo j
			ave1_squared_dis[j] = ave1_dis[j]*ave1_dis[j];

			sum2_dis[j] += ave1_dis[j];// c'è += così sommo M volte
			ave2_dis[j] = sum2_dis[j] / double(i+1);//calcolo la media al variare del numero di blocchi (c'è i+1 altrimenti se fosse solo i dividerei per N-1)

			sum3_dis[j] += ave1_squared_dis[j];
			ave2_squared_dis[j] = sum3_dis[j] / double(i+1);

			StDevMean_dis[j] = sqrt( (ave2_squared_dis[j] - ave2_dis[j]*ave2_dis[j]) / double(i+1) );

						
			//continuous walker
			ave1_con[j]         = sum1_con[j] / double(L);//calcolo le medie nei blocchi per ogni passo j
			ave1_squared_con[j] = ave1_con[j]*ave1_con[j];
			
			sum2_con[j] += ave1_con[j];// c'è += così sommo M volte
			ave2_con[j] = sum2_con[j] / double(i+1);//calcolo la media al variare del numero di blocchi (c'è i+1 altrimenti se fosse solo i dividerei per N-1)
			
			sum3_con[j] += ave1_squared_con[j];
			ave2_squared_con[j] = sum3_con[j] / double(i+1);

			StDevMean_con[j] = sqrt( (ave2_squared_con[j] - ave2_con[j]*ave2_con[j] ) / double(i+1) );
			

		output[0] << ave2_dis[i] <<" "<< StDevMean_dis[i] << endl;
		output[1] << ave2_con[i] <<" "<< StDevMean_con[i] << endl;
	}

	output[0].close();
	output[1].close();
	
	return 0;
}*/
