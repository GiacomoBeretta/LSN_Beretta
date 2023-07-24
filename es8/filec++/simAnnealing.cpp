#include "simAnnealing.hpp"

using namespace std;

SimulatedAnnealing :: SimulatedAnnealing() : Metropolis()
{
	T_min_ = 0;
	T_max_ = 0;
	n_temp_ = 0;
	temp_ = 0;
}

SimulatedAnnealing :: SimulatedAnnealing(double delta, int equilibration_time, double T_min, double T_max, int n_temp) : Metropolis(delta, equilibration_time)
{
	T_min_ = T_min;
	T_max_ = T_max;
	n_temp_ = n_temp;
	temp_ = T_max;
}

void SimulatedAnnealing :: setTMin(double T_min)
{
	T_min_ = T_min;
}

void SimulatedAnnealing :: setTMax(double T_max)
{
	T_max_ = T_max;	
}

void SimulatedAnnealing :: setNTemp(double n_temp)
{
	n_temp_ = n_temp;	
}

void SimulatedAnnealing :: setTemp(double temp)
{
	temp_ = temp;
}

void SimulatedAnnealing :: setSA_nstep(int SA_nstep)
{
	SA_nstep_ = SA_nstep;
}

double SimulatedAnnealing :: probability(vector<double> x)
{
	double costo = cost(x);
	//cout << "cost="<<costo<<endl;
	return exp(-costo/temp_);
	//return exp( cost(x)/(-temp_) );
}

void SimulatedAnnealing :: moveUniform()
{	
	vector<double> x_old = getX();
	vector<double> x_new( x_old.size() );
	int accepted = getAccepted();
	int attempted = getAttempted();
	
/*	for(unsigned int i=0; i<x_old.size(); i++)
	{
		x_new[i] = x_old[i] + ( Rannyu() - 0.5 )*getDelta()*pow(temp_, 3);
	}
	double p_old = probability(x_old);
	double cost_old = getCostValue();
	
	double p_new = probability(x_new);
	double cost_new = getCostValue();

	cost_ = cost_old;
	if( min( 1.0, p_new/p_old ) >= 1 )
	{
		cost_ = cost_new;
		x_old = x_new;
		accepted++;
	}else if( Rannyu() < p_new/p_old )
	{
		cost_ = cost_new;
		x_old = x_new;
		accepted++;
	}
	attempted++;

	setX(x_old);
	setAccepted(accepted);
	setAttempted(attempted);
	*/
	
	for(unsigned int i=0; i<x_old.size(); i++)
	{
		x_new[i] = x_old[i] + ( Rannyu() - 0.5 )*getDelta();

		double p_old = probability(x_old);
		double cost_old = getCostValue();
		
		double p_new = probability(x_new);
		double cost_new = getCostValue();
	
		if( p_new/p_old >= 1 )
		{
			cost_old = cost_new;
			x_old = x_new;
			accepted++;
		}else if( Rannyu() < p_new/p_old )
		{
			cost_old = cost_new;
			x_old = x_new;
			accepted++;
		}
		attempted++;		
		cost_ = cost_old;
	}	
	setX(x_old);
	setAccepted(accepted);
	setAttempted(attempted);
}

void SimulatedAnnealing :: PrintCostAndParameters(const string Filename)
{
	ofstream output( Filename );

	temp_ = T_max_;
	double delta_temp = (T_max_ - T_min_)/double(n_temp_-1);

	for(int i=0; i<n_temp_; i++)
	{		
		for(int j=0; j<3; j++) //moves three times the parameters before changing the temperature
		{
			//move parameters ( the cost is evaluated during the move )
			moveUniform();

			//print parameters values
			output << i << " " << temp_ << " ";
			for(unsigned int j=0; j<getX().size(); j++)
			{
				output << getX()[j] << " ";
			}
			output << cost_ << endl;
		}
		temp_ += delta_temp;
	}
	output.close();
}
