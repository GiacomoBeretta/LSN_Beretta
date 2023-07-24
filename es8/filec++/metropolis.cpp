#include "metropolis.hpp"

using namespace std;

Metropolis :: Metropolis() : Random()
{
	vector<double> v(1);
	x_old_ = v;
		
	delta_ = 0;
	equilibration_time_ = 0;

	accepted_ = 0;
	attempted_ = 0;
}

Metropolis :: Metropolis(double delta, int equilibration_time) : Random()
{
	vector<double> v(1);
	x_old_ = v;
		
	delta_ = delta;
	equilibration_time_ = equilibration_time;

	accepted_ = 0;
	attempted_ = 0;
}

void Metropolis :: setDelta(double delta)
{
	delta_ = delta;
}

void Metropolis :: setEquilibrationTime(int equilibration_time)
{
	equilibration_time_ = equilibration_time;
}

void Metropolis :: setX(vector<double> x)
{
	x_old_ = x;
}

void Metropolis :: setAccepted(int accepted)
{
	accepted_ = accepted;
}

void Metropolis :: setAttempted(int attempted)
{
	attempted_ = attempted;
}

void Metropolis :: moveUniform()
{
	//cout<<"move uniform metropolis"<< endl;
	vector<double> x_new( x_old_.size() );
	for(unsigned int i=0; i<x_old_.size(); i++)
	{
		x_new[i] = x_old_[i] + ( Rannyu() - 0.5 )*delta_;
	}
	
	double p_old = probability(x_old_);
	double p_new = probability(x_new);
/*
	cout <<"x_old = "<< x_old_[0] << endl;
	cout <<"x_new = "<< x_new[0] << endl;
	cout <<"p_old = "<< p_old << endl;
	cout <<"p_new = "<< p_new << endl;*/

	//if( min( 1.0, p_new/p_old ) >= 1 )
	if( p_old < p_new )
	{
		//cout <<"p_new/p_old >= 1 " << endl;
		x_old_ = x_new;
		accepted_++;
	}else if( Rannyu() < p_new/p_old )
	{
		x_old_ = x_new;
		//cout <<"Rannyu() < p_new/p_old " << endl;
		accepted_++;
	}else
	{
		//cout << "rejected" << endl;
	}
	attempted_++;
}

void Metropolis :: equilibrate()
{
	for(int i=0; i<equilibration_time_; i++)
	{
		moveUniform();
	}
}

double Metropolis :: computeAccRate(int nstep)
{
	setAccepted(0);
	setAttempted(0);

	for(int i=0; i<nstep; i++)
	{

		moveUniform();
	}

	return getAccRate();
}

void Metropolis :: calibrateDelta(int nstep)
{
	//cout << "recalibrate delta"<< endl;
	double delta_old, acc_rate_old, diff_old, delta_new, acc_rate_new, diff_new;
	//int max_attempt = 30;
	
	vector<double> x = getX();
	for(unsigned int i=0; i<x.size(); i++) x[i] = 0;
	
	delta_old = getDelta();
	//cout<<"ok1"<<endl;
	
	acc_rate_old = computeAccRate(nstep);
	//cout<<"ok2"<<endl;
	delta_new = delta_old + 1.0;
	setDelta(delta_new);
	acc_rate_new = computeAccRate(nstep);

	int k = 0;
	while( acc_rate_new < 0.45 || acc_rate_new > 0.55 )
	{
		/*if(acc_rate_new == 0)
		{
			cout << "k=" << k << endl;
			cout << "delta_old=" << delta_old << endl;
			cout << "delta_new=" << delta_new << endl;
			cout << "acc_rate_old=" << acc_rate_old << endl;
			cout << "acc_rate_new=" << acc_rate_new << endl;
		}*/
		
		if( acc_rate_old < 0.5 && acc_rate_new < 0.5 )
		{
			delta_old = delta_new;
			acc_rate_old = acc_rate_new;
			delta_new /= 2.0;
			
		}else if( acc_rate_old > 0.5 && acc_rate_new > 0.5 )
		{
			delta_old = delta_new;
			acc_rate_old = acc_rate_new;
			delta_new *= 2.0;		
			
		}else
		{
			diff_old = fabs(acc_rate_old - 0.5);
			diff_new = fabs(acc_rate_new - 0.5);
			if(diff_old < diff_new)
			{
				delta_new = delta_old + (delta_new - delta_old)/2.0;
			}else
			{
				delta_old = delta_new;
				acc_rate_old = acc_rate_new;
				delta_new = delta_new + 0.1;
			}
		}

		setDelta(delta_new);
		acc_rate_new = computeAccRate(nstep);
		
		k++;
		/*
		if(acc_rate_new == 0)
		{
			cout << "delta_old=" << delta_old << endl;
			cout << "delta_new=" << delta_new << endl;
			cout << "acc_rate_old=" << acc_rate_old << endl;
			cout << "acc_rate_new=" << acc_rate_new << endl;
		}*/

		/*if(k>max_attempt)
		{
			cout << "error: tried to recalibrate delta but it took over " << max_attempt << " attempts" << endl;
			exit(-1);
		}*/
	}

	setAccepted(0);
	setAttempted(0);
}
