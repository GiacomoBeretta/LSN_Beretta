#include "energy.hpp"

using namespace std;

Energy :: Energy() : Metropolis()
{
	mu_ = 0;
	sigma_ = 0;
	nblk_ = 0;
	nstep_ = 0;
	energy_ = 0;
	st_dev_mean_ = 0;
}
Energy :: Energy(double delta, int equilibration_time) : Metropolis(delta, equilibration_time)
{
	mu_ = 0;
	sigma_ = 0;
	nblk_ = 0;
	nstep_ = 0;
	energy_ = 0;
	st_dev_mean_ = 0;
}

void Energy :: setMu(double mu)
{
	mu_ = mu;
}

void Energy :: setSigma(double sigma)
{
	sigma_ = sigma;
}

void Energy :: setNblk(int nblk)
{
	nblk_ = nblk;
}

void Energy :: setNstep(int nstep)
{
	nstep_ = nstep;
}

void Energy :: setX(std::vector<double> x)
{
	if(x.size() > 1)
	{
		cout << "errore: il vettore x deve avere una sola componente" << endl;
		exit(-1);
	}
	Metropolis::setX(x);
}

void Energy :: setEnergy(double energy)
{
	energy_ = energy;
}

void Energy :: setStDevMean(double StDev)
{
	st_dev_mean_ = StDev;
}

/*
void Energy :: moveUniform()
{
	double x_old = getX()[0];
	double x_new = x_old + ( Rannyu() - 0.5 )*getDelta();
	
	double p_old = probability(x_old);
	double p_new = probability(x_new);

	int accepted = getAccepted();
	int attempted = getAttempted();

	cout <<"x_old = "<< x_old_[0] << endl;
	cout <<"x_new = "<< x_new[0] << endl;
	cout <<"p_old = "<< p_old << endl;
	cout <<"p_new = "<< p_new << endl;

	//if( min( 1.0, p_new/p_old ) >= 1 )
	if( p_old < p_new )
	{
		//cout <<"p_new/p_old >= 1 " << endl;
		x_old = x_new;
		accepted++;
	}else if( Rannyu() < p_new/p_old )
	{
		x_old_ = x_new;
		//cout <<"Rannyu() < p_new/p_old " << endl;
		accepted++;
	}else
	{
		//cout << "rejected" << endl;
	}
	attempted++;

	setAccepted(accepted);
	setAttempted(attempted);
}
*/

double Energy :: probability(vector<double> x)
{
	//cout<<"energy probability"<<endl;
	//cout << "x size="<<x.size()<<endl;
	double q = x[0];
	//cout<<"q="<<q<<endl;
	//cout<<"mu="<<mu_<<endl;
	//cout<<"sigma="<<sigma_<<endl;
		
	double a = ( q - mu_ ) / sigma_;
	double b = ( q + mu_ ) / sigma_;
	double c = exp( -a*a/2.0 );
	double d = exp( -b*b/2.0 );
	//cout<<"q=pow( ( q - mu_ ) / sigma_ , 2.0 )="<<a<<endl;
	//cout<<"b=pow( ( q + mu_ ) / sigma_ , 2.0 )="<<b<<endl;
	//cout<<"prob=pow( exp( -a/2.0 ) + exp( -b/2.0 ), 2)="<<pow( exp( -a/2.0 ) + exp( -b/2.0 ), 2)<<endl;

	return (c+d)*(c+d);
}

double Energy :: measure() const
{

	//cout<<"energy measure"<<endl;
	double x = getX()[0];
	//cout << "x="<<x<<endl;
	
	double a = ( x - mu_ ) / sigma_;
	double b = ( x + mu_ ) / sigma_;
	double c = exp( -a*a/2.0 );
	double d = exp( -b*b/2.0 );

	//cout << "a=" <<a<<endl;
	//cout << "b=" <<b<<endl;
	//cout << "c=" <<c<<endl;
	//cout << "d=" <<d<<endl;
	
	double kin = ( c*(1-a*a) + d*(1-b*b) ) / ( 2.0*sigma_*sigma_*(c+d) );
	//double kin = ( c*( sigma_*sigma_ - (x-mu_)*(x-mu_) ) + d*( sigma_*sigma_ - (x+mu_)*(x+mu_) ) ) / ( 2.0*(c+d)*pow(sigma_, 4) );
	
	double pot = pow(x,4.0) - 5.0*x*x/2.0;
	//cout << "kin="<<kin<<endl;
	//cout << "pot="<<pot<<endl;
	//cout << "measured energy = " << kin + pot << endl;
	
	return kin + pot;
}

void Energy :: printIntegralStDevMeanPerBlock(const string Filename)
{
	equilibrate();
	calibrateDelta(nstep_);
		
	ofstream out;
	out.open(Filename);
	
	double ave1;
	double ave1_squared;
	double ave2 = 0;
	double ave2_squared;
	double StDevMean = 0;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;

	//cout << "compute averages" << endl;
	for(int i=0; i<nblk_; i++)
	{
		sum1 = 0;

		setAccepted(0);
		setAttempted(0);

		for(int j=0; j<nstep_; j++)
		{
			//cout <<"j="<<j<<endl;
			moveUniform();
			sum1 += measure();
		}
		//calcolo le medie dei blocchi
		ave1 = sum1 / double(nstep_);
		ave1_squared = ave1 * ave1;
		
		//calcolo la media sui valori delle medie dei blocchi
		sum2 += ave1;
		sum3 += ave1_squared;
	
		ave2 = sum2 / double(i+1); //i+1 altrimenti sarebbe N-1 al denominatore
		
		ave2_squared = sum3 / double(i+1);
			
		StDevMean = sqrt( (ave2_squared - ave2*ave2) / double(i+1) );

		//scrivo su file i valori ottenuti
		out << i << " " << ave2 << " " << StDevMean << endl;
		
		cout << "iblk = " << i << endl;
		cout << "energy = " << ave2 << ", StDevMean = " << StDevMean << endl;
 		cout << "acceptance rate = " << getAccRate() << endl;
		cout << "----------------------------" << endl << endl;
	}

	//cout << "acceptance rate = " << getAccRate() << endl;
	//cout << "----------------------------" << endl << endl;
	out.close();
	
	energy_ = ave2;
	st_dev_mean_ = StDevMean;
}

void Energy :: computeIntegral()
{
	equilibrate();
	calibrateDelta(nstep_);
	//cout << "delta = " << getDelta() << endl;
		
	double ave1;
	double ave1_squared;
	double ave2;
	double ave2_squared;
	double StDevMean;
	double sum1 = 0;
	double sum2 = 0;
	double sum3 = 0;
	
	for(int i=0; i<nblk_; i++)
	{
		sum1 = 0;

		//setAccepted(0);
		//setAttempted(0);
		
		for(int j=0; j<nstep_; j++)
		{
			moveUniform();
			sum1 += measure();
		}

		//calcolo le medie dei blocchi
		ave1 = sum1 / double(nstep_);
		ave1_squared = ave1 * ave1;
	
		//calcolo la media sui valori delle medie dei blocchi
		sum2 += ave1;
		sum3 += ave1_squared;

		//cout << "iblk = " << i << endl;
		//cout << "acceptance rate = " << getAccRate() << endl;
		//cout << "----------------------------" << endl << endl;
	}
	
	ave2 = sum2 / double(nblk_); //i+1 altrimenti sarebbe N-1 al denominatore
	ave2_squared = sum3 / double(nblk_);		
	StDevMean = sqrt( (ave2_squared - ave2*ave2) / double(nblk_) );
	
	//cout << "acceptance rate = " << getAccRate() << endl;
	//cout << "----------------------------" << endl << endl;
		
	energy_ = ave2;
	st_dev_mean_ = StDevMean;
}

void Energy :: printIstoPsi(const std::string Filename)
{
	ofstream out;
	out.open(Filename);
	calibrateDelta(nstep_);

	for(int i=0; i<nblk_; i++)
	{
		cout << "i(isto)="<< i << endl;
		for(int j=0; j<nstep_; j++)
		{
			moveUniform();
			out << getX()[0] << endl;
		}

	}
	out.close();
}

SA_Energy :: SA_Energy(double delta_SA, int equilibration_time_SA, double T_min, double T_max, int n_temp, Energy ene) : SimulatedAnnealing(delta_SA, equilibration_time_SA, T_min, T_max, n_temp)
{
	ene_ = ene;
}

void SA_Energy :: moveUniform()
{
	vector<double> x_old = getX();
	vector<double> x_new( x_old.size() );
	int accepted = getAccepted();
	int attempted = getAttempted();
	
	do
	{
		x_new[0] = x_old[0] + ( Rannyu() - 0.5 )*getDelta();
		x_new[1] = x_old[1] + ( Rannyu() - 0.5 )*getDelta();
	}
	while( x_new[0] < 0 || x_new[1] < 0 );
	
	double p_old = probability(x_old);
	double cost_old = getCostValue();
	double energy_old = ene_.getEnergy();
	double StDev_old = ene_.getStDevMean();
	
	double p_new = probability(x_new);
	double cost_new = getCostValue();
	double energy_new = ene_.getEnergy();
	double StDev_new = ene_.getStDevMean();

	cost_ = cost_old;

	if(cost_old < cost_new)
	{
		mu_min_ = x_new[0];
		sigma_min_ = x_new[1];
	}else
	{
		mu_min_ = x_old[0];
		sigma_min_ = x_old[1];
	}
	
	if( p_old < p_new )
	{
		cost_ = cost_new;
		x_old = x_new;
		energy_old = energy_new;
		StDev_old = StDev_new;
		accepted++;
	}else if( Rannyu() < p_new/p_old )
	{
		cost_ = cost_new;
		x_old = x_new;
		energy_old = energy_new;
		StDev_old = StDev_new;
		accepted++;
	}
	attempted++;

	setX(x_old);
	ene_.setEnergy(energy_old);
	ene_.setStDevMean(StDev_old);
	setAccepted(accepted);
	setAttempted(attempted);	
}

void SA_Energy :: PrintEnergy_StDev_and_Parameters(const std::string Filename)
{	
	int nstep = ene_.getNstep();
	ofstream output( Filename );

	temp_ = T_max_;
	double delta_temp = (T_max_ - T_min_)/double(n_temp_-1);
	
	double SA_delta_min = 0.01; //delta for the SA algorithm for the temperature T=0.5
	double SA_delta_max = 0.5; //delta for the SA algorithm for the temperature T=2
	double SA_delta_var = (SA_delta_max - SA_delta_min)/(double)(n_temp_-1);
	double SA_delta = SA_delta_max;
	setDelta(SA_delta);

	double ave_mu;
	double ave_sigma;
	double sum_mu = 0;
	double sum_sigma = 0;
	double ave2_mu;
	double ave2_sigma;
	for(int i=0; i<n_temp_; i++)
	{		
		cout << "------------------------- temperature = " << temp_ <<" -----------------------------------" << endl;
		cout << "SA_delta=" << getDelta() << endl;
		//cout << "calibrate delta SA algorithm"<<endl;
	//	calibrateDelta(SA_nstep_); //it takes too long maybe it's unnecessary
		setAccepted(0);
		setAttempted(0);

		ave_mu = 0;
		ave_sigma = 0;
		for(int j=0; j<SA_nstep_; j++)
		{	
			if(j%100 == 0) cout << "step(SA_algo) = " << j << endl;
			//move parameters ( the cost is evaluated during the move )
			moveUniform();

			//print parameters values
			output << i << " " << temp_ << " " << getX()[0] << " " << getX()[1] << " " << ene_.getEnergy() << " " << ene_.getStDevMean() << endl;
			ave_mu += getX()[0];
			ave_sigma += getX()[1];

			/*cout << "mu = " << getX()[0] << ", sigma = " << getX()[1] << endl;
			cout << "energy = " << ene_.getEnergy() << ", StDevMean = " << ene_.getStDevMean() << endl;
			cout << "Acceptance rate (energy) = " << ene_.getAccRate() << endl << endl;*/
		}

		ave_mu = ave_mu/(double)SA_nstep_;
		ave_sigma = ave_sigma/(double)SA_nstep_;

		sum_mu += ave_mu;
		sum_sigma += ave_sigma;

		ave2_mu = sum_mu/(double)(i+1);
		ave2_sigma = sum_sigma/(double)(i+1);

		ene_.setMu( ave2_mu );
		ene_.setSigma( ave2_sigma );
		ene_.calibrateDelta(nstep);
		ene_.computeIntegral();

		cout << "ave_mu = " << ave2_mu << ", ave_sigma = " << ave2_sigma << endl << endl;
		cout << "energy = " << ene_.getEnergy() << ", StDevMean = " << ene_.getStDevMean() << endl;
		//cout << "Acceptance rate (energy) = " << ene_.getAccRate() << endl << endl;
		cout << "Acceptance rate (SA algorithm) = " << getAccRate() << endl;
		//output << i << " " << temp_ << " " << ave2_mu << " " << ave2_sigma << " " << ene_.getEnergy() << " " << ene_.getStDevMean() << endl;
		
		//temp_ -= delta_temp;
		temp_ = temp_/2.0;
		SA_delta -= SA_delta_var;
		//setDelta(SA_delta);
	}

	//mu_min_ = sum_mu / (double)n_temp_;
	//sigma_min_ = sum_sigma / (double)n_temp_;
	output.close();
}
		
double SA_Energy :: cost(std::vector<double> x)
{
	ene_.setMu( x[0] );
	ene_.setSigma( x[1] );
	ene_.computeIntegral();
	cost_ = ene_.getEnergy();
	return cost_;
}
