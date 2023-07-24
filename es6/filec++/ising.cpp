/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "ising.hpp"

using namespace std;

int main()
{ 
	Input(); //Inizialization
	OpenFinalFiles();
	for(int ipnt=1; ipnt <= n_temp; ++ipnt)
	{
		ModifyTemp(ipnt);
		if(autocor_on) OpenBlockFiles();
		Equilibrate();    
		for(int iblk = 1; iblk <= nblk; ++iblk)
		{
			Reset(iblk);   //Reset block averages
			for(int istep=1; istep <= nstep; ++istep)
			{
				Move();
				Measure();
				Accumulate(); //Update block averages
			}
			Averages(iblk);   //Print results for current block
		}
		//ConfFinal(); //Write final configuration
    
		PrintFinalValues(nblk);
  	
		if(autocor_on)
		{
			CloseBlockFiles();
			PrintEnergyAutoCorrelation();
			PrintMagnetizationAutoCorrelation();
			PrintMagnetSusceptibilityAutoCorrelation();
		}
	}
	CloseFinalFiles();
	return 0;
}

void Input(void)
{
	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;

	//Read input informations
	ifstream ReadInput, Primes, Seed;
	ReadInput.open("/mnt/c/Users/giaco/esercizi_python/6/input.dat");
  
	ReadInput >> restart;
	ReadInput >> nspin;
	ReadInput >> J;
	ReadInput >> h;
	ReadInput >> metro;
	ReadInput >> nblk;
	ReadInput >> nstep;
	ReadInput >> temp_min;
	ReadInput >> temp_max;
	ReadInput >> n_temp;
	ReadInput >> nstep_equilibrium;
	ReadInput >> autocor_on;

	ReadInput.close();

	cout << "Number of spins                        = " << nspin             << endl;
	cout << "Exchange interaction                   = " << J                 << endl;
	cout << "External field                         = " << h                 << endl;
	cout << "Number of blocks                       = " << nblk              << endl;
	cout << "Number of steps in one block           = " << nstep             << endl;
	cout << "Minimum temperature                    = " << temp_min          << endl;
	cout << "Maximum temperature                    = " << temp_max          << endl;
	cout << "Number of temperatures used            = " << n_temp            << endl;  
	cout << "Number of steps needed for equilibrium = " << nstep_equilibrium << endl << endl;

	if(metro==1)
	{
		algo_name = "metro";
		cout << "The program performs Metropolis moves" << endl;
	}
	else
	{
		algo_name = "gibbs";
		cout << "The program performs Gibbs moves" << endl;
	}
  
	if(autocor_on)
	{
		cout << "The program computes the autocorrelation of energy, magnetization and magnetic susceptibility and print block values" << endl << endl;
	}else
	{
		cout << "The program doesn't compute the autocorrelation and doesn't print block values" << endl << endl;
	}
 
	if(restart)
	{
		cout << "The program starts from a new random configuration" << endl;
		Seed.open("seed.in");
	}else
	{
		cout << "the program starts from a previous configuration" << endl;
		Seed.open("seed.out");
	}

	//Read seed for random numbers
	int p1, p2;
	Primes.open("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	Seed.close();

	//initial configuration
	if(restart)
	{
		for (int i=0; i<nspin; ++i)
		{
			if(rnd.Rannyu() >= 0.5) s[i] = 1;
			else s[i] = -1;
		}		
	}else
	{
		ReadConf();
	}  
  
	//Evaluate energy etc. of the initial configuration
	Measure();
 
	//Prepare arrays for measurements
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
 
	n_props = 4; //Number of observables

	temp = temp_max;
	cout << "starting the simulation at temperature = " << temp << endl; 
}

void Move()
{
	int o;
	double p, delta_u;

	for(int i=0; i<nspin; ++i)
	{
		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rnd.Rannyu()*nspin);
		attempted++;
	
		if(metro==1) //Metropolis
		{		
			delta_u = 2*J*s[o]*( s[Pbc(o-1)] + s[Pbc(o+1)] ) + 2*h*s[o];	
			p =  exp(-beta*delta_u);
			
			if( rnd.Rannyu() < p )
			{
				s[o] *= -1;
				accepted++;
			}
		}
		else //Gibbs sampling
		{
			p = 1 / ( 1 + exp( 2*beta*J*( s[Pbc(o-1)] + s[Pbc(o+1)] ) + 2*beta*h ) );
			if(rnd.Rannyu() < p)
			{
				s[o] = -1;
			}else
			{
				s[o] = +1;
			}
			accepted++;
		}
	}
}

double Boltzmann(int sm, int ip)
{
	double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
	return ene;
}

void Measure()
{
	double u = 0.0, m = 0.0;

	for (int i=0; i<nspin; ++i)
	{
		u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);
		m += s[i];
	}
	walker[iu] = u;
	walker[ic] = u*u;
	walker[im] = m;
	walker[ix] = m*m;

	if(autocor_on)
	{
		u_t.push_back(u);
		m_t.push_back(m);
		m_t.push_back(m*m);
	}
}

void Reset(int iblk) //Reset block averages
{
	if(iblk == 1)
	{
		for(int i=0; i<n_props; ++i)
		{
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; ++i)
	{
		blk_av[i] = 0;
	}
	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}

void Accumulate(void) //Update block averages
{
	for(int i=0; i<n_props; ++i)
	{
		blk_av[i] = blk_av[i] + walker[i];
	}
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk)
{    
	stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy (per particle)
	stima_c = beta*beta*( blk_av[ic] - blk_av[iu]*blk_av[iu]/blk_norm )/blk_norm/(double)nspin;//heat capacity (per particle)
	stima_m = blk_av[im]/blk_norm/(double)nspin; //magnetization ( per particle)
	//stima_x = beta*( blk_av[ix] - blk_av[im]*blk_av[im]/blk_norm )/blk_norm; //magnetic susceptibility
	stima_x = beta*( blk_av[ix] - blk_av[im]*blk_av[im]/blk_norm )/blk_norm/(double)nspin; //magnetic susceptibility (per particle)

	//stima_x = beta*( blk_av[ix]/blk_norm/(double)nspin - stima_m );//com'era scritta prima
	//stima_c = beta*beta*( blk_av[ic]/blk_norm/(double)nspin - stima_u ); //com'era scritta prima

	glob_av[iu]  += stima_u;
	glob_av2[iu] += stima_u*stima_u;
	err_u         = Error(glob_av[iu], glob_av2[iu], iblk);  
    
	glob_av[ic]  += stima_c;
	glob_av2[ic] += stima_c*stima_c;
	err_c         = Error(glob_av[ic], glob_av2[ic], iblk);
	
	glob_av[im]  += stima_m;
	glob_av2[im] += stima_m*stima_m;
	err_m         = Error(glob_av[im], glob_av2[im], iblk);
	  
	glob_av[ix]  += stima_x;
	glob_av2[ix] += stima_x*stima_x;
	err_x         = Error(glob_av[ix], glob_av2[ix], iblk);
	
	if(autocor_on) PrintBlockValues(iblk);

	//cout << "iblk=" << iblk << ", acceptance rate=" << accepted/(double)attempted << endl;	
}

void ConfFinal(void)
{
	ofstream WriteConf;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("final_config.txt");
	for (int i=0; i<nspin; ++i)
	{
		WriteConf << s[i] << endl;
	}
	WriteConf.close();

	rnd.SaveSeed();
}

void ReadConf(void)
{
	ifstream ReadConf("config.in");

	int numLines = 0;
	string unused;
	while( getline(ReadConf, unused) )
	{
		numLines++;
	}

	if(numLines != nspin)
	{
		cout << "tried to open previous configuration but the number of spins in the file and the number of spins written on the input file don't match" << endl;
		exit(-1);
	}
	
	ReadConf.clear();
	ReadConf.seekg(0);
	
	for(int i=0; i<nspin; ++i)
	{
		ReadConf >> s[i];
	}
}

int Pbc(int i)  //Algorithm for periodic boundary conditions
{
	if(i >= nspin) i = i - nspin;
	else if(i < 0) i = i + nspin;
	return i;
}

double Error(double sum, double sum2, int iblk)
{
	if(iblk==1) return 0.0;
	else return sqrt( ( sum2/(double)iblk - pow(sum/(double)iblk,2) )/(double)(iblk-1) );
}

void ModifyTemp(int ipnt)
{	
	if(ipnt != 1)
	{
		temp -= (temp_max - temp_min)/(double)(n_temp-1);
		cout << "Restarting the simulation with Temperature = " << temp << endl;
	}
	beta = 1.0/temp;
}

void Equilibrate(void)
{
	for(int i=0; i<nstep_equilibrium; i++)
	{
		Move();
	}
	//reinitialize u_t
	vector<double> v;
	u_t = v;
	m_t = v;
}
	 
void PrintEnergyAutoCorrelation(void)
{
	cout << "computing autocorrelation of the energy" << endl;
	ofstream output("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/energy_" + algo_name + "_autocorrelation_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");

	int tmax = nstep*nblk;
	double autocor_u;

	double sum1, sum2, sum3, var_u;
	double  u_sqrd = 0;
	for(int i=0; i<tmax; i++)
	{
		u_sqrd += u_t[ i ] * u_t[ i ]; 
	}
	var_u = u_sqrd/(double)tmax - pow( glob_av[iu]/(double)nblk, 2 );
	
	for(int t=1; t < tmax - nstep; t++)
	{
		if( t % 10000 == 0) cout << "t=" << t << endl;
		
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;

		for(int t2=0; t2 < tmax - t; t2++)
		{
			sum1 += u_t[ t2 ] * u_t[ t2 + t ];
			sum2 += u_t[ t2 ];
			sum3 += u_t[ t2 + t ];
		}
		autocor_u = ( sum1/(double)(tmax - t) - sum2*sum3/(double)( pow(tmax - t, 2) ) ) / var_u;
		output << setw(wd) << t << setw(wd) << autocor_u << endl;
	}	
	output.close();
}

void PrintMagnetizationAutoCorrelation(void)
{
	cout << "computing autocorrelation of the magnetization" << endl;
	ofstream output("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/magnetization_" + algo_name + "_autocorrelation_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");
	
	int tmax = nstep*nblk;
	double autocor_m;

	double sum1, sum2, sum3, var_m;
	double  m_sqrd = 0;
	for(int i=0; i<tmax; i++)
	{
		m_sqrd += m_t[ i ] * m_t[ i ]; 
	}
	var_m = m_sqrd/(double)tmax - pow( glob_av[im]/(double)nblk, 2 );
	
	for(int t=1; t < tmax - nstep; t++)
	{
		if( t % 10000 == 0) cout << "t=" << t << endl;
		
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;

		for(int t2=0; t2 < tmax - t; t2++)
		{
			sum1 += m_t[ t2 ] * m_t[ t2 + t ];
			sum2 += m_t[ t2 ];
			sum3 += m_t[ t2 + t ];
		}
		autocor_m = ( sum1/(double)(tmax - t) - sum2*sum3/(double)( pow(tmax - t, 2) ) ) / var_m;
		output << setw(wd) << t << setw(wd) << autocor_m << endl;
	}	
	output.close();
}

void PrintMagnetSusceptibilityAutoCorrelation()
{
	cout << "computing autocorrelation of magnetic susceptibility" << endl;
	ofstream output("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/magnetic_susceptibility_" + algo_name + "_autocorrelation_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");
	
	int tmax = nstep*nblk;
	double autocor_x;

	double sum1, sum2, sum3, var_x;
	double  x_sqrd = 0;
	for(int i=0; i<tmax; i++)
	{
		x_sqrd += x_t[ i ] * x_t[ i ]; 
	}
	var_x = x_sqrd/(double)tmax - pow( glob_av[ix]/(double)nblk, 2 );
	
	for(int t=1; t < tmax - nstep; t++)
	{
		if( t % 10000 == 0) cout << "t=" << t << endl;
		
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;

		for(int t2=0; t2 < tmax - t; t2++)
		{
			sum1 += x_t[ t2 ] * x_t[ t2 + t ];
			sum2 += x_t[ t2 ];
			sum3 += x_t[ t2 + t ];
		}
		autocor_x = ( sum1/(double)(tmax - t) - sum2*sum3/(double)( pow(tmax - t, 2) ) ) / var_x;
		output << setw(wd) << t << setw(wd) << autocor_x << endl;
	}	
	output.close();
}

void OpenBlockFiles(void)
{
	Ene_T.open("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/energy_" + algo_name + "_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");
	Heat_T.open("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/Heat_capacity_" + algo_name + "_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");
	Chi_T.open("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/magnetic_susceptibility_" + algo_name + "_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");  
	Mag_T.open("/mnt/c/Users/giaco/esercizi_python/6/values_per_block/magnetization_" + algo_name + "_T_" + to_string(temp) + "(h=" + to_string(h) + ").txt");
	
	Ene_T << setw(wd) << temp << setw(wd) << nblk << setw(wd) << nstep << setw(wd) << 0 <<endl;
	Heat_T << setw(wd) << temp << setw(wd) << nblk << setw(wd) << nstep << setw(wd) << 0 <<endl;
	Mag_T << setw(wd) << temp << setw(wd) << nblk << setw(wd) << nstep << setw(wd) << 0 <<endl;
	Chi_T << setw(wd) << temp << setw(wd) << nblk << setw(wd) << nstep << setw(wd) << 0 <<endl;		
}

void OpenFinalFiles(void)
{
 	Ene.open("/mnt/c/Users/giaco/esercizi_python/6/energy_" + algo_name + "(h=" + to_string(h) + ").txt");
 	Heat.open("/mnt/c/Users/giaco/esercizi_python/6/heat_capacity_" + algo_name + "(h=" + to_string(h) + ").txt");
 	Chi.open("/mnt/c/Users/giaco/esercizi_python/6/magnetic_susceptibility_" + algo_name + "(h=" + to_string(h) + ").txt");
 	Mag.open("/mnt/c/Users/giaco/esercizi_python/6/magnetization_" + algo_name + "(h=" + to_string(h) + ").txt");
}

void CloseBlockFiles(void)
{		
	Ene_T.close();
	Heat_T.close();
	Mag_T.close();
	Chi_T.close();
}

void CloseFinalFiles(void)
{
	Ene.close();
	Heat.close();
	Mag.close();
	Chi.close();
}

void PrintBlockValues(int iblk)
{
	Ene_T << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	Heat_T << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
	Mag_T << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
	Chi_T << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
}

void PrintFinalValues(int iblk)
{
	Ene << setw(wd) << temp << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
	Heat << setw(wd) << temp << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
	Mag << setw(wd) << temp << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
	Chi << setw(wd) << temp << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;	

	CloseBlockFiles();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
