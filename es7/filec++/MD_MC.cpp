/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include "MD_MC.hpp"

using namespace std;

int main(int argc, char* argv[])
{  
	Input(argc, argv); //Inizialization 
	OpenAveragesFiles(); 
	//OpenInstantValuesFiles();
	Equilibrate();
	int nconf = 1;
	iblk=1;
	while(iblk <= nblk) //Simulation
	{
		Reset();   //Reset block averages
		for(int istep=1; istep <= nstep; istep++)
		{
			Move();
			Measure();
			//PrintInstantValues(istep);
			Accumulate(); //Update block averages
			if(istep%10 == 0){
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
				nconf += 1;
			}
		}
		Averages(); //Compute averages and deviations for current block
		PrintAverages(iblk); //Print results for current block
		iblk++;
	}

	if(autocor_on)
	{
		PrintEpotAutoCorrelation();
		PrintPresAutoCorrelation();
	}
	Print_g_r_final_ave();
	CloseAveragesFiles();
	//CloseInstantValuesFiles();
	ConfFinal(); //Write final configuration
	return 0;
}

void Input(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "error: write ':/MD_MC.exe state' where state can be 'solid' or 'liquid' or 'gas'" << endl;
		exit(-1);
	}
    
	state = argv[1];
  
	if(state != "solid" && state != "liquid" && state != "gas" )
	{
		cout << "error: you may have mispelled the state of matter" << endl;
		exit(-1);
	}

	string input_name = "/mnt/c/Users/giaco/esercizi_python/7/input." + state;
	
	ifstream ReadInput, ReadConf, ReadVelocity, Primes, Seed;

	//Read input informations
	ReadInput.open(input_name);

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "MD(NVE) / MC(NVT) simulation       " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "Boltzmann weight exp(- beta * sum_{i<j} v(r_ij) ), beta = 1/T " << endl;
	cout << "The program uses Lennard-Jones units " << endl << endl;

	ReadInput >> iNVET;
	ReadInput >> restart;
	ReadInput >> temp_0;
	ReadInput >> npart;
	ReadInput >> rho;
	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nblk;
	ReadInput >> nstep;
	ReadInput >> nstep_equilibrium;
	ReadInput >> nbins;
	ReadInput >> autocor_on;

	if(iNVET == 1)
	{
		cout << "The program performs Metropolis moves with uniform translations" << endl << endl;
		algo_name = "metro";
	}
	else
	{
		cout << "The program performs Verlet moves" << endl << endl;
		algo_name = "verlet";
	}

	if(restart == 1)
	{
		cout << "The simulation restarts from a new configuration" << endl;
	}else
	{
		cout << "The simulation starts from a previous configuration and doesn't overwrite the results" << endl;
	}

	cout << "Temperature                                  = " << temp_0              << endl;
	cout << "Number of particles                          = " << npart             << endl;
	cout << "Density of particles                         = " << rho               << endl;

	vol = (double)npart/rho;
	box = pow(vol,1.0/3.0);
	beta = 1.0/temp_0;

	cout << "Volume of the simulation box                 = " << vol               << endl;
	cout << "Edge of the simulation box                   = " << box               << endl;
	cout << "Cutoff of the interatomic potential          = " << rcut              << endl;
	cout << "Moves parameter                              = " << delta             << endl;
	cout << "Number of blocks                             = " << nblk              << endl;
	cout << "Number of steps in one block                 = " << nstep             << endl;
	cout << "Number of steps necessary before equilibrium = " << nstep_equilibrium << endl;
	cout << "number of bins for the computation of g(r)   = " << nbins             << endl << endl;

	if(autocor_on)
	{
		cout << "The program computes autocorrelation of the potential energy" << endl << endl;
	}else
	{
		cout << "The program doesn't compute the autocorrelation" << endl << endl;
	}

	ReadInput.close();

	//Prepare arrays for measurements
	iv = 0; //Potential energy
	ip = 1; //Pressure
	if(iNVET == 1)
	{
		ig = 2;
		n_props = 2; //Number of observables
	}else
	{
		it = 2;
		ig = 3;
		n_props = 3;
	}

	//g(r)
	n_props += nbins;
	bin_size = box/2.0/(double)nbins;
	double r;
	for(int i=0; i<nbins; i++)
	{
		r = i*bin_size;
		deltaVol[i] = 4*pi*( pow(r + bin_size, 3) - pow(r, 3) )/3;
	}

	vtail = 8*pi*rho*( 1.0/9.0/pow(rcut, 9) - 1.0/3.0/pow(rcut, 3) ); //v tail correction per particle (division is left to right associative)
	ptail = -32*pi*rho*rho*( 1.0/9.0/pow(rcut, 9) - 1.0/6.0/pow(rcut, 3) ); //p tail correction

	cout << "tail correction to the potential energy = " << vtail << endl;
	cout << "tail correction to the pressure         = " << ptail << endl << endl;

	//Read seed for random numbers
	int p1, p2;
	Primes.open("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	if(restart == 0) Seed.open("seed.out");
	else Seed.open("seed.in");
	Seed >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	Seed.close();  

	//Read initial configuration
	cout << "Read initial configuration" << endl;
	if(restart == 0)
	{
		ReadConf.open("config(" + state + ").out");
		if(iNVET == 0)
		{
			ReadVelocity.open("velocity(" + state + ").out");
			for (int i=0; i<npart; ++i) ReadVelocity >> vx[i] >> vy[i] >> vz[i];
			ReadVelocity.close();
		}
	}
	else 
	{
		ReadConf.open("config.in");
		cout << "Prepare velocities with center of mass velocity equal to zero " << endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for(int i=0; i<npart; ++i)
		{
			vx[i] = rnd.Gauss(0.,sqrt(temp_0));
			vy[i] = rnd.Gauss(0.,sqrt(temp_0));
			vz[i] = rnd.Gauss(0.,sqrt(temp_0));
			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for(int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		double sumv2 = 0.0, fs;
		for(int i=0; i<npart; ++i)
		{
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];
			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;
		fs = sqrt(3 * temp_0 / sumv2);   // fs = velocity scale factor 
		cout << "velocity scale factor: " << fs << endl << endl; //togliendo la velocità media la temperatura si è ridotta devo quindi riscalare le velocità
		for(int i=0; i<npart; ++i)
		{
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;
		}
	}

	for(int i=0; i<npart; ++i)
	{
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = Pbc( x[i] * box ); //moltiplico per box perchè le coordinate in config.in sono scritte in unità del lato della "scatola" così moltiplico per il lato per riottenere le misure giuste
		y[i] = Pbc( y[i] * box );
		z[i] = Pbc( z[i] * box );
	}
	ReadConf.close();

	for(int i=0; i<npart; ++i)
	{
		if(iNVET)
		{
			xold[i] = x[i];
			yold[i] = y[i];
			zold[i] = z[i];
		}
		else
		{
			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		}
	}

	Measure();
	//Print initial values for measured properties
	cout << "Initial potential energy per particle = " << walker[iv]/(double)npart << endl;
	cout << "Initial pressure                      = " << walker[ip]               << endl << endl;
}

void Move()
{
	int o;
	double p, energy_old, energy_new;
	double xnew, ynew, znew;

	if(iNVET) // Monte Carlo (NVT) move
	{
		for(int i=0; i<npart; ++i)
		{
			//Select randomly a particle (for C++ syntax, 0 <= o <= npart-1)
			o = (int)(rnd.Rannyu()*npart);

			//Old
			energy_old = Boltzmann(x[o],y[o],z[o],o);

			//New
			x[o] = Pbc( x[o] + delta*(rnd.Rannyu() - 0.5) );
			y[o] = Pbc( y[o] + delta*(rnd.Rannyu() - 0.5) );
			z[o] = Pbc( z[o] + delta*(rnd.Rannyu() - 0.5) );

			energy_new = Boltzmann(x[o],y[o],z[o],o);

			//Metropolis test
			p = exp(beta*(energy_old-energy_new));
			if(p >= rnd.Rannyu())
			{
				//Update
				xold[o] = x[o];
				yold[o] = y[o];
				zold[o] = z[o];
				accepted += 1.0;
			}else
			{
				x[o] = xold[o];
				y[o] = yold[o];
				z[o] = zold[o];
			}
			attempted += 1.0;
		}
	}else // Molecular Dynamics (NVE) move
	{
		double fx[m_part], fy[m_part], fz[m_part];

		for(int i=0; i<npart; ++i) //Force acting on particle i
		{
			fx[i] = Force(i,0);
			fy[i] = Force(i,1);
			fz[i] = Force(i,2);
		}

		for(int i=0; i<npart; ++i) //Verlet integration scheme
		{
			xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
			ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
			znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

			vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
			vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
			vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

			xold[i] = x[i];
			yold[i] = y[i];
			zold[i] = z[i];

			x[i] = xnew;
			y[i] = ynew;
			z[i] = znew;
		}
	}
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r) //forza sulla i-esima particella nella idir direzione
	double f=0.0;
	double dvec[3], dr;

	for(int i=0; i<npart; ++i)
	{
		if(i != ip)
		{
    
			dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut)
			{
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
			}
		}
	}
	return f;
}

void Measure() //Properties measurement
{
	double v = 0.0, kin= 0.0, w = 0.0;
	double dx, dy, dz, dr;
	int g_bin;

	//cycle over pairs of particles
	for(int i=0; i<npart-1; ++i)
	{
		for(int j=i+1; j<npart; ++j)
		{
			// distance i-j in pbc
			dx = Pbc(x[i] - x[j]);
			dy = Pbc(y[i] - y[j]);
			dz = Pbc(z[i] - z[j]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut)
			{
				w += 2.0/pow(dr,12) - 1.0/pow(dr,6);
				v += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}

			//g(r)
			if(dr < box/2.0)
			{
				g_bin = ig + dr*2*nbins/box;
				walker[g_bin] += 2; // couple of particles found
			}
		}          
	}

	walker[iv] = 4.0 * v / (double)npart + vtail; // Potential energy
	
	if(iNVET == 1)
	{
		walker[ip] = rho*temp_0 + 24.0 * w/(3.0*vol) + ptail;
	}else
	{
		for(int i=0; i<npart; ++i) kin += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
		walker[it] = (2.0 / 3.0) * kin/(double)npart; // Temperature
		walker[ip] = rho*walker[it] + 24.0*w/(3.0*vol) + ptail; //pressure
	}
	if(autocor_on)
	{
		epot_t.push_back(walker[iv]);
		pres_t.push_back(walker[ip]);
	}

	//g(r)
	for(int i=0; i<n_props; i++)
	{
		walker[ig + i] = walker[ig + i]/rho/npart/deltaVol[i];
	}
}

void Reset() //Reset block averages
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

void Averages() //Print results for current block
{		
	stima_epot    = blk_av[iv]/blk_norm; //Potential energy per particle
	glob_av[iv]  += stima_epot;
	glob_av2[iv] += stima_epot*stima_epot;
	err_epot      = Error(glob_av[iv],glob_av2[iv],iblk);

	stima_pres    = blk_av[ip]/blk_norm; //Pressure
	glob_av[ip]  += stima_pres;
	glob_av2[ip] += stima_pres*stima_pres;
	err_pres      = Error(glob_av[ip],glob_av2[ip],iblk);

	if(iNVET == 0)
	{
		stima_temp    = blk_av[it]/blk_norm; //Pressure
		glob_av[it]  += stima_temp;
		glob_av2[it] += stima_temp*stima_temp;
		err_temp      = Error(glob_av[it],glob_av2[it],iblk);
	}

	//g(r)
	for(int i=0; i<nbins; i++)
	{
		stima_g_r[i]     = blk_av[ig + i]/blk_norm;
		glob_av[ig + i]  += stima_g_r[i];
		glob_av2[ig + i] += stima_g_r[i]*stima_g_r[i];
		err_g_r[i]       = Error(glob_av[ig + i], glob_av2[ig + i], iblk); 
	}  
}

void ConfFinal(void)
{
	ofstream WriteConf, WriteVelocity, WriteSeed;

	cout << "Print final configuration to file config.out" << endl << endl;
	WriteConf.open("config.out");
  
	for(int i=0; i<npart; ++i)
	{
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
		WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
	}
	WriteConf.close();

	if(iNVET == 0)
	{
		WriteVelocity.open("velocity.out");
		for(int i=0; i<npart; ++i)
		{
			WriteVelocity << vx[i] << "   " <<  vy[i] << "   " << vz[i] << endl;
		}
		WriteVelocity.close();
	}
	rnd.SaveSeed();
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for(int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double Pbc(double r)  //Algorithm for periodic boundary conditions with side L=box
{
	return r - box * rint(r/box);
}

double Error(double sum, double sum2, int iblk)
{
	return sqrt(fabs(sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void OpenAveragesFiles(void)
{
	if(restart == 1)
	{
		Epot.open("/mnt/c/Users/giaco/esercizi_python/7/ave_epot_" + state + "_" + algo_name + ".dat");
		Pres.open("/mnt/c/Users/giaco/esercizi_python/7/ave_pres_" + state + "_" + algo_name + ".dat");
		//Isto_g_r_block_ave.open("/mnt/c/Users/giaco/esercizi_python/7/g_r_block_ave_" + state + "_" + algo_name + ".dat");

		if(iNVET == 0) Temp.open("/mnt/c/Users/giaco/esercizi_python/7/ave_temp_" + state + "_" + algo_name + ".dat");
	}else
	{
		Epot.open("/mnt/c/Users/giaco/esercizi_python/7/ave_epot_" + state + "_" + algo_name + ".dat", ios::app);
		Pres.open("/mnt/c/Users/giaco/esercizi_python/7/ave_pres_" + state + "_" + algo_name + ".dat", ios::app);
		//Isto_g_r_block_ave.open("/mnt/c/Users/giaco/esercizi_python/7/g_r_block_ave_" + state + "_" + algo_name + ".dat", ios::app);

		if(iNVET == 0) Temp.open("/mnt/c/Users/giaco/esercizi_python/7/ave_temp_" + state + "_" + algo_name + ".dat", ios::app);
	}
}

void PrintAverages(int iblk)
{
	cout << "Block number = " << iblk << endl;
	if(iNVET) cout << "Acceptance rate = " << accepted/(double)attempted << endl << endl;
	
	Epot << setw(wd) << iblk <<  setw(wd) << stima_epot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_epot << endl;
	Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[ip]/(double)iblk << setw(wd) << err_pres << endl;
	if(iNVET == 0) Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
}

void CloseAveragesFiles(void)
{		
	Epot.close();
	Pres.close();
	if(iNVET == 0) Temp.close();
}

void Print_g_r_final_ave(void)
{
	Isto_g_r_final_ave.open("/mnt/c/Users/giaco/esercizi_python/7/g_r_final_ave_" + state + "_" + algo_name + ".dat");
	for(int i=0; i<nbins; i++)
	{
		Isto_g_r_final_ave << setw(wd) << i*bin_size + bin_size/2.0 << setw(wd) << glob_av[ig + i]/(double)nblk << setw(wd) << err_g_r[i] << endl;
	}
	Isto_g_r_final_ave.close();
}

void OpenInstantValuesFiles(void)
{
	if(restart == 1)
	{
		Epot.open("/mnt/c/Users/giaco/esercizi_python/7/epot_" + state + "_" + algo_name + ".dat");
		Pres.open("/mnt/c/Users/giaco/esercizi_python/7/pres_" + state + "_" + algo_name + ".dat");
	}else
	{
		Epot.open("/mnt/c/Users/giaco/esercizi_python/7/epot_" + state + "_" + algo_name + ".dat", ios::app);
		Pres.open("/mnt/c/Users/giaco/esercizi_python/7/pres_" + state + "_" + algo_name + ".dat", ios::app);
	}
}

void PrintInstantValues(int istep)
{	    
	if(istep % 1000 == 0) cout << "istep = " << istep << endl;
	Epot << setw(wd) << walker[iv] << endl;
	Pres << setw(wd) << walker[ip] << endl;
}

void CloseInstantValuesFiles(void)
{		
	Epot.close();
	Pres.close();
}

void Equilibrate(void)
{
	for(int i=0; i<nstep_equilibrium; i++)
	{
		Move();
		if(i % 10000 == 0) cout << "istep_equilibrium = " << i << endl;
	}
	//reinitialize Epot_t
	vector<double> v;
	epot_t = v;
	pres_t = v;
}

void PrintEpotAutoCorrelation(void)
{
	cout << "computing autocorrelation of the potential energy" << endl;
	ofstream output("/mnt/c/Users/giaco/esercizi_python/7/epot_" + state + "_" + algo_name + "_autocorrelation" + ".txt");

	int tmax = nstep*nblk;
	double autocor_epot;

	double sum1, sum2, sum3, var_epot;
	double  epot_sqrd = 0;
	for(int i=0; i<tmax; i++)
	{
		epot_sqrd += epot_t[ i ] * epot_t[ i ]; 
	}
	var_epot = epot_sqrd/(double)tmax - pow( glob_av[iv]/(double)nblk, 2 );
	
	for(int t=1; t < tmax - tmax/100.0; t++)
	{
		if( t % 10000 == 0) cout << "t=" << t << endl;
		
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;

		for(int t2=0; t2 < tmax - t; t2++)
		{
			sum1 += epot_t[ t2 ] * epot_t[ t2 + t ];
			sum2 += epot_t[ t2 ];
			sum3 += epot_t[ t2 + t ];
		}
		autocor_epot = ( sum1/(double)(tmax - t) - sum2*sum3/(double)( pow(tmax - t, 2) ) ) / var_epot;
		output << setw(wd) << t << setw(wd) << autocor_epot << endl;
	}	
	output.close();
}

void PrintPresAutoCorrelation(void)
{
	cout << "computing autocorrelation of the pressure" << endl;
	ofstream output("/mnt/c/Users/giaco/esercizi_python/7/pres_" + state + "_" + algo_name + "_autocorrelation" + ".txt");

	int tmax = nstep*nblk;
	double autocor_pres;

	double sum1, sum2, sum3, var_pres;
	double  pres_sqrd = 0;
	for(int i=0; i<tmax; i++)
	{
		pres_sqrd += pres_t[ i ] * pres_t[ i ]; 
	}
	var_pres = pres_sqrd/(double)tmax - pow( glob_av[ip]/(double)nblk, 2);
	
	for(int t=1; t < tmax - tmax/100.0; t++)
	{
		if( t % 10000 == 0) cout << "t=" << t << endl;
		
		sum1 = 0;
		sum2 = 0;
		sum3 = 0;

		for(int t2=0; t2 < tmax - t; t2++)
		{
			sum1 += pres_t[ t2 ] * pres_t[ t2 + t ];
			sum2 += pres_t[ t2 ];
			sum3 += pres_t[ t2 + t ];
		}
		autocor_pres = ( sum1/(double)(tmax - t) - sum2*sum3/(double)( pow(tmax - t, 2) ) ) / var_pres;
		output << setw(wd) << t << setw(wd) << autocor_pres << endl;
	}	
	output.close();
}

double Boltzmann(double xx, double yy, double zz, int ip)
{
	double ene=0.0;
	double dx, dy, dz, dr;

	for(int i=0; i<npart; ++i)
	{
		if(i != ip)
		{
			// distance ip-i in pbc
			dx = Pbc(xx - x[i]);
			dy = Pbc(yy - y[i]);
			dz = Pbc(zz - z[i]);

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			if(dr < rcut)
			{
				ene += 1.0/pow(dr,12) - 1.0/pow(dr,6);
			}
		}
	}
	return 4.0*ene;
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

