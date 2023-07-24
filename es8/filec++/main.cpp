#include "energy.hpp"

using namespace std;

int main(){
	int nblk, nstep, nstep_equilibrium, SA_nstep;
	double d, delta_SA, mu, sigma, T_min, T_max, n_temp;
	
	ifstream ReadInput("/mnt/c/Users/giaco/esercizi_python/8/input.txt");
	
	ReadInput >> nblk;
	ReadInput >> nstep;	
	ReadInput >> d;	
	ReadInput >> delta_SA;
	ReadInput >> SA_nstep;
	ReadInput >> nstep_equilibrium;
	ReadInput >> mu;
	ReadInput >> sigma;
	ReadInput >> T_min;
	ReadInput >> T_max;
	ReadInput >> n_temp;

	//T_min = pow(10, -10);
	ReadInput.close();

	cout << "nblk              = " << nblk              << endl;
	cout << "nstep             = " << nstep             << endl;
	cout << "d                 = " << d                 << endl;
	cout << "delta_SA          = " << delta_SA          << endl;
	cout << "SA_nstep          = " << SA_nstep          << endl;
	cout << "nstep_equilibrium = " << nstep_equilibrium << endl;
	cout << "mu                = " << mu                << endl;
	cout << "sigma             = " << sigma             << endl;
	cout << "T_min             = " << T_min             << endl;
	cout << "T_max             = " << T_max             << endl;
	cout << "n_temp            = " << n_temp            << endl;
		
	Energy qmc(d, nstep_equilibrium);
	
	qmc.SetRandomPrimes();
	qmc.setMu(mu);
	qmc.setSigma(sigma);
	qmc.setNblk(nblk);
	qmc.setNstep(nstep);

	vector<double> q(1);
	q[0] = 0;
	qmc.setX(q);

	cout << "esercizio 8.1" << endl;
	/*cout <<"delta="<<qmc.getDelta() <<"; nstep_equi="<<qmc.getEquilibrationTime()<<"; accepted="<<qmc.getAccepted()<<endl;
	cout<<"attempted="<<qmc.getAttempted()<<"; x.size()="<<qmc.getX().size()<<"; x="<<qmc.getX()[0]<<endl;
	cout<<"mu="<<qmc.getMu()<<"; sigma="<<qmc.getSigma()<<"; nblk="<<qmc.getNblk()<<endl;
	cout<<"nstep="<<qmc.getNstep()<<"; energy="<<qmc.getEnergy()<<"; stDev="<<qmc.getStDevMean()<<endl;*/
	qmc.printIntegralStDevMeanPerBlock( "/mnt/c/Users/giaco/esercizi_python/8/H_mean.txt" );

/*
	cout << "prova"<<endl;
	ofstream output("/mnt/c/Users/giaco/esercizi_python/8/prova.txt");
	double e;
	cout << "delta, prima =" << qmc.getDelta() << endl;
	qmc.equilibrate();
	qmc.calibrateDelta(nstep);
	cout << "delta, dopo =" << qmc.getDelta() << endl;
	for(int i=0; i<nblk; i++)
	{
		cout << "i="<<i << endl;
		qmc.setAccepted(0);
		qmc.setAttempted(0);

		for(int j=0; j<nstep; j++)
		{
			qmc.moveUniform();
			e = qmc.measure();
			output << e << endl;
		}
	}
	output.close();*/	

	//cout << "esercizio8.2" << endl;
	SA_Energy sa_qmc(delta_SA, 0, T_min, T_max, n_temp, qmc);
	
	sa_qmc.SetRandomPrimes();
	
	vector<double> x(2);
	x[0] = mu;
	x[1] = sigma;
	sa_qmc.setX(x);
	sa_qmc.setSA_nstep(SA_nstep);

	sa_qmc.PrintEnergy_StDev_and_Parameters("/mnt/c/Users/giaco/esercizi_python/8/SA_H.txt");
	
	double mu_min = sa_qmc.getMuMin();
	double sigma_min = sa_qmc.getSigmaMin();

	cout << "mu_min = " << mu_min << endl;
	cout << "sigma_min = " << sigma_min << endl;
	
	Energy qmc_min(d, nstep_equilibrium);

	qmc_min.SetRandomPrimes();
	qmc_min.setMu(mu_min);
	qmc_min.setSigma(sigma_min);
	qmc_min.setNblk(nblk);
	qmc_min.setNstep(nstep);
	qmc_min.setX(q);
	cout << "q="<<q[0]<<endl;
	
	qmc_min.printIntegralStDevMeanPerBlock( "/mnt/c/Users/giaco/esercizi_python/8/H_min.txt" );
	
	qmc_min.printIstoPsi( "/mnt/c/Users/giaco/esercizi_python/8/isto_psi_min.txt" );
	
	return 0;
}
