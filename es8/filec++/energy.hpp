#ifndef __Energy__
#define __Energy__

#include "simAnnealing.hpp"

class Energy : public Metropolis
{
	public:
		Energy();
		Energy(double delta, int equilibration_time);

		double getMu() const{return mu_;}
		double getSigma() const{return sigma_;}
		double getEnergy() const{return energy_;}
		double getStDevMean() const{return st_dev_mean_;}
		int getNblk() const{return nblk_;}
		int getNstep() const{return nstep_;}
		
		void setMu(double mu);
		void setSigma(double sigma);
		void setNblk(int nblk);
		void setNstep(int nstep);
		void setX(std::vector<double> x);
		void setEnergy(double energy);
		void setStDevMean(double StDev);

		double probability(std::vector<double> x) override;

		//void moveUniform() override; //this function is not necessary
		
		double measure() const;
		void computeIntegral();
		void printIntegralStDevMeanPerBlock(const std::string Filename);
		void printIstoPsi(const std::string Filename);
	private:		
		double mu_;
		double sigma_;

		int nblk_;
		int nstep_;
				
		double energy_;
		double st_dev_mean_;
};

class SA_Energy : public SimulatedAnnealing
{
	public:
		SA_Energy(double delta_SA, int equilibration_time_SA, double T_min, double T_max, int n_temp, Energy ene);

		double getMuMin() const{return mu_min_;}
		double getSigmaMin() const{return sigma_min_;}
		
		void moveUniform() override;
		void PrintEnergy_StDev_and_Parameters(const std::string Filename);
		
		double cost(std::vector<double> x) override; 
	private:
		Energy ene_;
		double mu_min_;
		double sigma_min_;		
};

#endif // __Energy__
