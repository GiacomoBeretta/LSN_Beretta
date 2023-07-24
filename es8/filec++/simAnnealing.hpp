#ifndef __Simulated_Annealing__
#define __Simulated_Annealing__

#include "metropolis.hpp"

class SimulatedAnnealing : public Metropolis
{
	public:
		SimulatedAnnealing();
		SimulatedAnnealing(double delta, int equilibration_time, double T_min, double T_max, int n_temp);

		double getTMin() const{return T_min_;}
		double getTMax() const{return T_max_;}
		int getNTemp() const{return n_temp_;}
		double getTemp() const{return temp_;}
		double getCostValue() const{return cost_;}
		int getSA_nstep() const{return SA_nstep_;}

		void setTMin(double T_min);
		void setTMax(double T_max);
		void setNTemp(double n_temp);
		void setTemp(double temp);
		void setSA_nstep(int SA_nstep);
		
		double probability(std::vector<double> x) override;
		void moveUniform() override;
		virtual double cost(std::vector<double> x) = 0;
		void PrintCostAndParameters(const std::string Filename);
	protected:
		double T_min_;
		double T_max_;
		int n_temp_;
		double temp_;
		int SA_nstep_;

		double cost_;
};

#endif //__Simulated_Annealing__

