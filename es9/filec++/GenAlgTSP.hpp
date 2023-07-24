#ifndef __GenAlgTSP__
#define __GenAlgTSP__

#include "path.hpp"
#include <algorithm> //sort, min

class GenAlgTSP
{
public:
	GenAlgTSP(const int pop_size, const Path &p, const double p_mut, const double p_cros, const int norm_power, Random &rand);
	GenAlgTSP(const std::vector<Path> &pop, const double p_mut, const double p_cros, const int norm_power, Random &rand);
	GenAlgTSP(const GenAlgTSP &alg, Random &rand);

	int get_n_cities() const { return n_cities_; }
	int get_pop_size() const { return pop_size_; }
	std::vector<Path> get_pop_vector() const { return pop_; }
	double get_p_mut() const { return p_mut_; }
	double get_p_cros() const { return p_cros_; }
	int getNormPower() const { return norm_power_; }
	Random getRandomGen() const { return rand_; }

	void setMutationProb(const double p_mut);
	void setCrossoverProb(const double p_cros);
	void setNormPower(const int norm_power);
	void setFromFile(std::ifstream &in);

	void checkPopulation() const;
	void computeLengths(const std::string config);
	void orderPopulation(); // order in ascending order by L^k
	int select(const double p);
	void mutationSwap();
	// void mutationShift();
	void mutationPermutation();
	void mutationInversion();
	void crossingOver(Path &father, Path &mother);
	void popCrossingOver();

	void reproduce(const std::string config);
	void evolve(const int n_generations, std::ofstream &out1, std::ofstream &out2, std::string config);

	void PrintPathsFile(std::ofstream &out) const;

	GenAlgTSP &operator=(const GenAlgTSP &tsp);

private:
	std::vector<Path> pop_;
	int n_cities_;
	int pop_size_;
	double p_mut_;
	double p_cros_;
	int norm_power_;
	Random rand_;
};

#endif // __GenAlgTSP__
