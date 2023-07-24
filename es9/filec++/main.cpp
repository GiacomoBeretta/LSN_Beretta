#include "GenAlgTSP.hpp"

using namespace std;

int main()
{
	double prob_mutation, prob_crossover;
	int n_cities, pop_size, norm_power, n_generations, restart;
	string config;
	ofstream output_best, output_average, output_xy, output_paths;

	Random rand;
	rand.SetRandomPrimes();

	ifstream input("/mnt/c/Users/giaco/esercizi_python/9/input.txt");

	input >> prob_mutation;
	input >> prob_crossover;
	input >> n_cities;
	input >> pop_size;
	input >> norm_power;
	input >> n_generations;
	input >> restart;
	input >> config;

	input.close();

	if (config == "american_capitals")
	{
		n_cities = 50;
	}

	cout << "Genetic Algorithm for the Travel Salesman Problem with:" << endl;
	cout << "Mutation Probability = " << prob_mutation << endl;
	cout << "Crossing Over Probability = " << prob_crossover << endl;
	cout << "Number of cities = " << n_cities << endl;
	cout << "Number of individuals = " << pop_size << endl;
	cout << "power of the distance = " << norm_power << endl;
	cout << "Number of generations = " << n_generations << endl;
	cout << "configurations of cities = " << config << endl
		 << endl;

	Path best_path(n_cities);
	Path path(n_cities);
	GenAlgTSP ts(pop_size, path, prob_mutation, prob_crossover, norm_power, rand);

	if (restart == 1)
	{
		cout << "the simulation starts from a new random configuration" << endl;
		if (config == "circumference")
		{
			path.setCitiesRandomlyOnCircumference(rand);
		}
		else if (config == "square")
		{
			path.setCitiesRandomlyInSquare(rand);
		}
		else
		{
			path.setAmericanCapitals("American_capitals.dat");
		}
		path.checkPath();
		GenAlgTSP ts2(pop_size, path, prob_mutation, prob_crossover, norm_power, rand);
		ts = ts2;
	}
	else
	{
		cout << "the simulation starts from a previous configuration" << endl;
		input.open("/mnt/c/LabSim/es9/paths_" + config + ".txt");
		ts.setFromFile(input);
		input.close();
	}

	output_paths.open("/mnt/c/LabSim/es9/paths_" + config + ".txt");

	if (config == "american_capitals")
	{
		output_best.open("/mnt/c/Users/giaco/esercizi_python/10/best_path_" + config + "_" + to_string(norm_power) + "norm_genetic.txt");
		output_average.open("/mnt/c/Users/giaco/esercizi_python/10/average_length_" + config + "_" + to_string(norm_power) + "norm_genetic.txt");
		output_xy.open("/mnt/c/Users/giaco/esercizi_python/10/best_path_" + config + "_xy_" + to_string(norm_power) + "norm_genetic.txt");
	}
	else
	{
		output_best.open("/mnt/c/Users/giaco/esercizi_python/9/best_path_" + config + "_" + to_string(norm_power) + "norm.txt");
		output_average.open("/mnt/c/Users/giaco/esercizi_python/9/average_length_" + config + "_" + to_string(norm_power) + "norm.txt");
		output_xy.open("/mnt/c/Users/giaco/esercizi_python/9/best_path_" + config + "_xy_" + to_string(norm_power) + "norm.txt");
	}

	ts.checkPopulation();
	ts.evolve(n_generations, output_best, output_average, config);
	ts.checkPopulation();
	ts.PrintPathsFile(output_paths);

	best_path = ts.get_pop_vector()[0];
	best_path.printPathFile(output_xy);

	output_best.close();
	output_average.close();
	output_xy.close();
	output_paths.close();

	return 0;
}
