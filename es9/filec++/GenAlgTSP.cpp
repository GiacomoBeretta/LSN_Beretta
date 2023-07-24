#include "GenAlgTSP.hpp"
#include <sstream>

using namespace std;

GenAlgTSP ::GenAlgTSP(const int pop_size, const Path &p, const double p_mut, const double p_cros, const int norm_power, Random &rand)
{
	n_cities_ = p.get_n_cities();
	pop_size_ = pop_size;
	p_mut_ = p_mut;
	p_cros_ = p_cros;
	norm_power_ = norm_power;
	rand_ = rand;

	vector<Path> pop(pop_size_);
	pop[0] = p;

	for (int i = 1; i < pop_size; i++)
	{
		pop[i] = pop[0].randomOrder(rand_);
	}
	pop_ = pop;
}

GenAlgTSP ::GenAlgTSP(const vector<Path> &pop, const double p_mut, const double p_cros, const int norm_power, Random &rand)
{
	pop_ = pop;
	n_cities_ = pop[0].get_n_cities();
	pop_size_ = pop.size();
	p_mut_ = p_mut;
	p_cros_ = p_cros;
	norm_power_ = norm_power;
	rand_ = rand;
}

GenAlgTSP ::GenAlgTSP(const GenAlgTSP &alg, Random &rand)
{
	n_cities_ = alg.get_n_cities();
	pop_size_ = alg.get_pop_size();
	pop_ = alg.get_pop_vector();
	p_mut_ = alg.get_p_mut();
	p_cros_ = alg.get_p_cros();
	norm_power_ = alg.getNormPower();
	rand_ = rand;
}

void GenAlgTSP ::setMutationProb(const double p_mut)
{
	p_mut_ = p_mut;
}

void GenAlgTSP ::setCrossoverProb(const double p_cros)
{
	p_mut_ = p_cros;
}

void GenAlgTSP ::setNormPower(const int norm_power)
{
	norm_power_ = norm_power;
}

void GenAlgTSP ::setFromFile(ifstream &in)
{
	string line;
	int pop_size = 0;
	while (!in.eof())
	{
		getline(in, line);
		pop_size++;
	}
	if (pop_size - 1 != pop_size_)
	{
		cout << "error:trying to read paths from file, the sizes of the populations don't match, in the file there are " << pop_size - 1 << " individuals and there should be " << pop_size_ << endl;
		exit(-1);
	}

	in.clear();
	in.seekg(0);

	getline(in, line);
	istringstream iss(line);

	int n_cities = 0;
	double dummy;
	while (iss >> dummy)
		n_cities++;
	if (n_cities / 3 - 1 != n_cities_)
	{
		cout << "error:trying to read paths from file, the number of cities don't match, in the file there are " << n_cities / 3 - 1 << " cities and there should be " << n_cities_ << endl;
		exit(-1);
	}

	iss.clear();
	in.clear();
	in.seekg(0);

	Path path(n_cities_);
	vector<int> indexes(n_cities_ + 1);
	double x, y;

	for (int i = 0; i < pop_size_; i++)
	{
		getline(in, line);
		iss.str(line);
		for (int j = 0; j < n_cities_ + 1; j++)
		{
			iss >> indexes[j];
			iss >> x;
			iss >> y;
			path[j].setX(x);
			path[j].setY(y);
		}
		path.setIndexes(indexes);
		pop_[i] = path;
		iss.clear();
	}
}

void GenAlgTSP ::checkPopulation() const
{
	for (int i = 0; i < pop_size_; i++)
	{
		pop_[i].checkPath();
	}
}

void GenAlgTSP ::computeLengths(const string config)
{
	for (int i = 0; i < pop_size_; i++)
	{
		if (config == "american_capitals")
		{
			pop_[i].pathLengthSphere(norm_power_);
		}
		else
		{
			pop_[i].pathLength(norm_power_);
		}
	}
}

void GenAlgTSP ::orderPopulation()
{
	sort(pop_.begin(), pop_.end(), [](Path a, Path b)
		 { return a.getPathLength() < b.getPathLength(); });
}

int GenAlgTSP ::select(const double p)
{
	int j = (pop_size_ - 1) * pow(rand_.Rannyu(), p);
	return j;
}

void GenAlgTSP ::mutationSwap()
{
	Path best = pop_[0];
	int a, b, j;

	for (int i = 0; i < pop_size_; i++)
	{
		j = select(0.5);

		if (rand_.Rannyu() <= p_mut_)
		{
			a = rand_.Integer(1, n_cities_ - 1);
			b = rand_.Integer(1, n_cities_ - 1);
			pop_[j].swapCities(a, b);
		}
	}
	Path temp = pop_[0];
	pop_[0] = best;
	pop_[pop_size_ - 1] = temp;
}

/*
void GenAlgTSP :: mutationShift()
{
	Path best = pop_[0];
	int a,n,shift,j;

	for(int i=0; i<pop_size_; i++)
	{
		j = select(0.5);

		if(rand_.Rannyu() <= p_mut_)
		{
			a = rand_.Integer(1, n_cities_ - 1);
			n = rand_.Integer(1, n_cities_ - 1);
			shift = rand_.Integer(1, n_cities_ - 2);

			pop_[j].shiftCities(a, n, shift);
		}
	}
	Path temp = pop_[0];
	pop_[0] = best;
	pop_[ pop_size_ - 1 ] = temp;
}*/

void GenAlgTSP ::mutationPermutation()
{
	Path best = pop_[0];
	int a, b, d, n, j;

	for (int i = 0; i < pop_size_; i++)
	{
		j = select(0.5);

		if (rand_.Rannyu() <= p_mut_)
		{
			a = rand_.Integer(1, n_cities_ - 1);
			b = rand_.Integer(1, n_cities_ - 1);
			d = boundaryDistance(a, b, n_cities_);

			if (d != 0)
			{
				n = rand_.Integer(1, d);
				pop_[j].permutateCities(a, n, b);
			}
		}
	}
	Path temp = pop_[0];
	pop_[0] = best;
	pop_[pop_size_ - 1] = temp;
}

void GenAlgTSP ::mutationInversion()
{
	Path best = pop_[0];
	int a, n, j;

	for (int i = 0; i < pop_size_; i++)
	{
		j = select(0.5);

		if (rand_.Rannyu() <= p_mut_)
		{
			a = rand_.Integer(1, n_cities_ - 1);
			n = rand_.Integer(1, n_cities_ - 1);
			pop_[j].reverseCities(a, n);
		}
	}
	Path temp = pop_[0];
	pop_[0] = best;
	pop_[pop_size_ - 1] = temp;
}

void GenAlgTSP ::crossingOver(Path &father, Path &mother)
{
	int cut_point, cut_length, k;

	Path father_copy = father;

	vector<int> index_father = father.getIndexes();
	vector<int> index_mother = mother.getIndexes();
	vector<int> index_father_copy = index_father;

	cut_point = rand_.Integer(2, n_cities_ - 2);
	cut_length = n_cities_ - cut_point;

	vector<int> index_father_cut(cut_length);
	vector<int> index_mother_cut(cut_length);

	for (int i = 0; i < cut_length; i++)
	{
		index_father_cut[i] = index_father[cut_point + i];
		index_mother_cut[i] = index_mother[cut_point + i];
	}

	k = 0;
	for (int i = 1; i < n_cities_; i++)
	{
		for (int j = 0; j < cut_length; j++)
		{
			if (index_mother[i] == index_father_cut[j])
			{
				father[cut_point + k] = mother[i];
				index_father[cut_point + k] = index_mother[i];
				k++;
				break;
			}
		}
	}
	k = 0;
	for (int i = 1; i < n_cities_; i++)
	{
		for (int j = 0; j < cut_length; j++)
		{
			if (index_father_copy[i] == index_mother_cut[j])
			{
				mother[cut_point + k] = father_copy[i];
				index_mother[cut_point + k] = index_father_copy[i];
				k++;
				break;
			}
		}
	}

	father.setIndexes(index_father);
	mother.setIndexes(index_mother);
}

void GenAlgTSP ::popCrossingOver()
{
	Path best = pop_[0];

	for (int i = 0; i < pop_size_; i++)
	{
		if (rand_.Rannyu() <= p_cros_)
		{
			crossingOver(pop_[select(2)], pop_[select(2)]);
		}
	}
	Path temp = pop_[0];
	pop_[0] = best;
	pop_[pop_size_ - 2] = temp;
}

void GenAlgTSP ::reproduce(const string config)
{
	popCrossingOver();
	computeLengths(config);
	orderPopulation();

	mutationSwap();
	computeLengths(config);
	orderPopulation();

	// mutationShift();

	mutationPermutation();
	computeLengths(config);
	orderPopulation();

	mutationInversion();
	computeLengths(config);
	orderPopulation();
}

void GenAlgTSP ::evolve(const int n_generations, ofstream &out1, ofstream &out2, const string config)
{
	double R = 6373000; // earth radius in meters
	double sum;
	int half_n_pop = pop_size_ / 2;
	for (int i = 0; i < n_generations; i++)
	{
		sum = 0;
		if (i % 1000 == 0)
			cout << "generation #" << i << endl;
		reproduce(config);

		for (int j = 0; j < half_n_pop; j++)
		{
			sum += pop_[j].getPathLength();
		}
		if (config == "american_capitals")
		{
			out1 << i << " " << R * pop_[0].getPathLength() << endl;
		}
		else
		{
			out1 << i << " " << pop_[0].getPathLength() << endl;
			out2 << i << " " << sum / (double)half_n_pop << endl;
		}
	}
}

void GenAlgTSP ::PrintPathsFile(std::ofstream &out) const
{
	for (int i = 0; i < pop_size_; i++)
	{
		for (int j = 0; j < n_cities_; j++)
		{
			out << pop_[i].getIndexes()[j] << " " << pop_[i][j].getX() << " " << pop_[i][j].getY() << " ";
		}
		out << pop_[i].getIndexes()[n_cities_] << " " << pop_[i][n_cities_].getX() << " " << pop_[i][n_cities_].getY() << endl;
	}
}

GenAlgTSP &GenAlgTSP ::operator=(const GenAlgTSP &tsp)
{
	pop_ = tsp.get_pop_vector();
	n_cities_ = tsp.get_n_cities();
	pop_size_ = tsp.get_pop_size();
	p_mut_ = tsp.get_p_mut();
	p_cros_ = tsp.get_p_cros();
	norm_power_ = tsp.getNormPower();
	rand_ = tsp.getRandomGen();

	return *this;
}
