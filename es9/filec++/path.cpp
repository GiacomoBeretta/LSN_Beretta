#include "path.hpp"

using namespace std;

const double pi = 3.14159265358979323846;

Path ::Path()
{
	vector<Position> path;
	vector<int> indexes;
	path_ = path;
	indexes_ = indexes;
	path_length_ = 0;
}

Path ::Path(int n_cities)
{
	vector<Position> path(n_cities + 1);
	vector<int> indexes(n_cities + 1);
	for (int i = 0; i < n_cities + 1; i++)
	{
		indexes[i] = i;
		path[i].setX(0);
		path[i].setY(0);
	}
	indexes[n_cities] = 0;

	path_ = path;
	path_length_ = 0;
	indexes_ = indexes;
}

Path ::Path(const vector<vector<double>> &v, int norm_power)
{
	unsigned int n_cities = v.size();

	vector<Position> path(n_cities + 1);
	vector<int> indexes(n_cities + 1);

	for (unsigned int i = 0; i < n_cities; i++)
	{
		Position pos(v[i][0], v[i][1]);
		path[i] = pos;
		indexes[i] = i;
	}
	path[n_cities] = path[0];
	indexes[n_cities] = 0;

	path_ = path;
	pathLength(norm_power);
	indexes_ = indexes;
}

Path ::Path(const std::vector<Position> &path, int norm_power)
{
	path_ = path;
	pathLength(norm_power);
	unsigned int n_cities = get_n_cities();
	vector<int> indexes(n_cities + 1);
	for (unsigned int i = 0; i < n_cities; i++)
	{
		indexes[i] = i;
	}
	indexes[n_cities] = 0;
	indexes_ = indexes;
}

Path ::Path(const Path &path)
{
	path_ = path.getPath();
	path_length_ = path.getPathLength();
	indexes_ = path.getIndexes();
}

void Path ::setCity(int i, const Position &p)
{
	path_[i] = p;
}

void Path ::setPath(const std::vector<Position> &path)
{
	path_ = path;
	pathLength(1);
}

void Path ::setIndexes(const vector<int> &indexes)
{
	indexes_ = indexes;
}

Path Path ::randomOrder(Random &rand) const
{
	unsigned int n_cities = get_n_cities();
	vector<int> indexes(n_cities + 1);

	for (unsigned int i = 0; i < n_cities + 1; i++)
	{
		indexes[i] = 0;
	}

	unsigned int j = 1;
	int r;
	while (j < n_cities) // n_cities=34, j=0,...,33
	{
		r = rand.Integer(1, n_cities - 1); // r=1,...,33
		if (indexes[r] == 0)
		{
			indexes[r] = j;
			j++;
		}
	}

	Path path_copy = *this;
	for (unsigned int i = 1; i < n_cities; i++)
	{
		path_copy.setCity(i, path_[indexes[i]]);
	}

	path_copy.setIndexes(indexes);

	return path_copy;
}

void Path ::setCitiesRandomlyOnCircumference(Random &rand)
{
	double r;
	for (unsigned int i = 0; i < get_n_cities(); i++)
	{
		r = rand.Rannyu(0, 2 * pi);
		path_[i].setX(cos(r)); // x on the unit circle
		path_[i].setY(sin(r)); // y on the unit circle
	}
	path_[get_n_cities()] = path_[0];
}

void Path ::setCitiesRandomlyInSquare(Random &rand)
{
	for (unsigned int i = 0; i < get_n_cities(); i++)
	{
		path_[i].setX(rand.Rannyu()); // x,y in the square (0,1)  (1,1)
		path_[i].setY(rand.Rannyu()); //                  (0,0)  (1,0)
	}
	path_[get_n_cities()] = path_[0];
}

void Path::setAmericanCapitals(const string filename)
{
	unsigned int n_cities = get_n_cities();
	ifstream input(filename);
	string line;
	getline(input, line); // first line doesn't contain data
	string word[3];
	double phi, theta;
	for (unsigned int i = 0; i < n_cities; i++)
	{
		getline(input, line);
		istringstream iss(line);
		while (iss >> word[2])
		{
			word[0] = word[1];
			word[1] = word[2];
		}
		phi = stod(word[0]);   // -180°<phi<180° longitude
		theta = stod(word[1]); // latitude

		phi += 180.0; // 0°<phi<360°
		theta =
			90.0 -
			theta; // so that theta is 0° at the north pole and 90° at the equator

		phi = phi * pi / 180.0; // phi and theta in radians
		theta = theta * pi / 180.0;
		path_[i].setX(phi);
		path_[i].setY(theta);
	}
	path_[n_cities] = path_[0];
}

void Path ::swapCities(int a, int b)
{
	Position temp_pos = path_[a];
	int temp_index = indexes_[a];

	path_[a] = path_[b];
	indexes_[a] = indexes_[b];

	path_[b] = temp_pos;
	indexes_[b] = temp_index;
}
/*
void Path :: shiftCities(int a, int n, int shift)
{
	int n_cities = get_n_cities();
	int c, d;

	for(int i=0; i<n; i++)
	{
		if( i+a < n_cities )
		{
			if( i+a+shift < n_cities )
			{
				swapCities(i+a, i+a+shift);
			}else
			{
				c = (i+a+shift) % n_cities;
				swapCities(i+a, i+a+shift-n_cities*c+1);
			}
		}else
		{
			c = (i+a+shift) % n_cities;
			d = (i+a)%n_cities;
			swapCities(i+a-n_cities*d+1, i+a+shift-n_cities*c+1);
		}
	}
}*/

void Path ::permutateCities(const int a, const unsigned int n, const int b)
{
	int c, d;
	int n_cities = get_n_cities();
	for (unsigned int i = 0; i < n; i++)
	{
		c = boundaryCondition(a + i, n_cities);
		d = boundaryCondition(b + i, n_cities);
		swapCities(c, d);
	}
}

void Path ::reverseCities(const int a, const unsigned int n)
{
	int b, c;
	int n_cities = get_n_cities();
	for (unsigned int i = 0; i < n / 2; i++)
	{
		b = boundaryCondition(a + i, n_cities);
		c = boundaryCondition(a + n - 1 - i, n_cities);
		swapCities(b, c);
	}
}

void Path ::checkPath() const
{
	unsigned int n_cities = get_n_cities();

	// check by index
	if (n_cities != indexes_.size() - 1)
	{
		cout << "error: the indexes vector hasn't the right dimension" << endl;
		exit(-1);
	}

	if (indexes_[0] != indexes_[n_cities])
	{
		cout << "error: the first city's index doesn't match the last city's" << endl;
		exit(-1);
	}

	for (unsigned int i = 0; i < n_cities - 1; i++)
	{
		for (unsigned int j = i + 1; j < n_cities; j++)
		{
			if (indexes_[i] == indexes_[j])
			{
				cout << "error: there are two cities with same index, the " << i << " city and the " << j << " city" << endl;
				exit(-1);
			}
		}
	}

	// check by position
	if (path_[0] != path_[n_cities])
	{
		cout << "error: the first city doesn't match the last city" << endl;
		exit(-1);
	}

	for (unsigned int i = 0; i < n_cities - 1; i++)
	{
		for (unsigned int j = i + 1; j < n_cities; j++)
		{
			if (path_[i] == path_[j])
			{
				cout << "error: there are two cities with same position, the " << i << " city and the " << j << " city" << endl;
				exit(-1);
			}
		}
	}
}

void Path ::pathLength(int norm_power)
{
	unsigned int n_cities = get_n_cities();
	double sum = 0;
	double d;
	for (unsigned int i = 0; i < n_cities; i++)
	{
		d = path_[i].distance(path_[i + 1]);
		sum += pow(d, norm_power);
	}
	path_length_ = pow(sum, 1 / (double)norm_power);
}

double Path::pathLengthSphere(int norm_power)
{
	unsigned int n_cities = get_n_cities();
	double d;
	double sum = 0;
	for (unsigned int i = 0; i < n_cities; i++)
	{
		d = path_[i].distanceOnSphere(path_[i + 1]);
		sum += pow(d, norm_power);
	}
	path_length_ = pow(sum, 1 / (double)norm_power);
	return path_length_;
}

Path &Path ::operator=(const Path &path)
{
	path_ = path.getPath();
	path_length_ = path.getPathLength();
	indexes_ = path.getIndexes();
	return *this;
}

Position &Path ::operator[](const int i)
{
	return path_[i];
}

const Position &Path ::operator[](const int i) const
{
	return path_[i];
}

void swap(Path &a, Path &b)
{
	Path temp = a;
	a = b;
	b = temp;
}

void Path ::printPath() const
{
	for (unsigned int i = 0; i < path_.size(); i++)
	{
		cout << "i = " << i << ", x = " << getXcity(i) << ", y = " << getYcity(i) << endl;
	}
}

void Path ::printPathFile(ofstream &out) const
{
	for (unsigned int i = 0; i < path_.size(); i++)
	{
		out << i << " " << getXcity(i) << " " << getYcity(i) << endl;
	}
}

void Path ::printIndexes() const
{
	for (unsigned int i = 0; i < path_.size(); i++)
	{
		cout << "i = " << i << ", index = " << indexes_[i] << endl;
	}
}

int boundaryCondition(const int a, const int n_cities)
{
	int b = a % n_cities;
	if (b == 0)
		return 1;
	else
		return b;
}

int boundaryDistance(const int a, const int b, const int n_cities)
{
	int d1 = abs(a - b);
	int d2 = 0;
	if (a < b)
	{
		d2 = n_cities - b + a - 2; //(n_cities - 1) - b + a - 1;
	}
	else
	{
		d2 = n_cities - a + b - 2;
	}
	return min(d1, d2);
}