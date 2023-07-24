#ifndef __Path__
#define __Path__

#include "random.hpp"
#include "position.hpp"
#include <sstream>
#include <vector>

class Path
{
public:
	Path();
	Path(int n_cities);
	Path(const std::vector<std::vector<double>> &v, int norm_power);
	Path(const std::vector<Position> &path, int norm_power);
	Path(const Path &path);

	std::vector<Position> getPath() const { return path_; }
	unsigned int get_n_cities() const { return path_.size() - 1; }
	double getPathLength() const { return path_length_; }
	std::vector<int> getIndexes() const { return indexes_; }
	double getXcity(int i) const { return path_[i].getX(); }
	double getYcity(int i) const { return path_[i].getY(); }

	void setCity(int i, const Position &p);
	void setPath(const std::vector<Position> &path);
	void setIndexes(const std::vector<int> &indexes);

	Path randomOrder(Random &rand) const;
	void setCitiesRandomlyOnCircumference(Random &rand);
	void setCitiesRandomlyInSquare(Random &rand);
	void setAmericanCapitals(const std::string filename);

	void swapCities(const int a, const int b);
	// void shiftCities(int a, int n, int shift);
	void permutateCities(const int a, const unsigned int n, const int b);
	void reverseCities(const int a, const unsigned int n);
	void checkPath() const;
	void pathLength(int norm_power);
	double pathLengthSphere(int norm_power);

	Path &operator=(const Path &path);
	Position &operator[](const int i);
	const Position &operator[](const int i) const;

	void printPath() const;
	void printPathFile(std::ofstream &out) const;
	void printIndexes() const;

private:
	std::vector<Position> path_;
	double path_length_;
	std::vector<int> indexes_;
};

int boundaryCondition(const int a, const int n_cities);
int boundaryDistance(const int a, const int b, const int n_cities);

#endif // __Path__
