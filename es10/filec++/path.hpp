#ifndef __Path__
#define __Path__

#include "position.hpp"
#include "random.hpp"
#include <sstream>
#include <vector>

class Path {
public:
  Path();
  Path(const int n_cities);
  Path(const int n_cities, const int norm_power);
  Path(const int n_cities, const int norm_power, const double temp);
  Path(const std::vector<std::vector<double>> &v, const int norm_power,
       const double temp);
  Path(const std::vector<Position> &path, const int norm_power,
       const double temp);
  Path(const Path &path);

  std::vector<Position> getPath() const { return path_; }
  unsigned int getNcities() const { return n_cities_; }
  int getNormPower() const { return norm_power_; }
  bool getSphere() const { return sphere_; }
  double getPathLength() const { return path_length_; }
  int getIndex(const int i) { return indexes_[i]; }
  std::vector<int> getIndexes() const { return indexes_; }
  double getXcity(const int i) const { return path_[i].getX(); }
  double getYcity(const int i) const { return path_[i].getY(); }
  double getTemp() const { return temp_; }
  int getAccepted() const { return accepted_; }
  int getAttempted() const { return attempted_; }
  double getAccRate() const { return accepted_ / (double)attempted_; }

  void setXcity(const int i, const double x);
  void setYcity(const int i, const double y);
  void setCity(const int i, const Position &p);
  void setPath(const std::vector<Position> &path);
  void setNormPower(const int norm_power);
  void setSphere(const bool sphere);
  void setPathLength(const double L);
  void setIndex(const int i, const int index);
  void setIndexes(const std::vector<int> &indexes);
  void setTemp(const double temp);
  void setAccepted(const int accepted);
  void setAttempted(const int attempted);

  void randomOrder(Random &rand);
  void setCitiesRandomlyOnCircumference(Random &rand);
  void setCitiesRandomlyInSquare(Random &rand);
  void setAmericanCapitals();

  void swapCities(const int a, const int b);
  // void shiftCities(int a, int n, int shift);
  void permutateCities(const int a, const unsigned int n, const int b);
  void reverseCities(const int a, const unsigned int n);
  void checkPath() const;
  void computePathLength();

  double probability();
  void metropolis(Random &rand); // it does three mutations and for each of
                                 // them computes the acceptance prob

  Path &operator=(const Path &path);
  Position &operator[](const int i);
  const Position &operator[](const int i) const;

  void printPath() const;
  void printPathFile(std::ofstream &out) const;
  void printIndexes() const;

private:
  std::vector<Position> path_;
  unsigned int n_cities_;
  int norm_power_;
  bool sphere_;
  double path_length_;
  std::vector<int> indexes_;

  double temp_;

  int accepted_;
  int attempted_;
};

int boundaryCondition(const int a, const int n_cities);
int boundaryDistance(const int a, const int b, const int n_cities);

#endif // __Path__
