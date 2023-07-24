#ifndef __ParallelTempering__
#define __ParallelTempering__

#include "mpi.h"
#include "path.hpp"

class ParallelTempering {
public:
  ParallelTempering();
  ParallelTempering(const double t_min, const double t_max, const int n_temp,
                    const int n_cities, const int norm_power, Random &rand);
  ParallelTempering(const Path &path, const double t_min, const double t_max,
                    const int n_temp, Random &rand);
  std::vector<Path> getPaths() const { return paths_; }
  Path getPath(const int i) const { return paths_[i]; }
  double getTmax() const { return t_max_; }
  double getTmin() const { return t_min_; }
  int getNtemp() const { return n_temp_; }
  int getNcities() const { return n_cities_; }
  int getNormPower() const { return norm_power_; }
  Random &getRandom() { return rand_; }

  void setTmax(const double tmax);
  void setTmx(const double tmin);
  void setNtemp(const double n_temp);
  void randomOrder();
  void setCitiesRandomlyOnCircumference();
  void setCitiesRandomlyInSquare();
  void setAmericanCapitals();
  void setRandom(Random &rand);

  void checkPaths() const;
  void computeLengths();

  void metropolis_temp(); // extended ensemble metropolis between paths
                          // with different temperatures which are in the
                          // same process
  // double exchange_prob_mpi(const int node1, const int node2);
  void exchange_paths_mpi(const int node1, const int node2);
  void metropolis_temp_mpi(); // extended ensemble metropolis between paths with
                              // different temperatures which are in adjacent
                              // processes

  // void sortPaths();  // sort the paths in ascending order of path_length
  void putBestOnTop();
  void metropolis(); // metropolis for each path in each node (with
                     // mutations)
  Path bestPath() const;
  Path printBestValues(const std::string filename, const int nstep);

private:
  std::vector<Path> paths_;
  double t_max_;
  double t_min_;
  int n_temp_;
  int n_cities_;
  int norm_power_;
  Random rand_;
  bool sphere_;
};

double exchange_prob(const double L1, const double L2, const double temp1,
                     const double temp2);

#endif // __ParallelTempering__
