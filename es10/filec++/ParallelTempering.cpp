#include "ParallelTempering.hpp"

using namespace std;

ParallelTempering::ParallelTempering() {
  vector<Path> v;
  paths_ = v;
  n_temp_ = 0;
  n_cities_ = 0;
  norm_power_ = 1;
  Random rand;
  rand_ = rand;
  sphere_ = false;
}

ParallelTempering::ParallelTempering(const double t_min, const double t_max,
                                     const int n_temp, const int n_cities,
                                     const int norm_power, Random &rand) {
  t_min_ = t_min;
  t_max_ = t_max;
  n_temp_ = n_temp;
  n_cities_ = n_cities;
  norm_power_ = norm_power;

  vector<Path> v(n_temp_);
  Path v2(n_cities_, norm_power_);
  double delta_t;

  if (n_temp_ != 1) {
    delta_t = (t_max_ - t_min_) / (double)(n_temp_ - 1);
  } else {
    delta_t = 0;
  }

  double temp;
  for (int i = 0; i < n_temp_; i++) {
    temp = t_min_ + delta_t * i;
    v2.setTemp(temp);
    v[i] = v2;
  }

  paths_ = v;
  rand_ = rand;
  sphere_ = false;
}

ParallelTempering::ParallelTempering(const Path &path, const double t_min,
                                     const double t_max, const int n_temp,
                                     Random &rand) {
  t_min_ = t_min;
  t_max_ = t_max;
  n_temp_ = n_temp;
  n_cities_ = path.getNcities();
  norm_power_ = path.getNormPower();
  sphere_ = path.getSphere();

  double delta_t;
  if (n_temp_ != 1) {
    delta_t = (t_max_ - t_min_) / (double)(n_temp_ - 1);
  } else {
    delta_t = 0;
  }

  vector<Path> paths(n_temp_);
  Path p = path;
  double temp;
  for (int i = 0; i < n_temp_; i++) {
    temp = t_min_ + delta_t * i;
    p.setTemp(temp);
    paths[i] = p;
  }
  paths_ = paths;
  rand_ = rand;
  randomOrder();
  computeLengths();
}

void ParallelTempering::setTmax(const double tmax) { t_max_ = tmax; }
void ParallelTempering::setTmx(const double tmin) { t_min_ = tmin; }
void ParallelTempering::setNtemp(const double n_temp) { n_temp_ = n_temp; }
void ParallelTempering::setRandom(Random &rand) { rand_ = rand; }

void ParallelTempering::randomOrder() {
  for (int i = 0; i < n_temp_; i++) {
    paths_[i].randomOrder(rand_);
  }
}

void ParallelTempering::setCitiesRandomlyOnCircumference() {
  paths_[0].setCitiesRandomlyOnCircumference(rand_);
  vector<Position> p0 = paths_[0].getPath();
  vector<int> index0 = paths_[0].getIndexes();

  for (int i = 1; i < n_temp_; i++) {
    paths_[i].setPath(p0);
    paths_[i].setIndexes(index0);
  }
  randomOrder();
}

void ParallelTempering::setCitiesRandomlyInSquare() {
  paths_[0].setCitiesRandomlyInSquare(rand_);

  vector<Position> p0 = paths_[0].getPath();
  vector<int> index0 = paths_[0].getIndexes();

  for (int i = 1; i < n_temp_; i++) {
    paths_[i].setPath(p0);
    paths_[i].setIndexes(index0);
  }
  randomOrder();
}

void ParallelTempering::setAmericanCapitals() {
  sphere_ = true;
  paths_[0].setAmericanCapitals();

  vector<Position> p0 = paths_[0].getPath();
  vector<int> index0 = paths_[0].getIndexes();
  for (int i = 1; i < n_temp_; i++) {
    paths_[i].setSphere(sphere_);
    paths_[i].setPath(p0);
    paths_[i].setIndexes(index0);
  }
  randomOrder();
}

void ParallelTempering::checkPaths() const {
  for (int i = 0; i < n_temp_; i++) {
    paths_[i].checkPath();
  }
}

void ParallelTempering::computeLengths() {
  for (int i = 0; i < n_temp_; i++) {
    paths_[i].computePathLength();
  }
}

void ParallelTempering::metropolis_temp() {
  double L1, L2, temp1, temp2, prob;
  vector<Position> path_tmp;
  vector<int> index_tmp;

  for (int i = 0; i < n_temp_ - 1; i++) {
    L1 = paths_[i].getPathLength();
    L2 = paths_[i + 1].getPathLength();
    temp1 = paths_[i].getTemp();
    temp2 = paths_[i + 1].getTemp();

    prob = exchange_prob(L1, L2, temp1, temp2);
    if (rand_.Rannyu() < prob) {
      path_tmp = paths_[i].getPath();
      index_tmp = paths_[i].getIndexes();

      paths_[i].setPath(paths_[i + 1].getPath());
      paths_[i].setIndexes(paths_[i + 1].getIndexes());
      paths_[i].setPathLength(L2);
      paths_[i + 1].setPath(path_tmp);
      paths_[i + 1].setIndexes(index_tmp);
      paths_[i + 1].setPathLength(L1);
    }
  }
}

void ParallelTempering::exchange_paths_mpi(const int node1, const int node2) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double x1, y1, L1, x2, y2, L2;
  int index1, index2;
  MPI_Status stat;

  if (rank == node1) {
    for (int i = 1; i < n_cities_; i++) {
      x1 = paths_[n_temp_ - 1].getXcity(i);
      y1 = paths_[n_temp_ - 1].getYcity(i);
      index1 = paths_[n_temp_ - 1].getIndex(i);

      MPI_Send(&x1, 1, MPI_REAL8, node2, 0, MPI_COMM_WORLD);
      MPI_Recv(&x2, 1, MPI_REAL8, node2, 0, MPI_COMM_WORLD, &stat);

      MPI_Send(&y1, 1, MPI_REAL8, node2, 0, MPI_COMM_WORLD);
      MPI_Recv(&y2, 1, MPI_REAL8, node2, 0, MPI_COMM_WORLD, &stat);

      MPI_Send(&index1, 1, MPI_INTEGER2, node2, 0, MPI_COMM_WORLD);
      MPI_Recv(&index2, 1, MPI_INTEGER2, node2, 0, MPI_COMM_WORLD, &stat);

      paths_[n_temp_ - 1][i].setX(x2);
      paths_[n_temp_ - 1][i].setY(y2);
      paths_[n_temp_ - 1].setIndex(i, index2);
    }

    L1 = paths_[n_temp_ - 1].getPathLength();
    MPI_Send(&L1, 1, MPI_REAL8, node2, 0, MPI_COMM_WORLD);
    MPI_Recv(&L2, 1, MPI_REAL8, node2, 0, MPI_COMM_WORLD, &stat);
    paths_[n_temp_ - 1].setPathLength(L2);
  }
  if (rank == node2) {
    for (int i = 1; i < n_cities_; i++) {
      x2 = paths_[0].getXcity(i);
      y2 = paths_[0].getYcity(i);
      index2 = paths_[0].getIndex(i);

      MPI_Recv(&x1, 1, MPI_REAL8, node1, 0, MPI_COMM_WORLD, &stat);
      MPI_Send(&x2, 1, MPI_REAL8, node1, 0, MPI_COMM_WORLD);

      MPI_Recv(&y1, 1, MPI_REAL8, node1, 0, MPI_COMM_WORLD, &stat);
      MPI_Send(&y2, 1, MPI_REAL8, node1, 0, MPI_COMM_WORLD);

      MPI_Recv(&index1, 1, MPI_INTEGER2, node1, 0, MPI_COMM_WORLD, &stat);
      MPI_Send(&index2, 1, MPI_INTEGER2, node1, 0, MPI_COMM_WORLD);

      paths_[0][i].setX(x1);
      paths_[0][i].setY(y1);
      paths_[0].setIndex(i, index1);
    }

    L2 = paths_[0].getPathLength();
    MPI_Recv(&L1, 1, MPI_REAL8, node1, 0, MPI_COMM_WORLD, &stat);
    MPI_Send(&L2, 1, MPI_REAL8, node1, 0, MPI_COMM_WORLD);
    paths_[0].setPathLength(L1);
  }
}

void ParallelTempering::metropolis_temp_mpi() {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat;

  double accept_prob = 0;
  double L1, L2, temp1, temp2;
  Random rand;

  bool proceeds = false;
  bool exchange = false;
  if (rank == 0) {
    proceeds = true;
  }

  if (proceeds == false) {
    L2 = paths_[0].getPathLength();
    temp2 = paths_[0].getTemp();
    MPI_Send(&L2, 1, MPI_REAL8, rank - 1, 0, MPI_COMM_WORLD);
    MPI_Send(&temp2, 1, MPI_REAL8, rank - 1, 0, MPI_COMM_WORLD);
    MPI_Recv(&exchange, 1, MPI_LOGICAL, rank - 1, 0, MPI_COMM_WORLD, &stat);
    if (exchange) {
      exchange_paths_mpi(rank - 1, rank);
    }
    if (rank != size - 1) {
      MPI_Recv(&proceeds, 1, MPI_LOGICAL, rank - 1, 0, MPI_COMM_WORLD, &stat);
    }
  }

  if (proceeds == true && rank != size - 1) {
    L1 = paths_[n_temp_ - 1].getPathLength();
    temp1 = paths_[n_temp_ - 1].getTemp();
    MPI_Recv(&L2, 1, MPI_REAL8, rank + 1, 0, MPI_COMM_WORLD, &stat);
    MPI_Recv(&temp2, 1, MPI_REAL8, rank + 1, 0, MPI_COMM_WORLD, &stat);

    accept_prob = exchange_prob(L1, L2, temp1, temp2);
    if (rand_.Rannyu() < accept_prob) {
      exchange = true;
      MPI_Send(&exchange, 1, MPI_LOGICAL, rank + 1, 0, MPI_COMM_WORLD);
      exchange_paths_mpi(rank, rank + 1);
    } else {
      exchange = false;
      MPI_Send(&exchange, 1, MPI_LOGICAL, rank + 1, 0, MPI_COMM_WORLD);
    }

    if (rank != size - 2) {
      MPI_Send(&proceeds, 1, MPI_LOGICAL, rank + 1, 0, MPI_COMM_WORLD);
    }
  }
}

void ParallelTempering::putBestOnTop() {
  vector<Position> best_path = paths_[0].getPath();
  vector<int> best_indexes = paths_[0].getIndexes();
  double best_L = paths_[0].getPathLength();

  vector<Position> path_tmp;
  vector<int> index_tmp;
  double L_tmp;

  for (int i = 1; i < n_temp_; i++) {
    L_tmp = paths_[i].getPathLength();
    if (L_tmp < best_L) {
      path_tmp = paths_[i].getPath();
      index_tmp = paths_[i].getIndexes();

      paths_[i].setPath(best_path);
      paths_[i].setIndexes(best_indexes);
      paths_[i].setPathLength(best_L);
      best_path = path_tmp;
      best_indexes = index_tmp;
      best_L = L_tmp;
    }
  }
  paths_[0].setPath(best_path);
  paths_[0].setIndexes(best_indexes);
  paths_[0].setPathLength(best_L);
}

void ParallelTempering::metropolis() {
  for (int i = 0; i < n_temp_; i++) {
    paths_[i].metropolis(rand_);
  }
}

Path ParallelTempering::bestPath() const {
  Path best_path = paths_[0];
  double best_L, L_tmp;

  for (int i = 1; i < n_temp_; i++) {
    L_tmp = paths_[i].getPathLength();
    best_L = best_path.getPathLength();

    if (L_tmp < best_L) {
      best_path = paths_[i];
    }
  }
  return best_path;
}

Path ParallelTempering::printBestValues(const string filename,
                                        const int nstep) {
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  ofstream output_best(filename);

  Path dummy_path, best_path;
  double dummy_L, best_L;

  best_path = bestPath();
  best_L = best_path.getPathLength();

  double R = 6373000;
  for (int i = 0; i < nstep; i++) {
    if (i % 10000 == 0 && rank == 0) {
      cout << "rank = " << rank << ", i = " << i << endl;
    }

    metropolis();

    if (rank == 0) {
      dummy_path = bestPath();
      dummy_L = dummy_path.getPathLength();
      if (dummy_L < best_L) {
        best_path = dummy_path;
        best_L = dummy_L;
      }
    }

    // putBestOnTop();

    metropolis_temp();

    if (rank == 0) {
      dummy_path = bestPath();
      dummy_L = dummy_path.getPathLength();
      if (dummy_L < best_L) {
        best_path = dummy_path;
        best_L = dummy_L;
      }
    }

    // putBestOnTop();

    metropolis_temp_mpi();

    if (rank == 0) {
      dummy_path = bestPath();
      dummy_L = dummy_path.getPathLength();
      if (dummy_L < best_L) {
        best_path = dummy_path;
        best_L = dummy_L;
      }
    }

    // putBestOnTop();

    if (rank == 0) {
      if (sphere_ == false) {
        output_best << i << " " << best_L << endl;
      } else {
        output_best << i << " " << R * best_L << endl;
      }
    }
  }
  output_best.close();
  return best_path;
}

double exchange_prob(const double L1, const double L2, const double temp1,
                     const double temp2) {
  double beta_1 = 1.0 / temp1;
  double beta_2 = 1.0 / temp2;

  return exp(-(beta_2 - beta_1) * (L1 - L2));
}