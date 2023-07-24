#include "ParallelTempering.hpp"

using namespace std;

int main(int argc, char *argv[]) {
  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double tstart = MPI_Wtime();

  int n_cities, norm_power, restart, n_temp, n_temp_per_node, nstep;
  double temp_min, temp_max, temp_min_node, temp_max_node;
  string config, input_namefile, best_namefile, best_xy_namefile;
  ofstream output_xy;

  input_namefile = "input.txt"; // computer laboratorio
  // input_namefile = "/mnt/c/Users/giaco/esercizi_python/10/input.txt"; // mio
  // computer

  ifstream input(input_namefile);

  input >> n_cities;
  input >> norm_power;
  input >> restart;
  input >> temp_min;
  input >> temp_max;
  input >> n_temp;
  input >> config;
  input >> nstep;

  input.close();

  if (n_temp % size != 0) {
    cout << "Please insert a number of temperauters divisible by the number of "
            "processes"
         << endl;
    exit(-1);
  }
  n_temp_per_node = n_temp / size;

  if (config == "american_capitals") {
    n_cities = 50;
  }

  if (rank == 0) {
    cout << endl;
    cout << "Parallel tempering to solve the Travel Salesman Problem:" << endl;
    cout << "number of processes = " << size << endl;
    cout << "Number of cities = " << n_cities << endl;
    cout << "power of the distance = " << norm_power << endl;
    cout << "min temperature = " << temp_min << endl;
    cout << "max temperature = " << temp_max << endl;
    cout << "number of temperatures = " << n_temp << endl;
    cout << "number of temperatures per node = " << n_temp_per_node << endl;
    cout << "configurations of cities = " << config << endl;
    cout << "number of step montecarlo = " << nstep << endl << endl;

    best_namefile =
        "/home/studenti/giacomo.beretta/LabSim/es10/python/best_path_" +
        config + "_" + to_string(norm_power) + "norm.txt";
    best_xy_namefile =
        "/home/studenti/giacomo.beretta/LabSim/es10/python/best_path_" +
        config + "_xy_" + to_string(norm_power) + "norm.txt";

    // sul mio computer
    // best_namefile =
    // "/home/studenti/giacomo.beretta/LabSim/es10/python/best_path_" + config +
    // "_" + to_string(norm_power) + "norm.txt"; best_xy_namefile =
    // "/mnt/c/Users/giaco/esercizi_python/10/best_path_" + config + "_xy_" +
    // to_string(norm_power) + "norm.txt";
  }

  double delta_temp = (temp_max - temp_min) / (double)(n_temp - 1);
  vector<double> temp(n_temp);
  for (int i = 0; i < n_temp; i++) {
    temp[i] = temp_min + delta_temp * i;
  }

  temp_min_node = temp[rank * n_temp_per_node];
  temp_max_node = temp[(rank + 1) * n_temp_per_node - 1];

  Random rand;
  rand.SetRandomPrimes(rank + 1);

  ParallelTempering pop(temp_min_node, temp_max_node, n_temp_per_node, n_cities,
                        norm_power, rand);
  if (config == "circumference") {
    pop.setCitiesRandomlyOnCircumference();
  } else if (config == "square") {
    pop.setCitiesRandomlyInSquare();
  } else {
    pop.setAmericanCapitals();
  }
  pop.checkPaths();

  Path best_path = pop.printBestValues(best_namefile, nstep);

  if (rank == 0) {
    output_xy.open(best_xy_namefile);
    best_path.printPathFile(output_xy);
    output_xy.close();
  }
  double tend = MPI_Wtime();
  double time = tend - tstart;

  if (rank == 0)
    cout << "tempo impiegato = " << time << endl;
  MPI_Finalize();
  return 0;
}
