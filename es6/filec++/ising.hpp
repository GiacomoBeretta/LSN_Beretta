/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __ISING__
#define __ISING__

//Random numbers
#include "random.hpp"
//#include <ostream>
#include <iomanip> //setw
#include <vector>

int seed[4];
Random rnd;

//parameters, observables
const int m_props = 1000;
int n_props, iu, ic, im, ix, ig;
double nbins;
double walker[m_props];
std::vector<double> u_t; //vector for computing autocorrelation function 
std::vector<double> m_t;
std::vector<double> x_t;

// averages
double blk_av[m_props], blk_norm;
double glob_av[m_props], glob_av2[m_props];
double stima_u, stima_c, stima_m, stima_x, stima_g;
double err_u, err_c, err_m, err_x, err_g;
int accepted ,attempted;

//configuration
const int m_spin=100;
double s[m_spin];

// thermodynamical state
int nspin;
double beta,temp,J,h;

// simulation
int nstep, metro, nblk, n_temp, restart, nstep_equilibrium;
double temp_min, temp_max;
bool autocor_on;

//output files
const int wd=15;
std::string algo_name;
std::ofstream Ene, Heat, Mag, Chi;
std::ofstream Ene_T, Heat_T, Mag_T, Chi_T;

//functions
void Input(void);
void Reset(int);
void Accumulate(void);
void Averages(int);
void Move(void);
void ConfFinal(void);
void ReadConf(void);
void Measure(void);
double Boltzmann(int, int);
int Pbc(int);
double Error(double,double,int);
void ModifyTemp(int);
void Equilibrate(void);
void PrintEnergyAutoCorrelation(void);
void PrintMagnetizationAutoCorrelation(void);
void PrintMagnetSusceptibilityAutoCorrelation(void);
void OpenBlockFiles(void);
void OpenFinalFiles(void);
void CloseBlockFiles(void);
void CloseFinalFiles(void);
void PrintBlockValues(int);
void PrintFinalValues(int);
#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
