/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __fluid__
#define __fluid__

//Random numbers
#include "random.hpp"
#include <ostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm> //transform
#include <sstream> //string stream
#include <stdio.h> //remove

int seed[4];
Random rnd;

//parameters, observables
const int m_props=1000; //per molecular dynamics m_props=n_props
int n_props, iv, ip, it, ig; //n_props è il numero di osservabili, iw è il viriale
double vtail, ptail, bin_size, nbins, deltaVol[m_props], sd;
double walker[m_props];
std::vector<double> epot_t;
std::vector<double> pres_t;

// average
double blk_av[m_props], blk_norm, accepted, attempted;
double glob_av[m_props], glob_av2[m_props];
double stima_epot, stima_pres, stima_temp, stima_g_r[m_props];
double err_epot, err_pres, err_temp, err_g_r[m_props], err_gdir;

//configuration
const int m_part=108;
double x[m_part],    y[m_part],    z[m_part];
double xold[m_part], yold[m_part], zold[m_part];
double vx[m_part],  vy[m_part],   vz[m_part];

// thermodynamical state
int npart;
double beta, temp_0 ,energy, vol, rho, box, rcut;

// simulation
int iNVET, restart, nblk, nstep, nstep_equilibrium, autocor_on, iblk; //0=MD(NVE) 1=MC(NVT), 0=restart, 1=start from seed.out, nstep=number of step in one block
double delta; //time step for integration
std::string state; //tells the program from which file extract the initial configuration(solid/liquid/gas)
std::string algo_name;
//output files
const int wd=15;
std::ofstream Epot, Pres, Temp, Isto_g_r_final_ave; //Isto_g_r,

//pigreco
const double pi=3.1415927;

//functions
void Input(int, char**);
void Reset(void);
void Accumulate(void);
void Averages(void);
void PrintAverages(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void);
double Boltzmann(double, double, double, int);
double Pbc(double);
double Error(double,double,int);
double Force(int, int);
void Equilibrate(void);
void PrintEnergyAutoCorrelation(void);
void OpenAveragesFiles(void);
void PrintAverages(int);
void CloseAveragesFiles(void);
void Print_g_r_final_ave(void);
void OpenInstantValuesFiles(void);
void PrintInstantValues(int);
void CloseInstantValuesFiles(void);
void PrintEpotAutoCorrelation(void);
void PrintPresAutoCorrelation(void);

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
