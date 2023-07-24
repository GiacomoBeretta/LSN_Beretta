#ifndef __Metropolis__
#define __Metropolis__

#include "random.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm> //transform
#include <functional> //plus

class Metropolis : public Random
{
	public:
		Metropolis();
		Metropolis(double delta, int equilibration_time);

		double getDelta() const{return delta_;}
		int getEquilibrationTime() const{return equilibration_time_;}
		std::vector<double> getX() const{return x_old_;}
		int getAccepted() const{return accepted_;}
		int getAttempted() const{return attempted_;}
		double getAccRate() const{return accepted_/(double)attempted_;}
	
		void setDelta(double delta);
		void setEquilibrationTime(int equilibration_time);
		virtual void setX(std::vector<double> x); // perchè virtual?
		void setAccepted(int accepted);
		void setAttempted(int attempted);
		
		virtual double probability(std::vector<double> x) = 0;//=0 makes this an abstract class which cannot be instantiated
		virtual void moveUniform();
		void equilibrate();

		double computeAccRate(int nstep);
		void calibrateDelta(int nstep);
	private:
		std::vector<double> x_old_; // vedere se è più veloce se metto un x_new_
		
		double delta_; //monte carlo step   
		int equilibration_time_;
		
		int accepted_; //number of accepted moves
		int attempted_; //number of attempted moves
};

#endif // __Metropolis__
