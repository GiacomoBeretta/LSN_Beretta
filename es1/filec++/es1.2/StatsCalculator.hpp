#ifndef __Stats_Calculator__
#define __Stats_Calculator__

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm> //transform
#include <functional> //plus


class StatsCalculator{

	public:

		StatsCalculator();
		StatsCalculator(int N, int M);
		
		~StatsCalculator();

		void SetNM(int N, int M);

		int GetN(){return N_;}
		int GetM(){return M_;}
		int GetL(){return L_;}
		int GetMean(){return mean_;}
		int GetStDevMean(){return StDevMean_;}

		virtual double GetValue() = 0;
		virtual std::vector<double> GetVectorValues() = 0;
		
		void Print_Mean_StDevMean(const std::string Filename);
		void Print_Mean_StDevMean(const std::vector<std::string> Filenames);
		virtual void Print_Parameters(const std::string Filename);
		
	protected:
		int N_; //number of throws
		int M_; //number of blocks
		int L_; //dimension of a block

		double mean_;
		double StDevMean_;
		
};


#endif // __Stats_Calculator__
