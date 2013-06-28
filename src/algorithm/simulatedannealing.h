/**
Author: Jialu Hu
Date: Jun. 27, 2013
File name: algorithm/simulatedannealing.h
Description: search an optimal solution for internal nodes.
**/

#ifndef SIMULATEDANNEALING_H_
#define SIMULATEDANNEALING_H_

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <deque>
#include <stack>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include "verbose.h"
#include "string.h"

template<typename PH, typename OP>
class SimulatedAnnealing
{
private:
	typedef PH Phylogeny;
	typedef OP Option;
public:
	float _tmax;
    float _tmin;
    unsigned _Kmax;
    unsigned _Nmax;
    float _k;
	SimulatedAnnealing();
	~SimulatedAnnealing(){};
	void run(Phylogeny&);
};

template<typename PH, typename OP>
SimulatedAnnealing<PH,OP>::SimulatedAnnealing()
{
}

template<typename PH, typename OP>
void
SimulatedAnnealing<PH,OP>::run(Phylogeny& phylogeny)
{
	unsigned k=0;
	float t=_tmax;
	float step=(_tmax-_tmin)/_Kmax;
	unsigned seed =std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::uniform_real_distribution<float> distribution(0.0,1.0);
	while(k<=_Kmax)
    {
		t = t-step;
		float beta = 1.0/(_k*t);
		for(unsigned n=0;n<_Nmax;++n)
		{
			float delta = 0.0;
			if(!phylogeny.interfere())continue;
			if(delta>0 || distribution(generator)<exp(beta*delta))
			{
				// update current state to the neighbor state.
			}
		}
	}
}

#endif
