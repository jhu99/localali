/**
Author: Jialu Hu
Date: Jun. 27, 2013
File name: algorithm/simulatedannealing.h
Description: search an optimal solution for internal nodes.
**/

#pragma once
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
	typedef typename Phylogeny::DeltaStructure DeltaStructure;
	typedef typename Phylogeny::Graph Graph;
	typedef typename Phylogeny::GraphData GraphData;
	TEMPLATE_GRAPH_TYPEDEFS(Graph);
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
SimulatedAnnealing<PH,OP>::SimulatedAnnealing():
_tmax(50),
_tmin(10),
_Kmax(20),
_Nmax(50),
_k(0.05)
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
	while(k++<=_Kmax)
	{
		t = t-step;// assert(t>_tmin);
		float beta = -1.0/(_k*t);
		for(unsigned n=0;n<_Nmax;++n)
		{
			DeltaStructure deltaData;
			if(!phylogeny.interfere(deltaData))continue;
			float sampledata=distribution(generator);
			if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
				std::cout <<deltaData.delta <<"\t" << sampledata<<"\t"<< exp(beta*deltaData.delta) <<"\n";
			if(deltaData.delta<0 || sampledata < exp(beta*deltaData.delta))
			{
				// update current state to the neighbor state and its interaction evolutionary score.
				int i=0;
				for(IncEdgeIt it(phylogeny._tree.g,deltaData.treenode);it!=lemon::INVALID;++it,++i)
				{
					phylogeny._tree.scoremap[it]=deltaData.updatedScores[i];
				}
			}
			else
			{
				Node treenode=deltaData.treenode;
				int nodeid=phylogeny._tree.g.id(treenode);
				GraphData* graphdata=phylogeny.node2graph[nodeid];
				if(graphdata->label2edge->find(deltaData.edgelabel)!=graphdata->label2edge->end())
				{
					Edge myedge=graphdata->label2edge->find(deltaData.edgelabel)->second;
					graphdata->deleteEdge(myedge,deltaData.edgelabel);
				}
				else
				{
					graphdata->addEdge(deltaData.nodeA,deltaData.nodeB);
				}
			}
		}
	}
	phylogeny.computeDist();
}

#endif
