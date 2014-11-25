/**
Author: Jialu Hu
Date: Jun. 10, 2013
File name: main.cpp
Description: Major body
**/

#include <iostream>
#include <fstream>
#include <vector>
#include "verbose.h"
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/time_measure.h>
#include "input/networkpool.h"
#include "input/layer_graphs.h"
#include "algorithm/subnet.h"
#include "algorithm/search.h"
#include "input/format.h"
#include "input/output_html.h"

using namespace std;

typedef struct _Option
{
  vector<std::string> speciesfiles;
  vector<std::string> networkfiles;
  std::string layerfile;
  std::string resultfolder;
  std::string profile;
  std::string treefile;
  std::string alignmentfile;
  std::string formatfile;
  double beta;
  double alpha;
  double score_threshold;
  int method;
  int task;
  int numspecies;
  int seedsize;
  int seedtries;
  int minext;
  int maxext;
  int extdist1;
  int extdist2;
  int numseeds;
  int numconnected;
  int numthreads;
  int numspinetries;
  int verbose;
  bool parallel;
  int seedReplication;
  _Option()
  {
    profile="./profile.txt";
	method=1;
	task=0;
    numspecies=3;
    seedsize=2;
    seedtries=2;
	minext=3;
	maxext=13;
	extdist1=2;
	extdist2=2;
    numseeds=2000;
    numconnected=3;
    numthreads=1;
	score_threshold=0.3;
	numspinetries=20;
	verbose=1;
	beta=2.0;
	alpha=0.2;
	parallel=false;
	seedReplication=1;
  }
}Option;

typedef lemon::ListGraph Graph;
typedef lemon::SmartGraph BpGraph;
typedef Tree<Graph, Option> MyTree;
typedef NetworkPool<Graph,BpGraph> InputGraph;
typedef Layer_graphs<BpGraph,InputGraph> LayerGraph;
typedef SubNet<InputGraph,LayerGraph> MySubNet;
typedef Search<InputGraph, MySubNet, LayerGraph, Option> MySearch;
typedef Format<InputGraph, Option> MyFormat;
typedef Output_html<Option> OutputHtml;

int main(int argc, char** argv)
{
   Option myoption;
   OutputHtml webpage;

  	g_verbosity=(VerbosityLevel)myoption.verbose;
  	webpage.set_header();
  	webpage.get_data(myoption);
  	
  	InputGraph networks;
    LayerGraph layergraph;
    MySearch localAlignment(myoption);
    Timer t(false);
	t.start();
	networks.initNetworkPool(myoption.networkfiles);
	layergraph.read(myoption.layerfile,networks);
	localAlignment.run(layergraph,networks);
	t.stop();
	std::cout <<"# Elapsed time: "<< t <<std::endl;
   return 1;
}
