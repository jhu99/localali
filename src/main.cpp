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
#include <lemon/arg_parser.h>
#include <lemon/time_measure.h>
#include "input/networkpool.h"
#include "input/processprofile.h"
#include "input/layer_graphs.h"
#include "algorithm/subnet.h"
#include "algorithm/search.h"

using namespace std;

typedef struct _Option
{
  vector<std::string> speciesfiles;
  vector<std::string> networkfiles;
  std::string layerfile;
  std::string resultfolder;
  std::string profile;
  std::string treefile;
  float beta;
  int numspecies;
  int seedsize;
  int seedtries;
  int numextention;
	int minext;
	int maxext;
  int numsamples;
  int numconnected;
  _Option()
  {
    profile="./test_profile.txt";
    numspecies=3;
    seedsize=5;
    seedtries=10;
    numextention=2;
		minext=5;
		maxext=10;
    numsamples=1000;
    numconnected=2;
  }
}Option;

typedef lemon::ListGraph Graph;
typedef lemon::SmartGraph BpGraph;
typedef Tree<Graph, Option> MyTree;
typedef NetworkPool<Graph,BpGraph> InputGraph;
typedef Layer_graphs<BpGraph,InputGraph> LayerGraph;
typedef SubNet<InputGraph,LayerGraph> MySubNet;
typedef Search<InputGraph, MySubNet, LayerGraph, Option> MySearch;

bool setParser(ArgParser& parser, Option& myoption)
{
	parser
	.refOption("profile","Configuration of various input parameters. Default is \"./profile.txt\".", myoption.profile)
	.refOption("numspecies","Number of the species compared. Default is 3.", myoption.numspecies)
	.refOption("seedtries","Number of tries for each refined seeds. Default is 100.", myoption.seedtries)
	.refOption("seedsize","Size of the seeds. Default is 5.", myoption.seedsize)
	.refOption("numextention","Number of the extention . Default is 10.", myoption.numextention)
	.refOption("numconnected","Number of connected subnetwork. Default is 2.", myoption.numconnected)
	.refOption("numsamples","Number of sampled seeds. Default is 1000.", myoption.numsamples);
	return true;
}

bool runParser(ArgParser& myparser, Option& myoption)
{
  std::string filename;
  ProcessProfile<Option> myprofile(myoption.profile);
  myprofile.getOption(myoption);
  myparser.run();
  filename.append(myoption.resultfolder);
  return true;
}

int main(int argc, char** argv)
{
  Option myoption;
  ArgParser myparser(argc,argv);
  setParser(myparser,myoption);
  runParser(myparser,myoption);
  
  InputGraph networks;
  LayerGraph layergraph;
  MySearch isearch(myoption);
  Timer t(false);

  g_verbosity=VERBOSE_ESSENTIAL;

  t.start();
  //Test read interface for PPI networks;
  networks.initNetworkPool(myoption.networkfiles);
  layergraph.read(myoption.layerfile,networks);

  //Test search implementation.
  isearch.test(layergraph,networks);
  t.stop();
  if(g_verbosity>=VERBOSE_ESSENTIAL)
  std::cerr <<"Elapsed time: "<< t <<std::endl;

  return 1;
}
