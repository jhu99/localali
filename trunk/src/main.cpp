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
#include "omp.h"
#include "input/format.h"
#include "analyse/analyse.h"

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
  double beta;
	double score_threshold;
	int task;
  int numspecies;
  int seedsize;
  int seedtries;
  int minext;
  int maxext;
  int numsamples;
  int numconnected;
  int numthreads;
  bool parallel;
  _Option()
  {
    profile="./profile.txt";
		task=0;
    numspecies=3;
    seedsize=3;
    seedtries=1;
		minext=6;
		maxext=12;
    numsamples=1000;
    numconnected=3;
    numthreads=1;
		score_threshold=0.4;
		parallel=false;
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

bool setParser(ArgParser& parser, Option& myoption)
{
	parser
	.boolOption("version","Show the version number.")
  .boolOption("alignment","Execute the alignment algorithm.")
  .boolOption("analyse","Make analysis on alignmenr result.")
  .boolOption("format","Process input or output file into proper format.")
	.optionGroup("method","version")
  .optionGroup("method","alignment")
  .optionGroup("method","analyse")
  .optionGroup("method","format")
  .onlyOneGroup("method")
  .mandatoryGroup("method")
	.refOption("task","Specify the task of each method. Default is 0.", myoption.task)
	.refOption("profile","Configuration of various input parameters. Default is \"./profile.txt\".", myoption.profile)
	.refOption("resultfolder","Configuration of various input parameters. Default is \"./result/dip/3-way/localali/", myoption.resultfolder)
	.refOption("numspecies","Number of the species compared. Default is 3.", myoption.numspecies)
	.refOption("seedtries","Number of tries for each refined seeds. Default is 3.", myoption.seedtries)
	.refOption("seedsize","Size of the seeds. Default is 3.", myoption.seedsize)
	.refOption("minext","Minimal number of the extension . Default is 1.", myoption.minext)
	.refOption("maxext","Maximal number of the extension . Default is 2.", myoption.maxext)
	.refOption("numconnected","Number of connected subnetwork. Default is 2.", myoption.numconnected)
	.refOption("numsamples","Number of sampled seeds. Default is 40000.", myoption.numsamples)
	.refOption("numthreads","Number of threads. Default is 1.", myoption.numthreads)
	.refOption("score_threshold","Score threshold of subnets which are qualified. Default is 0.4.", myoption.score_threshold)
	.refOption("parallel","Run LocalAli in parallel. Default is false.", myoption.parallel);
	return true;
}

bool runParser(ArgParser& myparser, Option& myoption)
{
  ProcessProfile<Option> myprofile;
  myparser.run();
  myprofile.getOption(myoption);
  if(myoption.parallel)
  {
	  std::cout << "This program will run with "<< myoption.numthreads <<" multiple threads." << std::endl;
  }
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
  MySearch localAlignment(myoption);
  Timer t(false);

  g_verbosity=VERBOSE_ESSENTIAL;//VERBOSE_ESSENTIAL;
	t.start();

	if(myparser.given("alignment"))
	{
  // Read interface for PPI networks;
		networks.initNetworkPool(myoption.networkfiles);
		layergraph.read(myoption.layerfile,networks);
		localAlignment.run(layergraph,networks);
	}
	else if(myparser.given("format"))
	{
		MyFormat myformat(myoption);
		if(myoption.task==0)
			// Convert Celeg20130707.txt.data etc.. -> Celeg20130707-int.txt (interactorA	interactorB 0.9)
		{
			myformat.extractInteractions();
		}
		else if(myoption.task==1)
			// Convert homology-list-20130826.b6.data -> homology-list-20130826.evals and homology-list-20130826.bscore
		{
			myformat.extractHomology();
		}
		else if(myoption.task==2)
			// Convert homology-list-20130826.evals -> Celeg20130707-Celeg20130707.evals Celeg20130707-Dmela20130707.evals...
		{
			networks.initNetworkPool(myoption.networkfiles);
			myformat.extractDatasetHomology(networks);
		}
		else if(myoption.task==3)
		{
			networks.initNetworkPool(myoption.networkfiles);
			myformat.extractDipAc(networks);
		}
		else if(myoption.task==4)
			// Generate a subtree for a set of species.
		{
			MyTree gtree;
			gtree.generateSubTree(myoption.treefile,myoption.speciesfiles);
		}
		else if(myoption.task==5)
			// write networkblast-m complexes with alignment score
		{
			myformat.writeSubnetworks(myoption);
		}else
		{
		}
	}
	else if(myparser.given("analyse"))
	{
		Analyse myanalyse;
		if(myoption.task==0)
		// translate DIP subnetworks to Uniprot subnetworks and remove redundant subnetworks;
		{
			myanalyse.readIdMap();
			myanalyse.translate(myoption.resultfolder,myoption.numspecies);
			// use gotermfinder-local.sh to calculate p-value of each subnetwork;
		}
		else if(myoption.task==1)
		// Assess the quality of subnetworks;
		{
			networks.initNetworkPool(myoption.networkfiles);
			myanalyse.assessQuality(myoption.resultfolder,myoption.numspecies);
			std::cout << "The percentage of covered proteins: " << myanalyse.numCoveredProtein/static_cast<float>(networks.allNodeNum) << std::endl;
		}
		else if(myoption.task==2)
		// translate DIP subnetworks to Uniprot subnetworks and remove redundant alignments;
		{
			myanalyse.readIdMap();
			myanalyse.translate_alignment(myoption.resultfolder,myoption.numspecies);
		}else
		{
		}
	}
	t.stop();
	if(g_verbosity>=VERBOSE_ESSENTIAL)
		std::cout <<"Elapsed time: "<< t <<std::endl;
  return 1;
}
