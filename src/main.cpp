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
    seedtries=1;
	minext=6;
	maxext=13;
	extdist1=1;
	extdist2=2;
    numseeds=400;
    numconnected=3;
    numthreads=1;
	score_threshold=0.2;
	numspinetries=5;
	verbose=0;
	beta=2.0;
	alpha=0.2;
	parallel=false;
	seedReplication=0;
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
typedef Analyse<InputGraph> MyAnalyse;

bool setParser(ArgParser& parser, Option& myoption)
{
	parser
	.boolOption("version","Show the version number.")
	.boolOption("alignment","Execute the alignment algorithm.")
	.boolOption("analyse","Make analysis on the alignment results.")
	.boolOption("format","Process input or output file into proper format.")
	.optionGroup("method","version")
	.optionGroup("method","alignment")
	.optionGroup("method","analyse")
	.optionGroup("method","format")
	.onlyOneGroup("method")
	.mandatoryGroup("method")
	.refOption("task","Specify the task of each method. Default is 0.", myoption.task)
	.refOption("method","Specify the method used for verification. LocalAli 1, NetworkBlastM 2. Default is 1.", myoption.method)
	.refOption("profile","Configuration of various input parameters. Default is \"./profile.txt\".", myoption.profile)
	.refOption("resultfolder","Configuration of various input parameters.", myoption.resultfolder)
	.refOption("formatfile","Input file which is used to analyse the quality of alignments.",myoption.formatfile)
	.refOption("numspecies","Number of the species compared. Default is 3.", myoption.numspecies)
	.refOption("seedtries","Number of tries for each refined seeds. Default is 1.", myoption.seedtries)
	.refOption("seedsize","Size of the seeds. Default is 2.", myoption.seedsize)
	.refOption("minext","Minimal number of the extension. Default is 11.", myoption.minext)
	.refOption("maxext","Maximal number of the extension. Default is 12.", myoption.maxext)
	.refOption("extdist1","Distance of neighbors in the search for seeds. Default is 1.", myoption.extdist1)
	.refOption("extdist2","Distance of neighbors in the search for subnets. Default is 2.", myoption.extdist2)
	.refOption("numconnected","Number of connected subnetwork. Default is 3.", myoption.numconnected)
	.refOption("numseeds","Number of refined seeds. Default is 200.", myoption.numseeds)
	.refOption("numthreads","Number of threads. Default is 1.", myoption.numthreads)
	.refOption("numspinetries","Number of tries for strongly connected spines. Default is 5.", myoption.numspinetries)
	.refOption("score_threshold","Score threshold of subnets which are qualified. Default is 0.2.", myoption.score_threshold)
	.refOption("alpha","Impact factor of the evolutionary rate. Default is 0.2", myoption.alpha)
	.refOption("beta","The second impact factor of the evolutionary rate of interactions. Default is 2.0.", myoption.beta)
	.refOption("verbose","Display standard output levle:0-3. Default is 0.", myoption.verbose)
	.refOption("parallel","Run LocalAli in parallel if it is true. Default is false.", myoption.parallel)
	.refOption("seedrep","Allow protein replicatioin in refined seeds. Default is false.", myoption.seedReplication);
	return true;
}

bool runParser(ArgParser& myparser, Option& myoption)
{
  ProcessProfile<Option> myprofile;
  myparser.run();
  myprofile.getOption(myoption);
  if(myoption.parallel)
  {
	  std::cout << "# This program will run with "<< myoption.numthreads <<" multiple threads." << std::endl;
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

  	g_verbosity=(VerbosityLevel)myoption.verbose;//VERBOSE_NON_ESSENTIAL;//VERBOSE_ESSENTIAL;
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
			myformat.extractIntActInteractions();
			//myformat.extractInteractions();
		}
		else if(myoption.task==1)
			// Convert homology-list-20130826.b6.data -> homology-list-20130826.evals and homology-list-20130826.bscore
		{
			networks.initNetworkPool(myoption.networkfiles);
			myformat.extractIntActHomology(networks);
			//myformat.extractHomology();
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
		}else if(myoption.task==6)
		{
			myformat.partitionGOA(myoption.formatfile);
		}
		else if(myoption.task==7)
		{
			myformat.generateAlignNemoPPI(myoption.formatfile);
		}
		else if(myoption.task==8)
		{
			myformat.generateAlignNemoSim(myoption.formatfile);
		}
		else if(myoption.task==9)
		{
			networks.initNetworkPool(myoption.networkfiles);
			myformat.convertAlignNemoNif(myoption.resultfolder,networks);
		}
		else if(myoption.task==10)
		{
			networks.initNetworkPool(myoption.networkfiles);
			myformat.convertNetBlastProp(myoption.resultfolder,networks);
		}
	}
	else if(myparser.given("analyse"))
	{
		MyAnalyse myanalyse(myoption.resultfolder);
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
			myanalyse.assessQuality(myoption.resultfolder,myoption.numspecies);
		}
		else if(myoption.task==2)
		// translate DIP subnetworks to Uniprot subnetworks and remove redundant alignments;
		{
			networks.initNetworkPool(myoption.networkfiles);
			//myanalyse.readIdMap();
			//myanalyse.translate_alignment(networks, myoption.resultfolder,myoption.numspecies);
			myanalyse.reduceRedundancy(networks, myoption.resultfolder,myoption.numspecies);
			std::cout << "The percentage of covered proteins: " << myanalyse.numCoveredProtein/static_cast<float>(networks.allNodeNum) << std::endl;
		}else if(myoption.task==3)
		{
			myanalyse.predictFunction(myoption.resultfolder,myoption.numspecies);
			myanalyse.countPrediction(myoption.formatfile);
		}
		else if(myoption.task==4)
		{
			myanalyse.countPrediction(myoption.formatfile);
		}
		else if(myoption.task==5)
		{
			myanalyse.countCrossVerification(myoption.resultfolder, myoption.method);
		}
		else if(myoption.task==6)
		{
			myanalyse.verifyPrediction(myoption.resultfolder);
		}
		else if(myoption.task==7)
		{
			myanalyse.verifyComplexes(myoption.formatfile);
		}
		else if(myoption.task==8)
		{
			myanalyse.timecheck(myoption.resultfolder);
		}
		else if(myoption.task==9)
		{
			myanalyse.ppvcheck(myoption.resultfolder);
		}
	}
	t.stop();
	if(g_verbosity>=VERBOSE_ESSENTIAL)
		std::cout <<"# Elapsed time: "<< t <<std::endl;
  return 1;
}
