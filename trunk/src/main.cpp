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

typedef lemon::ListGraph Graph;
typedef lemon::SmartGraph BpGraph;
typedef NetworkPool<Graph,BpGraph> InputGraph;
typedef Layer_graphs<BpGraph,InputGraph> LayerGraph;
typedef SubNet<InputGraph,LayerGraph> MySubNet;
typedef Search<InputGraph, MySubNet, LayerGraph> MySearch;

typedef struct _Option
{
  vector<std::string> networkfiles;
  std::string layerfile;
  std::string resultfolder;
  std::string profile;
  float beta;
  _Option()
  {
    profile="./profile.txt";
  }
}Option;

bool setParser(ArgParser& parser, Option& myoption)
{
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
  InputGraph networks;
  LayerGraph layergraph;
  MySubNet isubnet(5);
  MySearch isearch;
  g_verbosity=VERBOSE_NON_ESSENTIAL;

  runParser(myparser,myoption);

  //Test read interface for PPI networks;
  networks.initNetworkPool(myoption.networkfiles);
  layergraph.read(myoption.layerfile,networks);

  //Test search implementation.
  isearch.test(layergraph,networks);

  return 1;
}
