/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/subnet.h
Description: Data structure of a conserved subnetwork.
**/

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include "macro.h"

template<typename NP, typename LG>
class SubNet
{
private:
  typedef NP NetworkPool;
  typedef LG LayerGraph;
  typedef typename LayerGraph::Graph Graph;
public:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Labels of the nodes
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
  /// Mapping from labels to original nodes
  typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;

  Graph* data;
  OrigLabelNodeMap *node2label;
  InvOrigLabelNodeMap label2node;

  unsigned _numSpecies;// The number of observed species.
  struct _K_Spine
  {
	  std::array<Node,RESERVED_SPECIES> data;
	  std::array<bool,RESERVED_SPECIES> states;
	  _K_Spine():data()
	  {
		  states.fill(false);
	  }
	  ~_K_Spine(){};
  };
  typedef struct _K_Spine K_Spine;
  std::vector<K_Spine> net_spines;

  SubNet(unsigned);
  ~SubNet(){};
};

template<typename NP, typename LG>
SubNet<NP,LG>::SubNet(unsigned k=5):
net_spines()
{
	_numSpecies=k;
	//initSubNet();
}
