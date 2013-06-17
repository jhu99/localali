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

  unsigned _numSpecies;// The number of observed species.
  unsigned _seedSize;
  typedef struct _K_Spine
  {
	  std::array<Node,RESERVED_SPECIES> data;
	  std::array<bool,RESERVED_SPECIES> states;
	  _K_Spine():data()
	  {
		  states.fill(false);
	  }
	  ~_K_Spine(){};
  }K_Spine;
  typedef struct _GraphData
  {
	  Graph *g;
	  _GraphData()
	  {
		  g=new Graph();
	  }
	  ~_GraphData(){delete g;}
  }GraphData;
  std::vector<K_Spine> net_spines;
  std::vector<GraphData*> subgraphs;

  SubNet(unsigned);
  ~SubNet(){};
  bool induceSubgraphs(NetworkPool&, LayerGraph&);
};

template<typename NP, typename LG>
SubNet<NP,LG>::SubNet(unsigned k=5):
net_spines()
{
	_numSpecies=k;
	_seedSize=3;
	//initSubNet();
}

template<typename NP, typename LG>
bool
SubNet<NP,LG>::induceSubgraphs(NetworkPool& networks, LayerGraph& layergraph)
{
	//for(unsigned i=0;i<_numSpecies;++i)
	//{
		//Graph* graph= new Graph();
		//std::vector<std::string> nodeset;
		//for(unsigned j=0;j<_seedSize;++j)
		//{
			//std::string element=layergraph.node2label[subnet.net_spines[j].data[i]];
			//if(find(nodeset.begin(),nodeset.end(),element)!=nodeset.end())continue;
			//nodeset.push_back(element);
		//}
		//for(int p1=0;p1<nodeset.size();p1++)
		//{
			//std::string protein1=nodeset[p1];
			//for(int p2=p1+1;p2<nodeset.size();p2++)
			//{
				//std::string protein2=nodeset[p2];
			//}
		//}
	//}

	return true;
}
