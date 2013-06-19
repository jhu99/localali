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
	  OrigLabelNodeMap *node2label;
	  InvOrigLabelNodeMap *label2node;
	  _GraphData()
	  {
		  g=new Graph();
		  node2label=new OrigLabelNodeMap(*g);
		  label2node=new InvOrigLabelNodeMap();
	  }
	  ~_GraphData()
	  {
		  delete g;
		  delete node2label;
		  delete label2node;
	  }
  }GraphData;
  std::vector<K_Spine> net_spines;
  std::vector<GraphData*> subgraphs;

  SubNet(unsigned,unsigned);
  ~SubNet();
  bool induceSubgraphs(NetworkPool&, LayerGraph&);
  bool clearStructure();
};

template<typename NP, typename LG>
SubNet<NP,LG>::SubNet(unsigned num1, unsigned num2):
net_spines()
{
	_numSpecies=num1;
	_seedSize=num2;
	//initSubNet();
}

template<typename NP, typename LG>
SubNet<NP,LG>::~SubNet()
{
	for(unsigned i=0;i<subgraphs.size();i++)
		delete subgraphs[i];
}

template<typename NP, typename LG>
bool
SubNet<NP,LG>::clearStructure()
{
	for(unsigned i=0;i<_numSpecies;++i)
		delete subgraphs[i];
	return true;
}

template<typename NP, typename LG>
bool
SubNet<NP,LG>::induceSubgraphs(NetworkPool& networks, LayerGraph& layergraph)
{
	for(unsigned i=0;i<_numSpecies;++i)
	{
		GraphData* graphdata = new GraphData();
		std::vector<std::string> nodeset;
		for(unsigned j=0;j<_seedSize;++j)
		{
			std::string element=layergraph.node2label[net_spines[j].data[i]];
			if(find(nodeset.begin(),nodeset.end(),element)!=nodeset.end())continue;
			nodeset.push_back(element);
			Node node=graphdata->g->addNode();
			graphdata->node2label->set(node,element);
			(*graphdata->label2node)[element]=node;
			assert((*networks.getGraph(i)->invIdNodeMap).find(element)!=(*networks.getGraph(i)->invIdNodeMap).end());
		}
		for(unsigned p1=0;p1<nodeset.size();p1++)
		{
			std::string protein1=nodeset[p1];
			Node node1=(*graphdata->label2node)[protein1];
			for(unsigned p2=p1+1;p2<nodeset.size();p2++)
			{
				std::string protein2=nodeset[p2];
				std::string keystr;
				Node node2=(*graphdata->label2node)[protein2];
				if(protein1.compare(protein2)>0)
			     {
					 std::string tempstr = protein1;
					 protein1 = protein2;
					 protein2 = tempstr;
				 }
				 keystr.append(protein1);
				 keystr.append(protein2);
				 
				 if(networks.getGraph(i)->interactionmap.find(keystr)
				 ==networks.getGraph(i)->interactionmap.end())continue;
				 graphdata->g->addEdge(node1,node2);
			}
		}
		subgraphs.push_back(graphdata);
	}

	return true;
}
