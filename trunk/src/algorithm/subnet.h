/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/subnet.h
Description: Data structure of a conserved subnetwork.
**/
#ifndef SUBNET_H_
#define SUBNET_H_

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include "macro.h"
#include "algorithm/function.h"

template<typename NP, typename LG>
class SubNet
{
private:
  typedef NP NetworkPool;
public:
  typedef LG LayerGraph;
  typedef typename NetworkPool::Graph Graph;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Labels of the nodes
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
  /// Mapping from labels to original nodes
  typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;
  typedef std::unordered_map<std::string, typename Graph::Edge> InvLabelEdgeMap;

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
	  int nodeNum;
	  int edgeNum;
	  OrigLabelNodeMap *node2label;
	  InvOrigLabelNodeMap *label2node;
	  InvLabelEdgeMap *label2edge;
	  
	  std::vector<unsigned> offsprings;
	  _GraphData()
	  :nodeNum(0)
	  ,edgeNum(0)
	  {
		  g=new Graph();
		  node2label=new OrigLabelNodeMap(*g);
		  label2node=new InvOrigLabelNodeMap();
		  label2edge=new InvLabelEdgeMap();
	  }
	  ~_GraphData()
	  {
		  delete g;
		  delete node2label;
		  delete label2node;
		  delete label2edge;
	  }
	  std::string formEdgeLabel(Node& node1,Node& node2)
	  {
		  std::string label;
		  if(node2 < node1)
		  {
			  label.append(convert_num2str(g->id(node2)));
			  label.append(convert_num2str(g->id(node1)));
		  }else
		  {
			  label.append(convert_num2str(g->id(node1)));
			  label.append(convert_num2str(g->id(node2)));
		  }
		  return label;
	  }
  }GraphData;
    
  std::vector<K_Spine> net_spines;
  std::vector<GraphData*> subgraphs;

  SubNet(unsigned,unsigned);
  ~SubNet();
  bool induceSubgraphs(NetworkPool&, LayerGraph&);
  bool clearStructure();
  void output(LayerGraph&);
};

template<typename NP, typename LG>
SubNet<NP,LG>::SubNet(unsigned num1, unsigned num2):
net_spines(),
subgraphs()
{
	_numSpecies=num1;
	_seedSize=num2;
}

template<typename NP, typename LG>
SubNet<NP,LG>::~SubNet()
{
	for(unsigned i=0;i<subgraphs.size();i++)
		delete subgraphs[i];
}

template<typename NP, typename LG>
void
SubNet<NP,LG>::output(LayerGraph& layergraph)
{
	for(unsigned i=0;i<net_spines.size();i++)
	{
		for(unsigned j=0;j<_numSpecies;j++)
		{
			if(g_verbosity>=VERBOSE_ESSENTIAL)
			std::cout << layergraph.node2label[net_spines[i].data[j]]<<"\t";
		}
		std::cout << std::endl;
	}
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
	if(subgraphs.size()>0)subgraphs.clear();
	for(unsigned i=0;i<_numSpecies;++i)
	{
		GraphData* graphdata = new GraphData();
		graphdata->offsprings.push_back(i);
		std::vector<std::string> nodeset;
		for(unsigned j=0;j<_seedSize;++j)
		{
			std::string element=layergraph.node2label[net_spines[j].data[i]];
			if(find(nodeset.begin(),nodeset.end(),element)!=nodeset.end())continue;
			nodeset.push_back(element);
			Node node=graphdata->g->addNode();
			graphdata->nodeNum++;
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
				 graphdata->edgeNum++;
			}
		}
		subgraphs.push_back(graphdata);
	}

	return true;
}
#endif
