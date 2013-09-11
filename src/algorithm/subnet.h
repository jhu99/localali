/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/subnet.h
Description: Data structure of a conserved subnetwork.
**/
#pragma once
#ifndef SUBNET_H_
#define SUBNET_H_

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include "algorithm/function.h"
#include "macro.h"
#include "function.h"
#include "score.h"

template<typename NP, typename LG>
class SubNet
{
private:
  typedef NP NetworkPool;
public:
  typedef LG LayerGraph;
  typedef typename NetworkPool::Graph Graph;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Labels of the nodes.
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
  /// Degree of the nodes.
  typedef typename Graph::template NodeMap<int> DegreeMap;
  /// Mapping from labels to original nodes.
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
		void clear()
		{
		}
  }K_Spine;
  
  typedef struct _GraphData
  {
	  Graph *g;
	  int nodeNum;
	  int edgeNum;
	  int maxDegree;
	  DegreeMap *node2degree;
	  OrigLabelNodeMap *node2label;
	  InvOrigLabelNodeMap *label2node;
	  InvLabelEdgeMap *label2edge;
	  
	  std::vector<unsigned> offsprings;
	  _GraphData()
	  :nodeNum(0)
	  ,edgeNum(0)
		,maxDegree(0)
	  {
		  g=new Graph();
		  node2degree=new DegreeMap(*g);
		  node2label=new OrigLabelNodeMap(*g);
		  label2node=new InvOrigLabelNodeMap();
		  label2edge=new InvLabelEdgeMap();
	  }
	  ~_GraphData()
	  {
		  delete node2degree;
		  delete node2label;
		  delete label2node;
		  delete label2edge;
			delete g;
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
	  void deleteEdge(Edge& myedge,std::string& edgelabel)
	  {
		  Node node1,node2;
		  node1=g->u(myedge);
		  node2=g->v(myedge);
		  int degree1,degree2;
		  degree1=(*node2degree)[node1]--;
		  degree2=(*node2degree)[node2]--;
		  edgeNum--;
		  g->erase(myedge);
		  label2edge->erase(edgelabel);
		  if(degree1==maxDegree || degree2==maxDegree)
		  {
			  maxDegree--;
			  for(NodeIt node(*g);node!=lemon::INVALID;++node)
			  {
				  if((*node2degree)[node]>maxDegree)maxDegree++;
			  }
		  }
	  }	
	  void addEdge(Node& node1,Node& node2)
	  {
		  int degree1,degree2;
		  std::string edgelabel;
		  degree1=++(*node2degree)[node1];
		  degree2=++(*node2degree)[node2];
		  edgeNum++;
		  Edge newedge=g->addEdge(node1,node2);
		  if(degree1>maxDegree || degree2>maxDegree)maxDegree++;
		  edgelabel=formEdgeLabel(node1,node2);
		  (*label2edge)[edgelabel]=newedge;
	  }
  }GraphData;
    
  std::vector<K_Spine> net_spines;
  std::vector<GraphData*> subgraphs;

  SubNet();
  SubNet(unsigned,unsigned);
  ~SubNet();
  bool induceSubgraphs(NetworkPool&, LayerGraph&);
  bool clearStructure();
	void output(LayerGraph&,std::ofstream&);
	void outputSubgraphs(LayerGraph&,std::string&,int,int,int);
	void outputAlignment(float,LayerGraph&,std::string&,int,int,int);
};

template<typename NP, typename LG>
SubNet<NP,LG>::SubNet():
net_spines(),
subgraphs()
{
}

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
SubNet<NP,LG>::outputSubgraphs(LayerGraph& layergraph,std::string& folder,int seednum,int extnum,int triesnum)
{
	std::string filename;
	std::ofstream fout;
	GraphData* graphdata;
	for(unsigned i=0;i<subgraphs.size();i++)
	{
		filename=folder;
		filename.append("species_");
		filename.append(convert_num2str(i));
		filename.append("/complex_s");
		filename.append(convert_num2str(seednum));
		filename.append("_e");
		filename.append(convert_num2str(extnum));
		filename.append("_t");
		filename.append(convert_num2str(triesnum));
		filename.append(".txt");
		fout.open(filename.c_str(),std::ofstream::out | std::ofstream::trunc);
		graphdata=subgraphs[i];
		for(NodeIt it(*graphdata->g);it!=lemon::INVALID;++it)
		{
			fout << (*graphdata->node2label)[it] << std::endl;
		}
		fout.close();
	}
}

template<typename NP, typename LG>
void
SubNet<NP,LG>::output(LayerGraph& layergraph,std::ofstream& fout)
{
	for(unsigned i=0;i<net_spines.size();i++)
	{
		for(unsigned j=0;j<_numSpecies;j++)
		{
			if(g_verbosity>=VERBOSE_ESSENTIAL)
			fout << layergraph.node2label[net_spines[i].data[j]]<<"\t";
		}
		fout << std::endl;
	}
}

template<typename NP, typename LG>
void
SubNet<NP,LG>::outputAlignment(float overallscore,LayerGraph& layergraph,std::string& folder,int seednum,int extnum,int triesnum)
{
	std::string filename(folder);
	filename.append("alignments/ucomplex_s");
	filename.append(convert_num2str(seednum));
	filename.append("_e");
	filename.append(convert_num2str(extnum));
	filename.append("_t");
	filename.append(convert_num2str(triesnum));
	filename.append(".txt");
	std::ofstream fout(filename.c_str());
	fout <<"#Score:"<<overallscore<<std::endl;
	output(layergraph,fout);
	outputSubgraphs(layergraph,folder,seednum,extnum,triesnum);
	fout.close();
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
	if(subgraphs.size()>0)
	{
		subgraphs.clear();
	}
	for(unsigned i=0;i<_numSpecies;++i)
	{
		GraphData* graphdata = new GraphData();
		graphdata->offsprings.push_back(i);
		std::vector<std::string> nodeset;
		for(unsigned j=0;j<net_spines.size();++j)
		{
			std::string element=layergraph.node2label[net_spines[j].data[i]];
			if(find(nodeset.begin(),nodeset.end(),element)!=nodeset.end())continue;
			nodeset.push_back(element);
			Node node=graphdata->g->addNode();
			graphdata->nodeNum++;
			graphdata->node2label->set(node,element);
			graphdata->node2degree->set(node,0);
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
					keystr.append(protein2);
					keystr.append(protein1);
				 }
				else
				{
				 keystr.append(protein1);
				 keystr.append(protein2);
				}
				 
				if(networks.getGraph(i)->interactionmap.find(keystr)==networks.getGraph(i)->interactionmap.end())continue;
				 graphdata->g->addEdge(node1,node2);
				 (*graphdata->node2degree)[node1]++;
				 (*graphdata->node2degree)[node2]++;
				 if((*graphdata->node2degree)[node1]>graphdata->maxDegree)
					 graphdata->maxDegree=(*graphdata->node2degree)[node1];
				 if((*graphdata->node2degree)[node2]>graphdata->maxDegree)
					 graphdata->maxDegree=(*graphdata->node2degree)[node2];
				 graphdata->edgeNum++;
			}
		}
		subgraphs.push_back(graphdata);
	}
	return true;
}
#endif
