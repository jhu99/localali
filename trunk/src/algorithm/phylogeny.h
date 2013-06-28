/**
Author: Jialu Hu
Date: Jun. 25, 2013
File name: algorithm/phylogeny.h
Description: Data structure of the phylogeny of the observed functional modules.
**/

#ifndef PHYLOGENY_H_
#define PHYLOGENY_H_

#include <vector>
#include <array>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include "macro.h"

template<typename SN, typename TR>
class Phylogeny
{
private:
	typedef SN SubNet;
	typedef TR Tree;
	typedef typename SubNet::GraphData GraphData;
	typedef typename SubNet::LayerGraph LayerGraph;
	typedef typename Tree::Graph Graph;
	TEMPLATE_GRAPH_TYPEDEFS(Graph);
	//typedef struct _Distance
	//{
	//}Distance;
public:
	std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution;
	Tree _tree;
	std::vector<std::string> _speciesfiles;
	std::vector<Node> internalNode;
	std::vector<Node> externalNode;
	std::unordered_map<int,GraphData*> node2graph;
    
	Phylogeny(std::string&);
	~Phylogeny(){};
	bool initial(SubNet&,LayerGraph&);
	bool initialExternalNodes(SubNet&);
	bool interfere();
	bool existNode(std::vector<std::string>&,GraphData*);
	GraphData* constructInternalNodes(Node&,SubNet&,LayerGraph&);
};

template<typename SN, typename TR>
Phylogeny<SN,TR>::Phylogeny(std::string& filename)
:generator(std::chrono::system_clock::now().time_since_epoch().count())
,distribution(0,10000)
,_tree()
,node2graph()
{
	_tree.readTree(filename);
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::interfere()
{
	//float delta=0.0;
	int dice_roll=distribution(generator)%internalNode.size();
	int upperId,dice_1,dice_2;
	Node intNode = internalNode[dice_roll];
	GraphData* graphdata=node2graph[_tree.g.id(intNode)];
	upperId=graphdata->nodeNum;
	if(upperId<2) return false;
	Node nodeA,nodeB;
	dice_1 = distribution(generator)%upperId;
	nodeA = graphdata->g->nodeFromId(dice_roll);
	dice_2 = distribution(generator)%upperId;
	while(dice_1==dice_2)
	{
		dice_2 = distribution(generator)%upperId;
	}
	nodeB = graphdata->g->nodeFromId(dice_2);
	std::string edgelabel=graphdata->formEdgeLabel(nodeA,nodeB);
	if(graphdata->label2edge->find(edgelabel)!=graphdata->label2edge->end())
	{
		//delete edge
		graphdata->g->erase(graphdata->label2edge->find(edgelabel)->second);
		graphdata->label2edge->erase(edgelabel);
	}
	else{
		// add edge
		Edge newedge=graphdata->g->addEdge(nodeA,nodeB);
		(*graphdata->label2edge)[edgelabel]=newedge;
	}
	
	return true;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::initialExternalNodes(SubNet& subnet)
{
	for(unsigned i=0;i<_speciesfiles.size();i++)
	{
		std::string species=_speciesfiles[i];
		Node node=_tree.label2node[species];
		node2graph[_tree.g.id(node)]=subnet.subgraphs[i];
		externalNode.push_back(node);
	}
	return true;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::existNode(std::vector<std::string>& xspine,GraphData* graphdata)
{
	for(unsigned i=0;i<xspine.size();i++)
	{
		if(graphdata->label2node->find(xspine[i])!=graphdata->label2node->end())
			return true;
	}
	return false;			
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::initial(SubNet& subnet,LayerGraph& layergraph)
{
	initialExternalNodes(subnet);
	GraphData* graphdata=constructInternalNodes(_tree.root,subnet,layergraph);
	node2graph[_tree.g.id(_tree.root)]=graphdata;
	internalNode.push_back(_tree.root);
	return true;
}

template<typename SN, typename TR>
typename Phylogeny<SN,TR>::GraphData*
Phylogeny<SN,TR>::constructInternalNodes(Node& ancestor,SubNet& subnet,LayerGraph& layergraph)
{
	GraphData* graphdata=new GraphData();
	std::vector<Node> sonnodes;
	for(IncEdgeIt e(_tree.g,ancestor);e!=lemon::INVALID;++e)
	{
		GraphData* sondata;
		Node rnode=_tree.g.runningNode(e);
		if(ancestor<rnode)continue;
		sonnodes.push_back(rnode);
		if(node2graph.find(_tree.g.id(rnode))==node2graph.end())
		{
			sondata=constructInternalNodes(rnode,subnet,layergraph);
			node2graph[_tree.g.id(rnode)]=sondata;
			internalNode.push_back(rnode);
		}
		for(unsigned i=0;i<sondata->offsprings.size();++i)
			graphdata->offsprings.push_back(sondata->offsprings[i]);
	}
	/// construct internal nodes.
	typedef std::list<std::vector<std::string> > SpineList;
	SpineList incSpines;
	for(unsigned i=0;i<subnet.net_spines.size();++i)
	{
		std::vector<std::string> xspine;
		for(unsigned j=0;j<graphdata->offsprings.size();++j)
		{
			unsigned nspecies=graphdata->offsprings[j];
			Node leavenode=subnet.net_spines[i].data[nspecies];
			xspine.push_back(layergraph.node2label[leavenode]);
		}
		incSpines.push_back(xspine);
	}
	while(!incSpines.empty())
	{
		Node gnode = graphdata->g->addNode();
		std::vector<std::string> xspine=incSpines.front();
		incSpines.pop_front();
		for(std::vector<std::string>::iterator it=xspine.begin();it!=xspine.end();++it)
		{
			(*graphdata->label2node)[*it]=gnode;
		}
		for(SpineList::iterator it=incSpines.begin();it!=incSpines.end();++it)
		{
			xspine=*it;
			if(!existNode(xspine,graphdata))continue;
			for(unsigned i=0;i<xspine.size();i++)
			{
				(*graphdata->label2node)[xspine[i]]=gnode;
			}
			incSpines.erase(it);
		}
	}
	return graphdata;
}

#endif
