/**
Author: Jialu Hu
Date: Jun. 25, 2013
File name: algorithm/phylogeny.h
Description: Data structure of the phylogeny of the observed functional modules.
**/
#pragma once

#ifndef PHYLOGENY_H_
#define PHYLOGENY_H_

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <iostream>
#include <unordered_map>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include "macro.h"
#include "input/tree.h"

template<typename SN, typename TR>
class Phylogeny
{
private:
	typedef SN SubNet;
	typedef TR Tree;
	typedef typename SubNet::LayerGraph LayerGraph;
	typedef typename Tree::MatchingNodeMap MatchingNodeMap;
public:
	typedef typename SubNet::GraphData GraphData;
	typedef typename Tree::Graph Graph;
	TEMPLATE_GRAPH_TYPEDEFS(Graph);
	typedef struct _DeltaStructure
	{
		Node treenode;
		Node nodeA;
		Node nodeB;
		std::string edgelabel;
		float delta;
		Score updatedScores[3];// The neighbors of a node in a binary tree is upto 3.
		/// maxDegree and node degree will change along with the inference.
	}DeltaStructure;

	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution;
	Tree _tree;
	int _dsize;
	std::string _treefile;
	std::vector<std::string> _speciesfiles;
	std::vector<Node> internalNode;
	std::vector<Node> externalNode;
	std::unordered_map<int,GraphData*> node2graph;

	//Phylogeny();
	Phylogeny();//dsize is the number of k-spines.
	~Phylogeny();
	bool initial(std::string&,std::vector<std::string>&,SubNet&,LayerGraph&);
	bool initialExternalNodes(SubNet&);
	bool interfere(DeltaStructure&);
	bool initialBranchWeight();
	bool computeBranchWeight(EdgeIt&,Score&);
	bool computeScore(Score&,MatchingNodeMap*, GraphData*, GraphData*,float&);
	bool existNode(std::vector<std::string>&,GraphData*);
	bool computeDist();
	GraphData* constructInternalNodes(Node&,SubNet&,LayerGraph&);
	// Output the optimal graphs of internal nodes.
	void outputInternalGraphs(void);
	bool clearStructure(void);
	void setDsize(int);
};

template<typename SN, typename TR>
Phylogeny<SN,TR>::Phylogeny()
:generator(std::chrono::system_clock::now().time_since_epoch().count())
,distribution(0,10000)
,_tree()
,_dsize(0)
,_speciesfiles()
,internalNode()
,externalNode()
,node2graph()
{
}

template<typename SN, typename TR>
Phylogeny<SN,TR>::~Phylogeny()
{
}

template<typename SN, typename TR>
void
Phylogeny<SN,TR>::setDsize(int num)
{
	_dsize=num;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::clearStructure()
{
	internalNode.clear();
	externalNode.clear();
	node2graph.clear();
	return true;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::interfere(DeltaStructure& deltaStr)
{
	int dice_roll=distribution(generator)%internalNode.size();
	int upperId,dice_1,dice_2;
	Node intNode = internalNode[dice_roll];
	GraphData* graphdata=node2graph[_tree.g.id(intNode)];
	upperId=graphdata->nodeNum;
	if(upperId<2) return false;
	Node nodeA,nodeB;
	dice_1 = distribution(generator)%upperId;
	nodeA = graphdata->g->nodeFromId(dice_1);
	dice_2 = distribution(generator)%upperId;
	while(dice_1==dice_2)
	{
		dice_2 = distribution(generator)%upperId;
	}
	nodeB = graphdata->g->nodeFromId(dice_2);
	std::string edgelabel=graphdata->formEdgeLabel(nodeA,nodeB);
	deltaStr.treenode=intNode;
	deltaStr.nodeA=nodeA;
	deltaStr.nodeB=nodeB;
	deltaStr.edgelabel=edgelabel;
	if(graphdata->label2edge->find(edgelabel)!=graphdata->label2edge->end())
	{
		//delete edge
		Edge myedge=graphdata->label2edge->find(edgelabel)->second;
		graphdata->deleteEdge(myedge,edgelabel);
	}
	else{
		// add edge
		graphdata->addEdge(nodeA,nodeB);
	}

	// Compute the gap between the two scores, but without changes of their original score attritubtion.
	GraphData *ancestor, *descendant;
	Score deltaScore;
	int i=0;
	for(IncEdgeIt it(_tree.g,intNode);it!=lemon::INVALID;++it,++i)
	{
		Score updatedScore;
		Node neighbor=_tree.g.runningNode(it);
		descendant=node2graph[_tree.g.id(neighbor)];
		if(neighbor<intNode)
		{
			ancestor=graphdata;
		}else
		{
			ancestor=descendant;
			descendant=graphdata;
		}
		computeScore(updatedScore,_tree.matchingedgemap[it],ancestor,descendant,_tree.branchmap[it]);
		deltaScore+=updatedScore;
		deltaScore-=_tree.scoremap[it];
		deltaStr.updatedScores[i]=updatedScore;
	}
	deltaStr.delta=deltaScore.sumup();

	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
	std::cout << deltaStr.delta << std::endl;
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
Phylogeny<SN,TR>::computeDist()
{
	_tree.computeDistEvolution(_dsize);
	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
		std::cout << _tree.overallScore << std::endl;
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
Phylogeny<SN,TR>::initial(std::string& mytreefile,
						  std::vector<std::string>& myspeciesfiles,
						  SubNet& subnet,
						  LayerGraph& layergraph)
{
	_treefile=mytreefile;
	_speciesfiles=myspeciesfiles;
	_tree.readTree(_treefile);
	clearStructure();
	initialExternalNodes(subnet);
	GraphData* graphdata=constructInternalNodes(_tree.root,subnet,layergraph);
	node2graph[_tree.g.id(_tree.root)]=graphdata;
	internalNode.push_back(_tree.root);
	initialBranchWeight();
	return true;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::initialBranchWeight()
{
	//for(NodeIt in(_tree.g);in!=lemon::INVALID;++in)
	//{
	//	GraphData* graphdata=node2graph[_tree.g.id(in)];
	//	for(NodeIt node(*graphdata->g);node!=lemon::INVALID;++node)
	//	{
	//		graphdata->node2degree->set(node,0);
	//	}
	//	graphdata->maxDegree=0;
	//}
	for(EdgeIt ie(_tree.g);ie!=lemon::INVALID;++ie)
	{
		Score myscore;
		computeBranchWeight(ie,myscore);
		_tree.scoremap[ie]=myscore;
	}
	return true;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::computeScore(Score& score,MatchingNodeMap* matchingmap, GraphData* ancestor, GraphData* descedant,float& branchweight)
{
	int mykey;
	int maxdegree=ancestor->maxDegree;
	if(0==maxdegree)maxdegree=1;
	for(NodeIt anode(*ancestor->g);anode!=lemon::INVALID;++anode)
	{
		mykey=ancestor->g->id(anode);
		int numElement=matchingmap->count(mykey);
		if(numElement==1)
		/// Protein mutation.
		{
			score.fscore[0]+=(1-static_cast<float>((*ancestor->node2degree)[anode])/maxdegree)*branchweight;
		}else
		/// Protein duplication.
		{
			score.fscore[1]+=(1-static_cast<float>((*ancestor->node2degree)[anode])/maxdegree)*branchweight;
		}
	}
	/// The number of conserved edges.
	int conNum=0;
	for(EdgeIt edge(*ancestor->g);edge!=lemon::INVALID;++edge)
	{
		bool finderFlag(false);
		int id1,id2;
		Node node1,node2,node3,node4;
		node1=ancestor->g->u(edge);
		node2=ancestor->g->v(edge);
		id1=ancestor->g->id(node1);
		id2=ancestor->g->id(node2);
		auto range1=matchingmap->equal_range(id1);
		auto range2=matchingmap->equal_range(id2);
		for(typename MatchingNodeMap::iterator it1=range1.first;it1!=range1.second;++it1)
		{
			node3=it1->second;
			for(typename MatchingNodeMap::iterator it2=range2.first;it2!=range2.second;++it2)
			{
				node4=it2->second;
				std::string edgelabel=descedant->formEdgeLabel(node3,node4);
				if(descedant->label2edge->find(edgelabel)!=descedant->label2edge->end())
				{
					conNum++;
					finderFlag=true;
					break;
				}
			}
			if(finderFlag)break;
		}
	}
	// Interaction deletion.
	score.fscore[2]=_tree._beta*branchweight*(ancestor->edgeNum-conNum);
	// Interaction insertion.
	score.fscore[3]=_tree._beta*branchweight*(descedant->edgeNum-conNum);
	return true;
}

template<typename SN, typename TR>
bool
Phylogeny<SN,TR>::computeBranchWeight(EdgeIt& ie,Score& score)
{
	typedef typename SubNet::InvOrigLabelNodeMap::iterator ItInvOrigLabelNodeMap;
	typedef typename Tree::MatchingNodeMap MatchingNodeMap;
	typedef typename Tree::MatchingNodeMap::iterator MatchingNodeMapIt;
	
	Node node1,node2,anode,dnode;//node1, node2 are tree nodes, anode and dnode are graph nodes.
	int mykey,id1,id2;
	MatchingNodeMap* matchingmap = new MatchingNodeMap();
	GraphData *descedant, *ancestor;
	node1=_tree.g.u(ie);
	id1=_tree.g.id(node1);
	node2=_tree.g.v(ie);
	id2=_tree.g.id(node2);
	if(node2<node1)
	{
		descedant=node2graph[id2];
		ancestor=node2graph[id1];
	}
	else
	{
		descedant=node2graph[id1];
		ancestor=node2graph[id2];
	}
	
	for(ItInvOrigLabelNodeMap it=descedant->label2node->begin();it!=descedant->label2node->end();++it)
	{
		std::string protein=it->first;
		dnode=it->second;
		anode=(*ancestor->label2node)[protein];
		mykey=ancestor->g->id(anode);
		auto range=matchingmap->equal_range(mykey);
		bool isExist=false;
		for(MatchingNodeMapIt it=range.first;it!=range.second;++it)
		{
			if(it->second==dnode){isExist=true;break;}
		}
		if(isExist)continue;
		matchingmap->insert(std::make_pair(mykey,dnode));// A general error may look like: std::make_pair<int, Node>(mykey,dnode)
	}
	_tree.matchingedgemap[ie]=matchingmap;
	computeScore(score, matchingmap, ancestor, descedant,_tree.branchmap[ie]);
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
		}else{
			sondata=node2graph[_tree.g.id(rnode)];
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
		graphdata->nodeNum++;
		graphdata->node2degree->set(gnode,0);
		for(std::vector<std::string>::iterator it=xspine.begin();it!=xspine.end();++it)
		{
			(*graphdata->label2node)[*it]=gnode;
			///There are many labels for each internal node.
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
			it=incSpines.begin();
		}
	}
	return graphdata;
}

template<typename SN, typename TR>
void
Phylogeny<SN,TR>::outputInternalGraphs()
{
	for(unsigned i=0;i<internalNode.size();i++)
	{
		Node node=internalNode[i];
		Node node1,node2;
		int nodeid=_tree.g.id(node);
		GraphData* graphdata=node2graph[nodeid];
		std::cout << "The graph structure of internal node " << nodeid <<" is:" << std::endl;
		for(EdgeIt ie(*graphdata->g);ie!=lemon::INVALID;++ie)
		{
			node1=graphdata->g->u(ie);
			node2=graphdata->g->v(ie);
			std::cout << graphdata->g->id(node1) <<"\t" << graphdata->g->id(node2) << std::endl;
		}
	}
}

#endif
