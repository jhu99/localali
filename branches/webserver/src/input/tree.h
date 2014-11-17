/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/tree.h
Description: Searching high-scoring subnetworks.
**/
#pragma once
#ifndef TREE_H
#define TREE_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <deque>
#include <stack>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include "verbose.h"
#include "string.h"
#include "algorithm/score.h"
#include "algorithm/function.h"
#include <lemon/dfs.h>

template<typename GR, typename OP>
class Tree
{
private:
	typedef OP Option;
public:
	typedef GR Graph;
	TEMPLATE_GRAPH_TYPEDEFS(Graph);
	
	/// Labels of the nodes.
	typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
	/// Mapping from labels to original nodes.
	typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;
	/// Weights on original edges.
	typedef typename Graph::template EdgeMap<float> WeightEdgeMap;
	/// Evolutionary score for each branch.
	typedef typename Graph::template EdgeMap<Score> ScoreEdgeMap;
	/// Nodes mapping for a pair of functional modules.
	typedef std::unordered_multimap<int, typename Graph::Node> MatchingNodeMap;
	/// Nodes matching map for each branch.
	typedef typename Graph::template EdgeMap<MatchingNodeMap*> MatchingEdgeMap;
	typedef std::unordered_map<std::string, EdgeIt> LabelEdgeMap;
	typedef lemon::Bfs<Graph> BfsType;
	
	/// Hold the topology of this tree.
	Graph g;
	/// Root of the evolutionary tree.
	Node root;
	/// The distance score for the evolutionary tree.
	Score distEvolution;
	/// The evolutionary rate for interaction deletion and insertion events.
	float _beta;
	float overallScore;
	/// The labels for the tree nodes.
	OrigLabelNodeMap node2label;
	/// The nodes for the labels.
	InvOrigLabelNodeMap label2node;
	/// The labels for each edges.
	LabelEdgeMap	label2edge;
	/// The lengths for the branches in the tree.
	WeightEdgeMap branchmap;
	/// The local alignment score for each branch in the tree.
	ScoreEdgeMap  scoremap;
	/// Nodes matching map for each branch.
	MatchingEdgeMap matchingedgemap; 

	Tree();
	~Tree(){};
	bool readTree(std::string);
	bool constructTree(std::string);
	void compressNodes(std::stack<std::string>&);
	bool computeDistEvolution(int);
	void extractSubtree(std::string sentence, std::vector<std::string>&);
	bool generateSubTree(std::string,std::vector<std::string>&);
};

template<typename GR, typename OP>
Tree<GR,OP>::Tree()
:g()
,distEvolution()
,_beta(1.5)
,node2label(g)
,label2node()
,branchmap(g)
,scoremap(g)
,matchingedgemap(g)
{
}

template<typename GR, typename OP>
bool
Tree<GR,OP>::readTree(std::string filename)
{
	std::ifstream input(filename);
	if(!input.is_open())
	{
		std::cerr <<"Can't find the path of the evolutionary tree!"<<std::endl;
		return false;
	}
	std::string line;
	std::getline(input,line);// read the nested-parentheses line.
	constructTree(line);
	return true;
}

template<typename GR, typename OP>
bool
Tree<GR,OP>::constructTree(std::string sentence)
{
	std::deque<std::string> wordpipe;
	std::deque<char> parenthesepipe;
	std::stack<std::string> constructor;
	char * pSentence = new char [sentence.length()+1];
	strcpy (pSentence, sentence.c_str());
	char * pch = strtok(pSentence,"(),:;");
	if(g_verbosity>VERBOSE_NON_ESSENTIAL)
	std::cout << sentence << std::endl;
	while(pch!=0)
	{
		std::string word(pch);
		if(g_verbosity>VERBOSE_NON_ESSENTIAL)
		std::cout <<word<<std::endl;
		wordpipe.push_back(word);
		pch=strtok(NULL,"(),:;");
	}	
	delete pSentence;
	std::size_t found = sentence.find_first_of("(,)");
	while(found!=std::string::npos)
	{
		parenthesepipe.push_back(sentence[found]);
		if(g_verbosity>VERBOSE_NON_ESSENTIAL)
		std::cout << sentence[found];
		found=sentence.find_first_of("(,);",found+1);
	}
	if(g_verbosity>VERBOSE_NON_ESSENTIAL)
	std::cout << std::endl;
	
	size_t sz=parenthesepipe.size();
	std::string species_name, branch_weight;
	for(unsigned i=0;i<sz;i++)
	{
		if(parenthesepipe[i]=='(')
		{
			if(parenthesepipe[i+1]==',')
			{
				species_name=wordpipe.front();
				wordpipe.pop_front();
				branch_weight=wordpipe.front();
				wordpipe.pop_front();
				Node node=g.addNode();
				label2node[species_name]=node;
				
				constructor.push(species_name);
				constructor.push(convert_num2str(g.id(node)));
				constructor.push(branch_weight);
			}
		}
		else if(parenthesepipe[i]==')')
		{
			if(parenthesepipe[i-1]==',')
			{
				species_name=wordpipe.front();
				wordpipe.pop_front();
				branch_weight=wordpipe.front();
				wordpipe.pop_front();
				Node node=g.addNode();
				label2node[species_name]=node;
				
				constructor.push(species_name);				
				constructor.push(convert_num2str(g.id(node)));
				constructor.push(branch_weight);
			}
			
			compressNodes(constructor);
			if(parenthesepipe[i+1]!=';')
			{
				branch_weight=wordpipe.front();
				wordpipe.pop_front();
				constructor.push(branch_weight);
			}
			else{
				root=g.nodeFromId(std::stoi(constructor.top()));
				constructor.pop();
				constructor.pop();
			}			
		}
	}
	for(EdgeIt ie(g);ie!=lemon::INVALID;++ie)
	{
		int nodeid1,nodeid2;
		Node node1,node2;
		std::string edgelabel;
		node1=g.u(ie);
		node2=g.v(ie);
		nodeid1=g.id(node1);
		nodeid2=g.id(node2);
		if(nodeid1<nodeid2)
		{
			edgelabel.append(convert_num2str(nodeid1));
			edgelabel.append(convert_num2str(nodeid2));
		}else
		{
			edgelabel.append(convert_num2str(nodeid2));
			edgelabel.append(convert_num2str(nodeid1));
		}
		label2edge[edgelabel]=ie;
	}
	return true;
}

template<typename GR, typename OP>
void
Tree<GR,OP>::compressNodes(std::stack<std::string>& constructor)
{
	std::string species1, species2;
	float weight1,weight2;
	Node node1, node2, node3;
	Edge e1,e2;
	weight2=std::stof(constructor.top());
	constructor.pop();
	node2=g.nodeFromId(std::stoi(constructor.top()));
	constructor.pop();
	species2=constructor.top();
	constructor.pop();

	weight1=std::stof(constructor.top());
	constructor.pop();
	node1=g.nodeFromId(std::stoi(constructor.top()));
	constructor.pop();
	species1=constructor.top();
	constructor.pop();

	node3=g.addNode();
	constructor.push("-");
	constructor.push(convert_num2str(g.id(node3)));
	e1=g.addEdge(node1,node3);
	e2=g.addEdge(node2,node3);
	branchmap.set(e1,weight1);
	branchmap.set(e2,weight2);
	
	return;
}

template<typename GR, typename OP>
bool
Tree<GR,OP>::computeDistEvolution(int dsize)
{
	for(EdgeIt ie(g);ie!=lemon::INVALID;++ie)
	{
		distEvolution+=scoremap[ie];
	}
	float sumDist=distEvolution.sumup();
	overallScore=dsize/sumDist;
	return true;
}

template<typename GR, typename OP>
bool
Tree<GR,OP>::generateSubTree(std::string filename,std::vector<std::string>& speciesfiles)
{
	readTree(filename);
	BfsType mybfs(g);
	typename BfsType::PredMap pm(g);
	typename BfsType::ReachedMap rm(g);
	typename BfsType::DistMap dm(g);
	mybfs.predMap(pm).reachedMap(rm).distMap(dm);//.run(root);
	Node interestnode;
	mybfs.run(root);
	for(unsigned i=0;i<speciesfiles.size();i++)
	{
		if(label2node.find(speciesfiles[i])!=label2node.end())
		{
			std::cout <<speciesfiles[i]<<": ";
			interestnode=label2node[speciesfiles[i]];
			if(!mybfs.reached(interestnode))continue;
			Node tnode=interestnode;
			typename Graph::Arc arc;
			//typename BfsType::Path ph=mybfs.path(interestnode);
			std::cout << g.id(tnode) << " -> ";
			while ((arc=pm[tnode]) != lemon::INVALID)
			{
				std::cout << " (" <<branchmap[arc] <<") ";
				tnode =mybfs.predNode(tnode);
				std::cout << g.id(tnode) << " -> ";
			}
			std::cout << "end.\n";
		}
		else
		{
			std::cerr << "Species " << speciesfiles[i] <<" can't be found!";
			return false;
		}
	}
	return true;
}

template<typename GR, typename OP>
void
Tree<GR,OP>::extractSubtree(std::string sentence, std::vector<std::string>& speciesfiles)
{
}
#endif
