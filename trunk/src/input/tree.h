/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/tree.h
Description: Searching high-scoring subnetworks.
**/

#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <queue>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include "verbose.h"
#include "string.h"

template<typename GR, typename OP>
class Tree
{
	private:
	typedef GR Graph;
	typedef OP Option;
	TEMPLATE_GRAPH_TYPEDEFS(Graph);
	public:
	
	/// Labels of the nodes.
	typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
	/// Mapping from labels to original nodes.
	typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;
	/// Weights on original edges.
	typedef typename Graph::template EdgeMap<int> WeightEdgeMap;

	/// Hold the topology of this tree.
	Graph g;
	/// Root of the evolutionary tree.
	Node root;
	/// The labels for the tree nodes.
	OrigLabelNodeMap node2label;
	/// The lengths for the branches in this tree.
	WeightEdgeMap length;

	Tree();
	~Tree(){};
	bool readTree(std::string);
	bool constructTree(std::string);
};

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
	std::queue<std::string> fifopipe;
	char* pSentence=sentence.c_str();
	char* pch;
	pch = strtok(pSentence,"(),:;");
	while(pch!=NULL)
	{
		fifopipe.push(*pch);
		std::cout << *pch << std::endl;
	}
	return true;
}

