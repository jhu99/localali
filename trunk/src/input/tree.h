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
	WeightEdgeMap branchmap;

	Tree();
	~Tree(){};
	bool readTree(std::string);
	bool constructTree(std::string);
};

template<typename GR, typename OP>
Tree<GR,OP>::Tree()
:g()
,node2label(g)
,branchmap(g)
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
	std::queue<std::string> wordpipe;
	std::queue<char> parenthesepipe;
	char * pSentence = new char [sentence.length()+1];
	strcpy (pSentence, sentence.c_str());
	char * pch = strtok(pSentence,"(),:;");
	std::cout << sentence << std::endl;
	while(pch!=0)
	{
		std::string word(pch);
		if(g_verbosity>VERBOSE_ESSENTIAL)
		std::cout <<word<<std::endl;
		wordpipe.push(word);
		pch=strtok(NULL,"(),:;");
	}	
	delete pSentence;
	std::size_t found = sentence.find_first_of("()");
	while(found!=std::string::npos)
	{
		parenthesepipe.push(sentence[found]);
		if(g_verbosity>VERBOSE_ESSENTIAL)
		std::cout << sentence[found];
		found=sentence.find_first_of("()",found+1);
	}
	if(g_verbosity>VERBOSE_ESSENTIAL)
	std::cout << std::endl;
	//char element=parenthesepipe.front();
	while(!parenthesepipe.empty())
	{
		std::cout << parenthesepipe.front()<<" ";
		parenthesepipe.pop();
	}
	std::cout << std::endl;
	return true;
}

