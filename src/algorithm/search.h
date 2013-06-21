/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/subnet.h
Description: Searching high-scoring subnetworks.
**/

#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <chrono>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/connectivity.h>
#include "verbose.h"
//#include "input/tree.h"

template<typename NP, typename SN, typename LG, typename OP>
class Search
{
private:
	typedef NP NetworkPool;
	typedef SN SubNet;
	typedef LG LayerGraph;
	typedef OP Option;
	typedef typename LayerGraph::Graph Graph;
	typedef typename SubNet::K_Spine K_Spine;
public:
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Labels of the nodes.
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
  /// Mapping from labels to original nodes.
  typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;

    unsigned _numSpecies;
    int _seedSize;
    int _seedTries;
    int _numSamples;
    int _numExtention;
    int _numConnected;
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution;
    typename std::vector<SubNet*> refinedSeeds;
	Search(Option&);
	~Search(){};
	void test(LayerGraph&,NetworkPool&);
	void searchSeeds(LayerGraph&,NetworkPool&);
	void searchCandidates(std::vector<Node>&,K_Spine&,LayerGraph&,NetworkPool&);
	bool sampleSeed(SubNet&,LayerGraph&,NetworkPool&,std::discrete_distribution<int>&);
	int sampleStringElement(int);
	bool sampleKSpine(Node&,K_Spine&,LayerGraph&,NetworkPool&);
	bool expandspine(K_Spine,K_Spine&,std::vector<Node>,Node,LayerGraph&,NetworkPool&,unsigned);
	void verifyspine(LayerGraph&,NetworkPool&);
	bool checkConnection(SubNet&,LayerGraph&,NetworkPool&);
	void expandRefinedSeeds(SubNet&,LayerGraph&,NetworkPool&);
	bool heuristicSearch(SubNet&,std::vector<Node>&,LayerGraph&,NetworkPool&);
	
};

template<typename NP, typename SN, typename LG, typename OP>
Search<NP,SN,LG,OP>::Search(Option& myoption)
:generator(std::chrono::system_clock::now().time_since_epoch().count())
,distribution(0,10000)
,refinedSeeds()
{
	_numSpecies=myoption.numspecies;
	_seedSize=myoption.seedsize;
	_seedTries=myoption.seedtries;
	_numSamples=myoption.numsamples;
	_numConnected=myoption.numconnected;
	_numExtention=myoption.numextention;
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::test(LayerGraph& layergraph,NetworkPool& networks)
{
	searchSeeds(layergraph,networks);
	int numAll=0;
	for(unsigned i=0;i<refinedSeeds.size();i++)
	{
		int numQualified=0;
		for(int j=0; j<_seedTries;j++)
		{
			SubNet subnet=*refinedSeeds[i];
			expandRefinedSeeds(subnet,layergraph,networks);
			subnet.induceSubgraphs(networks,layergraph);
			if(checkConnection(subnet,layergraph,networks))
			{
				subnet.output(layergraph);
				numQualified++;
			}
		}
		numAll+=numQualified;
		if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
			std::cout <<numQualified<<" subnets have been found for seed "<< i <<"."<<std::endl;
	}
	std::cout <<numAll<<" subnets have been found for all seeds."<<std::endl;
}
template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::searchCandidates(std::vector<Node>& candidates,
									  K_Spine& pspine,
									  LayerGraph& layergraph,
									  NetworkPool& networks)
{
	std::vector<std::string> innerproteins;
	std::vector<std::string> neighborproteins;
	for(unsigned j=0;j<_numSpecies;j++)
	{
		innerproteins.push_back(layergraph.node2label[pspine.data[j]]);
	}
	networks.getNeighbors(innerproteins,neighborproteins);
	for(unsigned j=0;j<neighborproteins.size();++j)
	{
		std::string protein=neighborproteins[j];
		Node neighbor=layergraph.label2node[protein];
		int neighborid=layergraph.graph.id(neighbor);
		if(layergraph.validnodemap.find(neighborid)==layergraph.validnodemap.end())continue;
		candidates.push_back(neighbor);
	}
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::expandRefinedSeeds(SubNet& subnet,LayerGraph& layergraph,NetworkPool& networks)
{
	int num=0;
	std::vector<Node> candidates;
	// searching neighbors of subnet.
	for(unsigned i=0;i<subnet.net_spines.size();i++)
		searchCandidates(candidates, subnet.net_spines[i], layergraph, networks);
	while(num++<_numExtention)
	{
		heuristicSearch(subnet, candidates, layergraph,networks);
	}	
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::heuristicSearch(SubNet& subnet,std::vector<Node>& candidates,LayerGraph& layergraph,NetworkPool& networks)
{
	std::vector<Node> optimalcandidates;
	std::unordered_map<int,int> candidatesmap;
	int maxnum=1;
	for(unsigned i=0;i<candidates.size();i++)
	{
		int nodeid=layergraph.graph.id(candidates[i]);
		if(candidatesmap.find(nodeid)!=candidatesmap.end())
		{
			candidatesmap[nodeid]++;
			if(maxnum<candidatesmap[nodeid])
				maxnum=candidatesmap[nodeid];
		}
		else
			candidatesmap[nodeid]=1;
	}
	for(std::unordered_map<int,int>::iterator it=candidatesmap.begin();it!=candidatesmap.end();it++)
	{
		if(it->second==maxnum)
			optimalcandidates.push_back(layergraph.graph.nodeFromId(it->first));
	}
	///
	int dice_roll = distribution(generator)%optimalcandidates.size();
	Node node=optimalcandidates[dice_roll];
	K_Spine pspine;
    if(!sampleKSpine(node,pspine,layergraph,networks))
    {
		std::cerr <<"Invalid sample node!"<<std::endl;
		return false;
	}
	subnet.net_spines.push_back(pspine);
	searchCandidates(candidates,pspine,layergraph,networks);
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::searchSeeds(LayerGraph& layergraph,NetworkPool& networks)
{
	verifyspine(layergraph,networks);
	std::discrete_distribution<int> discrete(layergraph.density.begin(),layergraph.density.end());
	int num=0;
	while(num++<_numSamples)
	{
		SubNet* mysubnet= new SubNet(_numSpecies,_seedSize);
		if(!sampleSeed(*mysubnet,layergraph,networks,discrete))
		{
			delete mysubnet;
			if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
			std::cerr << "Failed to sample a refined seed!" << std::endl;
		}else
		{
			refinedSeeds.push_back(mysubnet);
		}
    }
    std::cerr << "There are a total of "<< refinedSeeds.size() <<" refined seeds."<< std::endl;
}

/// Randomly sample a d subnet from G_h
template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::sampleSeed(SubNet& subnet, LayerGraph& layergraph, NetworkPool& networks,std::discrete_distribution<int>& discrete)
{
	/// Randomly sample a k-spine containing node.
	int dice_roll = discrete(generator);
	Node firstnode = layergraph.validnodes[dice_roll];
	std::vector<Node> candidates;
	for(int i=0;i<_seedSize;++i)
	{
		Node node;
		if(0==i) node = firstnode;
		else if(candidates.size()==0)
		{
			dice_roll = discrete(generator);
			node=layergraph.validnodes[dice_roll];
		}
		else
		{
			//std::vector<Node> optimalcandidates;
			//std::unordered_map<int,int> candidatesmap;
			//int maxnum=1;
			//for(unsigned i=0;i<candidates.size();i++)
			//{
				//int nodeid=layergraph.graph.id(candidates[i]);
				//if(candidatesmap.find(nodeid)!=candidatesmap.end())
				//{
					//candidatesmap[nodeid]++;
					//if(maxnum<candidatesmap[nodeid])
						//maxnum=candidatesmap[nodeid];
				//}
				//else
					//candidatesmap[nodeid]=1;
			//}
			//for(std::unordered_map<int,int>::iterator it=candidatesmap.begin();it!=candidatesmap.end();it++)
			//{
				//if(it->second==maxnum)
					//optimalcandidates.push_back(layergraph.graph.nodeFromId(it->first));
			//}
			dice_roll = distribution(generator)%candidates.size();
			node=candidates[dice_roll];
		}
		K_Spine pspine;
	    if(!sampleKSpine(node,pspine,layergraph,networks))
	    {
			std::cerr <<"Invalid sample node!"<<std::endl;
			return false;
		}

		/// Output this sample.
		if(g_verbosity >= VERBOSE_NON_ESSENTIAL)
		{
			for(unsigned j=0;j<_numSpecies;++j)
				std::cout << layergraph.node2label[pspine.data[j]]<<" ";
			std::cout << std::endl;
		}
		searchCandidates(candidates, pspine, layergraph, networks);
		//std::vector<std::string> innerproteins;
		//std::vector<std::string> neighborproteins;
		//for(unsigned j=0;j<_numSpecies;j++)
		//{
			//innerproteins.push_back(layergraph.node2label[pspine.data[j]]);
		//}
		//networks.getNeighbors(innerproteins,neighborproteins);
		//for(unsigned j=0;j<neighborproteins.size();++j)
		//{
			//std::string protein=neighborproteins[j];
			//Node neighbor=layergraph.label2node[protein];
			//int neighborid=layergraph.graph.id(neighbor);
			//if(layergraph.validnodemap.find(neighborid)==layergraph.validnodemap.end())continue;
			//candidates.push_back(neighbor);
		//}
	    subnet.net_spines.push_back(pspine);
	}
	subnet.induceSubgraphs(networks,layergraph);
	return checkConnection(subnet,layergraph,networks);
}

/// Collect these nodes which can conduct a successful sample of k-spine.
template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::verifyspine(LayerGraph& layergraph,NetworkPool& networks)
{
	for(NodeIt node(layergraph.graph); node!=lemon::INVALID; ++node)
	{
		K_Spine pspine;
		if(sampleKSpine(node,pspine,layergraph,networks))
		{
			layergraph.setConfiguration(node);
		}
	}
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::checkConnection(SubNet& subnet,LayerGraph& layergraph,NetworkPool& networks)
{
	int numConnected=0;
	for(unsigned i=0;i<_numSpecies;++i)
	{
		if(subnet.subgraphs[i]->g->edgeNum()<(subnet.subgraphs[i]->g->nodeNum()-1))continue;
		numConnected++;
	}
	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
	std::cout << "There are "<<numConnected << " out of " <<_numSpecies <<" subnetwork(s) are connected!" << std::endl;
	return numConnected>=_numConnected;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::sampleKSpine(Node& node,K_Spine& pspine,LayerGraph& layergraph,NetworkPool& networks)
{
	std::vector<Node> candidates;//
	K_Spine spine;
	unsigned k=_numSpecies;
	unsigned host=networks.getHost(layergraph.node2label[node]);
	spine.data[host]=node;
	spine.states[host]=true;
	if(g_verbosity==VERBOSE_DEBUG)
	std::cerr << layergraph.node2label[node] << std::endl;
	
	for(IncEdgeIt it(layergraph.graph,node);it!=lemon::INVALID;++it)
	{
		Node rnode=layergraph.graph.runningNode(it);
		host=networks.getHost(layergraph.node2label[rnode]);
		if(g_verbosity==VERBOSE_DEBUG)
		std::cerr << layergraph.node2label[rnode] << std::endl;
		if(spine.states[host])continue;
		if(node < rnode) candidates.push_back(rnode);
	}
	
	return expandspine(spine,pspine,candidates, node,layergraph,networks,k);
}

template<typename NP, typename SN, typename LG, typename OP>
int
Search<NP,SN,LG,OP>::sampleStringElement(int up)
{
    int dice_roll = distribution(generator)%(up);
    return dice_roll;
}
template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::expandspine(	K_Spine spine
								, K_Spine& pspine
								, std::vector<Node> candidates
								, Node node					
								, LayerGraph& layergraph
								, NetworkPool& networks
								, unsigned k)
{

	if(0==--k)
	{
		pspine=spine;
		return true;
	}

	while(!candidates.empty())
	{
		unsigned host;
		int size_Candidates=candidates.size();
		typename std::vector<Node>::iterator in=candidates.begin()+sampleStringElement(size_Candidates);
		Node w = *in;
		candidates.erase(in);
		std::vector<Node> secondCandidates;//V'extention
		std::vector<Node> exclCandidates;//ExclCandidates  V'UN(V')
		/// Assign V'UN(V') to exclCandidates.
		for(unsigned i=0;i<_numSpecies;i++)
		{
			if(!spine.states[i])continue;
			Node spinenode=spine.data[i];
			exclCandidates.push_back(spinenode);
	
			for(IncEdgeIt it(layergraph.graph,spinenode);it!=lemon::INVALID;++it)
			{
				Node rnode=layergraph.graph.runningNode(it);
				host=networks.getHost(layergraph.node2label[rnode]);
				if(spine.states[host])continue;/// If yes, it's not a neighbor.			
				if( std::find(exclCandidates.begin(), exclCandidates.end(), rnode)
				    != exclCandidates.end())continue;
				exclCandidates.push_back(rnode);			
			}
		}
		/// Assign new value to secondCandidates.
		for(IncEdgeIt it(layergraph.graph,w);it!=lemon::INVALID;++it)
		{
			Node rnode=layergraph.graph.runningNode(it);
			host=networks.getHost(layergraph.node2label[rnode]);
			if(spine.states[host])continue;
			if( std::find(exclCandidates.begin(), exclCandidates.end(), rnode)
				    != exclCandidates.end())continue;
			if (node < rnode) secondCandidates.push_back(rnode);				
		}
		K_Spine newspine(spine);
		host=networks.getHost(layergraph.node2label[w]);
		if(g_verbosity==VERBOSE_DEBUG)
		std::cerr << layergraph.node2label[w] << std::endl;
		newspine.data[host]=w;
		newspine.states[host]=true;
		for(unsigned i=0;i<candidates.size();++i)
		{
			Node rnode=candidates[i];
			host=networks.getHost(layergraph.node2label[rnode]);
			if(newspine.states[host])continue;
			secondCandidates.push_back(rnode);
		}
		if(expandspine(newspine, pspine, secondCandidates, node, layergraph, networks, k))return true;
	}
	return false;
}
