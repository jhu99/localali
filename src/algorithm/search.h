/**
Author: Jialu Hu
Date: Jun. 11, 2013
File name: algorithm/search.h
Description: Searching high-scoring subnetworks.
**/
#pragma once
#ifndef SEARCH_H_
#define SEARCH_H_

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
#include "input/tree.h"
#include "algorithm/phylogeny.h"
#include "algorithm/simulatedannealing.h"
#include "function.h"
#include <omp.h>
#include "algorithm/score.h"
#include <unordered_map>
//#include <assert.h>

template<typename NP, typename SN, typename LG, typename OP>
class Search
{
public:
	typedef NP NetworkPool;
	typedef SN SubNet;
	typedef LG LayerGraph;
	typedef OP Option;
	typedef typename LayerGraph::Graph Graph;
	TEMPLATE_GRAPH_TYPEDEFS(Graph);
	typedef typename SubNet::K_Spine K_Spine;
	typedef typename SubNet::GraphData GraphData;
	typedef Tree<Graph, Option> MyTree;
	typedef Phylogeny<SubNet,MyTree> MyPhylogeny;
	typedef typename MyPhylogeny::DeltaStructure DeltaStructure;
	typedef SimulatedAnnealing<MyPhylogeny,Option> MySimulatedAnnealing;
	typedef typename SubNet::InvOrigLabelNodeMap::iterator ItInvOrigLabelNodeMap;
	typedef typename MyTree::MatchingNodeMap MatchingNodeMap;
	typedef typename MatchingNodeMap::iterator MatchingNodeMapIt;
	typedef typename Graph::template EdgeMap<Score> ScoreEdgeMap;
	typedef typename Graph::template EdgeMap<MatchingNodeMap*> MatchingEdgeMap;
  
	
  /// Labels of the nodes.
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
  /// Mapping from labels to original nodes.
  typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;
	typedef std::list<std::vector<std::string> > SpineList;

    unsigned _numSpecies;
    unsigned _seedSize;
    int _seedTries;
    int _numSamples;
	int _minExt;
	int _maxExt;
    int _numExtension;
    int _numConnected;
    int _numthreads;
    double _score_threshold;
	std::string _resultfolder;
	std::string _treefile;
	std::vector<std::string> _speciesfiles;
    
    std::default_random_engine generator;
    std::uniform_int_distribution<int> distribution;
    typename std::vector<SubNet*> refinedSeeds;

    typedef struct _PrivateVariable
		{
		SubNet subnet;
		SubNet *mysubnet;
		int k;
		int j;
		unsigned ei;
		unsigned ej;
		unsigned sj;
		int host;
		int num;
		int numExtension;
		int numConnected;
		int nodeid;
		int neighborid;
		int dice_roll;
		int sk;
		std::string protein,element,protein1,protein2;
		std::vector<Node> m_candidates;
		std::vector<Node> candidates;
		std::vector<Node> secondCandidates;//V'extention
		std::vector<Node> exclCandidates;//ExclCandidates  V'UN(V')
		std::unordered_map<int,bool> usedValidNodes;
		std::vector<std::string> innerproteins;
	  std::vector<std::string> neighborproteins;
		std::vector<std::string> nodeset,xspine;
	  std::uniform_int_distribution<int> discrete;
		GraphData *graphdata;
		Node neighbor,firstnode,node,rnode,wnode,spinenode,node1,node2;
		K_Spine pspine, spine, newspine;
		IncEdgeIt inc;
		std::default_random_engine generator;
		typename std::vector<Node>::iterator it;
		SpineList::iterator sit;
		_PrivateVariable(unsigned dm)
		:subnet(),k(0),j(0),ei(0),ej(0),host(0),num(0),numExtension(0),numConnected(0),m_candidates(),candidates(),secondCandidates(),exclCandidates(),usedValidNodes()
		,innerproteins(),neighborproteins(),discrete(0,dm),generator(std::chrono::system_clock::now().time_since_epoch().count())
		{
		}
	} PrivateVariable;
		typedef struct _PrivateVariablePlus
		{
			float step,beta,sampledata,sumDist,overallScore;
			std::default_random_engine generator;
			std::uniform_real_distribution<float> distribution;
			unsigned seed,si,sj,sk,st,mykey,id1,id2,numElement,maxdegree,conNum,dice_roll,upperId,dice_1,dice_2;
			int degree1,degree2;
			SubNet *subnet;
			GraphData *graphdata, *sondata, *descedant, *ancestor;
			Node node,node1,node2,node3,node4,rnode,anode,dnode,firstnode;
			NodeIt nit1,nit2;
			MySimulatedAnnealing simulatedannealing;
			IncEdgeIt incE;
			EdgeIt ie,ie1;
			Edge myedge;
			std::pair<MatchingNodeMapIt, MatchingNodeMapIt> range,range1,range2;
			bool isExist,finderFlag;
			std::ifstream *input;
			std::string line,edgelabel;
			std::string element,protein1,protein2,protein;
			std::vector<std::string> nodeset,xspine;
			std::unordered_map<int,bool> usedValidNodes;
			std::deque<std::string> wordpipe;
			std::deque<char> parenthesepipe;
			std::deque<Node> processnode;
			std::deque<IncEdgeIt> processinc;
			std::vector<std::string>::iterator it;
			std::vector<Node> candidates;
			std::unordered_map<std::string,bool> proteinmap;
			ItInvOrigLabelNodeMap cit;
			SpineList::iterator sit;
			MatchingNodeMapIt mit,mit1,mit2;
			SpineList incSpines;
			Score score, deltaScore, updatedScore, distEvolution;;
			ScoreEdgeMap *scoremap;
			MatchingNodeMap *matchingmap;
			MatchingEdgeMap *matchingedgemap; 
			float branchweight;
			DeltaStructure deltaData;
			MyPhylogeny phylogeny;
			std::ofstream *fout;
			_PrivateVariablePlus():
				generator(std::chrono::system_clock::now().time_since_epoch().count()),distribution(0.0,1.0),simulatedannealing(),deltaData()
			{
			}
			void clear()
			{
				step=0;beta=0;sampledata=0;sumDist=0;overallScore=0;
				seed=0;si=0;sj=0;sk=0;st=0;mykey=0;id1=0;id2=0;numElement=0;maxdegree=0;conNum=0;dice_roll=0;upperId=0;dice_1=0;dice_2=0;
				candidates.clear();usedValidNodes.clear();
			}
		}PrivateVariablePlus;
	Search(Option&);
	~Search(){};
	void run(LayerGraph&,NetworkPool&);
	void setExtension(PrivateVariable&);
	void searchSeeds(LayerGraph&,NetworkPool&);
	void searchSeedsParallel(LayerGraph&, NetworkPool&);
	void searchCandidates(std::unordered_map<int,bool>&,std::vector<Node>&,K_Spine&,LayerGraph&,NetworkPool&);
	void searchCandidatesParallel(PrivateVariable&,LayerGraph&,NetworkPool&);
	bool sampleSeed(SubNet*,LayerGraph&,NetworkPool&,std::uniform_int_distribution<int>&);
	bool sampleSeedParallel(PrivateVariable&,LayerGraph&,NetworkPool&);
	int sampleStringElement(int);
	bool sampleKSpineParallel(PrivateVariable&,LayerGraph&,NetworkPool&);
	bool sampleKSpine(Node&,K_Spine&,LayerGraph&,NetworkPool&);
	bool expandspineParallel(PrivateVariable&, K_Spine spine, std::vector<Node> candidates, LayerGraph&,NetworkPool&);
	bool expandspine(K_Spine,K_Spine&,std::vector<Node>,Node,LayerGraph&,NetworkPool&,unsigned);
	void verifyspine(LayerGraph&,NetworkPool&);
	void verifyspineParallel(LayerGraph&,NetworkPool&);
	bool checkConnection(SubNet*,LayerGraph&,NetworkPool&);
	bool checkConnectionParallel(PrivateVariable&,LayerGraph&,NetworkPool&);
	void expandRefinedSeeds(PrivateVariable&,LayerGraph&,NetworkPool&);
	bool heuristicSearch(PrivateVariable&,LayerGraph&,NetworkPool&);
	bool induceSubgraphs(PrivateVariablePlus&,LayerGraph&,NetworkPool&);
	bool induceSubgraphsPrivateVar(PrivateVariable&,LayerGraph&,NetworkPool&);
	bool initialPhylogy(PrivateVariablePlus&,LayerGraph&, const MyTree&);
	bool initialExternalNodes(PrivateVariablePlus&, const MyTree&);
	bool existNode(PrivateVariablePlus&);
	bool initialBranchWeight(PrivateVariablePlus&,const MyTree&);
	bool computeBranchWeight(PrivateVariablePlus&,const MyTree&);
	bool computeScore(PrivateVariablePlus&,const MyTree&);
	void computeDist(PrivateVariablePlus&,const MyTree&);
	bool clearScore(PrivateVariablePlus&);
	bool clearStructure(PrivateVariablePlus&,const MyTree&);
	void simulatedAnnealingMethod(PrivateVariablePlus&,const MyTree&);
	void sumupScore(PrivateVariablePlus&,Score&);
	GraphData* constructInternalNodes(Node,PrivateVariablePlus&,LayerGraph&,const MyTree&);
	bool interfere(PrivateVariablePlus&,const MyTree&);
	void deleteEdge(PrivateVariablePlus&);
	void addEdge(PrivateVariablePlus&);
	void formEdgeLabel(PrivateVariablePlus&);
	void outsubgraphs(LayerGraph&,PrivateVariablePlus&,std::ofstream& fout);
	void output(LayerGraph&,PrivateVariablePlus&,std::ofstream& fout);
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
	_numExtension=myoption.minext;
	_minExt=myoption.minext;
	_maxExt=myoption.maxext;
	_resultfolder=myoption.resultfolder;
	_numthreads=myoption.numthreads;
	_treefile=myoption.treefile;
	_speciesfiles=myoption.speciesfiles;
	_score_threshold=myoption.score_threshold;
}

template<typename NP, typename SN, typename LG, typename OP>
bool Search<NP,SN,LG,OP>::clearScore(PrivateVariablePlus& myPrivateVariablePlus)
{
	myPrivateVariablePlus.score.fscore.fill(0.0);
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool Search<NP,SN,LG,OP>::clearStructure(PrivateVariablePlus& myPrivateVariablePlus,const MyTree& localtree)
{
	myPrivateVariablePlus.clear();
	for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<myPrivateVariablePlus.phylogeny.node2graph.size();myPrivateVariablePlus.si++)
	{
		delete myPrivateVariablePlus.phylogeny.node2graph[myPrivateVariablePlus.si];
	}
	for(myPrivateVariablePlus.ie=EdgeIt(localtree.g);myPrivateVariablePlus.ie!=lemon::INVALID;++myPrivateVariablePlus.ie)
	{
		delete (*myPrivateVariablePlus.matchingedgemap)[myPrivateVariablePlus.ie];
	}
	delete myPrivateVariablePlus.scoremap;
	delete myPrivateVariablePlus.matchingedgemap;
	delete myPrivateVariablePlus.subnet;
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::initialExternalNodes(PrivateVariablePlus& myPrivateVariablePlus,const MyTree& localtree)
{
	for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<_speciesfiles.size();myPrivateVariablePlus.si++)
	{
		myPrivateVariablePlus.element=_speciesfiles[myPrivateVariablePlus.si];
		myPrivateVariablePlus.node=localtree.label2node.find(myPrivateVariablePlus.element)->second;
		myPrivateVariablePlus.phylogeny.node2graph[localtree.g.id(myPrivateVariablePlus.node)]=myPrivateVariablePlus.subnet->subgraphs[myPrivateVariablePlus.si];
		myPrivateVariablePlus.phylogeny.externalNode.push_back(myPrivateVariablePlus.node);
	}
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::existNode(PrivateVariablePlus& myPrivateVariablePlus)
{
	for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<myPrivateVariablePlus.xspine.size();myPrivateVariablePlus.si++)
	{
		if(myPrivateVariablePlus.graphdata->label2node->find(myPrivateVariablePlus.xspine[myPrivateVariablePlus.si])!=myPrivateVariablePlus.graphdata->label2node->end())
			return true;
	}
	return false;
}

template<typename NP, typename SN, typename LG, typename OP>
typename Search<NP,SN,LG,OP>::GraphData*
Search<NP,SN,LG,OP>::constructInternalNodes(Node ancestor, PrivateVariablePlus& myPrivateVariablePlus,LayerGraph& layergraph,const MyTree& localtree)
{
	myPrivateVariablePlus.graphdata=new GraphData();
	myPrivateVariablePlus.phylogeny.node2graph[localtree.g.id(ancestor)]=myPrivateVariablePlus.graphdata;
	for(myPrivateVariablePlus.incE=IncEdgeIt(localtree.g,ancestor);myPrivateVariablePlus.incE!=lemon::INVALID;++myPrivateVariablePlus.incE)
	{
		myPrivateVariablePlus.rnode=localtree.g.runningNode(myPrivateVariablePlus.incE);
		myPrivateVariablePlus.processinc.push_back(myPrivateVariablePlus.incE);
		if(ancestor<myPrivateVariablePlus.rnode)
		{
			myPrivateVariablePlus.processinc.pop_back();
			continue;
		}
		//sonnodes.push_back(myPrivateVariablePlus.rnode);
		if(myPrivateVariablePlus.phylogeny.node2graph.find(localtree.g.id(myPrivateVariablePlus.rnode))==myPrivateVariablePlus.phylogeny.node2graph.end())
		{
			myPrivateVariablePlus.processnode.push_back(myPrivateVariablePlus.rnode);
			myPrivateVariablePlus.sondata=constructInternalNodes(myPrivateVariablePlus.rnode,myPrivateVariablePlus,layergraph,localtree);
			myPrivateVariablePlus.rnode=myPrivateVariablePlus.processnode.back();
			myPrivateVariablePlus.processnode.pop_back();
			myPrivateVariablePlus.phylogeny.internalNode.push_back(myPrivateVariablePlus.rnode);
		}else{
			myPrivateVariablePlus.sondata=myPrivateVariablePlus.phylogeny.node2graph[localtree.g.id(myPrivateVariablePlus.rnode)];
		}
		myPrivateVariablePlus.graphdata=myPrivateVariablePlus.phylogeny.node2graph[localtree.g.id(ancestor)];
		myPrivateVariablePlus.incE=myPrivateVariablePlus.processinc.back();
		myPrivateVariablePlus.processinc.pop_back();
		for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<myPrivateVariablePlus.sondata->offsprings.size();++myPrivateVariablePlus.si)
			myPrivateVariablePlus.graphdata->offsprings.push_back(myPrivateVariablePlus.sondata->offsprings[myPrivateVariablePlus.si]);
	}
		/// construct internal nodes.
	myPrivateVariablePlus.incSpines.clear();
	for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<myPrivateVariablePlus.subnet->net_spines.size();++myPrivateVariablePlus.si)
	{
		myPrivateVariablePlus.xspine.clear();
		for(myPrivateVariablePlus.sj=0;myPrivateVariablePlus.sj<myPrivateVariablePlus.graphdata->offsprings.size();++myPrivateVariablePlus.sj)
		{
			myPrivateVariablePlus.sk=myPrivateVariablePlus.graphdata->offsprings[myPrivateVariablePlus.sj];
			myPrivateVariablePlus.node=myPrivateVariablePlus.subnet->net_spines[myPrivateVariablePlus.si].data[myPrivateVariablePlus.sk];
			myPrivateVariablePlus.xspine.push_back(layergraph.node2label[myPrivateVariablePlus.node]);
		}
		myPrivateVariablePlus.incSpines.push_back(myPrivateVariablePlus.xspine);
	}


	while(!myPrivateVariablePlus.incSpines.empty())
	{
		myPrivateVariablePlus.node = myPrivateVariablePlus.graphdata->g->addNode();
		myPrivateVariablePlus.xspine=myPrivateVariablePlus.incSpines.front();
		myPrivateVariablePlus.incSpines.pop_front();
		myPrivateVariablePlus.graphdata->nodeNum++;
		myPrivateVariablePlus.graphdata->node2degree->set(myPrivateVariablePlus.node,0);
		for(myPrivateVariablePlus.it=myPrivateVariablePlus.xspine.begin();myPrivateVariablePlus.it!=myPrivateVariablePlus.xspine.end();++myPrivateVariablePlus.it)
		{
			(*myPrivateVariablePlus.graphdata->label2node)[*myPrivateVariablePlus.it]=myPrivateVariablePlus.node;
			///There are many labels for each internal node.
		}
		myPrivateVariablePlus.sit=myPrivateVariablePlus.incSpines.begin();
		while(myPrivateVariablePlus.sit!=myPrivateVariablePlus.incSpines.end())
		{
				myPrivateVariablePlus.xspine=*myPrivateVariablePlus.sit;
				if(!existNode(myPrivateVariablePlus))
				{
					++myPrivateVariablePlus.sit;
					continue;
				}
				for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<myPrivateVariablePlus.xspine.size();myPrivateVariablePlus.si++)
				{
					(*myPrivateVariablePlus.graphdata->label2node)[myPrivateVariablePlus.xspine[myPrivateVariablePlus.si]]=myPrivateVariablePlus.node;
				}
				myPrivateVariablePlus.incSpines.erase(myPrivateVariablePlus.sit);
				myPrivateVariablePlus.sit=myPrivateVariablePlus.incSpines.begin();
		}
	}
	return myPrivateVariablePlus.graphdata;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::computeBranchWeight(PrivateVariablePlus& myPrivateVariablePlus,const MyTree& localtree)
{
	myPrivateVariablePlus.matchingmap = new MatchingNodeMap();
	myPrivateVariablePlus.node1=localtree.g.u(myPrivateVariablePlus.ie);
	myPrivateVariablePlus.id1=localtree.g.id(myPrivateVariablePlus.node1);
	myPrivateVariablePlus.node2=localtree.g.v(myPrivateVariablePlus.ie);
	myPrivateVariablePlus.id2=localtree.g.id(myPrivateVariablePlus.node2);
	if(myPrivateVariablePlus.node2<myPrivateVariablePlus.node1)
	{
		myPrivateVariablePlus.descedant=myPrivateVariablePlus.phylogeny.node2graph[myPrivateVariablePlus.id2];
		myPrivateVariablePlus.ancestor=myPrivateVariablePlus.phylogeny.node2graph[myPrivateVariablePlus.id1];
	}
	else
	{
		myPrivateVariablePlus.descedant=myPrivateVariablePlus.phylogeny.node2graph[myPrivateVariablePlus.id1];
		myPrivateVariablePlus.ancestor=myPrivateVariablePlus.phylogeny.node2graph[myPrivateVariablePlus.id2];
	}
	
	for( myPrivateVariablePlus.cit=myPrivateVariablePlus.descedant->label2node->begin();myPrivateVariablePlus.cit!=myPrivateVariablePlus.descedant->label2node->end();++myPrivateVariablePlus.cit)
	{
	    myPrivateVariablePlus.protein=myPrivateVariablePlus.cit->first;
			myPrivateVariablePlus.dnode=myPrivateVariablePlus.cit->second;
	  	myPrivateVariablePlus.anode=(*myPrivateVariablePlus.ancestor->label2node)[myPrivateVariablePlus.protein];
	  	myPrivateVariablePlus.mykey=myPrivateVariablePlus.ancestor->g->id(myPrivateVariablePlus.anode);
	  	myPrivateVariablePlus.range=myPrivateVariablePlus.matchingmap->equal_range(myPrivateVariablePlus.mykey);
	  	myPrivateVariablePlus.isExist=false;
			for(myPrivateVariablePlus.mit=myPrivateVariablePlus.range.first;myPrivateVariablePlus.mit!=myPrivateVariablePlus.range.second;++myPrivateVariablePlus.mit)
			{
				if(myPrivateVariablePlus.mit->second==myPrivateVariablePlus.dnode){myPrivateVariablePlus.isExist=true;break;}
			}
			if(myPrivateVariablePlus.isExist)continue;
		  myPrivateVariablePlus.matchingmap->insert(std::make_pair(myPrivateVariablePlus.mykey,myPrivateVariablePlus.dnode));// A general error may look like: std::make_pair<int, Node>(mykey,myPrivateVariablePlus.dnode)
	}
	myPrivateVariablePlus.matchingedgemap->set(myPrivateVariablePlus.ie,myPrivateVariablePlus.matchingmap);//[myPrivateVariablePlus.myedge]=myPrivateVariablePlus.matchingmap;
	computeScore(myPrivateVariablePlus,localtree);
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::computeScore(PrivateVariablePlus& myPrivateVariablePlus,const MyTree& localtree)
{
	myPrivateVariablePlus.maxdegree=myPrivateVariablePlus.ancestor->maxDegree;
	myPrivateVariablePlus.branchweight=localtree.branchmap[myPrivateVariablePlus.ie];
	if(0==myPrivateVariablePlus.maxdegree)myPrivateVariablePlus.maxdegree=1;
	myPrivateVariablePlus.score.fscore.fill(0.0);
	for(myPrivateVariablePlus.nit1=NodeIt(*myPrivateVariablePlus.ancestor->g);myPrivateVariablePlus.nit1!=lemon::INVALID;++myPrivateVariablePlus.nit1)
	{
		myPrivateVariablePlus.mykey=myPrivateVariablePlus.ancestor->g->id(myPrivateVariablePlus.nit1);
		myPrivateVariablePlus.numElement=myPrivateVariablePlus.matchingmap->count(myPrivateVariablePlus.mykey);
		if(myPrivateVariablePlus.numElement==1)
		/// Protein mutation.
		{
			myPrivateVariablePlus.score.fscore[0]+=(1-static_cast<float>((*myPrivateVariablePlus.ancestor->node2degree)[myPrivateVariablePlus.nit1])/myPrivateVariablePlus.maxdegree)*myPrivateVariablePlus.branchweight;
		}else
		/// Protein duplication.
		{
			myPrivateVariablePlus.score.fscore[1]+=(1-static_cast<float>((*myPrivateVariablePlus.ancestor->node2degree)[myPrivateVariablePlus.nit1])/myPrivateVariablePlus.maxdegree)*myPrivateVariablePlus.branchweight;
		}
	}
	/// The number of conserved edges.
	myPrivateVariablePlus.conNum=0;
	for(myPrivateVariablePlus.ie1=EdgeIt(*myPrivateVariablePlus.ancestor->g);myPrivateVariablePlus.ie1!=lemon::INVALID;++myPrivateVariablePlus.ie1)
	{
		myPrivateVariablePlus.finderFlag=false;
		myPrivateVariablePlus.node1=myPrivateVariablePlus.ancestor->g->u(myPrivateVariablePlus.ie1);
		myPrivateVariablePlus.node2=myPrivateVariablePlus.ancestor->g->v(myPrivateVariablePlus.ie1);
		myPrivateVariablePlus.id1=myPrivateVariablePlus.ancestor->g->id(myPrivateVariablePlus.node1);
		myPrivateVariablePlus.id2=myPrivateVariablePlus.ancestor->g->id(myPrivateVariablePlus.node2);
		myPrivateVariablePlus.range1=myPrivateVariablePlus.matchingmap->equal_range(myPrivateVariablePlus.id1);
		myPrivateVariablePlus.range2=myPrivateVariablePlus.matchingmap->equal_range(myPrivateVariablePlus.id2);
		for(myPrivateVariablePlus.mit1=myPrivateVariablePlus.range1.first;myPrivateVariablePlus.mit1!=myPrivateVariablePlus.range1.second;++myPrivateVariablePlus.mit1)
		{
			myPrivateVariablePlus.node1=myPrivateVariablePlus.mit1->second;
			for(myPrivateVariablePlus.mit2=myPrivateVariablePlus.range2.first;myPrivateVariablePlus.mit2!=myPrivateVariablePlus.range2.second;++myPrivateVariablePlus.mit2)
			{
				myPrivateVariablePlus.node2=myPrivateVariablePlus.mit2->second;
				//myPrivateVariablePlus.edgelabel=myPrivateVariablePlus.descedant->formEdgeLabel(myPrivateVariablePlus.node3,myPrivateVariablePlus.node4);
				formEdgeLabel(myPrivateVariablePlus);
				if(myPrivateVariablePlus.descedant->label2edge->find(myPrivateVariablePlus.edgelabel)!=myPrivateVariablePlus.descedant->label2edge->end())
				{
					myPrivateVariablePlus.conNum++;
					myPrivateVariablePlus.finderFlag=true;
					break;
				}
			}
			if(myPrivateVariablePlus.finderFlag)break;
		}
	}
	// Interaction deletion.
	myPrivateVariablePlus.score.fscore[2]=localtree._beta*myPrivateVariablePlus.branchweight*(myPrivateVariablePlus.ancestor->edgeNum-myPrivateVariablePlus.conNum);
	// Interaction insertion.
	myPrivateVariablePlus.score.fscore[3]=localtree._beta*myPrivateVariablePlus.branchweight*(myPrivateVariablePlus.descedant->edgeNum-myPrivateVariablePlus.conNum);
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::initialBranchWeight(PrivateVariablePlus& myPrivateVariablePlus,const MyTree& localtree)
{
	for(myPrivateVariablePlus.ie=EdgeIt(localtree.g);myPrivateVariablePlus.ie!=lemon::INVALID;++myPrivateVariablePlus.ie)
	{
		myPrivateVariablePlus.score.fscore.fill(0.0);
		computeBranchWeight(myPrivateVariablePlus,localtree);
		(*myPrivateVariablePlus.scoremap)[myPrivateVariablePlus.ie]=myPrivateVariablePlus.score;
	}
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::initialPhylogy(PrivateVariablePlus& myPrivateVariablePlus,LayerGraph& layergraph,const MyTree& localtree)
{
	myPrivateVariablePlus.phylogeny.internalNode.clear();
	myPrivateVariablePlus.phylogeny.externalNode.clear();
	myPrivateVariablePlus.phylogeny.node2graph.clear();
	initialExternalNodes(myPrivateVariablePlus,localtree);
	constructInternalNodes(localtree.root, myPrivateVariablePlus,layergraph,localtree);//shrink to 3 parameters
	myPrivateVariablePlus.phylogeny.internalNode.push_back(localtree.root);
	initialBranchWeight(myPrivateVariablePlus,localtree);
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
void 
Search<NP,SN,LG,OP>::deleteEdge(PrivateVariablePlus& myPrivateVariablePlus)
{
	myPrivateVariablePlus.node1=myPrivateVariablePlus.graphdata->g->u(myPrivateVariablePlus.myedge);
	myPrivateVariablePlus.node2=myPrivateVariablePlus.graphdata->g->v(myPrivateVariablePlus.myedge);
	myPrivateVariablePlus.degree1=(*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node1]--;
	myPrivateVariablePlus.degree2=(*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node2]--;
	myPrivateVariablePlus.graphdata->edgeNum--;
	myPrivateVariablePlus.graphdata->g->erase(myPrivateVariablePlus.myedge);
	myPrivateVariablePlus.graphdata->label2edge->erase(myPrivateVariablePlus.edgelabel);
	if(myPrivateVariablePlus.degree1==myPrivateVariablePlus.graphdata->maxDegree || myPrivateVariablePlus.degree2==myPrivateVariablePlus.graphdata->maxDegree)
	{
		myPrivateVariablePlus.graphdata->maxDegree--;
		for(myPrivateVariablePlus.nit1=NodeIt(*myPrivateVariablePlus.graphdata->g);myPrivateVariablePlus.nit1!=lemon::INVALID;++myPrivateVariablePlus.nit1)
		{
			if((*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.nit1]>myPrivateVariablePlus.graphdata->maxDegree)myPrivateVariablePlus.graphdata->maxDegree++;
		}
	}
}	

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::formEdgeLabel(PrivateVariablePlus& myPrivateVariablePlus)
{
	myPrivateVariablePlus.edgelabel.clear();
	if(myPrivateVariablePlus.node2 < myPrivateVariablePlus.node1)
	{
		myPrivateVariablePlus.edgelabel.append(convert_num2str(myPrivateVariablePlus.graphdata->g->id(myPrivateVariablePlus.node2)));
		myPrivateVariablePlus.edgelabel.append(convert_num2str(myPrivateVariablePlus.graphdata->g->id(myPrivateVariablePlus.node1)));
	}else
	{
		myPrivateVariablePlus.edgelabel.append(convert_num2str(myPrivateVariablePlus.graphdata->g->id(myPrivateVariablePlus.node1)));
		myPrivateVariablePlus.edgelabel.append(convert_num2str(myPrivateVariablePlus.graphdata->g->id(myPrivateVariablePlus.node2)));
	}
}

template<typename NP, typename SN, typename LG, typename OP>
void 
Search<NP,SN,LG,OP>::addEdge(PrivateVariablePlus& myPrivateVariablePlus)
{
	myPrivateVariablePlus.degree1=++(*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node1];
	myPrivateVariablePlus.degree2=++(*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node2];
	myPrivateVariablePlus.graphdata->edgeNum++;
	myPrivateVariablePlus.myedge=myPrivateVariablePlus.graphdata->g->addEdge(myPrivateVariablePlus.node1,myPrivateVariablePlus.node2);
	if(myPrivateVariablePlus.degree1>myPrivateVariablePlus.graphdata->maxDegree || myPrivateVariablePlus.degree2>myPrivateVariablePlus.graphdata->maxDegree)myPrivateVariablePlus.graphdata->maxDegree++;
	formEdgeLabel(myPrivateVariablePlus);
	(*myPrivateVariablePlus.graphdata->label2edge)[myPrivateVariablePlus.edgelabel]=myPrivateVariablePlus.myedge;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::interfere(PrivateVariablePlus& myPrivateVariablePlus,const MyTree& localtree)
{
	myPrivateVariablePlus.dice_roll=myPrivateVariablePlus.phylogeny.distribution(generator)%myPrivateVariablePlus.phylogeny.internalNode.size();
	myPrivateVariablePlus.node = myPrivateVariablePlus.phylogeny.internalNode[myPrivateVariablePlus.dice_roll];
	myPrivateVariablePlus.graphdata=myPrivateVariablePlus.phylogeny.node2graph[localtree.g.id(myPrivateVariablePlus.node)];
	myPrivateVariablePlus.upperId=myPrivateVariablePlus.graphdata->nodeNum;
	if(myPrivateVariablePlus.upperId<2) return false;
	myPrivateVariablePlus.dice_1 = myPrivateVariablePlus.phylogeny.distribution(generator)%myPrivateVariablePlus.upperId;
	myPrivateVariablePlus.node1 = myPrivateVariablePlus.graphdata->g->nodeFromId(myPrivateVariablePlus.dice_1);
	myPrivateVariablePlus.dice_2 = myPrivateVariablePlus.phylogeny.distribution(generator)%myPrivateVariablePlus.upperId;
	while(myPrivateVariablePlus.dice_1==myPrivateVariablePlus.dice_2)
	{
		myPrivateVariablePlus.dice_2 = myPrivateVariablePlus.phylogeny.distribution(generator)%myPrivateVariablePlus.upperId;
	}
	myPrivateVariablePlus.node2 = myPrivateVariablePlus.graphdata->g->nodeFromId(myPrivateVariablePlus.dice_2);
	formEdgeLabel(myPrivateVariablePlus);
	myPrivateVariablePlus.deltaData.treenode=myPrivateVariablePlus.node;
	myPrivateVariablePlus.deltaData.nodeA=myPrivateVariablePlus.node1;
	myPrivateVariablePlus.deltaData.nodeB=myPrivateVariablePlus.node2;
	myPrivateVariablePlus.deltaData.edgelabel=myPrivateVariablePlus.edgelabel;
	if(myPrivateVariablePlus.graphdata->label2edge->find(myPrivateVariablePlus.edgelabel)!=myPrivateVariablePlus.graphdata->label2edge->end())// label2edge empty?
	{
		////delete edge
		myPrivateVariablePlus.myedge=myPrivateVariablePlus.graphdata->label2edge->find(myPrivateVariablePlus.edgelabel)->second;
		deleteEdge(myPrivateVariablePlus);
	}
	else{
		//// add edge
		addEdge(myPrivateVariablePlus);
	}

	myPrivateVariablePlus.sj=0;
	myPrivateVariablePlus.deltaScore.fscore.fill(0);
	//// Compute the gap between the two scores, but without changes of their original score attritubtion.
	for(myPrivateVariablePlus.incE=IncEdgeIt(localtree.g,myPrivateVariablePlus.node);myPrivateVariablePlus.incE!=lemon::INVALID;++myPrivateVariablePlus.incE,++myPrivateVariablePlus.sj)
	{
		myPrivateVariablePlus.rnode=localtree.g.runningNode(myPrivateVariablePlus.incE);
		myPrivateVariablePlus.descedant=myPrivateVariablePlus.phylogeny.node2graph[localtree.g.id(myPrivateVariablePlus.rnode)];
		if(myPrivateVariablePlus.rnode<myPrivateVariablePlus.node)
		{
			myPrivateVariablePlus.ancestor=myPrivateVariablePlus.graphdata;
		}else
		{
			myPrivateVariablePlus.ancestor=myPrivateVariablePlus.descedant;
			myPrivateVariablePlus.descedant=myPrivateVariablePlus.graphdata;
		}
		myPrivateVariablePlus.edgelabel.clear();
		if(myPrivateVariablePlus.node < myPrivateVariablePlus.rnode)
		{
			myPrivateVariablePlus.edgelabel.append(convert_num2str(localtree.g.id(myPrivateVariablePlus.node)));
			myPrivateVariablePlus.edgelabel.append(convert_num2str(localtree.g.id(myPrivateVariablePlus.rnode)));
		}else
		{
			myPrivateVariablePlus.edgelabel.append(convert_num2str(localtree.g.id(myPrivateVariablePlus.rnode)));
			myPrivateVariablePlus.edgelabel.append(convert_num2str(localtree.g.id(myPrivateVariablePlus.node)));
		}
		myPrivateVariablePlus.ie=localtree.label2edge.at(myPrivateVariablePlus.edgelabel);
		computeScore(myPrivateVariablePlus,localtree);
		
		myPrivateVariablePlus.deltaScore+=myPrivateVariablePlus.score;
		myPrivateVariablePlus.deltaScore-=(*myPrivateVariablePlus.scoremap)[myPrivateVariablePlus.incE];
		myPrivateVariablePlus.deltaData.updatedScores[myPrivateVariablePlus.sj]=myPrivateVariablePlus.score;
	}
	//myPrivateVariablePlus.deltaData.delta=myPrivateVariablePlus.deltaScore.sumup();
	sumupScore(myPrivateVariablePlus,myPrivateVariablePlus.deltaScore);
	myPrivateVariablePlus.deltaData.delta=myPrivateVariablePlus.sumDist;
	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
	std::cout << myPrivateVariablePlus.deltaData.delta << std::endl;
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::induceSubgraphsPrivateVar(PrivateVariable& myPrivateVariable,LayerGraph& layergraph,NetworkPool& networks)
{
	myPrivateVariable.mysubnet->subgraphs.clear();
	for(myPrivateVariable.ei=0;myPrivateVariable.ei<_numSpecies;++myPrivateVariable.ei)
	{
		myPrivateVariable.graphdata = new GraphData();
		myPrivateVariable.graphdata->offsprings.push_back(myPrivateVariable.ei);
		myPrivateVariable.nodeset.clear();
		for(myPrivateVariable.ej=0;myPrivateVariable.ej<myPrivateVariable.mysubnet->net_spines.size();++myPrivateVariable.ej)
		{
			myPrivateVariable.element=layergraph.node2label[myPrivateVariable.mysubnet->net_spines[myPrivateVariable.ej].data[myPrivateVariable.ei]];
			if(find(myPrivateVariable.nodeset.begin(),myPrivateVariable.nodeset.end(),myPrivateVariable.element)!=myPrivateVariable.nodeset.end())continue;
			myPrivateVariable.nodeset.push_back(myPrivateVariable.element);
			myPrivateVariable.node=myPrivateVariable.graphdata->g->addNode();
			myPrivateVariable.graphdata->nodeNum++;
			myPrivateVariable.graphdata->node2label->set(myPrivateVariable.node,myPrivateVariable.element);
			myPrivateVariable.graphdata->node2degree->set(myPrivateVariable.node,0);
			(*myPrivateVariable.graphdata->label2node)[myPrivateVariable.element]=myPrivateVariable.node;
			//assert((*networks.getGraph(i)->invIdNodeMap).find(element)!=(*networks.getGraph(i)->invIdNodeMap).end());
		}
		for(myPrivateVariable.ej=0;myPrivateVariable.ej<myPrivateVariable.nodeset.size();myPrivateVariable.ej++)
		{
			myPrivateVariable.protein1=myPrivateVariable.nodeset[myPrivateVariable.ej];
			myPrivateVariable.node1=(*myPrivateVariable.graphdata->label2node)[myPrivateVariable.protein1];
			for(myPrivateVariable.sj=myPrivateVariable.ej+1;myPrivateVariable.sj<myPrivateVariable.nodeset.size();myPrivateVariable.sj++)
			{
				myPrivateVariable.protein2=myPrivateVariable.nodeset[myPrivateVariable.sj];
				myPrivateVariable.element.clear();
			  myPrivateVariable.node2=(*myPrivateVariable.graphdata->label2node)[myPrivateVariable.protein2];
				if(myPrivateVariable.protein1.compare(myPrivateVariable.protein2)>0)
				{
					myPrivateVariable.element.append(myPrivateVariable.protein2);
					myPrivateVariable.element.append(myPrivateVariable.protein1);
				 }
				else
				{
				 myPrivateVariable.element.append(myPrivateVariable.protein1);
				 myPrivateVariable.element.append(myPrivateVariable.protein2);
			  }
				 
				if(networks.getGraph(myPrivateVariable.ei)->interactionmap.find(myPrivateVariable.element)==networks.getGraph(myPrivateVariable.ei)->interactionmap.end())continue;
				 myPrivateVariable.graphdata->g->addEdge(myPrivateVariable.node1,myPrivateVariable.node2);
				 (*(myPrivateVariable.graphdata->node2degree))[myPrivateVariable.node1]++;
				 (*(myPrivateVariable.graphdata->node2degree))[myPrivateVariable.node2]++;
				 if((*myPrivateVariable.graphdata->node2degree)[myPrivateVariable.node1]>myPrivateVariable.graphdata->maxDegree)
					 myPrivateVariable.graphdata->maxDegree=(*myPrivateVariable.graphdata->node2degree)[myPrivateVariable.node1];
				 if((*myPrivateVariable.graphdata->node2degree)[myPrivateVariable.node2]>myPrivateVariable.graphdata->maxDegree)
					 myPrivateVariable.graphdata->maxDegree=(*myPrivateVariable.graphdata->node2degree)[myPrivateVariable.node2];
				 myPrivateVariable.graphdata->edgeNum++;
			}
		}
		myPrivateVariable.mysubnet->subgraphs.push_back(myPrivateVariable.graphdata);
	}
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::induceSubgraphs(PrivateVariablePlus& myPrivateVariablePlus,LayerGraph& layergraph,NetworkPool& networks)
{
	myPrivateVariablePlus.subnet->subgraphs.clear();
	for(myPrivateVariablePlus.si=0;myPrivateVariablePlus.si<_numSpecies;++myPrivateVariablePlus.si)
	{
		myPrivateVariablePlus.graphdata = new GraphData();
		myPrivateVariablePlus.graphdata->offsprings.push_back(myPrivateVariablePlus.si);
		myPrivateVariablePlus.nodeset.clear();
		for(myPrivateVariablePlus.sj=0;myPrivateVariablePlus.sj<myPrivateVariablePlus.subnet->net_spines.size();++myPrivateVariablePlus.sj)
		{
			myPrivateVariablePlus.element=layergraph.node2label[myPrivateVariablePlus.subnet->net_spines[myPrivateVariablePlus.sj].data[myPrivateVariablePlus.si]];
			if(find(myPrivateVariablePlus.nodeset.begin(),myPrivateVariablePlus.nodeset.end(),myPrivateVariablePlus.element)!=myPrivateVariablePlus.nodeset.end())continue;
			myPrivateVariablePlus.nodeset.push_back(myPrivateVariablePlus.element);
			myPrivateVariablePlus.node=myPrivateVariablePlus.graphdata->g->addNode();
			myPrivateVariablePlus.graphdata->nodeNum++;
			myPrivateVariablePlus.graphdata->node2label->set(myPrivateVariablePlus.node,myPrivateVariablePlus.element);
			myPrivateVariablePlus.graphdata->node2degree->set(myPrivateVariablePlus.node,0);
			(*myPrivateVariablePlus.graphdata->label2node)[myPrivateVariablePlus.element]=myPrivateVariablePlus.node;
			//assert((*networks.getGraph(i)->invIdNodeMap).find(element)!=(*networks.getGraph(i)->invIdNodeMap).end());
		}
		for(myPrivateVariablePlus.sj=0;myPrivateVariablePlus.sj<myPrivateVariablePlus.nodeset.size();myPrivateVariablePlus.sj++)
		{
			myPrivateVariablePlus.protein1=myPrivateVariablePlus.nodeset[myPrivateVariablePlus.sj];
			myPrivateVariablePlus.node1=(*myPrivateVariablePlus.graphdata->label2node)[myPrivateVariablePlus.protein1];
			for(myPrivateVariablePlus.sk=myPrivateVariablePlus.sj+1;myPrivateVariablePlus.sk<myPrivateVariablePlus.nodeset.size();myPrivateVariablePlus.sk++)
			{
				myPrivateVariablePlus.protein2=myPrivateVariablePlus.nodeset[myPrivateVariablePlus.sk];
				myPrivateVariablePlus.element.clear();
			  myPrivateVariablePlus.node2=(*myPrivateVariablePlus.graphdata->label2node)[myPrivateVariablePlus.protein2];
				if(myPrivateVariablePlus.protein1.compare(myPrivateVariablePlus.protein2)>0)
				{
					myPrivateVariablePlus.element.append(myPrivateVariablePlus.protein2);
					myPrivateVariablePlus.element.append(myPrivateVariablePlus.protein1);
				 }
				else
				{
				 myPrivateVariablePlus.element.append(myPrivateVariablePlus.protein1);
				 myPrivateVariablePlus.element.append(myPrivateVariablePlus.protein2);
				}
				 
				if(networks.getGraph(myPrivateVariablePlus.si)->interactionmap.find(myPrivateVariablePlus.element)==networks.getGraph(myPrivateVariablePlus.si)->interactionmap.end())continue;
				 myPrivateVariablePlus.graphdata->g->addEdge(myPrivateVariablePlus.node1,myPrivateVariablePlus.node2);
				 (*(myPrivateVariablePlus.graphdata->node2degree))[myPrivateVariablePlus.node1]++;
				 (*(myPrivateVariablePlus.graphdata->node2degree))[myPrivateVariablePlus.node2]++;
				 if((*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node1]>myPrivateVariablePlus.graphdata->maxDegree)
					 myPrivateVariablePlus.graphdata->maxDegree=(*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node1];
				 if((*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node2]>myPrivateVariablePlus.graphdata->maxDegree)
					 myPrivateVariablePlus.graphdata->maxDegree=(*myPrivateVariablePlus.graphdata->node2degree)[myPrivateVariablePlus.node2];
				 myPrivateVariablePlus.graphdata->edgeNum++;

			}
		}
		myPrivateVariablePlus.subnet->subgraphs.push_back(myPrivateVariablePlus.graphdata);
	}
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
void Search<NP, SN, LG, OP>::verifyspineParallel(LayerGraph& layergraph, NetworkPool& networks) {
	int nodenum=layergraph.nodeNum;
	PrivateVariable myPrivateVariable(1);
	std::vector<NodeIt> mynodeit;
	for(NodeIt myit=NodeIt(layergraph.graph);myit!=lemon::INVALID;++myit)
	{
		mynodeit.push_back(myit);
	}
	#pragma omp parallel for num_threads(_numthreads) schedule(dynamic,1) shared(layergraph,networks,mynodeit,nodenum) firstprivate(myPrivateVariable)
	for(int i=0;i<nodenum;i++)
	{
		myPrivateVariable.node=mynodeit[i];
		if(sampleKSpineParallel(myPrivateVariable,layergraph,networks))
		{
			#pragma omp critical
			{
				layergraph.setConfiguration(myPrivateVariable.node);
			}
		}
	}
}

template<typename NP, typename SN, typename LG, typename OP>
void Search<NP,SN,LG,OP>::output(LayerGraph& layergraph,PrivateVariablePlus& myPrivateVariablePlus,std::ofstream& fout)
{
	myPrivateVariablePlus.element.clear();
	myPrivateVariablePlus.element.append(_resultfolder);
	myPrivateVariablePlus.element.append("alignments/ucomplex_");
	myPrivateVariablePlus.element.append(convert_num2str(myPrivateVariablePlus.si));
	myPrivateVariablePlus.element.append(".txt");
	fout.open(myPrivateVariablePlus.element.c_str());
	outsubgraphs(layergraph,myPrivateVariablePlus,fout);
	fout.close();
	/*for(myPrivateVariablePlus.st=0;myPrivateVariablePlus.st<_numSpecies;++myPrivateVariablePlus.st)
	{
		myPrivateVariablePlus.element.clear();
		myPrivateVariablePlus.element.append(_resultfolder);
		myPrivateVariablePlus.element.append("species_");
		myPrivateVariablePlus.element.append(convert_num2str(myPrivateVariablePlus.st));
		myPrivateVariablePlus.element.append("/complex_");
		myPrivateVariablePlus.element.append(convert_num2str(myPrivateVariablePlus.si));
		myPrivateVariablePlus.element.append(".txt");
		fout.open(myPrivateVariablePlus.element.c_str());
		myPrivateVariablePlus.proteinmap.clear();
		fout <<"# score: " << myPrivateVariablePlus.overallScore << "\n";
		for(myPrivateVariablePlus.sk=0;myPrivateVariablePlus.sk<myPrivateVariablePlus.subnet->net_spines.size();++myPrivateVariablePlus.sk)
		{
			myPrivateVariablePlus.protein=layergraph.node2label[myPrivateVariablePlus.subnet->net_spines[myPrivateVariablePlus.sk].data[myPrivateVariablePlus.st]];
			if(myPrivateVariablePlus.proteinmap.find(myPrivateVariablePlus.protein)!=myPrivateVariablePlus.proteinmap.end())
				continue;
			else
				myPrivateVariablePlus.proteinmap[myPrivateVariablePlus.protein]=true;
			fout << myPrivateVariablePlus.protein <<"\n";
		}
		fout.close();
	}*/
}

template<typename NP, typename SN, typename LG, typename OP>
void Search<NP,SN,LG,OP>::outsubgraphs(LayerGraph& layergraph,PrivateVariablePlus& myPrivateVariablePlus,std::ofstream& fout)
{
	fout <<"# score: " << myPrivateVariablePlus.overallScore << "\n";
	for(myPrivateVariablePlus.st=0;myPrivateVariablePlus.st<myPrivateVariablePlus.subnet->net_spines.size();++myPrivateVariablePlus.st)
	{
		for(myPrivateVariablePlus.sj=0;myPrivateVariablePlus.sj<_numSpecies;++myPrivateVariablePlus.sj)
		{
			fout << layergraph.node2label[myPrivateVariablePlus.subnet->net_spines[myPrivateVariablePlus.st].data[myPrivateVariablePlus.sj]]<<"\t";
		}
		fout <<"\n";
	}
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::run(LayerGraph& layergraph,NetworkPool& networks)
{
	searchSeedsParallel(layergraph,networks);
	int csize=refinedSeeds.size();
	int outnum=0;
	//std::cout <<"Seed size: " << _seedSize << std::endl;
	//std::cout <<"Refind seeds: " << csize <<std::endl;
	//std::cout <<"Seeds tries: " << _seedTries <<std::endl;
	//std::cout << "Min subnet:" << _seedSize+_minExt << std::endl;
	//std::cout << "Max subnet:" << _seedSize+_maxExt << std::endl;
	MyTree localtree;
	PrivateVariable myPrivateVariable(layergraph.validnodes.size()-1);
	std::vector<SubNet*> mySubNetList;
	#pragma omp parallel for num_threads(_numthreads) schedule(dynamic,1) shared(layergraph,networks,csize,mySubNetList) firstprivate(myPrivateVariable)
	for(int i=0;i<csize;i++)
	{
		for(myPrivateVariable.k=_minExt;myPrivateVariable.k<=_maxExt;myPrivateVariable.k++)
		{	
			setExtension(myPrivateVariable);
			for(myPrivateVariable.j=0; myPrivateVariable.j<_seedTries;myPrivateVariable.j++)//_seedTries
			{
				myPrivateVariable.subnet=*refinedSeeds[i];// a copy of refinedSeeds;new SubNet();
			 	expandRefinedSeeds(myPrivateVariable,layergraph,networks);
#pragma omp critical
				{
				mySubNetList.push_back(new SubNet(myPrivateVariable.subnet));
				}
			}
		}
	}

	PrivateVariablePlus myPrivateVariablePlus;
	std::ofstream fout;
	csize=mySubNetList.size();
	localtree.readTree(_treefile);
	//std::cout <<"Subnets:" << csize << std::endl;
#pragma omp parallel for num_threads(_numthreads) schedule(dynamic,1) shared(layergraph,networks,mySubNetList,localtree,fout,outnum) firstprivate(myPrivateVariablePlus)
	for(int i=0;i<csize;i++)
	{
		#pragma omp critical
		{
			//std::cout << i+1 <<"/"<<csize << std::endl;
			myPrivateVariablePlus.scoremap=new ScoreEdgeMap(localtree.g);
			myPrivateVariablePlus.matchingedgemap=new MatchingEdgeMap(localtree.g);
		}
		myPrivateVariablePlus.subnet=mySubNetList[i];
		induceSubgraphs(myPrivateVariablePlus,layergraph,networks);
		myPrivateVariablePlus.phylogeny._dsize=myPrivateVariablePlus.subnet->net_spines.size();
		initialPhylogy(myPrivateVariablePlus,layergraph,localtree);
		simulatedAnnealingMethod(myPrivateVariablePlus,localtree);
		myPrivateVariablePlus.si=i;
		if(myPrivateVariablePlus.overallScore > _score_threshold)
		{
			#pragma omp critical
			{
				outnum++;
				output(layergraph,myPrivateVariablePlus,fout);
			}
		}
		clearStructure(myPrivateVariablePlus,localtree);
	}
	//std::cout << "Threshold: " << _score_threshold << std::endl;
	//std::cout << "Alignments: " << outnum << std::endl;
}

template<typename NP, typename SN, typename LG, typename OP>
void
	Search<NP,SN,LG,OP>::simulatedAnnealingMethod(PrivateVariablePlus& myPrivateVariable,const MyTree& localtree)
{
	myPrivateVariable.sk=0;
	myPrivateVariable.st=myPrivateVariable.simulatedannealing._tmax;
	myPrivateVariable.step=(myPrivateVariable.simulatedannealing._tmax-myPrivateVariable.simulatedannealing._tmin)/myPrivateVariable.simulatedannealing._Kmax;
	myPrivateVariable.seed =std::chrono::system_clock::now().time_since_epoch().count();
	myPrivateVariable.generator=std::default_random_engine(myPrivateVariable.seed);
	while(myPrivateVariable.sk++<=myPrivateVariable.simulatedannealing._Kmax)
	{
		myPrivateVariable.st = myPrivateVariable.st-myPrivateVariable.step;//
		myPrivateVariable.beta = -1.0/(myPrivateVariable.simulatedannealing._k*myPrivateVariable.st);
		for(myPrivateVariable.si=0;myPrivateVariable.si<myPrivateVariable.simulatedannealing._Nmax;++myPrivateVariable.si)
		{
			if(!interfere(myPrivateVariable,localtree))continue;// unsafe for example sumup();
			myPrivateVariable.sampledata=myPrivateVariable.distribution(myPrivateVariable.generator);
			if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
				std::cout <<myPrivateVariable.deltaData.delta <<"\t" << myPrivateVariable.sampledata<<"\t"<< exp(myPrivateVariable.beta*myPrivateVariable.deltaData.delta) <<"\n";
			if(myPrivateVariable.deltaData.delta<0 || myPrivateVariable.sampledata < exp(myPrivateVariable.beta*myPrivateVariable.deltaData.delta))
			{
				// update current state to the neighbor state and its interaction evolutionary score.
				myPrivateVariable.sj=0;
				for(myPrivateVariable.incE=IncEdgeIt(localtree.g,myPrivateVariable.deltaData.treenode);myPrivateVariable.incE!=lemon::INVALID;++myPrivateVariable.incE,++myPrivateVariable.sj)
				{
					(*myPrivateVariable.scoremap)[myPrivateVariable.incE]=myPrivateVariable.deltaData.updatedScores[myPrivateVariable.sj];
				}
			}
			else
			{
				myPrivateVariable.node=myPrivateVariable.deltaData.treenode;
				myPrivateVariable.id1=localtree.g.id(myPrivateVariable.node);
				myPrivateVariable.graphdata=myPrivateVariable.phylogeny.node2graph[myPrivateVariable.id1];
				if(myPrivateVariable.graphdata->label2edge->find(myPrivateVariable.deltaData.edgelabel)!=myPrivateVariable.graphdata->label2edge->end())
				{
					 myPrivateVariable.myedge=myPrivateVariable.graphdata->label2edge->find(myPrivateVariable.deltaData.edgelabel)->second;
					 myPrivateVariable.edgelabel=myPrivateVariable.deltaData.edgelabel;
					 deleteEdge(myPrivateVariable);
				}
				else
				{
					myPrivateVariable.node1=myPrivateVariable.deltaData.nodeA;
					myPrivateVariable.node2=myPrivateVariable.deltaData.nodeB;
					addEdge(myPrivateVariable);
				}
			}
		}
	}
	computeDist(myPrivateVariable,localtree);
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::computeDist(PrivateVariablePlus& myPrivateVariable,const MyTree& localtree)
{
	myPrivateVariable.distEvolution.clear();
	for(myPrivateVariable.ie=EdgeIt(localtree.g);myPrivateVariable.ie!=lemon::INVALID;++myPrivateVariable.ie)
	{
		myPrivateVariable.distEvolution+=(*myPrivateVariable.scoremap)[myPrivateVariable.ie];
	}
	sumupScore(myPrivateVariable,myPrivateVariable.distEvolution);
	myPrivateVariable.overallScore=myPrivateVariable.phylogeny._dsize/myPrivateVariable.sumDist;
	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
		std::cout << localtree.overallScore << std::endl;
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::sumupScore(PrivateVariablePlus& myPrivateVariable,Score& myscore)
{
	myPrivateVariable.sumDist=0;
	for(myPrivateVariable.sj=0;myPrivateVariable.sj<4;myPrivateVariable.sj++)
	{
		myPrivateVariable.sumDist+=myscore.fscore[myPrivateVariable.sj];
	}
}


template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::setExtension(PrivateVariable& myPrivateVariable)
{
	myPrivateVariable.numExtension=myPrivateVariable.k;
}

template<typename NP, typename SN, typename LG, typename OP> 
void 
Search<NP,SN,LG,OP>::searchCandidates(std::unordered_map<int,bool>& usedValidNodes,
																			std::vector<Node>& candidates,
																			K_Spine& pspine,
																			LayerGraph& layergraph,
																			NetworkPool& networks) 
{
	std::vector<std::string> innerproteins;
	std::vector<std::string> neighborproteins;
	for(unsigned j=0;j<_numSpecies;j++)
	{
		innerproteins.push_back(layergraph.node2label[pspine.data[j]]);
	}         networks.getNeighbors(innerproteins,neighborproteins);
	for(unsigned j=0;j<neighborproteins.size();++j)
	{
		std::string protein=neighborproteins[j];
		Node neighbor=layergraph.label2node[protein];
		int neighborid=layergraph.graph.id(neighbor);
		if(layergraph.validnodemap.find(neighborid)==layergraph.validnodemap.end())continue;
		if(usedValidNodes.find(neighborid)!=usedValidNodes.end())continue;
		if(find(candidates.begin(),candidates.end(),neighbor)!=candidates.end())continue;
		candidates.push_back(neighbor);
	} 
} 

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::searchCandidatesParallel(PrivateVariable& myprivateVariable,
																			LayerGraph& layergraph,
																			NetworkPool& networks)
{
	myprivateVariable.innerproteins.clear();
	myprivateVariable.neighborproteins.clear();
	for(myprivateVariable.sj=0;myprivateVariable.sj<_numSpecies;myprivateVariable.sj++)
	{
		myprivateVariable.innerproteins.push_back(layergraph.node2label[myprivateVariable.pspine.data[myprivateVariable.sj]]);
	}
	networks.getNeighbors(myprivateVariable.innerproteins,myprivateVariable.neighborproteins);
	for(myprivateVariable.sj=0;myprivateVariable.sj<myprivateVariable.neighborproteins.size();++myprivateVariable.sj)
	{
		myprivateVariable.protein=myprivateVariable.neighborproteins[myprivateVariable.sj];
		myprivateVariable.neighbor=layergraph.label2node[myprivateVariable.protein];
		myprivateVariable.neighborid=layergraph.graph.id(myprivateVariable.neighbor);
		if(layergraph.validnodemap.find(myprivateVariable.neighborid)==layergraph.validnodemap.end())continue;
		if(myprivateVariable.usedValidNodes.find(myprivateVariable.neighborid)!=myprivateVariable.usedValidNodes.end())continue;
		if(find(myprivateVariable.candidates.begin(),myprivateVariable.candidates.end(),myprivateVariable.neighbor)!=myprivateVariable.candidates.end())continue;
		myprivateVariable.candidates.push_back(myprivateVariable.neighbor);
	}
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::expandRefinedSeeds(PrivateVariable& myprivateVariable,
										LayerGraph& layergraph,
										NetworkPool& networks)
{
	
	//// searching neighbors of subnet.
	myprivateVariable.usedValidNodes.clear();
	myprivateVariable.candidates.clear();
	for(myprivateVariable.ei=0;myprivateVariable.ei<myprivateVariable.subnet.net_spines.size();myprivateVariable.ei++)
	{
		myprivateVariable.pspine=myprivateVariable.subnet.net_spines[myprivateVariable.ei];
		for(myprivateVariable.ej=0;myprivateVariable.ej<_numSpecies;myprivateVariable.ej++)
		{
			myprivateVariable.node=myprivateVariable.pspine.data[myprivateVariable.ej];
			myprivateVariable.nodeid=layergraph.graph.id(myprivateVariable.node);
			if(layergraph.validnodemap.find(myprivateVariable.nodeid)==layergraph.validnodemap.end())continue;
			myprivateVariable.usedValidNodes[layergraph.graph.id(myprivateVariable.node)]=true;
		}
		
		searchCandidatesParallel(myprivateVariable, 
						layergraph, 
						networks);
	}
	myprivateVariable.it=myprivateVariable.candidates.begin();
	while(myprivateVariable.it!=myprivateVariable.candidates.end())
	{
		if(myprivateVariable.usedValidNodes.find(layergraph.graph.id(*myprivateVariable.it))
		  !=myprivateVariable.usedValidNodes.end())
			myprivateVariable.it=myprivateVariable.candidates.erase(myprivateVariable.it);
		else
			++myprivateVariable.it;
	}
	myprivateVariable.num=0;
	while(myprivateVariable.num++<myprivateVariable.numExtension)
	{
		heuristicSearch(myprivateVariable,
						layergraph,
						networks);
	}
	//assert(myprivateVariable.subnet.net_spines.size()==static_cast<unsigned>(myprivateVariable.numExtension+_seedSize));
	
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::heuristicSearch(PrivateVariable& myprivateVariable,LayerGraph& layergraph,NetworkPool& networks)
{
	if(myprivateVariable.candidates.size()== 0)
	{
		myprivateVariable.dice_roll = myprivateVariable.discrete(generator);
		myprivateVariable.node=layergraph.validnodes[myprivateVariable.dice_roll];
		while(myprivateVariable.usedValidNodes.find(layergraph.graph.id(myprivateVariable.node))!=myprivateVariable.usedValidNodes.end())
		{
			myprivateVariable.dice_roll = myprivateVariable.discrete(generator);
			myprivateVariable.node=layergraph.validnodes[myprivateVariable.dice_roll];
		}
	}else
	{
		myprivateVariable.dice_roll = distribution(generator)%myprivateVariable.candidates.size();
		myprivateVariable.node=myprivateVariable.candidates[myprivateVariable.dice_roll];
		myprivateVariable.candidates.erase(myprivateVariable.candidates.begin()+myprivateVariable.dice_roll);
	}
	
	myprivateVariable.usedValidNodes[layergraph.graph.id(myprivateVariable.node)]=true;
	//if(!sampleKSpine(myprivateVariable.node,myprivateVariable.pspine,layergraph,networks))
	myprivateVariable.spine.clear();
	if(!sampleKSpineParallel(myprivateVariable,layergraph,networks))// parallize sampleSpine
	{
		std::cerr <<"Invalid sample node!"<<std::endl;
		return false;
	}
	myprivateVariable.subnet.net_spines.push_back(myprivateVariable.pspine);
	searchCandidatesParallel(myprivateVariable,layergraph,networks);
	
	return true;
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::searchSeeds(LayerGraph& layergraph,NetworkPool& networks)
{
	verifyspine(layergraph,networks);
	std::uniform_int_distribution<int> discrete(0,layergraph.validnodes.size()-1);
	int num=0;
	while(num++<_numSamples)
	{
		SubNet* mysubnet= new SubNet(_numSpecies,_seedSize);
	 if(!sampleSeed(mysubnet,layergraph,networks,discrete))
		{
			delete mysubnet;
			if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
			std::cerr << "Failed to sample a refined seed!" << std::endl;
		}else
		{
			refinedSeeds.push_back(mysubnet);
		}
	}
	std::cout << "There are a total of "<< refinedSeeds.size() <<" refined seeds."<< std::endl;
}

template<typename NP, typename SN, typename LG, typename OP>
void
Search<NP,SN,LG,OP>::searchSeedsParallel(LayerGraph& layergraph, NetworkPool& networks)
{
	verifyspineParallel(layergraph, networks);
	PrivateVariable myPrivateVariable(layergraph.validnodes.size()-1);
#pragma omp parallel for num_threads(_numthreads) schedule(dynamic,1) shared(layergraph,networks) firstprivate(myPrivateVariable)
	for(int i=0; i<_numSamples;i++)
	{
		//myPrivateVariablePlus.clear();
		myPrivateVariable.mysubnet=new SubNet(_numSpecies,_seedSize);// used local variable i, mysubnet 
		if(!sampleSeedParallel(myPrivateVariable,layergraph,networks))
		{
			delete myPrivateVariable.mysubnet;
		}
		else
		{
#pragma omp critical
			{
				refinedSeeds.push_back(myPrivateVariable.mysubnet);
			}
		}
	}
}

/// Randomly sample a d subnet from G_h
template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::sampleSeed(SubNet* subnet, LayerGraph& layergraph, NetworkPool& networks,std::uniform_int_distribution<int>& discrete)
{
	/// Randomly sample a k-spine containing node.
	int dice_roll = discrete(generator);
	Node firstnode = layergraph.validnodes[dice_roll];
	std::vector<Node> candidates;
	std::unordered_map<int,bool> usedValidNodes;
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
			/// randomly sample a neighbor in candidates.
			dice_roll = distribution(generator)%candidates.size();
			node=candidates[dice_roll];
			candidates.erase(candidates.begin()+dice_roll);
		}
		K_Spine pspine;
		usedValidNodes[layergraph.graph.id(node)]=true;
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
		searchCandidates(usedValidNodes,candidates, pspine, layergraph, networks);
		subnet->net_spines.push_back(pspine);
	}
	subnet->induceSubgraphs(networks,layergraph);
	return checkConnection(subnet,layergraph,networks);
}

/// Randomly sample a d subnet from G_h
template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::sampleSeedParallel(PrivateVariable& myPrivateVariable, LayerGraph& layergraph, NetworkPool& networks)
// used local variable i, mysubnet 
{
	myPrivateVariable.dice_roll = myPrivateVariable.discrete(myPrivateVariable.generator);
	myPrivateVariable.firstnode = layergraph.validnodes[myPrivateVariable.dice_roll];
	////std::vector<Node> candidates;
	myPrivateVariable.candidates.clear();
	myPrivateVariable.usedValidNodes.clear();
	////std::unordered_map<int,bool> usedValidNodes;
	for(myPrivateVariable.ei=0;myPrivateVariable.ei<_seedSize;++myPrivateVariable.ei)
	{
		if(0==myPrivateVariable.ei) myPrivateVariable.node = myPrivateVariable.firstnode;// firstnode, node
		else if(myPrivateVariable.candidates.size()==0)
		{
			myPrivateVariable.dice_roll = myPrivateVariable.discrete(myPrivateVariable.generator);
			myPrivateVariable.node=layergraph.validnodes[myPrivateVariable.dice_roll];
		}
	else
	{
	//		/// randomly sample a neighbor in candidates.
			myPrivateVariable.dice_roll = distribution(myPrivateVariable.generator)%myPrivateVariable.candidates.size();
			myPrivateVariable.node=myPrivateVariable.candidates[myPrivateVariable.dice_roll];
			myPrivateVariable.candidates.erase(myPrivateVariable.candidates.begin()+myPrivateVariable.dice_roll);
		}
	////	//K_Spine pspine;
		myPrivateVariable.usedValidNodes[layergraph.graph.id(myPrivateVariable.node)]=true;
		if(!sampleKSpineParallel(myPrivateVariable,layergraph,networks))
		{
			std::cerr <<"Invalid sample node!"<<std::endl;
			return false;
		}

		/// Output this sample.
		if(g_verbosity >= VERBOSE_NON_ESSENTIAL)
		{
			for(myPrivateVariable.ej=0;myPrivateVariable.ej<_numSpecies;++myPrivateVariable.ej)
				std::cout << layergraph.node2label[myPrivateVariable.pspine.data[myPrivateVariable.ej]]<<" ";
			std::cout << std::endl;
		}

		searchCandidatesParallel(myPrivateVariable,layergraph,networks);
		myPrivateVariable.mysubnet->net_spines.push_back(myPrivateVariable.pspine);
	}
	//myPrivateVariable.mysubnet->induceSubgraphs(networks,layergraph);
	induceSubgraphsPrivateVar(myPrivateVariable,layergraph,networks);
	return checkConnectionParallel(myPrivateVariable,layergraph,networks);
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
	Search<NP,SN,LG,OP>::checkConnection(SubNet* subnet,LayerGraph& layergraph,NetworkPool& networks)
{
	int numConnected=0;
	for(unsigned i=0;i<_numSpecies;++i)
	{
		if(subnet->subgraphs[i]->edgeNum< (subnet->subgraphs[i]->nodeNum-NUM_GAP_EDGE) )continue;/// To fulfills the minimal number of edges. 
		numConnected++;
	}
	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
	std::cout << "There are "<< numConnected << " out of " <<_numSpecies <<" subnetwork(s) are connected!" << std::endl;
	return numConnected >= _numConnected;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::checkConnectionParallel(PrivateVariable& myprivateVariable,LayerGraph& layergraph,NetworkPool& networks)
{
	myprivateVariable.numConnected=0;
	for( myprivateVariable.ei=0;myprivateVariable.ei<_numSpecies;++myprivateVariable.ei)
	{
		if(myprivateVariable.mysubnet->subgraphs[myprivateVariable.ei]->edgeNum< (myprivateVariable.mysubnet->subgraphs[myprivateVariable.ei]->nodeNum-NUM_GAP_EDGE) )continue;/// To fulfills the minimal number of edges. 
		myprivateVariable.numConnected++;
	}
	if(g_verbosity>=VERBOSE_NON_ESSENTIAL)
	std::cout << "There are "<<myprivateVariable.numConnected << " out of " <<_numSpecies <<" subnetwork(s) are connected!" << std::endl;
	return myprivateVariable.numConnected >= _numConnected;
}

template<typename NP, typename SN, typename LG, typename OP>
bool
Search<NP,SN,LG,OP>::sampleKSpineParallel(PrivateVariable& myprivateVariable,LayerGraph& layergraph,NetworkPool& networks)
{
	
	myprivateVariable.m_candidates.clear();//m_candidates -> candidates for expand to a kspine;
	myprivateVariable.host=networks.getHost(layergraph.node2label[myprivateVariable.node]);
	myprivateVariable.spine.clear();
	myprivateVariable.spine.data[myprivateVariable.host]=myprivateVariable.node;
	myprivateVariable.spine.states[myprivateVariable.host]=true;
	if(g_verbosity==VERBOSE_DEBUG)
	std::cerr << layergraph.node2label[myprivateVariable.node] << std::endl;
	
	for(myprivateVariable.inc=IncEdgeIt(layergraph.graph,myprivateVariable.node);myprivateVariable.inc!=lemon::INVALID;++myprivateVariable.inc)
	{
		myprivateVariable.rnode=layergraph.graph.runningNode(myprivateVariable.inc);
		myprivateVariable.host=networks.getHost(layergraph.node2label[myprivateVariable.rnode]);
		if(g_verbosity==VERBOSE_DEBUG)
		std::cerr << layergraph.node2label[myprivateVariable.rnode] << std::endl;
		if(myprivateVariable.spine.states[myprivateVariable.host])continue;
		if(myprivateVariable.node < myprivateVariable.rnode) myprivateVariable.m_candidates.push_back(myprivateVariable.rnode);
	}
	myprivateVariable.sk=_numSpecies;
	return expandspineParallel(myprivateVariable,myprivateVariable.spine,myprivateVariable.m_candidates, layergraph,networks);
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
Search<NP,SN,LG,OP>::expandspineParallel(
									PrivateVariable& myprivateVariable
								, K_Spine spine
								, std::vector<Node> candidates
								, LayerGraph& layergraph
								, NetworkPool& networks)
{
	if(0==--myprivateVariable.sk)
	{
		myprivateVariable.pspine=spine;
		return true;
	}
	while(!candidates.empty())
	{
		//unsigned host;
		myprivateVariable.dice_roll = distribution(generator)%(candidates.size());
		myprivateVariable.wnode = candidates[myprivateVariable.dice_roll];
		candidates.erase(candidates.begin()+myprivateVariable.dice_roll);
		myprivateVariable.secondCandidates.clear();
		myprivateVariable.exclCandidates.clear();
		/// Assign V'UN(V') to exclCandidates.
		for(myprivateVariable.sj=0;myprivateVariable.sj<_numSpecies;myprivateVariable.sj++)
		{
			if(!spine.states[myprivateVariable.sj])continue;
			myprivateVariable.spinenode=spine.data[myprivateVariable.sj];
			myprivateVariable.exclCandidates.push_back(myprivateVariable.spinenode);
	
			for(myprivateVariable.inc=IncEdgeIt(layergraph.graph,myprivateVariable.spinenode);myprivateVariable.inc!=lemon::INVALID;++myprivateVariable.inc)
			{
				myprivateVariable.rnode=layergraph.graph.runningNode(myprivateVariable.inc);
				myprivateVariable.host=networks.getHost(layergraph.node2label[myprivateVariable.rnode]);
				if(spine.states[myprivateVariable.host])continue;/// If yes, it's not a neighbor.			
				if( std::find(myprivateVariable.exclCandidates.begin(), myprivateVariable.exclCandidates.end(), myprivateVariable.rnode)
				    != myprivateVariable.exclCandidates.end())continue;
				myprivateVariable.exclCandidates.push_back(myprivateVariable.rnode);			
			}
		}
		/// Assign new value to secondCandidates.
		for(myprivateVariable.inc=IncEdgeIt(layergraph.graph,myprivateVariable.wnode);myprivateVariable.inc!=lemon::INVALID;++myprivateVariable.inc)
		{
			myprivateVariable.rnode=layergraph.graph.runningNode(myprivateVariable.inc);
			myprivateVariable.host=networks.getHost(layergraph.node2label[myprivateVariable.rnode]);
			if(spine.states[myprivateVariable.host])continue;
			if( std::find(myprivateVariable.exclCandidates.begin(), myprivateVariable.exclCandidates.end(), myprivateVariable.rnode)
				    != myprivateVariable.exclCandidates.end())continue;
			if (myprivateVariable.node < myprivateVariable.rnode) myprivateVariable.secondCandidates.push_back(myprivateVariable.rnode);				
		}
		myprivateVariable.newspine=spine;
		myprivateVariable.host=networks.getHost(layergraph.node2label[myprivateVariable.wnode]);
		if(g_verbosity==VERBOSE_DEBUG)
		std::cerr << layergraph.node2label[myprivateVariable.wnode] << std::endl;
		myprivateVariable.newspine.data[myprivateVariable.host]=myprivateVariable.wnode;
		myprivateVariable.newspine.states[myprivateVariable.host]=true;
		for(myprivateVariable.sj=0;myprivateVariable.sj<candidates.size();++myprivateVariable.sj)
		{
			myprivateVariable.rnode=candidates[myprivateVariable.sj];
		  myprivateVariable.host=networks.getHost(layergraph.node2label[myprivateVariable.rnode]);
			if(myprivateVariable.newspine.states[myprivateVariable.host])continue;
			myprivateVariable.secondCandidates.push_back(myprivateVariable.rnode);
		}
		if(expandspineParallel(myprivateVariable, myprivateVariable.newspine, myprivateVariable.secondCandidates, layergraph, networks))return true;
	}
	myprivateVariable.sk++;
	return false;
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
#endif
