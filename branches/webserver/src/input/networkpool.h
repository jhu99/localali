/* networkpool.h
Author: Jialu Hu
Date: 18.06.2012*/

#ifndef NETWORKPOOL_H_
#define NETWORKPOOL_H_

#include <vector>
#include <array>
#include <lemon/core.h>
#include <lemon/bits/graph_extender.h>
#include <lemon/connectivity.h>
#include <unordered_map>
#include "verbose.h"
#include <assert.h>
#include "macro.h"
#include <time.h>
#include <stdlib.h>
#include "algorithm/function.h"

using namespace lemon;

template<typename GR, typename BP>
class NetworkPool
{
private:
	struct EdgeElement
	{
		std::string startNode;
		std::string endNode;
	};

public:
  typedef GR Graph;
  typedef BP BpGraph;
  TEMPLATE_GRAPH_TYPEDEFS(Graph);

  /// Labels of the nodes
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
  /// Mapping from labels to original nodes
  typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;
  /// Weights on original edges
  typedef typename Graph::template EdgeMap<int> WeightEdgeMap;
  /// Host of proteins for all input networks
  std::unordered_map<std::string,short> proteinHost;
  typedef typename Graph::template NodeMap<unsigned> DegreeNodeMap;
  int allNodeNum;
  typedef std::array<EdgeElement,MAX_NETWORK_EDGE> EdgeArray;
  typedef std::array<std::string,MAX_NETWORK_NODE> NodeArray;

  struct GraphData
  {
    Graph *g;
    OrigLabelNodeMap *label;
    InvOrigLabelNodeMap *invIdNodeMap;
    DegreeNodeMap *degreeMap;
    WeightEdgeMap *weight;
    std::unordered_map<std::string,int> interactionmap;
    lemon::ArcLookUp<Graph> *arcLookUpG;
    int nodeNum;
    int edgeNum;
    GraphData()
    {
      nodeNum=0;
      edgeNum=0;
      g = new Graph();
      label = new OrigLabelNodeMap(*g);
      invIdNodeMap = new InvOrigLabelNodeMap();
      degreeMap = new DegreeNodeMap(*g);
      weight = new WeightEdgeMap(*g);
      arcLookUpG = new lemon::ArcLookUp<Graph>(*g);
    }
    ~GraphData()
    {
      delete g;
      delete label;
      delete invIdNodeMap;
      delete degreeMap;
      delete weight;
      delete arcLookUpG;
    }
  };
  std::vector<GraphData*> _graphSet;
  
  NetworkPool();
  ~NetworkPool();
  bool initNetworkPool(std::vector<std::string>&);
  bool readNetwork(std::string&,short);
  GraphData* getGraph(int);
  unsigned getHost(std::string);
  bool existNode(std::string);
  bool readSubgraph(GraphData*,std::vector<std::string>&,Graph&);
  bool getNeighbors(std::vector<std::string>&, std::vector<std::string>&);
  bool outputNetworks();
  bool generateRandNetworks(EdgeArray&,int,std::string);
  bool generateUniNetworks(NodeArray&,int,int,std::string);
  bool outputMaWIShNetworks();
  void outputMaWIShSimilarity(std::string);
};

template<typename GR, typename BP>
NetworkPool<GR,BP>::NetworkPool()
:allNodeNum(0)
{
}

template<typename GR, typename BP>
NetworkPool<GR,BP>::~NetworkPool()
{
  for(unsigned i=0;i<_graphSet.size();i++)
    delete _graphSet[i];
}

template<typename GR, typename BP>
typename NetworkPool<GR,BP>::GraphData* NetworkPool<GR,BP>::getGraph(int i)
{
  //std::cout << _graphSet.size() << std::endl;
  return _graphSet[i];
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::getNeighbors(
std::vector<std::string>& innerproteins,
std::vector<std::string>& neighborproteins)
{
	for(unsigned i=0;i<innerproteins.size();++i)
	{
		std::string protein=innerproteins[i];
		GraphData *graphdata=getGraph(i);
		Node node=(*graphdata->invIdNodeMap)[protein];
		for (IncEdgeIt e(*graphdata->g,node);e!=lemon::INVALID;++e)
		{
			Node rnode=graphdata->g->runningNode(e);
			neighborproteins.push_back((*graphdata->label)[rnode]);
		}
	}
	return true;
}

template<typename GR, typename BP>
unsigned NetworkPool<GR,BP>::getHost(std::string protein)
{
  if(proteinHost.find(protein)==proteinHost.end()) return 100;// 100 represents protein doesn't exist in input networks.
  return proteinHost[protein];
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::existNode(std::string protein)
{
  if(proteinHost.find(protein)!=proteinHost.end())
  return true;
  else
  return false;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::readSubgraph(GraphData* data, std::vector<std::string>& nodeset, Graph&)
{
	return true;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::initNetworkPool(std::vector<std::string> &filelist)
{
  std::vector<std::string>::iterator it;
  short i=0;
  for(it=filelist.begin();it!=filelist.end();++i,++it)
  {
    readNetwork(*it,i);
  }
  return 1;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::readNetwork(std::string &filename,short i)
{
  GraphData *data = new GraphData();
  std::string line;
  
  std::ifstream input(filename.c_str());
  if(!input.good())
  {
    std::cerr << filename <<" cannot be opened!"<<std::endl;
    return 0;
  }
  _graphSet.push_back(data);
  // std::getline(input,line);/// Skip header line: INTERACTOR A INTERACTOR B
  while(std::getline(input,line))
  {
    std::string protein1, protein2, tempstr, keystr;
    std::stringstream lineStream(line);
    lineStream >> protein1 >> protein2;
    if(protein1.compare(protein2)==0)
      continue;
    else if(protein1.compare(protein2)>0)
     {
		 tempstr = protein1;
		 protein1 = protein2;
		 protein2 = tempstr;
	 }
	 keystr.append(protein1);
	 keystr.append(protein2);
	 if(data->interactionmap.find(keystr)==data->interactionmap.end())
		data->interactionmap[keystr]=1;
	 else
	 {
		 data->interactionmap[keystr]++;
		 continue;
	 }
    Node node1,node2;
    if(data->invIdNodeMap->find(protein1) == data->invIdNodeMap->end())
    {
      node1 = data->g->addNode();
      data->label->set(node1,protein1);
      (*(data->invIdNodeMap))[protein1] = node1;
      proteinHost[protein1] = i;
      data->nodeNum++;
    }else
    {
      node1 = (*(data->invIdNodeMap))[protein1];
    }
    if(data->invIdNodeMap->find(protein2) == data->invIdNodeMap->end())
    {
      node2 = data->g->addNode();
      data->label->set(node2,protein2);
      (*data->invIdNodeMap)[protein2] = node2;
      proteinHost[protein2] = i;
      data->nodeNum++;
    }else
    {
      node2 = (*data->invIdNodeMap)[protein2];
    }

    data->g->addEdge(node1, node2);
    data->edgeNum++;
  }
  data->arcLookUpG->refresh();
  unsigned maxNode=0;
  for(NodeIt mynode(*(data->g));mynode!=lemon::INVALID;++mynode)
  {
	  for(IncEdgeIt myedge(*(data->g),mynode);myedge!=lemon::INVALID;++myedge)
	  {
		  (*data->degreeMap)[mynode]++;
	  }
	  if(maxNode < (*data->degreeMap)[mynode])
		maxNode=(*data->degreeMap)[mynode];
  }
  if(g_verbosity>=VERBOSE_NONE)
  {
    std::cout <<"# " <<filename <<" has been read successfully!<br>"<<std::endl;
    std::cout <<"# number of proteins:"<< data->nodeNum <<"<br>"<<std::endl;
    std::cout <<"# number of interactions:"<<data->edgeNum<<"<br>"<<std::endl;
    std::cout <<"# the largest degree:"<< maxNode <<"<br>"<< std::endl;
  }
  allNodeNum+=data->nodeNum;
  return 1;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::outputNetworks()
{
	unsigned tSize=_graphSet.size();
	std::vector<std::string> prefixName;
	prefixName.push_back("1-human9606");prefixName.push_back("2-worm6239");
	prefixName.push_back("3-fly7227");prefixName.push_back("4-yeast4932");prefixName.push_back("7-ecoli562");
	for(unsigned i=0;i<tSize;i++)
	{
		GraphData* data=getGraph(i);
		//EdgeArray edgelist;
		NodeArray nodelist;
		int numEdge=0;
		int numNode=0;
		for(NodeIt inode(*(data->g));inode!=lemon::INVALID;++inode)
		{
			nodelist[numNode]=(*(data->label))[inode];
			numNode++;
		}
		numEdge=data->edgeNum;
		generateUniNetworks(nodelist,numNode,numEdge,prefixName[i]);
	}
	return true;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::generateRandNetworks(EdgeArray& edgelist,int num,std::string prename)
{
	std::string tmp;
	for(int i=1;i<=NUM_RANDOM_NETWORK;i++)
	{
		for(int j=1;j<=NUM_ITERATIONS;j++)
		{
			 srand (time(NULL));
			for(int m=0;m<num;m++)
			{
				int k=rand()%num;
				tmp = edgelist[m].endNode;
				edgelist[m].endNode=edgelist[k].endNode;
				edgelist[k].endNode=tmp;				
			}
		}
		std::string filename("../crosslink/dataset/rawdata/random/");
		filename.append(prename);filename.append("_r");filename.append(convert_num2str(i));filename.append(".txt");
		std::ofstream output(filename.c_str());
		for(int m=0;m<num;m++)
		{
			output << edgelist[m].startNode <<"\t" << edgelist[m].endNode <<"\t" << "0.9\n";
		}
		output.close();
	}
	return true;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::generateUniNetworks(NodeArray& nodelist,int numNode,int numEdge,std::string prename)
{
	for(int i=1;i<=NUM_RANDOM_NETWORK;i++)
	{
		EdgeArray edgelist;
		std::unordered_map<std::string,bool> uniedgemap;
		 srand (time(NULL));
		for(int m=0;m<numEdge;m++)
		{
			std::string edgecode;
			int k=rand()%numNode;
			int g=rand()%numNode;
			if(g<k){ int h=k;k=g;g=h;}
			edgecode.append(convert_num2str(k));edgecode.append(convert_num2str(g));
			if((g==k) || uniedgemap.find(edgecode)!=uniedgemap.end() )
			{
				m--;continue;
			}
			edgelist[m].startNode=nodelist[k];
			edgelist[m].endNode=nodelist[g];				
		}
		std::string filename("../crosslink/dataset/rawdata/random/");
		filename.append(prename);filename.append("_r");filename.append(convert_num2str(i));filename.append(".txt");
		std::ofstream output(filename.c_str());
		for(int m=0;m<numEdge;m++)
		{
			output << edgelist[m].startNode <<"\t" << edgelist[m].endNode <<"\t" << "0.9\n";
		}
		output.close();
	}
	return true;
}

template<typename GR, typename BP>
bool NetworkPool<GR,BP>::outputMaWIShNetworks()
{
	unsigned tSize=_graphSet.size();
	std::vector<std::string> prefixName;
	prefixName.push_back("1-human9606.txt.si1.mawish");
	prefixName.push_back("2-worm6239.txt.si1.mawish");
	prefixName.push_back("3-fly7227.txt.si1.mawish");
	prefixName.push_back("4-yeast4932.txt.si1.mawish");
	prefixName.push_back("7-ecoli562.txt.si1.mawish");
	for(unsigned i=0;i<tSize;i++)
	{
		GraphData* data=getGraph(i);
		std::string filename("../crosslink/dataset/rawdata/10022014/");
		filename.append(prefixName[i]);
		std::ofstream output(filename.c_str());
		for(EdgeIt e(*(data->g));e!=lemon::INVALID;++e)
		{
			Node startnode=data->g->u(e);
			Node endnode=data->g->v(e);
			output <<  (*(data->label))[startnode]<<","<< (*(data->label))[endnode]<<std::endl;
		}
		output.close();
	}
	return true;
}

template<typename GR, typename BP>
void NetworkPool<GR,BP>::outputMaWIShSimilarity(std::string folder)
{
	std::string inputname,outputname;
	inputname.append(folder);inputname.append("input_blast_evals.txt");
	outputname.append(folder);outputname.append("input_blast_mawish.txt");
	std::ifstream input(inputname.c_str());
	std::ofstream output(outputname.c_str());
	std::string line,protein1, protein2;
	double evalue,weight;
	while(std::getline(input,line))
	{
		std::stringstream lineStream(line);
		lineStream >> protein1 >> protein2 >> evalue;
		if(evalue<1e-125)
			weight=1.0;
		else
			weight=-log10(evalue)/125.0;
		output << protein1 <<"\t"<< protein2 <<"\t"<<weight<<std::endl;
	}
	input.close();
	output.close();
}


#endif //NETWORKPOOL_H_
