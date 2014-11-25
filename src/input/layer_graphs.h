/* layer_graphs.h
Author: Jialu Hu
Date: 10.06.2012*/

#pragma once
#ifndef LAYER_GRAPHS_H
#define LAYER_GRAPHS_H

#include <vector>
#include <fstream>
#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <unordered_map>

template<typename GR, typename NP>
class Layer_graphs
{
private:
  
  typedef NP NetworkPool;
  typedef typename NetworkPool::EdgeArray EdgeArray;
  typedef typename NetworkPool::NodeArray NodeArray;
public:
  typedef typename NetworkPool::Graph Graph;
  
  TEMPLATE_GRAPH_TYPEDEFS(Graph);
  /// Labels of the nodes
  typedef typename Graph::template NodeMap<std::string> OrigLabelNodeMap;
    /// Verified Nodes in the search table. 
  typedef std::unordered_map<int, bool> VerifiedNodeMap;
  typedef std::vector<typename Graph::Node> ValidNodeSeq;
  typedef std::vector<int> DensitySeq;
  /// Mapping from labels to original nodes
  typedef std::unordered_map<std::string, typename Graph::Node> InvOrigLabelNodeMap;
  typedef std::unordered_map<std::string, int> EdgeNumMap; 

  Graph graph;
  OrigLabelNodeMap node2label;
  InvOrigLabelNodeMap label2node;
  ValidNodeSeq validnodes;
  VerifiedNodeMap validnodemap;
  DensitySeq density;
  EdgeNumMap edgenum;
  int nodeNum;
  int edgeNum;

  Layer_graphs();
  ~Layer_graphs(){};
  bool read(std::string&,NetworkPool&);
  void setConfiguration(Node&);
  void generateRandKlayer(NetworkPool&,std::string&);
};

template<typename GR, typename NP>
Layer_graphs<GR,NP>::Layer_graphs()
  :graph()
  ,node2label(graph)
  ,label2node()
  ,validnodes()
  ,density()
  ,edgenum()
  ,nodeNum(0)
  ,edgeNum(0)
{
}

template<typename GR, typename NP>
void
Layer_graphs<GR,NP>::setConfiguration(Node& node)
{
	validnodes.push_back(node);
	validnodemap[graph.id(node)]=true;
	density.push_back(1);// assign a probability density function for validnodes.
}

template<typename GR, typename NP>
bool
Layer_graphs<GR,NP>::read(std::string& filename,NetworkPool& networks)
{
  std::ifstream input(filename.c_str());
  if(!input.good())
  {
    std::cerr <<"# Can't open "<< filename <<"!"<<std::endl;
    return 0;
  }
  std::string line;
  std::string protein1, protein2;
  double evalue;
  graph.reserveNode(networks.allNodeNum);
  while(std::getline(input,line))
  {
    std::stringstream lineStream(line);
    lineStream >> protein1 >> protein2>> evalue;
    Node node1,node2;
    /// Proteins require to be available in PPI networks.
    if(!(networks.existNode(protein1) && networks.existNode(protein2)))
      continue;
    /// There are no edges for paralogous proteins.
    if(networks.getHost(protein1)==networks.getHost(protein2))
	  continue;
    if(protein1.compare(protein2)==0)continue;
    if(protein1.compare(protein2)>0)
    {
		std::string temp=protein1;
		protein1=protein2;
		protein2=temp;
	}
	std::string edgelabel;
	edgelabel.append(protein1);
	edgelabel.append(protein2);
	if(edgenum.find(edgelabel)!=edgenum.end())continue;
    if(label2node.find(protein1)!=label2node.end())
    {
      node1 = label2node[protein1];
      //exist node
    }else
    {
      node1 = graph.addNode();
      node2label.set(node1,protein1);
      label2node[protein1] = node1;
      //insert new node
      nodeNum++;
    }
    if(label2node.find(protein2)!=label2node.end())
    {
      node2=label2node[protein2];
      //exist node
    }else
    {
      node2=graph.addNode();
      node2label.set(node2,protein2);
      label2node[protein2]=node2;
      //insert new node
      nodeNum++;
    }
    graph.addEdge(node1,node2);
    edgenum[edgelabel]=1;
    edgeNum++;
  }
  if(g_verbosity>=VERBOSE_ESSENTIAL)
  {
	std::cout <<"# "<<filename <<" has been read successfully!<br>"<<std::endl;
	std::cout <<"# of proteins:"<< nodeNum<<"\t<br>"<<std::endl;
	std::cout <<"# of homologous pairs:"<< edgeNum <<"<br>"<<std::endl;
  }
  return true;
}

template<typename GR, typename NP>
void
Layer_graphs<GR,NP>::generateRandKlayer(NetworkPool& networks,std::string& randfilename)
{
	NodeArray nodelist;
	EdgeArray edgelist;
	int i=0;
	for(NodeIt nodeid(graph);nodeid!=lemon::INVALID;++nodeid)
	{
		nodelist[i]=node2label[nodeid];
		i++;
	}
	std::unordered_map<std::string,bool> uniedgemap;
	 srand (time(NULL));
	for(i=0;i<edgeNum;i++)
	{
		std::string edgecode;
		int k=rand()%nodeNum;
		int g=rand()%nodeNum;
		if(g<k){ int h=k;k=g;g=h;}
		edgecode.append(convert_num2str(k));edgecode.append(convert_num2str(g));
		if((g==k) || uniedgemap.find(edgecode)!=uniedgemap.end() || networks.getHost(nodelist[k])==networks.getHost(nodelist[g]))
		{
			i--;continue;
		}
		
		edgelist[i].startNode=nodelist[k];
		edgelist[i].endNode=nodelist[g];
	}
	std::string filename(randfilename);
	std::ofstream output(filename.c_str());
	for(int m=0;m<edgeNum;m++)
	{
		output << edgelist[m].startNode <<"\t" << edgelist[m].endNode <<"\t" << "1.0E-50\n";
	}
	output.close();
}

#endif
