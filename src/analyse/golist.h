/* golist.h
Author: Jialu Hu
Data: 23.07.2012*/
#ifndef GO_LIST_H_
#define GO_LIST_H_
#include <iostream>
#include <iomanip> 
#include <fstream>
#include <string>
#include <set>
#include <unordered_map>
#include <cassert>
#include <unordered_set>
#include "macro.h"
#include "stdlib.h"


class GoList
{
private:
  typedef std::unordered_set<std::string> GoTerms;
  typedef struct _GOntology
  {
    GoTerms MF;
    GoTerms BP;
    GoTerms CC;
  }GOntology;
  typedef std::unordered_map<std::string,GOntology> GOntologyMap;
public:
  GOntologyMap go_map;
  std::vector<std::string> associationfiles;
 
  int numMF,numBP,numCC;/// Keep the number of annotated proteins in alignment graph
  GoList()
  {
    numMF=0;numBP=0;numCC=0;
  }
  ~GoList()
  {
  }
  bool readGeneOntology(const char*);/// Read gene ontology file for one species.
  bool outputGOterm(std::ofstream&,std::string);
  bool goInitial();
  GoTerms& getGoTerm(std::ofstream&,GOntology&,short);

};

bool GoList::readGeneOntology(const char* filename)
{
  /// Format of input file should be like this: Acc GO_term Inferrence Ontologies 
  std::ifstream input(filename);
  std::string line;
  while(std::getline(input,line))
  {
    if(line[0]=='!')continue;
    std::array<std::string,4> gorecord;
    std::stringstream streamline(line);
    for(int i=0;i<4;i++)
    {
      if(streamline.good())
      {
        streamline >> gorecord[i];
        //std::cout << gorecord[i] <<"\t";
      }
      else break;
    }

    /// If GO_term was inferred from IEA or ISS, ignore it.
    if(gorecord[2].compare("IEA")==0 ||gorecord[2].compare("ISS")==0)
      continue; 
    GOntology& record=go_map[gorecord[0]];
    if(gorecord[3].compare("C")==0)
    {
       // record.CC.push_back(gorecord[1]); We don't write these proteins that only have Cell Component annotations.
      /// Cellular Component annotation is not needed in our analysis
    }
    else if(gorecord[3].compare("F")==0)
    {
      record.MF.insert(gorecord[1]);
    }
    else if(gorecord[3].compare("P")==0)
    {
      record.BP.insert(gorecord[1]);
    }
    else{}
  }
  input.close();
  return true;
}

bool GoList::goInitial()
{
  for(unsigned i=0; i<associationfiles.size(); i++)
  {
    readGeneOntology(associationfiles[i].c_str());
  }
  return true;
}



bool GoList::outputGOterm(std::ofstream& output,std::string protein)
{
  if(go_map.find(protein)==go_map.end())
    return false;
  output << ">"<<protein << std::endl;
  GOntology& go=go_map[protein];
  for(short i=0;i<2;i++)
  {
    GoTerms& goterm=getGoTerm(output,go,i);
    for(GoTerms::iterator it=goterm.begin(); it!=goterm.end(); it++)
    {
      output << *it <<" ";
    }
    output << std::endl;
  }
  return true;
}


GoList::GoTerms& GoList::getGoTerm(std::ofstream& output,GOntology& go, short i)
{
  GoTerms* goterm=0;
  switch(i)
  {
    case 0:goterm=&go.MF;output <<"MF:";break;
    case 1:goterm=&go.BP;output <<"BP:";break;
    case 2:goterm=&go.CC;output <<"CC:";break;
  }
  return *goterm;
}

#endif /// GO_LIST_H_
