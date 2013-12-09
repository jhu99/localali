/* module.h
Author: Jialu Hu
Data: 27.11.2013*/

#pragma once
#ifndef MODULE_H_
#define MODULE_H_

#include <string>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <sstream>
#include "macro.h"

class Module
{
private:
	enum Category {P,C,F};
public:
	Module(std::string folder)
	{
		resultFolder=folder;
	};
	~Module(){};
	std::string resultFolder;
	std::unordered_map<std::string,bool> consideredProteins;
	std::unordered_map<std::string,bool> unconsideredProteins;

	typedef struct GeneOntologyInfo
	{
		Category type;
		std::string goid;
		std::string description;
		double pvalue;
		double cpvalue;
		bool isAnnotated;
		std::vector<std::string> enrichedProteins;
		GeneOntologyInfo()
		{
			clear();
		}
		void clear()
		{
			type=P;
			goid="";
			description="";
			pvalue=0.0;
			cpvalue=0.0;
			isAnnotated=false;
			enrichedProteins.clear();
		}
	}GOInfo;

	std::vector<GOInfo> goCandidates;
	GOInfo currentGO;
	void readConsideredProteins(std::string&);
	void readUnconsideredProteins(std::string&,int);
	void readEnrichedGeneOntology(std::string&,int,std::unordered_map<std::string,bool>&);
	bool isEnrichedGO(GOInfo&,std::unordered_map<std::string,bool>&);
	bool isHigherLevel(std::string,std::unordered_map<std::string,bool>&);
	void outputPredictedFunction(std::ofstream&);
};

void Module::readConsideredProteins(std::string& line)
{
	std::stringstream streamline(line);
	std::string proteinid, otherid;
	streamline >> proteinid >> otherid;
	if(consideredProteins.find(proteinid)==consideredProteins.end())
		consideredProteins[proteinid]=true;
}

void Module::readUnconsideredProteins(std::string& line,int linenum)
{
	if(linenum==0) return;
	std::stringstream streamline(line);
	std::string proteinid, otherid;
	streamline >> proteinid >> otherid;
	if(unconsideredProteins.find(proteinid)==unconsideredProteins.end())
		unconsideredProteins[proteinid]=true;
}

void Module::readEnrichedGeneOntology(std::string& line,int linenum, std::unordered_map<std::string,bool>& ancestormap)
{
	std::stringstream streamline(line);
	std::string s1,s2,s3,s4;
	std::size_t found;
	switch(linenum)
	{
	case 1:
		streamline>>s1>>s2;
		currentGO.goid=s2;
		if(s2.compare("GO:XXXXXXX")!=0)
			currentGO.isAnnotated=true;
		break;
	case 2:
		currentGO.description=line;
		break;
	case 3:
		streamline >> s1 >> s2 >> currentGO.pvalue;
		break;
	case 4:
		streamline >> s1 >> s2 >> currentGO.cpvalue;
		break;
	case 9:
		if(!currentGO.isAnnotated)break;
		while(streamline.good())
		{
			streamline >> s1;
			found = s1.find_first_of(",");
			if(found!=std::string::npos)
			{
				s2 = s1.substr(0,found);
				currentGO.enrichedProteins.push_back(s2);
			}else
			{
				currentGO.enrichedProteins.push_back(s1);
			}
		}
		if(isEnrichedGO(currentGO,ancestormap))
		goCandidates.push_back(currentGO);
		currentGO.clear();
		break;
	default:
		break;
	}
}

bool
	Module::isEnrichedGO(GOInfo& currentGO,std::unordered_map<std::string,bool>& ancestormap)
{
	unsigned num=currentGO.enrichedProteins.size();
	float rate=0.0;
	rate=num/(1.0*consideredProteins.size());
	if(num>=NUM_PROTEIN_ENRICHED && rate>0.5 && isHigherLevel(currentGO.goid,ancestormap))
		return true;
	return false;
}

bool
	Module::isHigherLevel(std::string goid, std::unordered_map<std::string,bool>& ancestormap)
{
	std::string filename,goterm,line,commandline("./bin/go_ancestor_finder.sh ");
	filename.append("./dataset/goancestors/");
	filename.append(goid);
	filename.append(".ancestors");
	if(ancestormap.find(filename)==ancestormap.end())
	{
		commandline.append(goid);
		system(commandline.c_str());
		ancestormap[filename]=true;
	}
	std::ifstream input(filename.c_str());
	bool switcher=false;
	std::vector<std::string> ancestor_go;
	while(std::getline(input,line))
	{
		std::stringstream streamline(line);
		streamline >> goterm;
		if(goterm.compare("GO:0003673")==0)// GO:0003673 is the root of go tree.
			switcher=true;
		if(switcher)
			ancestor_go.push_back(goterm);
		if(ancestor_go.size()>=NUM_GOTERM_LEVEL)
			return true;
		if(goterm.compare(goid)==0)
		{
			switcher=false;
			ancestor_go.clear();
		}
	}
	input.close();
	return false;
}

void
	Module::outputPredictedFunction(std::ofstream& output)
{
	GOInfo currentGO=goCandidates.front();
	for(std::vector<std::string>::iterator it=currentGO.enrichedProteins.begin();it!=currentGO.enrichedProteins.end();++it)
	{
		consideredProteins[*it]=false;
	}
	for(std::unordered_map<std::string,bool>::iterator it=consideredProteins.begin();it!=consideredProteins.end();++it)
	{
		if(it->second)
			output << it->first <<"\t" << currentGO.goid << std::endl;
	}
	for(std::unordered_map<std::string,bool>::iterator it=unconsideredProteins.begin();it!=unconsideredProteins.end();++it)
	{
		output << it->first <<"\t" << currentGO.goid << std::endl;
	}
}
#endif
