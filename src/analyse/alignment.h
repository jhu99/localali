/* alignment.h
Author: Jialu Hu
Data: 02.10.2013*/
#pragma once
#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <unordered_map>
#include "algorithm/function.h"

class Alignment
{
private:
	//typedef NetworkPoolType NetworkType;
public:
	typedef std::unordered_map<std::string,int> AlignmentMap;
	Alignment();
	~Alignment(){}
	void readAlignment(std::string, std::string, int);
	void writeAlignmentFile(//NetworkType&,
							std::string,
							int,
							std::unordered_map<std::string,std::string>&,
							std::unordered_map<std::string,int>&);
	AlignmentMap alignmap;
	std::string alignmentfile;
	float score;
};


Alignment::Alignment()
{
}


void Alignment::readAlignment(std::string folder, std::string filename, int numspecies)
{
	alignmentfile=filename;
	filename=folder;filename.append("alignments/");filename.append(alignmentfile);
	std::ifstream input(filename);
	if(!input.is_open())
	{
		std::cerr << filename <<" does not exist!" << std::endl;
	}
	std::string line, start, ss, protein;
	bool linenum=false;
	while(std::getline(input,line))
	{
		std::size_t found = line.find_first_of(",");
		while (found!=std::string::npos)
		{
			line[found]=' ';
			found=line.find_first_of(",",found+1);
		}
		std::stringstream streamline(line);
		if(!linenum)
		{
			streamline >> start >> ss >> score;
			linenum=true;
			continue;
		}
		for(int i=0; i < numspecies; i++)
		{
			streamline >> protein;
			if(alignmap.find(protein)==alignmap.end())
			{
				alignmap[protein] = i;
			}
		}
	}
	input.close();
}


void Alignment::writeAlignmentFile(//NetworkType& networks,
									std::string folder,
									int numspecies,
									std::unordered_map<std::string,std::string>& idmap,
									std::unordered_map<std::string,int>& coveredProteinMap)
{
	typedef std::unordered_map<std::string, bool> ProteinList;
	typedef ProteinList::iterator Iter;
	std::string line,start,ss,protein,outfilename,infilename;
	float score;
	std::vector<ProteinList*> subnetworks;
	infilename.append(folder);infilename.append("alignments/");infilename.append(alignmentfile);
	std::ifstream input(infilename);
	if(!input.is_open())
	{
		std::cerr << infilename <<" is not existed!" << std::endl;
	}
	for(int i=0;i<numspecies;i++)
	{
		subnetworks.push_back(new std::unordered_map<std::string, bool>());
	}
	bool linenum=false;
	while(getline(input,line))
	{
		std::stringstream linestream(line);
		if(!linenum)
		{
			linestream >> start >> ss >> score;
			linenum=true;
			continue;
		}
		std::size_t found = line.find_first_of(",");
		while (found!=std::string::npos)
		{
			line[found]=' ';
			found=line.find_first_of(",",found+1);
		}
		for(int i=0; i<numspecies; i++)
		{ 
			linestream >> protein;
			if(coveredProteinMap.find(protein)==coveredProteinMap.end())
				coveredProteinMap[protein]=1;
			//int si=
			if(subnetworks[i]->find(protein)==subnetworks[i]->end())
			{
				(*subnetworks[i])[protein]=true;
			}
		}
	}
	for(int i=0; i<numspecies; i++)
	{
		outfilename.clear();
		outfilename.append(folder);
		outfilename.append("species_");
		outfilename.append(convert_num2str(i));
		outfilename.append("/");
		outfilename.append(alignmentfile);
		std::ofstream output(outfilename);
		ProteinList *subnetwork;
		subnetwork=subnetworks[i];
		if(!output.is_open())
		{
			std::cerr << outfilename <<" is not existed!" << std::endl;
		}
		output <<"# score: " << score << std::endl; 
		for(Iter it=subnetwork->begin();it!=subnetwork->end();++it)
		{
			if(idmap.find(it->first)!=idmap.end())
			output << idmap[it->first] << std::endl;
		}
		output.close();
	}
	input.close();
	for(int i=0;i<numspecies;i++)
	{
		delete subnetworks[i];
	}
}
#endif //ALIGNMENT_H_
