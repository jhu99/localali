/* subnetwork.h
Author: Jialu Hu
Data: 02.10.2013*/
#pragma once
#ifndef SUBNETWORK_H_
#define SUBNETWORK_H_

#include <vector>
#include <fstream>
#include <sstream>
#include <unordered_map>

class Subnetwork
{
public:
	Subnetwork();
	~Subnetwork(){};
	float score;
	std::string uniprotname;
	std::vector<std::string> proteinlist;
	void readsubnetwork(std::string,std::string);
	void readsubnetwork2(std::string);
	void writegenelist(std::unordered_map<std::string,std::string>&);
};

Subnetwork::Subnetwork():score(0.0),proteinlist()
{
}

void Subnetwork::readsubnetwork(std::string filename,std::string translatename)
{
	std::ifstream input(filename.c_str());
	std::string line;
	uniprotname=translatename;
	bool firstline=true;
	while(std::getline(input,line))
	{
		std::stringstream linestream(line);
		std::string start, protein;
		if(firstline)
		{
			linestream >> start >> protein >> score;
			firstline=false;
			continue;
		}else
		{
			linestream >> protein;
			proteinlist.push_back(protein);
		}
	}
}

void Subnetwork::readsubnetwork2(std::string filename)
{
	std::ifstream input(filename.c_str());
	std::string line;
	bool firstline=true;
	while(std::getline(input,line))
	{
		std::stringstream linestream(line);
		std::string start, protein;
		if(firstline)
		{
			linestream >> start >> protein >> score;
			firstline=false;
		}else
		{
			linestream >> protein;
			proteinlist.push_back(protein);
		}
	}
}

void Subnetwork::writegenelist(std::unordered_map<std::string,std::string>& idmap)
{
	std::ofstream output(uniprotname.c_str());
	output << "# Score: " << score << std::endl;
	for(unsigned i=0;i<proteinlist.size();i++)
	{
		std::string protein=proteinlist[i];
		if(idmap.find(protein)!=idmap.end())
		{
			output << idmap[protein] << std::endl;
		}
	}
	output.close();
}

#endif
