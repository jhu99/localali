/* subnetwork.h
Author: Jialu Hu
Data: 02.10.2013*/
#ifndef SUBNETWORK_H_
#define SUBNETWORK_H_

#include <vector>
#include <fstream>

class Subnetwork
{
public:
	Subnetwork();
	~Subnetwork(){};
	float score;
	std::vector<std::string> proteinlist;
	void readsubnetwork(std::string);
	void writegenelist(std::string,std::unordered_map<std::string,std::string>&);
};

Subnetwork::Subnetwork():score(0.0),proteinlist()
{
}

void Subnetwork::readsubnetwork(std::string filename)
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
			continue;
		}else
		{
			linestream >> protein;
			proteinlist.push_back(protein);
		}
	}
}

void Subnetwork::writegenelist(std::string filename,std::unordered_map<std::string,std::string>& idmap)
{
	std::ofstream output(filename.c_str());
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