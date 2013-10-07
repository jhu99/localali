/* analyse.h
Author: Jialu Hu
Data: 02.10.2013*/
#ifndef ANALYSE_H_
#define ANALYSE_H_
#include "analyse/subnetwork.h"
#include <stdio.h>

class Analyse
{
public:
	Analyse();
	~Analyse(){};
	std::unordered_map<std::string,std::string> protein_dip_uniprot_map;
	void translate(std::string,int);
	void readIdMap();
};

Analyse::Analyse(){}

void Analyse::readIdMap()
{
	std::string filename;
	std::ifstream input;
	std::string line,dipId,uniprotId;
	for(int i=0;i<=10;i++)
	{
		filename.clear();
		filename.append("./dataset/dip/dip-rawdata/20130707/proteinlist_");
		filename.append(convert_num2str(i));
		filename.append("_map.txt");
		input.open(filename.c_str());
		std::getline(input,line);
		while (std::getline(input,line))
		{
			input >> dipId >> uniprotId;
			protein_dip_uniprot_map[dipId]=uniprotId;
		}
		input.close();
	}
}

void Analyse::translate(std::string folder, int speciesnum)
{
	std::vector<std::string> filelist;
	std::string filename1,filename2;
	filename1.append(folder);
	filename1.append("filenames.txt");
	std::ifstream input(filename1.c_str());
	std::ofstream output(filename2.c_str());
	std::string item;
	while (std::getline(input,item))
	{
		std::string filename;
		unsigned pos = item.find(".");
		filename=item.substr(0,pos);
		filelist.push_back(filename);
	}
	for(int k=0;k<speciesnum;k++)
	{
		for(unsigned i=0;i<filelist.size();i++)
		{
			filename1.clear();filename1.append(folder);filename1.append("species_");filename1.append(convert_num2str(k));filename1.append("/");
			filename1.append(filelist[i]);filename1.append(".txt");
			filename2.clear();filename2.append(folder);filename2.append("species_");filename2.append(convert_num2str(k));filename2.append("/");
			filename2.append(filelist[i]);filename2.append("_genelist.txt");
			Subnetwork sub;
			sub.readsubnetwork(filename1);
			sub.writegenelist(filename2,protein_dip_uniprot_map);
		}
	}
}
#endif