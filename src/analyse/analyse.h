/* analyse.h
Author: Jialu Hu
Data: 02.10.2013*/
#ifndef ANALYSE_H_
#define ANALYSE_H_
#include "analyse/subnetwork.h"
#include <stdio.h>
#include <vector>
#include <algorithm>

class Analyse
{
public:
	Analyse();
	~Analyse(){};
	std::unordered_map<std::string,std::string> protein_dip_uniprot_map;
	void translate(std::string,int);
	void readIdMap();
	void removeRedundant(std::string,int);
	bool checkRedundante(Subnetwork*, Subnetwork*);
	struct Compare_Sub
	{
		bool operator() (Subnetwork* sub1, Subnetwork* sub2) {return sub1->score > sub2->score;}
	} compare_obj;
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
		std::vector<Subnetwork*> sublist;
		for(unsigned i=0;i<filelist.size();i++)
		{
			filename1.clear();filename1.append(folder);filename1.append("species_");filename1.append(convert_num2str(k));filename1.append("/");
			filename1.append(filelist[i]);filename1.append(".txt");
			filename2.clear();filename2.append(folder);filename2.append("species_");filename2.append(convert_num2str(k));filename2.append("/");
			filename2.append(filelist[i]);filename2.append("_genelist.txt");
			Subnetwork* sub=new Subnetwork();
			sub->readsubnetwork(filename1,filename2);
			sublist.push_back(sub);
			//sub->writegenelist(filename2,protein_dip_uniprot_map);
		}
		// sort subnetworks in the list according to their score.
		std::stable_sort(sublist.begin(),sublist.end(),compare_obj);
		//remove redundant subnetworks
		for(std::vector<Subnetwork*>::iterator it1=sublist.begin();it1!=sublist.end();++it1)
		{
			std::vector<Subnetwork*>::iterator it2=it1+1;
			(*it1)->writegenelist(protein_dip_uniprot_map);
			while(it2!=sublist.end())
			{
				if(checkRedundante(*it1,*it2))
				{
					it2=sublist.erase(it2);
				}else
				{
					++it2;
				}
			}
		}
	}
}

bool Analyse::checkRedundante(Subnetwork* sub1, Subnetwork* sub2)
{
	std::unordered_map<std::string,bool> promap;
	unsigned i,j,intersectNum,minProteinNum;
	float rate=0.0;
	for(i=0;i<sub1->proteinlist.size();i++)
	{
		promap[sub1->proteinlist[i]]=true;
	}
	i=sub1->proteinlist.size();
	intersectNum=0;
	for(j=0;j<sub2->proteinlist.size();j++)
	{
		std::string protein=sub2->proteinlist[j];
		if(promap.find(protein)!=promap.end())
			intersectNum++;
	}
	j=sub2->proteinlist.size();
	if(i<j)minProteinNum=i;
	else minProteinNum=j;
	rate=intersectNum/static_cast<float>(minProteinNum);
	return rate>0.5;
}

void Analyse::removeRedundant(std::string folder, int speciesnum)
{
	std::vector<std::string> filelist;
	std::string filename1,filename2;
	filename1.append(folder);
	filename1.append("filenames.txt");
	std::ifstream input1(filename1.c_str());
	std::ifstream input2(filename2.c_str());
	std::string item;
	while (std::getline(input1,item))
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
			filename2.append(filelist[i]);filename2.append(".txt");
		}
	}
}
#endif
