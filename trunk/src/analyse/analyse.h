/* analyse.h
Author: Jialu Hu
Data: 02.10.2013*/

#pragma once
#ifndef ANALYSE_H_
#define ANALYSE_H_
#include "analyse/subnetwork.h"
#include "analyse/alignment.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <iomanip>

class Analyse
{
public:
	Analyse();
	~Analyse(){};
	int numCoveredProtein;
	std::unordered_map<std::string,std::string> protein_dip_uniprot_map;
	std::unordered_map<std::string,int> coveredProteinMap;
	std::vector<Subnetwork*> sublist;
	std::vector<Alignment*> alignmentlist;
	void translate(std::string,int);// subnetworks and remove redundant alignments.
	void translate_alignment(std::string, int);// translate alignment and remove redundant alignments.
	void readIdMap();
	void removeRedundant(std::string,int);
	bool checkRedundante(Subnetwork*, Subnetwork*);
	bool checkRedundanteAli(Alignment*,Alignment*);
	void assessQuality(std::string,int);
	struct Compare_Sub
	{
		bool operator() (Subnetwork* sub1, Subnetwork* sub2) {return sub1->score > sub2->score;}
	} compare_obj;
	struct Compare_Alignment
	{
		bool operator() (Alignment* align1, Alignment* align2) {return align1->score > align2->score;}
	} compare_ali;
};

Analyse::Analyse()
{}

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

void Analyse::translate_alignment(std::string folder, int speciesnum)
{
	std::vector<std::string> filelist;
	std::string filename1,filename2;
	filename1.append(folder);
	filename1.append("alignmentfiles.txt");
	std::ifstream input(filename1.c_str());
	std::string item;
	while (std::getline(input,item))
	{
		filelist.push_back(item);
	}
	input.close();
	for(unsigned i=0;i<filelist.size();i++)
	{
		filename2=filelist[i];
		Alignment *myalignment=new Alignment();
		myalignment->readAlignment(folder, filename2, speciesnum);
		alignmentlist.push_back(myalignment);
	}
	std::stable_sort(alignmentlist.begin(),alignmentlist.end(),compare_ali);
	for(std::vector<Alignment*>::iterator it1=alignmentlist.begin();it1!=alignmentlist.end();++it1)
	{
		std::vector<Alignment*>::iterator it2=it1+1;
		(*it1)->writeAlignmentFile(folder,speciesnum, protein_dip_uniprot_map,coveredProteinMap);
		// write subnetworks for it1;
		while(it2!=alignmentlist.end())
		{
			if(checkRedundanteAli(*it1,*it2))
			{
				it2=alignmentlist.erase(it2);
			}else
			{
				++it2;
			}
		}
	}
	numCoveredProtein=coveredProteinMap.size();
}

bool Analyse::checkRedundanteAli(Alignment* alig1, Alignment* alig2)
{
	int protein_num,num1,num2,conserved_num;
	std::string protein;
	float portion;
	num1=alig1->alignmap.size();num2=alig2->alignmap.size();
	if(num1<num2) protein_num=num1;
	else protein_num=num2;
	conserved_num=0;
	for(std::unordered_map<std::string,int>::iterator it=alig2->alignmap.begin();it!=alig2->alignmap.end();++it)
	{
		protein=it->first;
		if(alig1->alignmap.find(protein)!=alig1->alignmap.end())
			conserved_num++;
	}
	portion=conserved_num/static_cast<float>(protein_num);
	if(portion>0.5)
		return true;
	else
	{
		return false;
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

void Analyse::assessQuality(std::string folder,int speciesnum)
{
	std::string filename;
	int numDiscovered,numCoherent;
	float rate;
	std::vector<std::string> filelist;
	std::ifstream input1;
	std::string item,line,keystr,valstr;
	bool switcher1,switcher2,isCoherent;
	for(int k=0;k<speciesnum;k++)
	{
		std::unordered_map<std::string,bool> go_category;
		filename.clear();filename.append(folder);filename.append("species_");filename.append(convert_num2str(k));filename.append("/termfile.txt");
		input1.open(filename.c_str());
		filelist.clear();
		while (std::getline(input1,item))
		{
			filelist.push_back(item);
		}
		input1.close();
		numCoherent=0;
		numDiscovered=filelist.size();
		for(int i=0;i<numDiscovered;i++)
		{
			switcher1=false;
			switcher2=false;
			isCoherent=true;
			filename.clear();filename.append(filelist[i]);
			input1.open(filename.c_str());
			while(std::getline(input1,line))
			{
				if(!switcher1 && !switcher2 )
				{
					if(line.compare("Finding terms for P")==0)
					{
						switcher2=true;
						continue;
					}
					else if(line.compare("None of the gene names were recognized")==0)
					{
						isCoherent=false;
						break;
					}
					else
					{
						continue;
					}
				}
				else if(!switcher1 && switcher2)
				{
					if(line.compare("Finding terms for C")==0)break;
					else if(line.compare("No terms were found for this aspect with a corrected P-value <= 0.05.")==0)
					{
						isCoherent=false;
						break;
					}
					else
					{
						std::stringstream linestream(line);
						linestream >> keystr >> valstr;
						if(keystr.compare("GOID")==0 && go_category.find(valstr)==go_category.end())
						{
							go_category[valstr]=true;
						}
					}
				}
			}
			input1.close();
			if(isCoherent)numCoherent++;
		}
		rate=numCoherent/static_cast<float>(numDiscovered);
		std::cout <<"Species " << k <<":"<< std::endl;
		std::cout << "The number of coherent subnetworks they cover: " << numCoherent << std::endl;
		std::cout << "The number of distinct GO categories they cover: " << go_category.size() << std::endl;
		std::cout << std::setprecision(3) << "The percent of functionally coherent subnetworks discovered: "<< 100*rate <<"%" << std::endl;
	}
}
#endif
