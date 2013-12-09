/* analyse.h
Author: Jialu Hu
Data: 02.10.2013*/

#pragma once
#ifndef ANALYSE_H_
#define ANALYSE_H_
#include "analyse/subnetwork.h"
#include "analyse/alignment.h"
#include "analyse/module.h"
#include "analyse/golist.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include <iomanip>

class Analyse
{
	struct CrossVerification
	{
		std::vector<int> predictions;
		std::vector<int> corrections;
	};
public:
	Analyse();
	~Analyse(){};
	int numCoveredProtein;
	std::string resultFolder;
	std::unordered_map<std::string,std::string> protein_dip_uniprot_map;
	std::unordered_map<std::string,int> coveredProteinMap;
	std::vector<Subnetwork*> sublist;
	std::vector<Alignment*> alignmentlist;
	std::vector<Module> modulelist;
	void translate(std::string,int);// subnetworks and remove redundant alignments.
	void translate_alignment(std::string, int);// translate alignment and remove redundant alignments.
	void readIdMap();
	void removeRedundant(std::string,int);
	bool checkRedundante(Subnetwork*, Subnetwork*);
	bool checkRedundanteAli(Alignment*,Alignment*);
	void assessQuality(std::string,int);
	void predictFunction(std::string,int);
	void parseTermFile(std::ifstream&,std::unordered_map<std::string,bool>&);
	void countPrediction(std::string);
	void countVerification(std::string,CrossVerification&);
	void countCrossVerification(std::string,int);
	void verifyPrediction(std::string);
	
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

void Analyse::verifyPrediction(std::string folder)
{
	std::string filename1,filename2,filename3,commandline,parameter;
	GoList golist;
	std::unordered_map<std::string,bool> ancestormap;
	golist.readGeneOntology("./dataset/goa/gene_association.goa_target.fsst");
	filename1.append(folder);filename1.append("function_prediction.txt");
	filename2.append(folder);filename2.append("verify_prediction.txt");
	std::ifstream input(filename1),inputlist;
	std::ofstream output(filename2,std::ofstream::out);
	filename1.clear();
	filename1="./dataset/goancestors/list.txt";
	inputlist.open(filename1);
	std::string line;
	while(std::getline(inputlist,line))
	{
		std::stringstream streamline(line);
		ancestormap[line]=true;
	}
	inputlist.close();
	
	if(!input.is_open())
	{
		std::cerr << "Cann't open " <<filename1 <<".\n";return;
	}
	while(std::getline(input,line))
	{
		if(line[0] =='#')
		{
			output << line << std::endl;continue;
		}
		std::stringstream streamline(line);
		std::string protein, goid;
		streamline >> protein >> goid;
		bool exist=false;
		if(golist.go_map.find(protein)!=golist.go_map.end())
		{
			if(golist.go_map[protein].BP.find(goid)!=golist.go_map[protein].BP.end())
			{
				exist=true;
			}
			else if(!golist.go_map[protein].BP.empty())
			{
				parameter.clear();
				for(auto it=golist.go_map[protein].BP.begin();it!=golist.go_map[protein].BP.end();++it)
				{
					filename3.clear();filename3.append("./dataset/goancestors/");filename3.append(*it);filename3.append(".ancestors");
					if(ancestormap.find(filename3)==ancestormap.end())
					{
						commandline.clear();commandline.append("./bin/go_ancestor_finder.sh ");commandline.append(*it);system(commandline.c_str());ancestormap[filename3]=true;
					}
					parameter.append(filename3);parameter.append("  ");
				}
				parameter.append(" > ");
				parameter.append(folder); parameter.append("gotree.txt");
				commandline.clear();commandline.append("cat ");commandline.append(parameter);
				system(commandline.c_str());
				filename3.clear();filename3.append(folder);filename3.append("gotree.txt");
				inputlist.open(filename3);
				std::string line_gotree;
				while(std::getline(inputlist,line))
				{
					line_gotree.append(line);line_gotree.append(";");
				}
				inputlist.close();
				if(line_gotree.find(goid)!=std::string::npos)
					exist=true;
			}
		}
		if(exist)
		{
			output << "YES\t" << protein <<"\t" << goid << std::endl;
		}else
		{
			output << "NO\t" << protein << "\t" << goid << std::endl;
		}
	}
	output.close();
}

void Analyse::countCrossVerification(std::string folder,int methods)
{
	std::vector<CrossVerification> cvList;
	for(int i=0;i<3;i++)
	{
		std::string filename=folder;
		CrossVerification cv;
		if(methods==1) filename.append("localali_");
		else filename.append("netblastm_");
		filename.append(convert_num2str(i));
		if(methods==1) filename.append("/sample_5");
		else filename.append("/sample_0");
		filename.append("/verify_prediction.txt");
		countVerification(filename,cv);
		cvList.push_back(cv);
	}
	unsigned numspecies=0;
	numspecies=cvList[0].corrections.size();
	for(int i=0;i<3;i++)
	{
		CrossVerification cv=cvList[i];
		for(unsigned j=0;j<cv.corrections.size();j++)
		{
			std::cout << "species " <<j <<":" << cv.corrections[j] <<"\t" << cv.predictions[j] <<"\t" ;
			std::cout << std::setprecision(3) << 100*cv.corrections[j]/(1.0*cv.predictions[j])<< std::endl;
		}
	}
	std::cout <<" Summary: " << std::endl;
	for(unsigned j=0;j<numspecies;j++)
	{
		int corrections=0,predictions=0;
		for(int i=0;i<3;i++)
		{
			corrections+=cvList[i].corrections[j];
			predictions+=cvList[i].predictions[j];
		}
		std::cout << "species" << j << ":" << corrections <<"\t" << predictions << "\t" ;
		std::cout << std::setprecision(3) << 100*corrections/(1.0*predictions)<< std::endl;
	}
}

void Analyse::countVerification(std::string formatfile,CrossVerification& cv)
{
	std::ifstream input(formatfile.c_str());
	int species=-1,numCorrect=0,numPrediction=0;
	std::string yn,goid,protein,line;
	while(std::getline(input,line))
	{
		if(line[0]=='#')
		{
			if(species==-1)
			{
				species++;
				std::getline(input,line);
				continue;
			}
			std::getline(input,line);
			cv.predictions.push_back(numPrediction);
			cv.corrections.push_back(numCorrect);
			numCorrect=0;
			numPrediction=0;
			species++;
			continue;
		}
		std::stringstream streamline(line);
		streamline >> yn >> protein >> goid;
		if(yn.compare("YES")==0)
		{numCorrect++;numPrediction++;
		}
		else if(yn.compare("NO")==0)numPrediction++;
		else continue;
	}
	cv.predictions.push_back(numPrediction);
	cv.corrections.push_back(numCorrect);
}

void Analyse::countPrediction(std::string formatfile)
{
	std::ifstream input(formatfile.c_str());
	int numPredictions=0,numAllPredictions=0,species=-1;
	std::vector<int> numProteinsList,numPredictionsList;
	std::unordered_map<std::string,std::string> predictionMap,predictionAllMap;
	std::string line,protein,goid;
	while(std::getline(input,line))
	{
		if(line[0]=='#')
		{
			if(species==-1)
			{
				species++;
				std::getline(input,line);
				continue;
			}
			std::getline(input,line);
			std::cout << "Species "<< species <<": "<< numPredictions <<" predictions, " << predictionMap.size() <<" proteins" << std::endl;
			species++;
			numPredictions=0;
			predictionMap.clear();continue;
		}
		std::stringstream streamline(line);
		streamline >> protein >> goid;
		numPredictions++;numAllPredictions++;
		if(predictionMap.find(protein)==predictionMap.end())
			predictionMap[protein]=goid;
		if(predictionAllMap.find(protein)==predictionAllMap.end())
			predictionAllMap[protein]=goid;
	}
	std::cout << "Species "<< species <<": "<< numPredictions <<" predictions, " << predictionMap.size() <<" proteins" << std::endl;
	std::cout << "Alignment: "<< numAllPredictions <<" predictions, " << predictionAllMap.size() <<" proteins" <<std::endl;
	input.close();
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

void Analyse::predictFunction(std::string folder,int speciesnum)
{
	std::string filename,item;
	std::ifstream input,ingolist;
	std::vector<std::string> filelist;
	std::unordered_map<std::string,bool> ancestormap;
	resultFolder=folder;
	filename=(folder);filename.append("function_prediction.txt");
	std::ofstream output(filename.c_str());
	ingolist.open("./dataset/goancestors/list.txt");
	while(std::getline(ingolist,item))
	{
		ancestormap[item]=true;
	}
	ingolist.close();
	for(int k=0;k<speciesnum;k++)
	{
		filename.clear();filename.append(folder);filename.append("species_");filename.append(convert_num2str(k));filename.append("/termfile.txt");
		input.open(filename.c_str());
		filelist.clear();
		while (std::getline(input,item))
		{
			filelist.push_back(item);// number of reported hits
		}
		input.close();
		modulelist.clear();
		for(std::vector<std::string>::iterator it=filelist.begin();it!=filelist.end();++it)
		{
			filename=*it;
			input.open(filename.c_str());
			parseTermFile(input,ancestormap);
			input.close();
		}
		// output everything.
		output << "# Predicted function for species " << k <<":"<< std::endl;
		output << "# Functionally coherent modules: " << modulelist.size() << std::endl;
		for(std::vector<Module>::iterator it=modulelist.begin();it!=modulelist.end();++it)
		{
			it->outputPredictedFunction(output);
		}
	}
	output.close();
}

void Analyse::parseTermFile(std::ifstream& input, std::unordered_map<std::string,bool>& ancestormap)
{
	std::string line;
	int signal=0,myindex=0;
	Module mymodule(resultFolder);
	while(std::getline(input,line))
	{
		if(line.compare("The following gene(s) will be considered:")==0)
		{
			signal=1;myindex=0;continue;
		}
		else if(line.compare("The following gene(s) were not recognized, and will not be considered:")==0)
		{
			signal=2;myindex=0;continue;
		}
		else if(line.compare("Finding terms for P")==0)
		{
			signal=3;myindex=0;continue;
		}
		else if(line.compare("Finding terms for C")==0)
		{
			signal=4;myindex=0;continue;
		}
		else if(line.compare("Finding terms for F")==0)
		{
			signal=5;myindex=0;continue;
		}
		else if(line.compare("")==0)
		{
			continue;
		}
		switch (signal)
		{
		case 1:mymodule.readConsideredProteins(line);break;//readConsideredGenes(line,myindex);break;
		case 2:mymodule.readUnconsideredProteins(line,myindex);myindex++;break;
		case 3:mymodule.readEnrichedGeneOntology(line,myindex,ancestormap);myindex++;if(myindex==10)myindex=0;break;
		case 4:;break;
		default:
			break;
		}
	}
	if(!mymodule.goCandidates.empty())
		modulelist.push_back(mymodule);
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
					else if(line.compare("No terms were found for this aspect with a corrected P-value <= 0.01.")==0)
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
