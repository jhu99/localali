/* complex.h
Author: Jialu Hu
Data: 27.11.2013*/

#pragma once
#ifndef COMPLEX_H_
#define COMPLEX_H_

#include <iostream>
#include <fstream>
#include <string>
#include <unordered_set>
#include <cstddef> // std::size_t

class Complex
{
public:
	typedef struct _ComplexType
	{
		std::string complexId;
		std::string FunCatCategories;
		std::string organism;
		std::string uniprotId;
		std::string entrezId;
		std::string purificationMethod;
		std::string otherInfo;
		std::unordered_set<std::string> proteinlist;
	}ComplexType;
	std::unordered_map<std::string,ComplexType> complexMap;
	std::unordered_map<std::string,std::string> proteinKeyMap;
	Complex(){}
	~Complex(){}
	void readComplexes(std::string& filename);
};

void Complex::readComplexes(std::string& filename)
{
	std::ifstream input(filename.c_str());
	std::string line,item;
	std::getline(input,line);
	while(std::getline(input,line))
	{
		ComplexType complexInstance;
		std::size_t pos=line.find(";");
		int num=0;
		while(pos!=std::string::npos)
		{
			item=line.substr(0,pos);
			line=line.substr(pos+1);
			pos=line.find(";");
			switch(num++)
			{
				case 0:complexInstance.complexId=item;break;
				case 3:complexInstance.organism=item;break;
				case 4:complexInstance.uniprotId=item;break;
				case 5:complexInstance.entrezId=item;break;
				case 6:complexInstance.purificationMethod=item;break;
				case 8:complexInstance.FunCatCategories=item;break;
				default:complexInstance.otherInfo.append(item);break;
			}
		}
		if(complexMap.find(complexInstance.complexId)==complexMap.end())
		{
			line=complexInstance.uniprotId;
			std::size_t pos=line.find_first_of(",()");
			while(pos!=std::string::npos)
			{
				line[pos]=' ';
				pos=line.find_first_of(",()");
			}
			std::stringstream streamline(line);
			while(streamline.good())
			{
				streamline >> item;
				complexInstance.proteinlist.insert(item);
				proteinKeyMap[item]=complexInstance.complexId;
			}
			if(complexInstance.proteinlist.size()<3)continue;
			if(complexInstance.organism.compare("Human")==0 ||
			   complexInstance.organism.compare("Mouse")==0)
			complexMap[complexInstance.complexId]=complexInstance;
		}
	}
	std::cout << "Real protein complexes in CORUM: "<< complexMap.size() << std::endl;
}
#endif
