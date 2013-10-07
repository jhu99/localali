/* format.h
Author: Jialu Hu
Data: 30.09.2013*/
#ifndef FORMAT_H_
#define FORMAT_H_
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <unordered_map>

template<typename NetworksType,typename MyOption>
class Format
{
	public:
  Format(MyOption&);
  ~Format(){};
	std::string fromatfilename;
	void extractInteractions();
	void extractHomology();
	void extractDatasetHomology(NetworksType& networks);
	void extractDipAc(NetworksType& networks);
};

template<typename NetworksType,typename MyOption>
Format<NetworksType,MyOption>::Format(MyOption& myoption)
{
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::extractInteractions()
{
	std::string filename1,filename2;
	std::vector<std::string> filelist;			
	filelist.push_back("Celeg20130707");			filelist.push_back("Dmela20130707");			filelist.push_back("Ecoli20130707");			filelist.push_back("Hpylo20130707");
	filelist.push_back("Hsapi20130707");			filelist.push_back("Mmusc20130707");			filelist.push_back("Rnorv20130707");			filelist.push_back("Scere20130707");
	for(unsigned i=0; i<filelist.size(); ++i)
	{
		filename1.clear();
		filename1.append("./dataset/dip/dip-rawdata/20130707/");
		filename1.append(filelist[i]);
		filename1.append(".txt.data");
		filename2.clear();
		filename2.append("./dataset/dip/dip-rawdata/20130707/");
		filename2.append(filelist[i]);
		filename2.append("-int.txt");
		std::ifstream input(filename1.c_str());
		std::ofstream output(filename2.c_str());
		std::string line;
		if(!input.is_open())
		{
			std::cerr << "Cannot open " <<filename1 <<"!" << std::endl;
		}
		std::size_t found1,found2;
		while (std::getline(input,line))
		{
			found1 = line.find_first_of("(");
			found2 = line.find_last_of("\t");
			line.erase(line.begin()+found1,line.begin()+found2);
			found1 = line.find_first_of("(");
			line.erase(line.begin()+found1,line.end());

			std::string protein1,protein2,species1,species2,species3,species4;
			std::stringstream streamline(line);
			streamline >> protein1 >> protein2 >> species1 >> species2;
			if(protein1.compare("DIP-29165N")==0 && protein2.compare("DIP-29165N")==0)
				protein1="DIP-29165N";
			if(species1.compare(species2)==0 && protein1.compare(protein2)!=0)
			{
				output << protein1 <<"\t" << protein2 <<"\t" << 0.9 << "\n";
			}
		}
		input.close();
		output.close();
	}
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::extractHomology()
{
	std::string filename1,filename2,filename3;
	filename1.append("./dataset/dip/dip-rawdata/20130707/homology-list-20130826.b6.data");
	filename2.append("./dataset/dip/dip-rawdata/20130707/homology-list-20130826.evals");
	filename3.append("./dataset/dip/dip-rawdata/20130707/homology-list-20130826.bscore");
	std::ifstream input(filename1.c_str());
	std::ofstream output1(filename2.c_str());
	std::ofstream output2(filename3.c_str());
	std::string line,protein1,protein2;
	float bscore=0.0;
	double evalue=0.0;
	while (std::getline(input,line))
	{
		std::stringstream streamstr(line);
		std::string substr1,substr2;
		streamstr >> protein1 >> protein2 >> bscore >> evalue;
		if(protein1.compare(protein2)==0)continue;
		unsigned pos = protein1.find(":");
		substr1 = protein1.substr(pos+1);
		pos=protein1.find(":");
		substr2 = protein2.substr(pos+1);
		output1 << substr1 <<"\t"<< substr2 <<"\t"<< evalue<<"\n";
		output2 << substr1 <<"\t"<< substr2 <<"\t"<< bscore<<"\n";
	}
	output1.close();
	output2.close();
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::extractDatasetHomology(NetworksType& networks)
{
	std::ifstream input("./dataset/dip/dip-rawdata/20130707/homology-list-20130826.evals");
	std::string inputNet[]={"Celeg20130707","Dmela20130707","Ecoli20130707","Hpylo20130707","Hsapi20130707","Mmusc20130707","Rnorv20130707","Scere20130707"};
	std::string line;
	std::unordered_map<std::string,int> checklist;
	while(std::getline(input,line))
	{
		std::stringstream streamline(line);
		std::string protein1,protein2;
		double evalue;
		streamline >> protein1 >> protein2 >> evalue;
		unsigned i=networks.getHost(protein1);
		unsigned j=networks.getHost(protein2);
		if(i==100 || j==100) continue;
		if(i>j)
		{
			unsigned temp=i;
			std::string protein=protein1;
			i=j;
			j=temp;
			protein1=protein2;
			protein2=protein;
		}
		if(checklist.find(protein1)==checklist.end())
		{
		  checklist[protein1]=1;
		  std::cout << protein1 <<std::endl;
	    }
		if(checklist.find(protein2)==checklist.end())
		{
		  checklist[protein2]=1;
		  std::cout << protein2 <<std::endl;
	    }
		std::string outfilename("./dataset/dip/dip-rawdata/20130707/");
    outfilename.append(inputNet[i]);
		outfilename.append("-");
		outfilename.append(inputNet[j]);
		outfilename.append(".evals");
		std::ofstream output(outfilename.c_str(),std::ios_base::out|std::ios_base::app);
		output << protein1 <<"\t" << protein2 << "\t" << evalue << std::endl;
		output.close();		
	}
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::extractDipAc(NetworksType& networks)
{
	typedef std::unordered_map<std::string,short>::iterator TIterator;
	int i=0;
	std::ofstream output;
	for(TIterator it=networks.proteinHost.begin();it!=networks.proteinHost.end();++it,++i)
	{
		if((i%2000)==0)
		{
			if(i>0)output.close();
			std::string filename="./dataset/dip/dip-rawdata/20130707/proteinlist_";
			filename.append(convert_num2str(i/2000));
			filename.append(".txt");
			output.open(filename.c_str());
		}
		output << it->first << std::endl;
	}
}
#endif