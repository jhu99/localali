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
	void extractIntActInteractions();
	void extractHomology();
	void extractIntActHomology(NetworksType& networks);
	void extractDatasetHomology(NetworksType& networks);
	void extractDipAc(NetworksType& networks);
	void writeAlignmentFile(std::string,std::string);// usefulless
	void writeSubnetworks(MyOption&);
	void partitionGOA(std::string);
	void generateAlignNemoPPI(std::string);
	void generateAlignNemoSim(std::string);
	int numspecies;
};

template<typename NetworksType,typename MyOption>
Format<NetworksType,MyOption>::Format(MyOption& myoption)
{
	numspecies=myoption.numspecies;
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::generateAlignNemoPPI(std::string inputfilename)
{
	std::ifstream input(inputfilename.c_str());
	std::string outputfilename(inputfilename);outputfilename.append(".nemo");
	std::ofstream output(outputfilename.c_str());
	std::string line,protein1,protein2;
	float intscore;
	while (std::getline(input,line))
	{
		std::stringstream linestream(line);
		linestream >> protein1 >> protein2 >> intscore;
		output << protein1 <<"\t" << protein2 <<"\tUSELESS\t" << intscore << std::endl;
	}
	input.close();
	output.close();
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::generateAlignNemoSim(std::string inputfilename)
{
	std::ifstream input(inputfilename.c_str());
	std::string outputfilename(inputfilename);outputfilename.append(".nemo");
	std::ofstream output(outputfilename.c_str());
	std::string line,protein1,protein2;
	double escore;
	while (std::getline(input,line))
	{
		std::stringstream linestream(line);
		linestream >> protein1 >> protein2 >> escore;
		output << protein1 <<"\t" << protein2 <<"\t1.000" << std::endl;
	}
	input.close();
	output.close();
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
void Format<NetworksType,MyOption>::extractIntActInteractions()
{
	std::string filename1,filename2;
	std::vector<std::string> filelist;			
	filelist.push_back("1-human9606");			filelist.push_back("2-worm6239");			filelist.push_back("3-fly7227");			filelist.push_back("4-yeast4932");
	filelist.push_back("5-rat10116");			filelist.push_back("6-mouse10090");			filelist.push_back("7-ecoli562");
	for(unsigned i=0; i<filelist.size(); ++i)
	{
		std::unordered_map<std::string,int> interactionmap;
		filename1.clear();
		filename1.append("./dataset/rawdata/10022014/");
		filename1.append(filelist[i]);
		filename1.append(".txt.si");
		filename2.clear();
		filename2.append("./dataset/rawdata/10022014/");
		filename2.append(filelist[i]);
		filename2.append(".txt.si1");
		std::ifstream input(filename1.c_str());
		std::ofstream output(filename2.c_str());
		std::string line;
		if(!input.is_open())
		{
			std::cerr << "Cannot open " <<filename1 <<"!" << std::endl;
		}
		int numedge=0;
		while (std::getline(input,line))
		{
			std::string protein1,protein2,tempstr,keystr;
			std::stringstream streamline(line);
			streamline >> protein1 >> protein2;
			if(protein1.compare(protein2)==0)
				continue;
			else if(protein1.compare(protein2)>0)
		     {
				 tempstr = protein1;
				 protein1 = protein2;
				 protein2 = tempstr;
			 }
			 keystr.append(protein1);
			 keystr.append(protein2);
			 if(interactionmap.find(keystr)==interactionmap.end())
			 {
				output << line << std::endl;
				interactionmap[keystr]=1;
				numedge++;
			 }
		}
		std::cout << "network " << i <<" : " <<numedge << std::endl;
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
void Format<NetworksType,MyOption>::extractIntActHomology(NetworksType& networks)
{
	std::string filename1,filename2,filename3,filename4,line,protein1,protein2;
	float bscore;
	double evalue;
	std::vector<std::string> filelist;			
	filelist.push_back("1-human9606");			filelist.push_back("2-worm6239");			filelist.push_back("3-fly7227");			filelist.push_back("4-yeast4932");
	filelist.push_back("5-rat10116");			filelist.push_back("6-mouse10090");			filelist.push_back("7-ecoli562");
	for(int i=0;i<7;i++)
	{
		for(int j=i;j<7;j++)
		{
			std::unordered_map<std::string,int> edgenum;
			filename1.clear();filename1.append("./dataset/bldata/blastdata10022014/");filename1.append(filelist[i]);filename1.append("-");filename1.append(filelist[j]);filename1.append(".b6.data");
			filename2.clear();filename2.append("./dataset/bldata/blastdata10022014/");filename2.append(filelist[i]);filename2.append("-");filename2.append(filelist[j]);filename2.append(".b6.data.si");
			filename3.clear();filename3.append("./dataset/bldata/blastdata10022014/");filename3.append(filelist[i]);filename3.append("-");filename3.append(filelist[j]);filename3.append(".bscore");
			filename4.clear();filename4.append("./dataset/bldata/blastdata10022014/");filename4.append(filelist[i]);filename4.append("-");filename4.append(filelist[j]);filename4.append(".evals");
			std::ifstream input(filename1.c_str());
			if(!input.is_open())
			{
				std::cerr <<"Can't open "<< filename1 <<"!"<<std::endl;
				return;
			}
			std::ofstream output(filename2.c_str());
			std::ofstream output1(filename3.c_str());
			std::ofstream output2(filename4.c_str());
			while(std::getline(input,line))
			{
				std::stringstream lineStream(line);
			    lineStream >> protein1 >> protein2 >> bscore >>evalue;
			    /// Proteins require to be available in PPI networks.
			    if(!(networks.existNode(protein1) && networks.existNode(protein2)))
			      continue;
			    if(protein1.compare(protein2)==0)continue;
			    if(protein1.compare(protein2)>0)
			    {
					std::string temp=protein1;
					protein1=protein2;
					protein2=temp;
				}
				std::string edgelabel;
				edgelabel.append(protein1);
				edgelabel.append(protein2);
				if(edgenum.find(edgelabel)!=edgenum.end())continue;
				edgenum[edgelabel]=1;
				output << line << std::endl;
				output1 << protein1<<"\t"<<protein2<<"\t"<<bscore<<std::endl;
				output2 << protein1<<"\t"<<protein2<<"\t"<<evalue<<std::endl;
			}
			input.close();
			output.close();
			output1.close();
			output2.close();
		}
	}
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

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::writeSubnetworks(MyOption& myoption)
{
	std::string filename,filename1,filename2,line;
	std::vector<std::string> filelist;
	filename.append(myoption.resultfolder);
	filename.append("alignmentfiles.txt");
	std::ifstream input(filename);
	if(!input.is_open())
	{
		std::cerr << filename <<" is not existed!" << std::endl;
	}
	while(getline(input,line))
	{
		std::string ss;
		std::stringstream linestream(line);
		linestream >> ss;
		filelist.push_back(ss);
	}
	for(unsigned i=0;i<filelist.size();i++)
	{
		writeAlignmentFile(filelist[i],myoption.resultfolder);
	}
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::partitionGOA(std::string filename)
{
	std::ifstream input(filename);
	std::ofstream outputs[NUM_GOA_PARTS];
	int linenum=0;
	std::string outfilename,line;
	for(int i=0;i<NUM_GOA_PARTS;i++)
	{
		outfilename="./dataset/goa/partition/gene_association.goa_target_";
		outfilename.append(convert_num2str(i));
		outputs[i].open(outfilename);
	}
	int mod=0;
	while(std::getline(input,line))
	{
		if(linenum < 49493) mod=0;
		else if(linenum < 2*49493) mod=1;
		else mod=2;
		outputs[mod] << line << std::endl;
		linenum++;
	}
	for(int i=0;i<NUM_GOA_PARTS;i++)
	{
		outputs[i].close();
	}
	input.close();
}

template<typename NetworksType,typename MyOption>
void Format<NetworksType,MyOption>::writeAlignmentFile(std::string filename, std::string folder)
{
	typedef std::unordered_map<std::string, bool> ProteinList;
	typedef ProteinList::iterator Iter;
	std::string line,start,ss,protein,outfilename,infilename;
	float score;
	std::vector<ProteinList*> subnetworks;
	infilename.append(folder);
	infilename.append("alignments/");
	infilename.append(filename);
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
		std::size_t found = line.find_first_of(",");
		while (found!=std::string::npos)
		{
			line[found]=' ';
			found=line.find_first_of(",",found+1);
		}
		std::stringstream linestream(line);
		if(!linenum)
		{
			linestream >> start >> ss >> score;
			linenum=true;
			continue;
		}
		for(int i=0; i<numspecies; i++)
		{ 
			linestream >> protein;
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
		outfilename.append(filename);
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
			output << it->first << std::endl;
		}
	}
	for(int i=0;i<numspecies;i++)
	{
		delete subnetworks[i];
	}
}
#endif
