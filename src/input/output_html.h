/* format.h
Author: Jialu Hu
Data: 10.04.2014*/
#ifndef OUTPUT_HTML_H_
#define OUTPUT_HTML_H_

#include <exception>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include "getpost.h"
#include "function.h"

using namespace std;

template<typename Option>
class Output_html
{
public:
	std::string jobid;
	Output_html()
	{
	}
	void set_header();
	void set_footer();
	void get_data(Option&);
	void update_table();
};

template<typename Option>
void Output_html<Option>::set_header()
{
	std::cout <<"Content-Type: text/html;charset=us-ascii\n\n";
}

template<typename Option>
void Output_html<Option>::set_footer()
{
	std::cout <<"</body></html>\n";
}

template<typename Option>
void Output_html<Option>::get_data(Option& myoption)
{
	std::unordered_map<std::string,std::string> formData;
	std::string numspecies,seedsize,algorithm,numspinetries,numseeds,threshold,seedtries,minext,maxext,extdist1,extdist2,seedrep,verbose;
	if(!getPost(formData))return;
	jobid=formData["jobid"];
	numspecies=formData["numspecies"];
	seedsize=formData["seedsize"];
	algorithm=formData["algorithm"];
	numspinetries=formData["numspinetries"];
	numseeds=formData["numseeds"];
	threshold=formData["threshold"];
	seedtries=formData["seedtries"];
	minext=formData["minext"];
	maxext=formData["maxext"];
	extdist1=formData["extdist1"];
	extdist2=formData["extdist2"];
	seedrep=formData["seedrep"];
	verbose=formData["verbose"];
	std::cout <<"Job id:\t"<<jobid<<"<br>"<<std::endl;
	std::cout <<"numspecies:\t"<<numspecies<<"<br>"<<std::endl;
	std::cout <<"seedsize:\t"<<seedsize<<"<br>"<<std::endl;
	std::cout <<"numseeds:\t"<<numseeds<<"<br>"<<std::endl;
	std::cout <<"threshold:\t"<<threshold<<"<br>"<<std::endl;
	std::cout <<"seedtries:\t"<<seedtries<<"<br>"<<std::endl;
	std::cout <<"minext:\t"<<minext<<"<br>"<<std::endl;
	std::cout <<"maxext:\t"<<maxext<<"<br>"<<std::endl;
	std::cout <<"extdist1:\t"<<extdist1<<"<br>"<<std::endl;
	std::cout <<"extdist2:\t"<<extdist2<<"<br>"<<std::endl;
	std::cout <<"seedrep:\t"<<seedrep<<"<br>"<<std::endl;
	std::cout <<"algorithm:\t"<<algorithm<<"<br>"<<std::endl;
	int num=stoi(numspecies);
	myoption.numspecies=num;
	myoption.seedsize=stoi(seedsize);
	myoption.numspinetries=stoi(numspinetries);
	myoption.numseeds=stoi(numseeds);
	myoption.score_threshold=stof(threshold);
	myoption.seedtries=stoi(seedtries);
	myoption.minext=stoi(minext);
	myoption.maxext=stoi(maxext);
	myoption.extdist1=stoi(extdist1);
	myoption.extdist2=stoi(extdist2);
	
	for(int i=0;i<num;i++)
	{
		std::string ppifilename="/var/www/html/mnetali/data-rw/uploadfiles/";
		ppifilename.append(jobid);ppifilename.append("/ppi");
		ppifilename.append(convert_num2str(i));
		ppifilename.append(".txt");
		myoption.networkfiles.push_back(ppifilename);
	}
	myoption.layerfile="/var/www/html/mnetali/data-rw/uploadfiles/";
	myoption.layerfile.append(jobid);
	myoption.treefile=myoption.layerfile;
	myoption.resultfolder=myoption.layerfile;
	myoption.layerfile.append("/input_blast.txt");
	myoption.treefile.append("/tree.txt");
	myoption.resultfolder.append("/");
}

#endif //OUTPUT_HTML_H_
