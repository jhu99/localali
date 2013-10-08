#!/bin/bash
########################################################
# 1. Perform all-against-all pairwise sequence alignment of DIP sequences using BLASTP.
# 2. ### 2.1 Format the raw DIP interaction data to a simple framework. For example: DIP-1081N|uniprotkb:P04273	DIP-339N|refseq:NP_313151|uniprotkb:P0A6F5	 -	-	-	-	MI:0019(coimmunoprecipitation)	-	pubmed:8676499|pubmed:DIP-98S	taxid:10036(Mesocricetus auratus)	taxid:83333(Escherichia coli K12)	MI:0218(physical interaction)	MI:0465(dip)	DIP-58E	dip-quality-status:core	dip:0002(small scale)		- -> DIP-1081N DIP-339N ### 2.2 Format BLASTP output file. For example: gnl|dip|DIP-1N|ref|NP_113971||sp|P19527	gnl|dip|DIP-1N|ref|NP_113971||sp|P19527	1072	0.0  --> DIP-1N	DIP-1N 1072 0.0
# 3. Assign each interactions of DIP a confident value, e.g. 0.9.
# 4. Rename the names of DIP interaction files (input-int-*.txt) and extract homologyfile (input-blast.txt), so that it can run on networkblast, graemlin, cappi, and networkblast-m.
# 6. RUN NETWORKBLASTM on large data.
########################################################
OPEN_BLASTP_ORTHOLOGY=0 # step 1
OPEN_R_FORMAT=0 # step 2
OPEN_WEIGHT_INTERACTION=0 # step 3
OPEN_FILTER_HOMOLOGY=0 # step 4
OPEN_RENAME=0 # step 5
########################################################
RUN_NETWORKBLASTM=0 # step 6
RUN_NETWORKBLAST=1 # step 7
RUN_CAPPI=0 # step
########################################################
BLASTP=~/software/ncbi-blast-2.2.28+/bin/blastp
DIPDATA=~/RAID/localali/dataset/dip/dip-rawdata/20130707
CURRENTFOLDER=$(pwd)
SPECIES=(Celeg20130707 Dmela20130707 Ecoli20130707 Hpylo20130707 Hsapi20130707 Mmusc20130707 Rnorv20130707 Scere20130707)
########################################################
########################################################
########################################################
### step 1 Perform all-against-all pairwise sequence alignment of DIP sequences using BLASTP.
### fasta***.seq -> homology-list-***.b6
if [ ${OPEN_BLASTP_ORTHOLOGY} -eq 1 ] 
then
	#echo ${DIPDATA}/fasta20130901.seq
	${BLASTP} -query ${DIPDATA}/fasta20130826.seq -subject ${DIPDATA}/fasta20130826.seq -outfmt "6 qacc sacc bitscore evalue" -evalue 1.0e-7 -out ${DIPDATA}/homology-list-20130826.b6
fi
########################################################
### step 2 Format the raw DIP interaction data to a simple framework.
### Celeg20130707.txt -> Celeg20130707.txt.data
### homology-list-***.b6 -> homolog-list-***.b6.data
if [ ${OPEN_R_FORMAT} -eq 1 ]
then
	#echo ${CURRENTFOLDER}
	R -e "setwd(\"${CURRENTFOLDER}\");source(\"./bin/parser.R\"); filelist<-list.files(path=\"./dataset/dip/dip-rawdata/20130707\",pattern=\"*txt\",full.names=TRUE); format_dip(filelist);format_dip_blast(\"./dataset/dip/dip-rawdata/20130707/homology-list-20130826.b6\")";
fi
########################################################
### step 3  First, convert Celeg20130707.txt.data etc.. -> Celeg20130707-int.txt (interactorA	interactorB 0.9), discarding interactions belonged to different species.
### Convert homology-list-20130826.b6.data -> homology-list-20130826.evals and homology-list-20130826.bscore.
if [ ${OPEN_WEIGHT_INTERACTION} -eq 1 ]
then
	./bin/localali -format
	./bin/localali -format -task 1
fi
########################################################
### step 4 Filter homology list for different datasets. homology-list-20130826.evals -> celeg-celeg.txt celeg-dme.txt ....
if [ ${OPEN_FILTER_HOMOLOGY} -eq 1 ]
then
	./bin/localali -format -task 2
fi
########################################################
### step 5 Rename the names of DIP interaction files and homologyfile, so that it can run on networkblast, cappi, and networkblast-m.
if [ ${OPEN_RENAME} -eq 1 ]
then
	rm ${CURRENTFOLDER}/dataset/dip/3-way/input*
	j=0
	for i in 0 1 7
	do
		ln -s ${DIPDATA}/${SPECIES[$i]}-int.txt  ${CURRENTFOLDER}/dataset/dip/3-way/input_int-$j.txt
		j=$(($j+1))
		for k in 0 1 7
		do
			if [ $k -gt $i ]
			then
			cat ${DIPDATA}/${SPECIES[$i]}-${SPECIES[$k]}.evals >> ${CURRENTFOLDER}/dataset/dip/3-way/input_blast.txt
			fi
		done
	done
fi
########################################################
### Run NETWORKBLASTM ON TEST DATASET.
if [ ${RUN_NETWORKBLASTM} -eq 1 ]
then
	./bin/bnm-local.sh 3 3-way >> ./result/measure_time_local.txt
fi
########################################################
### Run NETWORKBLAST ON TEST DATASET.
if [ ${RUN_NETWORKBLAST} -eq 1 ]
then
	./benchmark/networkblast/networkblast -i ./dataset/dip/3-way/input -o ./result/dip/3-way/networkblast/output
fi
