#!/bin/sh
ANALYS_GO=1
GO_DIR="./bin/GO-TermFinder-0.86"
ANALYSEGO="${GO_DIR}/examples/analyze.pl ./dataset/goa/gene_association.goa_combination 10000 dataset/goa/gene_ontology.1_2.obo"
if [ ${ANALYS_GO} -eq 1 ]
then
perl $ANALYSEGO ./result/dip/3-way/netblastm/species_0/complex_*_genelist.txt
perl $ANALYSEGO ./result/dip/3-way/netblastm/species_1/complex_*_genelist.txt
perl $ANALYSEGO ./result/dip/3-way/netblastm/species_2/complex_*_genelist.txt
fi
