# Overview #

We developed a fast and scalable local network alignment tool so-called LocalAli for the identification of functionally conserved modules in multiple networks. In this algorithm, we firstly proposed a new framework to reconstruct the evolution history of conserved
modules based on a maximum-parsimony evolutionary model. By relying on this model, LocalAli facilitates interpretation of resulting local alignments in terms of conserved modules which have been evolved from a common ancestral module through a series of evolutionary events. A meta-heuristic method simulated annealing was used to search for the optimal or near-optimal inner nodes (i.e. ancestral modules) of the evolutionary tree. To evaluate the performance and the statistical significance, LocalAli were tested on a total of 26 real datasets and 1040 randomly generated datasets. The results suggest that LocalAli outperforms all existing algorithms in terms of coverage, consistency and scalability, meanwhile retains a high precision in the identification of functionally coherent subnetworks.

# Features #

  * fast
  * scalable
  * accurate
  * evolutionary model
  * support parallelism

## Results ##

Quality alignments can be used to predict protein function by the method of annotation transfer. To validate the prediction results, we used a 10-fold cross-validation to assess the precision of protein function prediction.

![https://localali.googlecode.com/svn/trunk/images/success_rate.png](https://localali.googlecode.com/svn/trunk/images/success_rate.png)

We executed LocalAli on random datasets of each possible
combination of input species to verify the statistical significance of
our results. As a result, we found all these data about hits, FCS and
precision are non-random and statistically significant.

![https://localali.googlecode.com/svn/trunk/images/significance.png](https://localali.googlecode.com/svn/trunk/images/significance.png)

# Download #
  * [README](http://ftp.mi.fu-berlin.de/jhu/LocalAli/README.txt) (It describe how to compile the source code and run an example under Linux, Mac Os X and Windows.)
  * [dataset](http://ftp.mi.fu-berlin.de/jhu/LocalAli/dataset.tar.gz)
  * [lemon-1.2.3](http://ftp.mi.fu-berlin.de/jhu/LocalAli/lemon-1.2.3.tar.gz)
  * [localali](http://ftp.mi.fu-berlin.de/jhu/LocalAli/localali_linux_x86_64.tar.gz) (It is the executable code for Linux\_x86\_64.)
  * [application example](http://ftp.mi.fu-berlin.de/jhu/LocalAli/prediction/worm_fruitfly.txt) (new protein functions predicted by running LocalAli 20 times with th=0.5)

# Corrections #
  * Please find corrections of our published paper by clicking [here](https://code.google.com/p/localali/wiki/Corrections)!

# Others #
  * [NetCoffee](https://code.google.com/p/netcoffee/)
  * [NetworkBlast-M](http://www.cs.tau.ac.il/~bnet/License-nbm.htm)

# Contact Us #

[Jialu Hu](https://www.researchgate.net/profile/Jialu_Hu)

Email: Jialu.Hu{at}fu-berlin.de

# Reference #
  * Jialu Hu and Knut Reinert, [LocalAli: An Evolutionary-based Local Alignment Approach to Identify Functionally Conserved Modules in Multiple Networks, Bioinformatics (2015) 31 (3): 363-372 first published online October 4, 2014 doi:10.1093/bioinformatics/btu652](http://bioinformatics.oxfordjournals.org/content/31/3/363.full)
