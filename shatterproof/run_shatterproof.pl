#!/usr/bin/perl
# Author: Gonzalo Lopez, PhD
# e-mail: gonzolgarcia@gmail.com; gonzalo.lopezgarcia@mssm.edu

use strict;

# change dir to local path
chdir "~/git/svpluscnv_doc_code/"; 

# obtain the list of samples
my @ls = `cat shatterproof/list_of_samples.txt`;

for(@ls){
	chop($_);
	print "$_\n";

	# cp single sample files to tmp directories and create output directorie
	`cp shatterproof/cn/$_.spc  shatterproof/cntmp/$_.spc`;
	`cp shatterproof/tra/$_.spt  shatterproof/tratmp/$_.spt`;
	`mkdir shatterproof/out/$_`;
	# run shatterproof
	`perl -w ~/install/Shatterproof-0.13/scripts/shatterproof.pl --cnv shatterproof/cntmp/ --trans shatterproof/tratmp/ --config shatterproof/config.pl --output  shatterproof/out/$_`;
	# remove temporary cn and tra files
	`rm shatterproof/cntmp/$_.spc`;
	`rm shatterproof/tratmp/$_.spt`;
	# remove big log files to save space
	`rm -R shatterproof/out/$_/mutation_clustering`; 
}


