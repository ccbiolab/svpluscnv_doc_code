#!/usr/bin/perl

use strict;

#my @ls = `ls ~/tmp/shatterproof/cn/`;
my @ls = `cat /Users/lopezg16/tmp/shatterproof/breast_samples.txt`;
for(@ls){
	chop($_);
	print "aqui $_\n";
	print "cp ~/tmp/shatterproof/cn/$_.spc  ~/tmp/shatterproof/cnbrca/$_.spc\n";
	`cp ~/tmp/shatterproof/cn/$_.spc  ~/tmp/shatterproof/cnbrca/$_.spc`;
	`cp ~/tmp/shatterproof/tra/$_.spt  ~/tmp/shatterproof/trabrca/$_.spt`;
	`mkdir ~/tmp/shatterproof/out/$_`;
	`perl -w /Users/lopezg16/install/Shatterproof-0.13/scripts/shatterproof.pl --cnv /Users/lopezg16/tmp/shatterproof/cnbrca/ --trans /Users/lopezg16/tmp/shatterproof/trabrca/ --config /Users/lopezg16/tmp/shatterproof/config_brca.pl --output  /Users/lopezg16/tmp/shatterproof/outbrca/$_`;
	print "perl -w /Users/lopezg16/install/Shatterproof-0.13/scripts/shatterproof.pl --cnv ~/tmp/shatterproof/cnbrca/ --trans ~/tmp/shatterproof/trabrca/ --config /Users/lopezg16/tmp/shatterproof/config_brca.pl --output  ~/tmp/shatterproof/outbrca/$_\n";
	`rm ~/tmp/shatterproof/cnbrca/$_.spc`;
	`rm ~/tmp/shatterproof/trabrca/$_.spt`;
	`rm -R ~/tmp/shatterproof/outbrca/$_/mutation_clustering`;
}


