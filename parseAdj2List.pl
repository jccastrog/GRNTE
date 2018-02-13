#!/usr/bin/perl
################################################################################
#Name:    parseAdj2List.pl
#Author:  Juan C. Castro <jcastro37@gatech.edu>
#         Diego M. Riaño P. <diriano@gmail.com>
#Update:  10-Feb-2018
#Version: 1.2
#License: GNU General Public License v3.0.
#===============================================================================
#
################################################################################
#================= 1.0 Load packages and initialize variables==================
use warnings;
use strict;
my $dataFile = shift @ARGV; #Initial file
my @values;
my @YEAST_ID;

#=================== 2.0 Open the file and run through it=====================#
open IN, "<",$dataFile or die "Cannot read $dataFile\n"; #Open initial file
while (<IN>) {
	chomp;
	if ($_ =~ /YPL230W\tYGL162W/) {
		@YEAST_ID = split('\t',$_);
	}
	else {
		my $line = $_;
		my (@fields) = split('\t',$line);
		my ($ID) = $fields[0];
		@values = split('\t',$line);
		@values = splice(@values,1);
		#@values = (1-@values );
		#===============================================================
		for (my $i=0; $i < ($#YEAST_ID+1); $i++) {
			if ($YEAST_ID[$i] eq $ID){
				print "";
			}else {
				print "$YEAST_ID[$i]\t$ID\t$values[$i]\n"; #print the values
			}
		}			
	}
}

#======================= 3.0 Close the file and exit =========================#
close IN; #Close file
#=============================================================================#
