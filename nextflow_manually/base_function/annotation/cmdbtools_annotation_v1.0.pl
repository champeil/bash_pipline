#!/usr/bin/env perl -w
use strict;
use warnings;
use File::Basename;

=head1 introduction
# this script is for annotate the vcf file with cmdbtools
# author: laojp
# time: 2023.11.17
# position: SYSUCC bioinformatics platform
# version: 1.0
# desription: input the germline vcf file (or vcf file with the necessary column 1-8: #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO)
# process: 
#   split the vcf file according to the alt
#   annotate each line with cmdbtools
#   merge the variants with the cmdbtools annotated
# usage
#   process_vcf_alt_split(vcf,output_vcf)
#   process_cmdb_annotate(vcf,output_vcf)
#   process_merge_cmdb_vcf(vcf,cmdb_annotated_vcf,output_vcf)
=cut

sub process_vcf_alt_split {
    # first read the parameter
    my $input = $@[0] or die "Please input the vcf file"
    my $output = $@[1] or die "Please input the output vcf file"
    open (IN, "<", $input) or die "$!";
    open (OUT, ">", $output) or die "$!";

    # then read the file, found #CHROM and get header, and read lines to judge the replicate, split and output
    my (%header, %main);
    while(<IN>){
        chomp $_;
        if ($_ ~= /^CHROM/) {
            my @F = split(/\t/,$_);
            for (my $i = 0; $i <= $#F; $i++) { # header hash {key: name, value: index}
                $header{$F[$i]} = $i;
            }
        }
        last;
    }
    while(<IN>){
        chomp $_;
        my @F = split(/\t/,$_);
        my @ALT = split(/,/,$F[$header{"ALT"}]);
        if
    }
}