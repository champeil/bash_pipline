#!/usr/bin/perl -w
use strict;
use File::Basename;
use warnings;

=head1
#====================================
# this script is for modify the peaks regions, cause the MACS will cross over the genome regions when extend regions
# usage: perl modify_peak_region.pl [input_narrowPeak_file] [genome_fai] [output_dir]
# author: laojp
# organization: sysucc
# time: 2024.06.28
#====================================
=cut

# to judge the existence of parameters
my $raw_peak = shift or die("please input the narrowpeak file\n$!");
my $genome = shift or die("please input the genome fai file\n$!");
my $outdir = shift or die("please input the out dir\n$!");

if( ! -e $outdir){
    `mkdir -p $outdir`;
}

my $filenames = basename($raw_peak);
$filenames =~ /(\S+)\.narrowPeak/;
open(MOD,">","$outdir/$1_modify.narrowPeak") or die $!;
open(GEN,"<",$genome) or die $!;
open(NAR,"<",$raw_peak) or die $!;


print "########## first input the fai file and obtain the length of genome ##########\n";
my %genome;
while(<GEN>){
    chomp $_;
    my @F = split(/\t/,$_);
    $genome{$F[0]} = $F[1];
}
close(GEN);

print "########## then read narrowPeak and judge the region ##########\n";
while(<NAR>){
    chomp $_;
    my @F = split(/\t/,$_);
    if ($F[1]<0){
        print "\t $F[0]-$F[1]:$F[2] outbound of the start of genome \n";
        $F[1]=0;
    }
    if($F[2]>$genome{$F[0]}){
        print "\t $F[0]-$F[1]:$F[2] outbound of the end of genome \n";
        $F[2]=$genome{$F[0]};
    }
    print MOD join("\t", @F)."\n";
}
close(NAR);
close(MOD);

print "########## finished ##########\n";