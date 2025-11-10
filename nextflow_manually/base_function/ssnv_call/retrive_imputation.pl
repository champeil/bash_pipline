#!/usr/bin/perl -w
use strict;
use File::Basename;
use warnings;

=head1
#====================================
# this script is for retrive the mis-filtered mutations by imputation
# usage: perl retrive_imputation.pl vep_annotated_maf_file imputed_maf_file outdir
# author: laojianpei
# organization: sysucc
# time: 2022.12.27
#====================================
=cut

# to judge the existence of parameters
my $raw_maf = shift or die("please input the raw_maf file\n$!");
my $imputated_maf = shift or die("please input the imputed_maf file\n$!");
my $outdir = shift or die("please input the outdir\n$!");

if( ! -e $outdir){
    `mkdir -p $outdir`;
}
my $filenames = basename($raw_maf);
$filenames =~ /(\S+)\.maf/;

my (%raw, %imputed, %retrive, %id);
open(RAW,"-|","sed -e 's/chr//g' '$raw_maf'") or die $!;
open(IMP,"<",$imputated_maf) or die $!;
open(RET,">","$outdir/$1_retrive.maf") or die $!;

my $header_raw = <RAW>;
$header_raw = <RAW> unless $header_raw !~ /^#/;
chomp $header_raw;
my $header_imputed = <IMP>;
chomp $header_imputed;
print RET $header_imputed."\tVAF\n";

my @ID = split(/\t/, $header_raw);
for(0..$#ID){
    $id{$ID[$_]} = $_ #dict: keys: header, values: order
}

while(<RAW>){
    chomp $_;
    my @F = split(/\t/,$_);
    $raw{$F[$id{"Tumor_Sample_Barcode"}]}{$F[$id{"Hugo_Symbol"}]}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
}
while(<IMP>){
    chomp $_;
    my @F = split(/\t/,$_);
    $imputed{$F[$id{"Tumor_Sample_Barcode"}]}{$F[$id{"Hugo_Symbol"}]}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
    print RET $_."\n";
}

my $n=0;
for my $k1(keys %raw){
    for my $k2(keys %{$raw{$k1}}){
        for my $k3(keys %{$raw{$k1}{$k2}}){
            for my $k4(keys %{$raw{$k1}{$k2}{$k3}}){
                if (exists($imputed{$k1}{$k2}{$k3}{$k4})) {
                    next;
                }
                else{
		    my @F = split(/\t/,$raw{$k1}{$k2}{$k3}{$k4});
		    my $vaf = int($F[$id{"t_alt_count"}] * 10000 / ($F[$id{"t_ref_count"}] + $F[$id{"t_alt_count"}])) / 10000;
                    print RET $raw{$k1}{$k2}{$k3}{$k4}."\t".$vaf."\n";
                    $n = $n+1;
                }
            }
        }
}
}

print "new added $n mutations\n";
close(RAW);
close(IMP);
close(RET);









