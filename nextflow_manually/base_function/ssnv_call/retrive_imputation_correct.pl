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
# time: 2024.06.28
# ps: new corrected: 
#   cause the logic of retrive_imputation.pl is wrong, it check each line in raw whether exist in imput, if yes next, no will output the raw
#   based on the idea of imputation, it should be, check the exist of the mutation in raw of imputation, if exist, output the raw, not then output the imput
#   the imputation is based on all variants of the patients, include the sample itself, so need to judge the exist of the raw
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

print "########## start read file for $filenames ##########\n";
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

my $raw=0;
my $imp=0;
my ($raw_name,$imp_name);
while(<RAW>){
    chomp $_;
    my @F = split(/\t/,$_);
    $raw{$F[$id{"Hugo_Symbol"}]}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
    my $vaf = int($F[$id{"t_alt_count"}] * 10000 / ($F[$id{"t_ref_count"}] + $F[$id{"t_alt_count"}])) / 10000;
    print RET $_."\t".$vaf."\n";
    $raw=$raw+1;
    $raw_name=$F[$id{"Tumor_Sample_Barcode"}];
}
while(<IMP>){
    chomp $_;
    my @F = split(/\t/,$_);
    $imputed{$F[$id{"Hugo_Symbol"}]}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
    $imp=$imp+1;
    $imp_name=$F[$id{"Tumor_Sample_Barcode"}];
}
print "\twe got $raw from annotation maf and $imp from imputation\n";

print "########## then judge the exist of imputation on raw and add imputation which isn't exist ##########\n";

my $n_imp=0;
my $n_raw=0;
for my $k1(keys %imputed){
    for my $k2(keys %{$imputed{$k1}}){
        for my $k3(keys %{$imputed{$k1}{$k2}}){
            if (exists($raw{$k1}{$k2}{$k3})) { # exist in raw, the use the raw
                $n_raw = $n_raw+1;
            }
            else{ # or else use the imputated
		        my @F = split(/\t/,$imputed{$k1}{$k2}{$k3});
                $F[$id{"Tumor_Sample_Barcode"}]=$raw_name;
		        my $vaf = int($F[$id{"t_alt_count"}] * 10000 / ($F[$id{"t_ref_count"}] + $F[$id{"t_alt_count"}])) / 10000;
                print RET join("\t", @F)."\t".$vaf."\n";
                $n_imp = $n_imp+1;
            }
        }
}
}

print "\tthe imputed variants have $n_raw mutations exist in raw, and new added $n_imp mutations\n";
close(RAW);
close(IMP);
close(RET);









