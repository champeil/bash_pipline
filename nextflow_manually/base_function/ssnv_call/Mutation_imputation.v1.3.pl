#!/usr/bin/perl -w
=head1 
#===============================================================================
#        USAGE: perl Mutation_imputation.pl <Maf_file> [bamfile] [outdir]
#
#  DESCRIPTION: Imputate mutation sites of samples which come from the same patient (Muti-sampling)
#
#  INPUT FILES: Merged mutation maf file
#
# REQUIREMENTS:
#        NOTES: None
#       AUTHOR: Xiaolong Zhang, zhangxiaol1@sysucc.org.cn
# ORGANIZATION: Bioinformatics Center, Sun Yat-sen University Cancer Center
#      VERSION: 1.3
#      CREATED: //2021-09-03
#     REVISION: ---
#
# PS for user laojianpei
# 	1. pattern of tumor_sample_barcode need to modify in line 58 and line 161
# 	2. the Chromosome column of maf file need to modify in line 76
# 	3. need to modify samtools in line 81
# 	4. maf_file include all filtered mutations of tissues in the same patients
# 	5. bam_file is the bqsr bam of tissues need to be imputation 
#	6. output imput_raw_maf file contains the raw records of each tested mutations
#	7. the next step is retrive imputed mutations with imput.maf file and combine with the raw maf files
#===============================================================================
=cut
#die `pod2text $0` unless @ARGV == 2;
use strict;
use File::Basename;

print "########## set min_depth as $min_depth and set min_alt_reads as $min_alt_reads ##########\n";
my $min_depth = 10;
my $min_alt_reads = 2;

print "########## read files and create the output file ##########\n";
my $maf = shift or die $!;
my $bamfile = shift or die $!;
my $outdir = shift or die $!;
my (%info, %mut, %type, %ref, %alt, %sample, %id);
my $file_name = basename($bamfile);
$file_name =~ /(\S+)_bqsr\.bam/;

open IN, $maf or die $!;
open OUT, ">$outdir/$1.imputed.maf" or die $!;
open OUT2, ">$outdir/$1.imputed_raw.maf" or die $!;
my $header = <IN>;
chomp $header;
print OUT $header."\tVAF\n";
print OUT2 $header."\tVAF\n";
chomp $header;
my @ID = split /\t/, $header; # read the header of maf file: key: header---value: position
for (0..$#ID){
	$id{$ID[$_]} = $_;
}

while(<IN>){ 
	chomp;
	my @F = split /\t/;
	$F[$id{"Tumor_Sample_Barcode"}] =~ /(LLM-.*)/; # the pattern of the patient-sample, get the Tumor_Sample_Barcode name according to the pattern and save to $1
	$info{$1}{$F[$id{"Tumor_Sample_Barcode"}]}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_; # info: key: Tumor_Sample_Barcode-Tumor_Sample_Barcode-Chromosome-Start_Position---value: $_
	$mut{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $_;
	$type{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $F[$id{"Variant_Type"}];
	$ref{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $F[$id{"Reference_Allele"}];
	$alt{$1}{$F[$id{"Chromosome"}]}{$F[$id{"Start_Position"}]} = $F[$id{"Tumor_Seq_Allele2"}];
	$sample{$1}{$F[$id{"Tumor_Sample_Barcode"}]} = (); # sample: key: Tumor_Sample_Barcode-Tumor_Sample_Barcode---value: empty
}
$file_name =~ /(\S+)_bqsr\.bam/;
`mkdir -p $outdir/Binary_mat_$1`;

print "########## start imputation ##########\n";
my (%binary, %count);
for my $k1(sort keys %sample){ #k1 stores each sample names
	my $m = keys %{$sample{$k1}};
	my $n = 0;
	for my $k2(sort keys %{$sample{$k1}}){ #k2 stores each sample names 
		$n ++;
		for my $k3(1..22,"X","Y"){ # need to modify the colunm of Chromosome in maf file
			next unless exists $mut{$k1}{$k3}; # choose the start_position exists variants
			for my $k4(sort {$a<=>$b} keys %{$mut{$k1}{$k3}}){ # with sample chr and position as key, $a and $b is temped for storing the compared position, and sort to compare the position and sort.
				print "\t########## $k1-$k2-$k3-$k4 started ##########\n";
				my @F = split /\t/, $mut{$k1}{$k3}{$k4}; # k4 stores start_position
				my $line = `samtools mpileup -A -x -B -d 10000 -q 20 -Q 20 -r chr$k3:$k4-$k4 -f /home/laojp/database/human/GATK/hg38/Homo_sapiens_assembly38.fasta $bamfile`; # need to change the reference file dir, samtools to return the result as line
				$binary{$k1}{"$k3:$k4-$F[0]"} .= "\tNA" and next unless $line;  # if no mapped reads, will return empty line. sample chr position - symbol, if no results in mpileup, then set sample chr:pos-symbol as NA and next loop
				my @G = split /\t/, $line; # mapping state from mpileup
				my ($mutant, $alt_read_count);
				if ($F[$id{"Variant_Type"}] eq "SNP"){ #line 75, maf line
					$mutant = $F[$id{"Tumor_Seq_Allele2"}];
					$alt_read_count = () = $G[4] =~ /($mutant)/ig; #mapping reads num
				}elsif($F[$id{"Variant_Type"}] eq "INS"){
					$mutant = $F[$id{"Tumor_Seq_Allele2"}];
					$mutant =~ s/^$F[$id{"Reference_Allele"}]//;
					$alt_read_count = () = $G[4] =~ /($mutant)/ig; #fitst overall and ignore the upper and lower write and search variants base in mpileup reads mapping states column, and then send mapping results to an array, then send the array to a variants , obtain the times the mutant appears.  
				}else{
					$alt_read_count = () = $G[4] =~ /\*/g; #ambigous base(two or more types of bases)
				}
				my $ref_read_count = () = $G[4] =~ /(\.|,)/ig; # match with + and - strand
				$binary{$k1}{"$k3:$k4-$F[0]"} .= "\tNA" and next if $alt_read_count + $ref_read_count == 0;  ## Special type SNV: DEL, DNP and so on
				my $vaf = int($alt_read_count * 10000 / ($alt_read_count + $ref_read_count)) / 10000;
				print "\t########## $k1-$k2-$k3-$k4\t$alt_read_count\t$ref_read_count\t$vaf ##########\n";
				$F[$id{"t_ref_count"}] = $ref_read_count;
				$F[$id{"t_alt_count"}] = $alt_read_count;
				$F[$id{"Tumor_Sample_Barcode"}] = $k2;
				if ($alt_read_count < $min_alt_reads and $G[3] <= $min_depth){  ## Reads depth is not enough or no difference between cancer and ref
					print "\t\t########## $k1-$k2-$k3-$k4\tabandon cause the not enough reads depth or not enough difference between ref and alt ##########\n";
					$binary{$k1}{"$k3:$k4-$F[0]"} .= "\tNA"; #.= connect string
					for my $i(@ID){
						print OUT2 "$F[$id{$i}]\t"; # do not change anything and output
					}
					print OUT2 "$vaf\n"; # output the new vaf
					next;
				}elsif($alt_read_count < $min_alt_reads and $G[3] > $min_depth){  ## No mutations
					print "\t\t########## $k1-$k2-$k3-$k4\tabandon cause the not enough count of alt, so no mutation ##########\n";
					$binary{$k1}{"$k3:$k4-$F[0]"} .= "\t0";
					for my $i(@ID){
						print OUT2 "$F[$id{$i}]\t";
					}
					print OUT2 "$vaf\n";
					next;
				}elsif($alt_read_count >= $min_alt_reads){  ## Inputed mutations
					print "\t\t########## $k1-$k2-$k3-$k4\taccepted as inputed mutations ##########\n";
					$binary{$k1}{"$k3:$k4-$F[0]"} .= "\t1";
					$count{$k1}{"$k3:$k4-$F[0]"} += 1 * 10**($m-$n);  ## Score this mut_site of this sample, to distinguish which sample has imputed mutations, 1 means gets, 0 means none, position means tissues name
				}
				for my $i(@ID){
					print OUT "$F[$id{$i}]\t"; 
					print OUT2 "$F[$id{$i}]\t";
				}
				print OUT "$vaf\n"; # OUT is the imputed mutations
				print OUT2 "$vaf\n"; # OUT2 is the record of all tested mutations
			}
		}
	}
}

print "########## end imputation ##########\n";

###############代码逻辑##################
#根据maf文件提取样本名称、突变染色体：位置：基因：突变碱基
#使用samtools mpileup对每一个样本的每个突变在bam文件中提取该位置的比对情况，并且对ref、alt reads数重新计算，并计算VAF值，用于添加在maf文件最后一列，mpileup里面的信息有：染色体、位置、参考碱基、测序深度、每条reads在该位点的碱基比对结果、每条reads在该位点的碱基质量
#根据maf文件中突变类型进行判断，只处理SNP与IND类型的突变，其他则计算模糊碱基的数量作为突变alt数（*）
#根据重新计算的ref、alt数进行判断，如果满足alt大于min_alt_reads且位点覆盖reads数（测序深度）大于min_depth，则定义该样本的该位置的突变为imputed mutation，否则分为非突变点与测序深度不足
#如果是imputed，则在binary中对应的样本处写作1，否则写作0（非突变）与NA（测序深度不足），此处记录方式为使用字典在样本、染色体-位置-基因名处指定字符，binary的值为字符串，里面记录了\t相隔的imputation的情况
#如果是imputed，则在count中记录为1，否则记录为0
#输出重新计算的ref、alt值以及输出VAF值到新的maf文件中，maf文件中所有类型的突变原封不动输出到新的maf文件（除了ref、alt count以及新增vaf列）
###############代码逻辑##################

###############输出结果##################
#OUT也就是imputed.maf为maf在bam中为imputation的那部分突变信息并重新计算ref、alt count、vaf
#OUT2也就是imputed_raw.maf为maf每条突变信息在bam文件中的位点情况，包括imputation、no mutation、no enough depth
#Binary里面记录了每一个样本查到的含有突变的区域（记录为1），其中每一个bam文件都会有对应的Binary目录，里面记载了maf文件中每一个样本在该bam文件的imputation的情况
###############输出结果##################

##  Binary matrix output
print "########## output imputation as Binary matrix ##########\n";
$file_name =~ /(\S+)_bqsr\.bam/;
my $sample_2 = $1;
for my $k1(sort keys %sample){
	print "\tsave sample in $sample_2 with $k1\n";
	open BM, ">$outdir/Binary_mat_$sample_2/$k1\_binary_mat.txt" or die $!;
	print BM "\t";
	my $line = join "\t", sort keys %{$sample{$k1}}; # output sample names as header
	$line =~ s/(LLM-.*)//g; #need to modify the pattern
	print BM "\t$1\tN\n"; #N set as 0 refer to no imputed mutations
	for my $k2(sort {$count{$k1}{$a} <=> $count{$k1}{$b}} keys %{$count{$k1}}){ #sort according to the values of dict, and return the order, then map to the keys of dict, then return the new ordered keys of dict (here is the chr:pos-symbol)
		print BM "$k2$binary{$k1}{$k2}\t0\n"; #print the mutation state, imputed position\timputetion state
	}
	close BM;
}

print "########## finished ##########\n";
###############代码逻辑##################
#count变量记录了样本在染色体-位置-基因处是否为imputed mutation
#首先根据样本count字典的值进行排序，得到一个新的序列map进样本count的key中，得到该样本新的顺序的keys
#该区域所有样本都没有imputation，则不输出该区域
#输出样本-chr-区域-基因名以及该样本的imputation情况，以及N=0
################代码逻辑#################

################理论基础#################
#
#甲病人测了ABC三个样本
#	1. imputation是把A样本中的突变提取，然后去看看BC样本中有没有跟A样本一样的突变，有就回补进去
#	2. 假如原本BC没有检测到TP53的突变，A检测到了，imputation的时候发现BC也有，那就回补回去，BC也有TP53
#	3. 可能BC的TP53突变丰度低一些，不做imputation，能找回的突变数目就少了
#
#imputation假设测序质量比较高以及测序深度足够深，以确保可靠性
#肿瘤具有异质性，一个肿瘤组织里面，有的肿瘤细胞有这个突变，有的没有，在这种情况下那些本来低丰度的突变或者是不占有优势的突变容易被忽略掉
#过滤突变是一个有一定主观性的而非简单客观的，尤其是那些处于标准线边缘的那些突变
#
################理论基础#################

################使用方法#################
# 1. 将同一个病人过滤后的突变的maf文件整理成一个maf文件，里面染色体的列去除chr
# 2. 输入病人maf文件，以及需要被imputation的bqsr bam文件（可以每一个样本都imputation一下）
# 3. 在bam文件中对maf文件中每一个突变位置都找一下，看看在cutoff前提下突变存不存在，并重新计算该突变的ref、alt count写入new imputed maf文件中（该突变在该样本的参考位点测序深度、突变位点测序深度以及vaf）
# 4. 根据binary的矩阵结果提取为1的区域的突变
################使用方法#################













