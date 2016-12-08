##########################################################################################
# VirusFinder.pl script identifies respiratory viruses in the human-unmapped reads from  #
# WTS RNA-seq data. It takes the .nomapping SAM file produced by GSNAP as input, then    #
# another round of mapping to all known human sequences from NCBI with SNAP to remove    #
# lower quality human reads. Then the remaining reads are BLASTed against the NCBI NT    #
# database, viral reads are identified and assembled into longer contigs. Then the 		 #
# closest viral sequence is identified with BLAST, the sequence is downloaded from NCBI  #
# and all the human unmapped reads are mapped to it to get the final virus coverage.     #
# Finally, a coverage plot of the identified viral genome is produced.					 #																			 # If <95% of viral reads can be mapped to the best identified viral sequence, and >50    #
# reads are still unmapped, then attempt another iteration of assembly and finding 		 #
# a potential secondary viral infection.												 #
#																						 #
# Run:																					 #
# The following modules have to be loaded:												 #
# module load perl5																		 #
# module load snap																		 #
# module load blast+																	 #
# module load samtools/0.1.19															 #
# module load R/3.2.0																	 #
# module load velvet																	 #
# module load bedtools																	 #
# 																						 #
# ./VirusFinder.pl foo.nomapping														 #
#																						 #
##########################################################################################

#!/usr/bin/perl
use File::Basename;
use List::Util qw(max);
use Bio::DB::GenBank;
use Bio::SeqIO;
use XML::Simple qw(:strict);
use Data::Dumper;
use LWP;

use strict;

########## Input file and directory ##########

my $file = $ARGV[0]; 	### .nomapping file from GSNAP
if ($file !~ m/nomapping$/){ die "$file has to be a SAM file with .nomapping extension.\n"; }
my $base = basename($file, ".nomapping");
unless (-d $base){mkdir $base;}
unless (-d "$base/Unmapped_reads"){mkdir "$base/Unmapped_reads";}
unless (-e "$base/Unmapped_reads/$base.nomapping"){
	system "cp $file $base/Unmapped_reads/$base.nomapping";
}

open(LOGOUT,'>',"$base/$base.log");

## define directory with taxonomy files, and other scripts
my $vir_dir = "/Seibold/proj/Nasal_microbiome/Virus_detection_pipeline";
## number of threads for mapping and blasting
my $n_threads = 12;

########## Convert .nomapping to .fq ##########
my $get_fq = "awk '{if (\$1 !~ /^\@/){print \"\@\" \$1 \"\\n\" \$10 \"\\n+\\n\" \$11}}' $file > $base/Unmapped_reads/$base.unmapped.fq";
system $get_fq;

########## SNAP mapping to all human sequences from NCBI ##########
### this is split into two parts, and has to run on the fat nodes ###
### because SNAP is loading the whole index into memory ###
unless (-e "$base/Unmapped_reads/$base.hg19.1.sam" && -e "$base/Unmapped_reads/$base.hg19.2.sam"){
	my $snap1 = "bsub -q fat -n $n_threads -J \"$base.snap1\" -e \"$base/Unmapped_reads/$base.snap1.err\" -o \"$base/Unmapped_reads/$base.snap1.out\" \"snap single /notbackedup2/reference/hg19_rRNA_mito_Hsapiens_rna.1/ $base/Unmapped_reads/$base.unmapped.fq -d 12 -t $n_threads -o $base/Unmapped_reads/$base.hg19.1.sam\"";
	my $snap1_job = `$snap1`;
	my $snap1_job_id = $1 if $snap1_job =~ m/<(\d+)>/;

	my $snap2 = "bsub -q fat -n $n_threads -J \"$base.snap2\" -e \"$base/Unmapped_reads/$base.snap2.err\" -o \"$base/Unmapped_reads/$base.snap2.out\" \"snap single /notbackedup2/reference/hg19_rRNA_mito_Hsapiens_rna.2/ $base/Unmapped_reads/$base.unmapped.fq -d 12 -t $n_threads -o $base/Unmapped_reads/$base.hg19.2.sam\"";
	my $snap2_job = `$snap2`;
	my $snap2_job_id = $1 if $snap2_job =~ m/<(\d+)>/;

	## check if the SNAP jobs are done:
	my $done = 0;
	while ($done == 0){
		if (`bjobs | grep $snap1_job_id` eq '' && `bjobs | grep $snap2_job_id` eq ''){
			$done = 1;
		} else {sleep 60;}
	}
}

### find reads unmapped in both SNAP SAM files ###
my @snap1_unmapped = `awk '{if (\$3 == \"*\"){print \$1 \"\\t\" \$10}}' $base/Unmapped_reads/$base.hg19.1.sam`;
my @snap2_unmapped = `awk '{if (\$3 == \"*\"){print \$1 \"\\t\" \$10}}' $base/Unmapped_reads/$base.hg19.2.sam`;

my %snap1_unmapped = map { $_ => 1 } @snap1_unmapped;
my @unmapped_both;
foreach my $read ( @snap2_unmapped) {
  if (exists $snap1_unmapped{$read}) {
    push @unmapped_both, $read;
  }
}
print LOGOUT "SNAP unmapped reads:\t", scalar(@unmapped_both),"\n";

########## Create a FASTA file with all SNAP unmapped reads longer than 50 nt ##########
my $count = 0;
open(OUT,'>',"$base/Unmapped_reads/$base.unmapped.fa") or die;
foreach my $read (@unmapped_both){
	my ($id, $seq) = split("\t",$read);
	if (length($seq)>=50){
		print OUT ">", $id,"\n",$seq;
		$count++;
	}
}
close OUT;
print LOGOUT ">=50nt long:\t$count\n";

##########  BLAST unmapped reads against the NCBI NT database ##########
unless (-d "$base/Blast_NT"){mkdir "$base/Blast_NT";}
my $blastn = "bsub -n $n_threads -J \"$base.blast\" -e \"$base/Blast_NT/$base.blast.err\" -o \"$base/Blast_NT/$base.blast.out\" \"blastn -query $base/Unmapped_reads/$base.unmapped.fa -out $base/Blast_NT/$base.nt_blast.txt -outfmt 6 -db /data/reference/BLAST/nt -num_threads $n_threads -max_target_seqs 1  -best_hit_score_edge 0.05\"";

unless (-e "$base/Blast_NT/$base.nt_blast.txt"){
	my $blast_job = `$blastn`;
	my $blast_job_id = $1 if $blast_job =~ m/<(\d+)>/;

	### check if BLAST finished
	my $blast_done = 0;
	while ($blast_done == 0){
		if (`bjobs | grep $blast_job_id` eq ''){
			$blast_done = 1;
		} else {sleep 60;}
	}
}

##########  get taxonomy annotation for BLAST results ##########
system "perl $vir_dir/annotate_blast.pl $base/Blast_NT/$base.nt_blast.txt > $base/Blast_NT/$base.nt_blast.annotated.txt";

##########  Count virus and respiratory virus hits ##########

### Respiratory viruses: rhinovirus, parainfluenza, RSV, influenza, adenovirus, 
### coronavirus, enterovirus, metapneumovirus
### exclude hits to: KF772945.1	Influenza A virus (A/wild bird/Chile/1805/2008(H5N9)) segment 1 polymerase PB2 (PB2) gene, complete cds	Viruses; ssRNA viruses; ssRNA negative-strand viruses; Orthomyxoviridae; Influenzavirus A
### this is not a correct viral sequence and has been fixed in the latest release of NCBI

my @virus_reads;
my @respiratory_virus_reads;
open(IN,'<',"$base/Blast_NT/$base.nt_blast.annotated.txt") or die;
open(OUT,'>',"$base/Blast_NT/$base.virus.nt_blast.txt") or die;
while (my $line = <IN>){
	chomp $line;
	my @line = split("\t",$line);
	if ($line[13] =~ m/^Viruses; /){
		unless ($line[1] =~ m/KF772945\.1/){ ### exclude hits to Influenza A, those are wrong
			push(@virus_reads, $line[0]);
			print OUT $line,"\n";
			### now check if respiratory
			if ($line[12] =~ m/rhinovirus|parainfluenza|PIV|RSV|syncytial|influenza|adenovirus|coronavirus|enterovirus|metapneumovirus/i){
				push(@respiratory_virus_reads, $line[0]);
			}
		}
	}
}
close IN;
close OUT;
my @unique_vr = do { my %seen; grep { !$seen{$_}++ } @virus_reads };
my @unique_resp_vr = do { my %seen; grep { !$seen{$_}++ } @respiratory_virus_reads };

print LOGOUT "Total virus reads:\t", scalar(@unique_vr),"\n";
print LOGOUT "Respiratory virus reads:\t",scalar(@unique_resp_vr),"\n";

##########  Check if there were any respiratory virus reads, if not - exit ##########

if (scalar(@unique_resp_vr) == 0){ exit;}

##########  Create a file with viral reads id's ##########
unless (-d "$base/Viral_reads"){mkdir "$base/Viral_reads";}
open(OUT,'>',"$base/Viral_reads/$base.viral_reads_ids.txt") or die;
print OUT join("\n",@unique_vr),"\n";
close OUT;

##########  Find the most common virus hit from BLAST ##########
my $top_result = `awk -F "\t" '{print \$2}' $base/Blast_NT/$base.virus.nt_blast.txt | sort | uniq -c | sort -n | tail -1  | awk '{print \$2}' | awk -F "\|" '{printf \$4}'`;

print LOGOUT "Top virus accession - BLAST:\t$top_result\n";
my $ann = `grep $top_result $vir_dir/all_taxonomy | awk -F "\t" '{print \$2}'`;
print LOGOUT "Top virus annotation - BLAST:\t$ann";

##########  Download viral sequence from NCBI, get length, estimate coverage ##########
my $gb = Bio::DB::GenBank->new;
my $seq = $gb->get_Seq_by_id($top_result);
my $sequence = $seq->primary_seq->seq;
my $virus_length = length($sequence);

### we estimate coverage by assuming 120nt read length
my $est_cov = 120 * scalar(@unique_vr) / $virus_length;
print LOGOUT "Virus length:\t$virus_length\n";
print LOGOUT "Estimated coverage:\t$est_cov\n";	

########## Create a FASTQ file with viral reads ##########
unless (-e "$base/Unmapped_reads/$base.nomapping.bam"){system "samtools view -Sb $file > $base/Unmapped_reads/$base.nomapping.bam";}

my $R_command = "source(\"$vir_dir/subset_BAM_to_FASTQ.R\");subset_bam(\"$base/Unmapped_reads/$base.nomapping.bam\",\"$base/Viral_reads/$base.viral_reads_ids.txt\",\"$base/Viral_reads/$base.viral_reads.fq\")";
system "Rscript --vanilla -e '$R_command' >/dev/null";

########## If >10 reads, assemble reads into longer contigs ##########

if (scalar(@unique_resp_vr) <= 10){ exit;}

system "velveth $base/Velvet 31 -fastq $base/Viral_reads/$base.viral_reads.fq >/dev/null";
my $min_cov = roundup($est_cov/100);
system "velvetg $base/Velvet -min_contig_lgth 200 -exp_cov $est_cov -cov_cutoff $min_cov >/dev/null";

########## Annotate resulting contigs with BLAST ##########
my $contig_no = `grep -c \">\" $base/Velvet/contigs.fa`;
chomp $contig_no;

system "blastn -query $base/Velvet/contigs.fa -out $base/Velvet/$base.contigs.blast.txt -outfmt 6 -db /data/reference/BLAST/nt -num_threads $n_threads -max_target_seqs 3  -best_hit_score_edge 0.05";

########## Report the best hit ##########
$top_result = `awk '{print \$2}' $base/Velvet/$base.contigs.blast.txt | sort | uniq -c | sort -n | tail -1 | awk -F \"\|\" '{printf \$4}'`;

print LOGOUT "Top virus accession - Velvet:\t$top_result\n";
my $ncbi_ann = `curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$top_result&retmode=xml\" | grep -E 'definition'`;
my $def = $1 if $ncbi_ann =~ m/<GBSeq_definition>(.+)<\/GBSeq_definition>/;

print LOGOUT "Top virus annotation - Velvet:\t$def\n";

########## Download the best hit sequence from NCBI ##########
unless (-d "$base/$top_result"){mkdir "$base/$top_result";}
my $gb = Bio::DB::GenBank->new;
my $seq = $gb->get_Seq_by_id($top_result);
my $sequence = $seq->primary_seq->seq;
open(OUT,'>',"$base/$top_result/$top_result.fasta") or die;
print OUT ">$top_result\n",$sequence,"\n";
close OUT;

########## Create SNAP index for the best hit and map all viral reads against it #########
system "snap index $base/$top_result/$top_result.fasta $base/$top_result/$top_result >/dev/null";
system "snap single $base/$top_result/$top_result/ $base/Viral_reads/$base.viral_reads.fq -o $base/$top_result/$base.$top_result.bam -t $n_threads -d 30 >/dev/null";

########## Calculate coverage along the viral genome #########
open(OUT,'>',"$base/$top_result/$top_result.bed") or die;
print OUT "$top_result\t1\t", length($sequence),"\n";
close OUT;

system "coverageBed -abam $base/$top_result/$base.$top_result.bam -b $base/$top_result/$top_result.bed -d > $base/$top_result/$base.$top_result.coverage";

my $dp = `awk '{sum += \$5} END {printf sum/NR}' $base/$top_result/$base.$top_result.coverage`;
print LOGOUT "Viral sequence depth:\t$dp\n";

########## Create coverage plot for the viral genome #########

##### Get gene structure from NCBI #####

system "curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$top_result&retmode=xml\" > $base/$top_result/$top_result.xml";
my $simple = XML::Simple->new( ); 
my $config = $simple->XMLin("$base/$top_result/$top_result.xml",ForceArray =>0, KeyAttr => []);

open(OUT,'>',"$base/$top_result/$top_result.gene_structure.txt") or die;
my @feature_array = @{$config->{ GBSeq }->{ 'GBSeq_feature-table' }->{'GBFeature'}};
foreach my $f (@feature_array){
	my $gene = 0;
	if ($f->{'GBFeature_key'} eq "CDS"){
		my $loc = $f->{'GBFeature_location'};
		if ($loc =~ m/^(\d+)\.\.(\d+)$/){print OUT $1,"\t",$2,"\t";}
		elsif ($loc =~ m/^(\d+)\.\..?(\d+)$/){print OUT $1,"\t",$2,"\t";}
		elsif ($loc =~ m/^join\((\d+)\.\..+(\d+)\)$/){print OUT $1,"\t",$2,"\t";}
		foreach my $q (@{$f->{'GBFeature_quals'}->{'GBQualifier'}}){
			if ( $q->{'GBQualifier_name'} eq "gene"){
				print OUT $q->{'GBQualifier_value'},"\n";
				$gene=1;
			}
		}
		if ($gene == 0){	
		foreach my $q (@{$f->{'GBFeature_quals'}->{'GBQualifier'}}){	
			if ( $q->{'GBQualifier_name'} eq "product"){
				print OUT $q->{'GBQualifier_value'},"\n";
			}
		}
		}
	}
}
close OUT;

##### Plot coverage #####
my $R_plot = "source(\"$vir_dir/cov_plot.R\");cov_plot(\"$base/$top_result/$base.$top_result.coverage\",\"$base/Velvet/$base.contigs.blast.txt\",\"$base/$top_result/$top_result.gene_structure.txt\",\"$base\",\"$top_result\")";
system "Rscript --vanilla -e '$R_plot' >/dev/null";

########## Check what percentage of the viral reads identified by BLAST mapped ###########
########## to the target virus sequence with SNAP, if <95% and >50 reads left  ###########
########## then try to assemble reads further to detect secondary infections     ###########

my $vir_mapped = `samtools view $base/$top_result/$base.$top_result.bam | awk '{if (\$3 == "$top_result"){sum += 1}} END {printf sum}'`;
my $vir_map_pc = $vir_mapped/scalar(@unique_resp_vr)*100;
my $vir_unmap = scalar(@unique_resp_vr) - $vir_mapped;

print LOGOUT "Detected virus accounts for $vir_map_pc % of the respiratory viral reads.\n$vir_unmap reads not accounted for\n";

my $iterate = 0;
my $iteration = 1;

### store all the viruses identified so far...
my @viruses;

if ($vir_map_pc < 95 && $vir_unmap >= 50){
	unless (-d "$base/iteration.$iteration"){mkdir "$base/iteration.$iteration";}
	unless (-d "$base/iteration.$iteration/$top_result"){mkdir "$base/iteration.$iteration/$top_result";}
	system "cp $base/$top_result/$base.$top_result.bam $base/iteration.$iteration/$top_result/$base.$top_result.bam";
	
	push(@viruses, $top_result);
	
	$iterate = 1;
	$iteration++;	
} else {$iterate = 0;}

while ($iterate == 1){
	print LOGOUT "Running iteration $iteration of the pipeline to find secondary virus infections\n";
	### create a directory for the next virus hit
	unless (-d "$base/iteration.$iteration"){mkdir "$base/iteration.$iteration";}
	### create a FQ file with reads unmapped to previous hit
	my $p_iteration = $iteration -1;
	system "samtools view $base/iteration.$p_iteration/$top_result/$base.$top_result.bam | awk '{if (\$3 != \"$top_result\"){print \"\@\" \$1 \"\\n\" \$10 \"\\n+\\n\" \$11}}' > $base/iteration.$iteration/$base.unmapped.fq";
	### now assemble reads
	my $est_cov = 120 * $vir_unmap / $virus_length;
	my $min_cov = roundup($est_cov/100);
	system "velveth $base/iteration.$iteration/Velvet 31 -fastq $base/iteration.$iteration/$base.unmapped.fq >/dev/null";
	system "velvetg $base/iteration.$iteration/Velvet -min_contig_lgth 200 -exp_cov $est_cov -cov_cutoff $min_cov >/dev/null";

	if (-s "$base/iteration.$iteration/Velvet/contigs.fa"){
	########## Annotate resulting contigs with BLAST ##########
	system "blastn -query $base/iteration.$iteration/Velvet/contigs.fa -out $base/iteration.$iteration/Velvet/$base.contigs.blast.txt -outfmt 6 -db /data/reference/BLAST/nt -num_threads $n_threads -max_target_seqs 3  -best_hit_score_edge 0.05";

	########## Report the best hit ##########
	my $res_so_far = join('\\|',@viruses);
	$top_result = `awk '{print \$2}' $base/iteration.$iteration/Velvet/$base.contigs.blast.txt | grep -v "$res_so_far" | sort | uniq -c | sort -n | tail -1 | awk -F \"\|\" '{printf \$4}'`;
	if ($top_result){
	push (@viruses, $top_result);
	
	print LOGOUT "Top virus accession - Velvet:\t$top_result\n";
	my $ncbi_ann = `curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$top_result&retmode=xml\" | grep -E 'definition'`;
	my $def = $1 if $ncbi_ann =~ m/<GBSeq_definition>(.+)<\/GBSeq_definition>/;
	print LOGOUT "Top virus annotation - Velvet:\t$def\n";

	########## Download the best hit sequence from NCBI ##########
	unless (-d "$base/iteration.$iteration/$top_result"){mkdir "$base/iteration.$iteration/$top_result";}
	my $gb = Bio::DB::GenBank->new;
	my $seq = $gb->get_Seq_by_id($top_result);
	my $sequence = $seq->primary_seq->seq;
	open(OUT,'>',"$base/iteration.$iteration/$top_result/$top_result.fasta") or die;
	print OUT ">$top_result\n",$sequence,"\n";
	close OUT;

	########## Create SNAP index for the best hit and map all viral reads against it #########
	system "snap index $base/iteration.$iteration/$top_result/$top_result.fasta $base/iteration.$iteration/$top_result/$top_result >/dev/null";
	system "snap single $base/iteration.$iteration/$top_result/$top_result/ $base/iteration.$iteration/$base.unmapped.fq -o $base/iteration.$iteration/$top_result/$base.$top_result.bam -t $n_threads -d 30 >/dev/null";

	########## Calculate coverage along the viral genome #########
	open(OUT,'>',"$base/iteration.$iteration/$top_result/$top_result.bed") or die;
	print OUT "$top_result\t1\t", length($sequence),"\n";
	close OUT;

	system "coverageBed -abam $base/iteration.$iteration/$top_result/$base.$top_result.bam -b $base/iteration.$iteration/$top_result/$top_result.bed -d > $base/iteration.$iteration/$top_result/$base.$top_result.coverage";

	my $dp = `awk '{sum += \$5} END {printf sum/NR}' $base/iteration.$iteration/$top_result/$base.$top_result.coverage`;
	print LOGOUT "Viral sequence depth:\t$dp\n";

	########## Create coverage plot for the viral genome #########

	##### Get gene structure from NCBI #####

	system "curl -s \"http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$top_result&retmode=xml\" > $base/iteration.$iteration/$top_result/$top_result.xml";
	my $simple = XML::Simple->new( ); 
	my $config = $simple->XMLin("$base/iteration.$iteration/$top_result/$top_result.xml",ForceArray =>0, KeyAttr => []);

	open(OUT,'>',"$base/iteration.$iteration/$top_result/$top_result.gene_structure.txt") or die;
	my @feature_array = @{$config->{ GBSeq }->{ 'GBSeq_feature-table' }->{'GBFeature'}};
	foreach my $f (@feature_array){
		my $gene = 0;
		if ($f->{'GBFeature_key'} eq "CDS"){
			my $loc = $f->{'GBFeature_location'};
			if ($loc =~ m/^(\d+)\.\.(\d+)$/){print OUT $1,"\t",$2,"\t";}
			elsif ($loc =~ m/^(\d+)\.\..?(\d+)$/){print OUT $1,"\t",$2,"\t";}
			elsif ($loc =~ m/^join\((\d+)\.\..+(\d+)\)$/){print OUT $1,"\t",$2,"\t";}
			foreach my $q (@{$f->{'GBFeature_quals'}->{'GBQualifier'}}){
				if ( $q->{'GBQualifier_name'} eq "gene"){
					print OUT $q->{'GBQualifier_value'},"\n";
					$gene=1;
				}
			}
			if ($gene == 0){	
			foreach my $q (@{$f->{'GBFeature_quals'}->{'GBQualifier'}}){	
				if ( $q->{'GBQualifier_name'} eq "product"){
					print OUT $q->{'GBQualifier_value'},"\n";
				}
			}
			}
		}
	}
	close OUT;

	##### Plot coverage #####
	my $R_plot = "source(\"$vir_dir/cov_plot.R\");cov_plot(\"$base/iteration.$iteration/$top_result/$base.$top_result.coverage\",\"$base/iteration.$iteration/Velvet/$base.contigs.blast.txt\",\"$base/iteration.$iteration/$top_result/$top_result.gene_structure.txt\",\"$base\",\"$top_result\")";
	system "Rscript --vanilla -e '$R_plot' >/dev/null";
	
	##### Check if next iteration should be run:
	my $vir_mapped = `samtools view $base/iteration.$iteration/$top_result/$base.$top_result.bam | awk '{if (\$3 == "$top_result"){sum += 1}} END {printf sum}'`;
	my $vir_map_pc = $vir_mapped/scalar(@unique_resp_vr)*100;
	$vir_unmap = $vir_unmap - $vir_mapped;
	print LOGOUT "Detected virus accounts for $vir_map_pc % of the respiratory viral reads.\n$vir_unmap reads not accounted for\n";

	if ($vir_map_pc < 95 && $vir_unmap >= 50){
		$iterate = 1;
		$iteration++;
	} else {$iterate = 0;}
	} else {$iterate = 0; print LOGOUT "No new viral sequences among remaining contigs\n";} ## if all the contigs match already identified sequences
	} else {$iterate = 0; print LOGOUT "Can't assemble any more contigs\n";} ## if can't assemble any more contigs
}




close LOGOUT;

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}
