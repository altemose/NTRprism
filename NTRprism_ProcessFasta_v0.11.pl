#Nicolas Altemose, April 2021
#creates a text file with k-mer interval information for a bed file of genomic regions

use warnings;
use strict;
my $tic = time;
print "\n\n";

my $usage = "usage: perl NTRprism_ProcessFasta_v0.11.pl region_fasta_file_path[REQUIRED] output_prefix[NTRspectrum] bin_size[1] total_span[30000] kmer_min_count[10] kmer_length[6]\n";

#initialize parameters and hardcode defaults
my $fastafile = "";#fasta file containing sequences to analyze
my $outbase = "NTRprism"; #prefix for output files
my $binsize = 1; #size of interval length bins (highest resolution =1)
my $span = 30000; #max interval size to interrogate, in bp. anything larger will be collapsed into span+1. minimum 1000.
my $countthresh = 10; #minimum number of times a k-mer is seen to count towards heatmap/spectrum. minimum 2
my $klen = 6; #kmer length


#read in and check command line arguments
if(defined $ARGV[0]){
	$fastafile = $ARGV[0];
	chomp($fastafile);
}else{
	die "$usage\n";
}
if(defined $ARGV[1]){
	$outbase = $ARGV[1];
	chomp($outbase);
}
if(defined $ARGV[2]){
	$binsize = $ARGV[2];
	chomp($binsize);
	if($binsize<1){
		$binsize=1;
	}
}
if(defined $ARGV[3]){
	$span = $ARGV[3];
	chomp($span);
	if($span<1000){
		$span=1000;
	}
}
if(defined $ARGV[4]){
	$countthresh = $ARGV[4];
	chomp($countthresh);
	if($countthresh<2){
		$countthresh=2;
	}
}
if(defined $ARGV[5]){
	$klen = $ARGV[5];
	chomp($klen);
	if($klen<1){
		$klen=1;
	}
}

my $tophitsfile = "";
if($fastafile=~/.+\/(.+)\..+?$/){
	$tophitsfile = "$1_NTRprism_TopHits.txt";
}elsif($fastafile=~/^(.+)\..+?$/){
	$tophitsfile = "$1_NTRprism_TopHits.txt";
}else{
	$tophitsfile = "$fastafile.NTRprism_TopHits.txt";
}

open(TOP,'>'.$tophitsfile);
#my $finalstring="";
#my $topscores = "";
#open fasta and loop through each sequence one at a time
if(1){
	local $/ = ">"; #change the record separator (locally)-i.e. split files by '>' instead of newlines
	open(FASTA, $fastafile) or die "ERROR: unable to open $fastafile\n"; #open input file for reading

	my $ct2=0; #initialise a counter to count the number of sequences read in
	my $junk = <FASTA>; #remove first > from beginning of file
	while(my $frecord = <FASTA>){ #read in input file one fasta record at a time
		chomp($frecord); #remove any newlines from the end of the record
		$ct2++;
		if($frecord=~/^(\S+?)\n(.*)$/s){ #check whether each record matches the expected format

			my $regionname = $1;
			my $fullseq=$2;
			
			$regionname=~s/[\s\:\-]/\./g;
			#$finalstring.=" $regionname";
			
			#initialize variables
			my $totcount=0;
			my %kmers; #hash to store the frequency of each kmer
			my %kmerslastpos; #hash to store the last seen position of each k-mer
			my %kmersdist; #hash to store the interval lengths for each k-mer
		
			my $outfile = $outbase.".region_".$regionname.".span".$span.".k".$klen.".mincount".$countthresh.".bin".$binsize.".txt";
			$fullseq=~s/[\n\s]//g;
			$fullseq=~s/[LRUN]+/N/g; #convert/compress special characters to N. could remove if eg using protein sequence
			my $seqlen = length($fullseq);

			#count kmers in sequence and record interval lengths
			for(my $i = 0; $i<$seqlen-$klen;$i++){
				my $sub = substr($fullseq,$i,$klen);
				$kmers{$sub}++;
				if(exists $kmerslastpos{$sub}){
					$kmersdist{$sub}.= ",".($i-$kmerslastpos{$sub}); #store list of interval lengths as comma-delimited string
				}
				$kmerslastpos{$sub}=$i;	
			}
			my @orderedkmers = sort {$kmers{$b}<=>$kmers{$a}} keys %kmers;#sort kmers by descending frequency

			my @colsums = (0) x (int($span/$binsize));
			my @nums = 1..(scalar(@colsums));
			#print interval distributions for each k-mer to text outfile
			open(OUT2,'>'.$outfile);
			foreach my $mer (@orderedkmers){
					my $kcount = $kmers{$mer};
					my $prop = $kcount/$seqlen;
					next unless $kcount>=$countthresh;
					my $freq  = int(1/$prop);
					print OUT2 "$mer\t$kcount\t$prop\t$freq";
					my @bins = (0) x (int($span/$binsize)+1); #initialize array with all histogram bins
					my $dists = $kmersdist{$mer}; #get comma-delimited string of interval lengths
					$dists =~ s/^\,//; #remove leading comma
					my @list = split(",",$dists); #split to array on commas
					foreach my $dist(@list){
						if($dist>$span){
							$bins[$span/$binsize]++;  #store and intervals larger than $span in the last bin
						}else{
							$bins[int($dist/$binsize)-1]++; #increment bin corresponding to each interval length
						}
					}
					my $printstr = join("\t",@bins); #convert array to tsv string
					for(my $i=0;$i<scalar(@colsums);$i++){
						$colsums[$i]+=$bins[$i]/($kcount-1);
					}
					print OUT2 "\t$printstr\n";
			}
			close OUT2;
			my @maxlist = sort {$b<=>$a} @colsums;
			my @whichmaxlist = sort {$colsums[$b-1]<=>$colsums[$a-1]} @nums;
			#my $max = $whichmaxlist[0];
			#$topscores.=','.$max;
			my $topten1 = join("\t",@whichmaxlist[0..9]);
			#my $topten2 = join("\t",@maxlist[0..9]);
			print TOP "$regionname\t$topten1\n";
			
			print "printed to $outfile\n";#"\t$topten1\n\t$topten2\n";

		}else{
			die "Error reading in header number $ct2\n";
		}#closes if	
	}#closes while
	close FASTA;
}

close TOP;
#$finalstring =~ s/^\s+//;
#print "\n$finalstring\n";
#$topscores =~ s/^\,//;
#print "\n$topscores\n\n";

#my @topscoresA = split(',',$topscores);
#my @finalstringA = split(/\s+/,$finalstring);
#for(my $i=0;$i<scalar(@topscoresA);$i++){
#	print "$finalstringA[$i]\t$topscoresA[$i]\n";
#}

#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
