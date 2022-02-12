#Nicolas Altemose, Feb 2022
#input: a fasta file with DNA sequences
#output 1: text file with top repeat period peaks and column sums
#output 2: text file with matrix of k-mer interval information (optional) 

use warnings;
use strict;
srand(0);
my $tic = time;
print "\n\n";

my $usage = "usage: perl NTRprism_ProcessFasta_v0.22.pl region_fasta_file_path[REQUIRED] output_prefix[NTRspectrum] bin_size[1] total_span[30000] kmer_min_count[10] kmer_length[6] suppress_matrix_output[0]\n";

#initialize parameters and hardcode defaults
my $fastafile = "";#fasta file containing sequences to analyze
my $outbase = "NTRprism"; #prefix for output files
my $binsize = 1; #size of interval length bins (highest resolution =1)
my $span = 30000; #max interval size to interrogate, in bp. anything larger will be collapsed into span+1. minimum 1000.
my $countthresh = 10; #minimum number of times a k-mer is seen to count towards heatmap/spectrum. minimum 2
my $klen = 6; #kmer length
my $suppress = 0; #suppress matrix output [set to any value except 0 to suppress]

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
if(defined $ARGV[6]){
	$suppress = $ARGV[6];
	chomp($suppress);
	if($suppress != 0){
		$suppress = 1;
	}
}

my $tophitsfile = "$outbase.bin"."$binsize.NTRprism_TopHits.txt";

my $maxbin = $span/$binsize+1;
open(TOP,'>'.$tophitsfile);
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
			
			#initialize variables
			my $totcount=0;
			my %kmers; #hash to store the frequency of each kmer
			my %kmerslastpos; #hash to store the last seen position of each k-mer
			my %kmersdist; #hash to store the interval lengths for each k-mer
		
			my $outfile = $outbase.".region_".$regionname.".span".$span.".k".$klen.".mincount".$countthresh.".bin".$binsize.".txt";
			$fullseq=~s/[\n\s]//g;
			my $seqlen = length($fullseq);

			#count kmers in sequence and record interval lengths
			for(my $i = 0; $i<$seqlen-$klen;$i++){
				my $sub = substr($fullseq,$i,$klen);
				$kmers{$sub}++;
				if(exists $kmerslastpos{$sub}){
					my $interval = $i-$kmerslastpos{$sub};
					my $bin1 = int($interval/$binsize);
					if($interval>$span){
						$bin1 = $maxbin; #set all interval lengths above the max span to the max span +1
					}
					$kmersdist{$sub}{$bin1}++;
				}
				$kmerslastpos{$sub}=$i;
			}
			my @orderedkmers = sort {$kmers{$b}<=>$kmers{$a}} keys %kmers;#sort kmers by descending frequency

			my @colsums = (0) x ($maxbin+1);
			my @nums = 0..($maxbin);
			my $totalk = 0;
			#print interval distributions for each k-mer to text outfile
			open(OUT2,'>'.$outfile) unless($suppress==1);
			foreach my $mer (@orderedkmers){
					my $kcount = $kmers{$mer};
					last if($kcount<$countthresh);
					my $prop = $kcount/$seqlen;
					$totalk++;
					my $freq  = int(1/$prop);
					print OUT2 "$mer\t$kcount\t$prop\t$freq" unless($suppress==1);
					my @bins = (0) x ($maxbin+1);
					for my $i (keys %{ $kmersdist{$mer} }){
						$colsums[$i]+=$kmersdist{$mer}{$i}/($kcount-1); #normalize to sum to 1
						$bins[$i]=$kmersdist{$mer}{$i};
					}
					my $printstr = join("\t",@bins); #convert array to tsv string
					print OUT2 "\t$printstr\n" unless($suppress==1);
			}
			close OUT2 unless($suppress==1);
			if($totalk>0){
				$colsums[$maxbin]=0; #exclude bins above span from maximum column sum computations
				for(my $n=0;$n<scalar(@colsums);$n++){
					$colsums[$n]=$colsums[$n]/$totalk; #normalize to total number of k-mers above count threshold
				}
				my @maxlist = sort {$b<=>$a} @colsums;
				my @whichmaxlist = sort {$colsums[$b]<=>$colsums[$a]} @nums;

				my $topten1 = join("\t",@whichmaxlist[0..4]);
				my $topten2 = join("\t",@maxlist[0..4]);
				print TOP "$regionname\t$topten1\t$topten2\n";
				print "$ct2: printed to $outfile; $topten1\t$topten2\n";
			}else{
				print TOP "$regionname\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
				print "$ct2: no k-mers satisfy kmer_min_count ($countthresh)\n";
			}
		}else{
			die "Error reading in header number $ct2\n";
		}#closes if	
	}#closes while
	close FASTA;
}

close TOP;

#calculate and display runtime
my $toc = time;
my $elapsed = $toc-$tic;
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($elapsed / 3600), int(($elapsed % 3600) / 60), int($elapsed % 60));
