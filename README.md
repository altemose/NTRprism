# README for NTRprism scripts v0.22
Nicolas Altemose, 2021
- purpose: creates a Nested Tandem Repeat (NTR) "spectrum" indicating the most abundant tandem repeat periodicities found in an input DNA sequence
- summary: takes a DNA sequence in fasta format as input and decomposes it into inter-k-mer intervals then combines interval length distributions across k-mers, and plots the results (see pseudocode below)
- version notes: uses a 2D hash to speed up counting, resulting in 2X faster runtime relative to v0.11, fixes small bugs relative to v0.21

## Required files
- NTRprism_ProcessFasta_v0.22.pl
- NTRprism_PlotSpectrum.r
- NTRprism_PlotHeatmap.r
- also requires r packages ggplot2 and reshape2
- a test input file from beta satellite arrays on human chr1 and chr9 is provided as test.fa

## Usage
### Step 1) Obtain a fasta file containing sequences from regions of interest. One spectrum will be produced per sequence in the fasta file. 
e.g. use samtools faidx to extract regions from a reference genome


### Step 2) run perl script to generate human readable tables
- usage: perl NTRprism_ProcessFasta_v0.21.pl region_fasta_file_path[REQUIRED] output_prefix[NTRspectrum] bin_size[1] total_span[30000] kmer_min_count[10] kmer_length[6] suppress_matrix_output[0]
- output 1: a *OUTPUT_PREFIX.binBINSIZE.NTRprism_TopHits.txt* file, with one line per fasta entry listing the entry's header and top 5 peaks (one per column, e.g. column 2 lists the largest peak), followed by the final normalized column sums for those top 5 peaks (e.g. column 7 in the output contains the normalized column sum for the top peak)
- output 2 (only if suppress_matrix_output == 0): produces one output file per region of the format: OUTPUT_PREFIX.regionREGION_NAME_FROM_FASTA_HEADER.spanTOTAL_SPAN.kKMER_LENGTH.mincountKMER_MINCOUNT.binBIN_SIZE.txt
- note 1: I usually run with both bin size 1 (better for spectrum) and bin size 100 (better for heatmap). for example:
- note 2: the method can be sensitive to the choice of kmer_min_count. a value should be chosen that allows for sufficiently granular interval length distribution inference, i.e. setting this value too low will flood the combined distribution with rare k-mers that may only have a few intervals in the entire sequence and setting it too high will discard k-mers that occur with a frequency proportional to the true predominant repeat unit length
- note 3: using k=6 is sufficient to detect most peri/centromeric satellite periodicities in the human genome. Increasing k may improve sensitivity to very long periodicities but increases memory complexity.
- note 4: increasing the maximum span increases runtime substantially. It can help to start by running with a large span and large bin size, then examining the heatmap and setting the max span to capture the largest substantial peaks.
```
perl NTRprism_ProcessFasta_v0.22.pl test.fa NTRprism_TEST 1 5000 30 6 0
perl NTRprism_ProcessFasta_v0.22.pl test.fa NTRprism_TEST 100 5000 30 6 0
```

### Step 3) make plots in R
- usage: Rscript NTRprism_PlotSpectrum.r --args perl_output_filepath region_name max_period output_prefix
- prints region name and top 10 peak lengths to STDOUT (and on plot)
- usage: Rscript NTRprism_PlotHeatmap.r --args perl_output_filepath region_name max_period bin_size output_prefix

make spectral plots for all regions using bin size 1
```
for i in *.bin1.txt;do
	j=$(perl -lae 'if(m/region_(.+)\.span/){print "$1"}' <(echo $i)) #automatically extract region name from filename
	Rscript NTRprism_PlotSpectrum.r --args $i $j 5000 NTRprism_Plot
done
```
make heatmaps for all regions using bin size 100
```
for i in *.bin100.txt;do
	j=$(perl -lae 'if(m/region_(.+)\.span/){print "$1"}' <(echo $i))
	Rscript NTRprism_PlotHeatmap.r --args $i $j 5000 100 NTRprism_Plot
done
```


## Full test run
using provided test.fa file (beta satellite from chr1 and chr9, with 68 bp and 2596 bp units)
```
perl NTRprism_ProcessFasta_v0.22.pl test.fa NTRprism_TEST 100 5000 30 6 0
perl NTRprism_ProcessFasta_v0.22.pl test.fa NTRprism_TEST 1 5000 30 6 0
for i in *.bin1.txt;do
	j=$(perl -lae 'if(m/region_(.+)\.span/){print "$1"}' <(echo $i)) #automatically extract region name from filename
	Rscript NTRprism_PlotSpectrum.r --args $i $j 5000 NTRprism_TEST
done
for i in *.bin100.txt;do
	j=$(perl -lae 'if(m/region_(.+)\.span/){print "$1"}' <(echo $i))
	Rscript NTRprism_PlotHeatmap.r --args $i $j 5000 100 NTRprism_TEST
done
```

## Pseudocode
```
#NTRprism pseudocode
#Let S be a string of [ACGT] of length s.
#Let S[i:i+k-1] be the substring of S starting at 0-based index i, 
#    and having length k.
#Let K be a dictionary of k-mers, used to store the last observed
#    position of each k-mer.
#Let j be a particular k-mer of length k.
#Let L be a dictionary of k-mers, with each entry L[j] as an array of 
#    0s of length m+1 [0-based indexing], with m equal to the longest 
#    allowed repeat periodicity.
#Let t be the threshold for minimum observed k-mer count [default 2].

i=0
while(i < s-k){ #loop through all substrings of length k in the sequence  
	j = S[i:i+k-1] 
	if(exists L[j]){ #if this particular k-mer j has been seen before, record the  
                    interval length between the last occurrence and this occurrence
		interval = i-K[j] 
		if(interval > m+1){ #if longer than the max span m, store it as m+1
			interval = m+1
		}
		L[j][interval-1]++ #increment the interval length???s count in j???s array
	}else{ #if this is the first occurrence of j, initialize its count array
		L[j]=array(0,m+1)
	}
	K[j]=i #store the current occurrence's position
	i++
}
sorted = sort(j in L, sum(L[j]), descending) #create an array of all observed k-mers       
                                         in descending order of their frequency in S

ColSums = array(0,m+1) #initialize an array of 0s to store normalized column sums
c=0
while(c < length(sorted)){ #loop through all observed k-mers
	j = sorted[c]
	count = sum(L[j])+1 #report the count of j as the number of intervals +1
	if(count >= t){
		L[j] = L[j]/(count-1) #normalize interval count array to sum to 1
	}
last if(count < t)
i=0 #increment column sums
	while(i<m+1){
		ColSums[i] = ColSums[i] + L[j][i]
            i++                           
	}
	c++
}
ColSums=ColSums/c #normalize col. sums by the tot. no. of k-mers above thresh. t
Lmat = L[sorted[0:c]] #create a matrix from L with each row corresponding to a k-mer 
#above the count threshold m, sorted by decreasing frequency, 
#and each column corresponding to an interval length
print(Lmat) #print the normalized, sorted, thresholded matrix
heatmap(Lmat) #plot the matrix as a heatmap
barplot(ColSums) #plot ColSums as a barplot to produce an NTR spectrum plot
TopPeak = argmax(ColSums[0:m-1])+1 #report the most common repeat periodicity (excluding m+1)
print(TopPeak)

```

