# README for NTRprism scripts v0.11
Nicolas Altemose, 2021
- purpose: to discover repeat periodicity in large tandem repeat arrays
- summary: takes a DNA sequence in fasta format as input and decomposes it into inter-k-mer intervals then combines interval length distributions across k-mers, and plots the results


## Required files
- NTRprism_ProcessFasta_v0.11.pl
- NTRprism_PlotSpectrum.r
- NTRprism_PlotHeatmap.r
- also requires r packages ggplot2 and reshape2
- a test input file from beta satellite arrays on human chr1 and chr9 is provided as test.fa

## Usage
### Step 1) Obtain a fasta file containing sequences from regions of interest. One spectrum will be produced per sequence in the fasta file. 
e.g. use samtools faidx to extract regions from a reference genome


### Step 2) run perl script to generate human readable tables
- usage: perl NTRprism_ProcessFasta_v0.11.pl region_fasta_file_path[REQUIRED] output_prefix[NTRspectrum] bin_size[1] total_span[30000] kmer_min_count[10] kmer_length[6]
- outputs: produces one output file per region of the format: OUTPUT_PREFIX.regionREGION_NAME_FROM_FASTA_HEADER.spanTOTAL_SPAN.kKMER_LENGTH.mincountKMER_MINCOUNT.binBIN_SIZE.txt
- also outputs a _NTRprism_TopHits.txt file, with one line per fasta entry listing the entry's header and top 10 peaks (one per column, e.g. column 2 lists the largest peak) 
- note 1: I usually run with both bin size 1 (better for spectrum) and bin size 100 (better for heatmap). for example:
- note 2: the method can be sensitive to the choice of kmer_min_count. a value should be chosen that allows for sufficiently granular interval length distribution inference, i.e. setting this value too low will flood the combined distribution with rare k-mers that may only have a few intervals in the entire sequence and setting it too high will discard k-mers that occur with a frequency proportional to the true predominant repeat unit length
- note 3: using k=6 is sufficient to detect most peri/centromeric satellite periodicities in the human genome. Increasing k may improve sensitivity to very long periodicities but increases memory complexity.
- note 4: increasing the maximum span increases runtime substantially. It can help to start by running with a large span and large bin size, then examining the heatmap and setting the max span to capture the largest substantial peaks.
```
perl NTRprism_ProcessFasta_v0.11.pl test.fa NTRprism_TEST 1 5000 30 6
perl NTRprism_ProcessFasta_v0.11.pl test.fa NTRprism_TEST 100 5000 30 6
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
perl NTRprism_ProcessFasta_v0.11.pl test.fa NTRprism_TEST 100 5000 30 6
perl NTRprism_ProcessFasta_v0.11.pl test.fa NTRprism_TEST 1 5000 30 6
for i in *.bin1.txt;do
	j=$(perl -lae 'if(m/region_(.+)\.span/){print "$1"}' <(echo $i)) #automatically extract region name from filename
	Rscript NTRprism_PlotSpectrum.r --args $i $j 5000 NTRprism_TEST
done
for i in *.bin100.txt;do
	j=$(perl -lae 'if(m/region_(.+)\.span/){print "$1"}' <(echo $i))
	Rscript NTRprism_PlotHeatmap.r --args $i $j 5000 100 NTRprism_TEST
done
```

