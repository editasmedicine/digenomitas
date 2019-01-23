# Digenomitas

This repository contains software for analyzing [digneome-sequencing](https://www.ncbi.nlm.nih.gov/pubmed/25664545) data.

## License

See the License file in this project. 

## Pre-requisites

1. A Java Development Kit v8 or higher available from [the Oracle Java site ](http://www.oracle.com/technetwork/java/javase/downloads/index.html).
2. A working installation of the Scala Built Tool aka _sbt_ v1.1 or higher. Instructions for installing sbt can be found on the [SBT website](https://www.scala-sbt.org/download.html).
3. A working internet connection which SBT will use to automatically download dependencies

## Building

The software has been built and tested on Linux and Mac OS X.  The software should build and run on any operating system that supports Java >= 8 and SBT, but other platforms have not been tested.

To compile and test the software run the following in the same directory as this README file:

```
sbt clean test
```

This process may take several minutes the first time it is run as sbt will download needed software libraries from the internet.  Repeated executions should take 15-30s.  A successful run of the tests will result in lines similar to the following being printed before sbt exits:

```
[info] All tests passed.
[info] Passed: Total 34, Failed 0, Errors 0, Passed 34
[success] Total time: 12 s, completed Oct 2, 2018 11:53:18 AM
```

The digenomitas software is compiled into a standalone [JAR](https://en.wikipedia.org/wiki/JAR_(file_format)) file for use.  To build the JAR file run:

```
sbt assembly
```

As with the tests, the last line printed should contain `[success]`, e.g.: 

```
[success] Total time: 17 s, completed Oct 2, 2018 11:54:27 AM
```

This will produce a JAR file located at `./digenome/target/scala-2.12/digenome.jar`.

## Running

The primary tool is run by invoking the following command with appropriate options:

```
java -Xmx8g -jar digenome/target/scala-2.12/digenome.jar IdentifyCutSites
```

### Input sequencing data

The tool expects as it's primary input a BAM file containing sorted, duplicate marked, aligned sequencing reads from a digenome experiment.  In our experience the following set of command yields results that work well with the software (though any reasonable WGS alignment pipeline should be compatible):

```
# Align the data and convert the output to BAM
bwa mem -p -R '<read-group-string>' -t 32 hg38.fa r1.fastq.gz r2.fastq.gz \
  | sambamba view --sam-input -l 1 -t 3 -f bam -o aligned.bam /dev/stdin

# Sort the data into coordinate order
sambamba sort -m 4G --tmpdir . -l 5 -t 16 -o sorted.bam aligned.bam

# Duplicate mark the data
picard -Duse_async_io_write_samtools=true -Xmx8g MarkDuplicates \
  I=sorted.bam O=deduped.bam M=dupe_metrics.txt CREATE_INDEX=true
```

The tools used in the above commands can be installed easily using [conda](https://conda.io/miniconda.html) and the included [conda-requirements.txt](conda-requirements.txt) as follows:

```
conda create -n digneome -c bioconda --file conda-requirements.txt
conda activate digenome
```

### IdentifyCutSites Options

|Name|Flag|Type|Description|Required?|Max Values|Default Value(s)|
|----|----|----|-----------|---------|----------|----------------|
|input|i|BAM file|Input BAM file of aligned, sorted and de-duped reads.|Y|Unlimited||
|output|o|Text file|Output tab-delimited file of cut site predictions.|Y|1||
|guide|g|String|The sequence of the guide used in the experiment.|N|1||
|pam-five-prime|p|String|The PAM sequence associated with the guide, if found at the 5' end of the guide.|N|1||
|pam-three-prime|P|String|The PAM sequence associated with the guide, if found at the 3' end of the guide|N|1||
|enzyme|e|String|The name of the cutting enzyme used. Only used to label outputs.|N|1||
|ref|r|FASTA file|The FASTA file for the genome to which the BAM is aligned. The reference must have been indexed, e.g. with `samtools faidx`.|Y|1||
|intervals|l|IntervalList File|An optional file of regions to restrict the analysis to, in Picard IntervalList format.|N|1||
|min-map-q|m|Integer|The minimum mapping quality for a read to be used in the analysis.|N|1|30|
|overhang||Integer|The expected overhang at the cut site (varies by enzyme). Positive for 5' overhang, negative for 3'.|N|1|0|
|max-offset|x|Integer|The maximum distance a read-start can be from the expected cut site and still be counted.|N|1|2|
|min-depth|N|Integer|Minimum sequencing depth at/near a cut site for the cut site to be emitted.|N|1|10|
|max-depth|X|Integer|Maximum sequencing depth at/near a cut site for the cut site to be emitted.|N|1|300|
|max-low-mapq-fraction|Q|Decimal|The maximum fraction of reads supporting a cut site with `mapq < min-map-q`.|N|1|0.3|
|min-forward-reads|F|Integer|Minimum number of forward strand reads supporting a cut site.|N|1|4|
|min-reverse-reads|R|Integer|Minimum number of reverse strand reads supporting a cut site.|N|1|4|
|min-supporting-reads|M|Integer|Minimum total number of reads supporting a cut site|.|.|1|8|
|min-supporting-fraction|S|Decimal|Minimum fraction of reads supporting a cut site [`cut reads / (cut reads + spanning reads)`].).|N|1|0.2|

### Example Invocation

An example invocation of the toolkit might look like:

```
java -Xmx8g -jar digenome/target/scala-2.12/digenome.jar IdentifyCutSites \
  --input=s1.aligned.sorted.deduped.bam \
  --output=s1.cut_sites.txt \
  --guide=ACGTTACGAAACGCCACGGGA \
  --pam-three-prime=NNGRRN \
  --enzyme=SauCas9 \
  --ref=/refs/hg38/hg38.fa  
```

### Output Format

The output file generated is a tab-delimited text file with the following columns:

|Column Name|Description|Example Value(s)|
|-----------|-----------|----------------|
|digenomitas\_version   |The version of the software that produced the results.|20181002-bd34ba2-dirty|
|sample                 |The name of the smaple used int he experiment (taken from the BAM file header).|EMX1-example-sample|
|guide                  |The guide sequenced used in the experiment.|GTTCTGTCCTCAGTAAAAGGTA|
|enzyme                 |The enzyme used in the experiment.|SauCas9|
|expected\_overhang     |The expected overhang of the enzyme (as passed to `--overhang` when running analysis.|0|
|window\_size           |The `window-size` option given to the tool when running the analysis.|5|
|mapq\_cutoff           |The mapq cutoff given to the tool when running the analysis.|30|
|chrom                  |The chromosome on which the putative cut site resides.|chr12|
|pos                    |The most likely position of the cut site.|88102160|
|strand                 |The strand of the cut site.|F|
|low\_mapq\_fraction    |The fraction of reads supporting the cut site that had low mapping quality.|0|
|forward\_starts        |The number of F strand reads supporting the cut site.|50|
|reverse\_starts        |The number of R strand reads supporting the cut site.|53|
|read\_depth            |The total number of reads at the cut site (including spanning reads).|104|
|read\_fraction\_cut    |The fraction of reads at the cut site that support the cut site (i.e. do not span it).|0.990385|
|read\_score            |A score representing confidence in the cut site. Scales with `read_depth` and `read_fraction_cut`.|132.129971|
|template\_depth        |The number of templates (aka inserts) at the cut site (including spanning templates).|105|
|template\_fraction\_cut|The fraction of templates that appear to be cut at the cut site.|0.980952|
|template\_score        |A template-based score representing confidence in the cut site. Scales with `template_depth` and `template_fraction_cut`.|145.65576|
|median\_overhang       |The median overhang of reads at the cut site.|0|
|overhang\_distribution |The distribution of read overhangs at the cut site. By default contains five semi-colon separated values that represent the count of reas with overhang = -2, -1, 0, 1, 2.|0;2;48;0;0|
|neighborhood\_ratio    |A measure of whether read-starts are suppressed in the neighborhood surrounding the cut site.  When there are no read starts in the area around a cut site this will be ~0, when there are the expected number of read starts for the observed coverage and random read start positions, the value will be ~1.  0 or close to 0 is expected for a true cut site.|0|
|aln\_start             |The 1-based start position of the best alignment of the guide+pam at or near the cut site.|88102142|
|aln\_end               |The 1-based end position of the best alignment of the guide+pam at or near the cut site.|88102170|
|aln\_strand            |The strand of the genome on which the best alignment of the guide+pam resides.|F|
|aln\_padded\_guide     |The padded guide sequence at the alignment position.|`GTTCTGTCCTCAGTAAAAGGTAnngrrn`|
|aln\_alignment\_string |A string representing the pairwise alignment of the guide+pam to the reference.|`||||||||||||||||||||||||||||`|
|aln\_padded\_target    |The padded sequence of the reference at the alignment position.|`GTTCTGTCCTCAGTAAAAGGTATAGAGT`|
|aln\_mismatches        |The number of mismatches in the alignment.|0|
|aln\_gap\_bases        |The number of gapped bases in the alignment (i.e. a gap of length 2 is counted as 2).|0|
|aln\_mm\_and\_gaps     |The sum of `aln_mismatches` and `aln_gap_bases`.|0|

The padded alignment strings are intended to be concatenated with line breaks to help visualize the alignment.  The following is an example of a more complicated alignment represented in the same format, with three mismatches and one bulge/indel.

```
GTTCTGTCCTCAGTAAAAGGTAnngrrn
|||||..||||||||||| |.|||||||
GTTCTACCCTCAGTAAAA-GAACTGGGT
```


## Example Dataset

Files are located at:

1. [digeome-emx1-example-data.bam](https://github.com/editasmedicine/digenome/blob/master/example-data/digeome-emx1-example-data.bam)
2. [digeome-emx1-example-data.bai](https://github.com/editasmedicine/digenome/blob/master/example-data/digeome-emx1-example-data.bai)


An example dataset is included in this repository.  The dataset is contains sequencing reads over nine 10kb wide regions from a sample treated with EMX1 guide `GTTCTGTCCTCAGTAAAAGGTA` and Staph aureus CAS9.  The regions include the expected on-target cut site and a number of putative off-target cut sites (including some likely false positives).

To analyze the data try running:

```
java -Xmx8g -jar digenome/target/scala-2.12/digenome.jar IdentifyCutSites \
  --input=digeome-emx1-example-data.bam \
  --output=cut_sites.txt \
  --guide=GTTCTGTCCTCAGTAAAAGGTA \
  --pam-three-prime=NNGRRN \
  --enzyme=SauCas9 \
  --ref=/refs/hg38/hg38.fa  
```
