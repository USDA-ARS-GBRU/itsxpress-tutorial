# Q2-ITSxpress: a QIIME 2 plugin to trim ITS sequences
_Adam R. Rivers - USDA Agricultural Research Service_

## Background

The internally transcribed spacer (ITS) region a widely used phylogenetic marker for fungi and other taxa. Previous work by [Nilsson et al. (2009)](https://doi.org/10.1111/j.1574-6968.2009.01618.x) showed that removing the conserved regions around the ITS results in more accurate taxonomic classification. An existing program, [ITSx](https://doi.org/10.1111/2041-210X.12073), can trim FASTA sequences by matching HMM profiles to the ends of the flanking conserved genes. ITSxpress is designed to extend this technique to trim the FASTQ files needed for the newer exact sequence variant methods used by in [QIIME 2](https://qiime2.org/): [Dada2](https://doi.org/10.1038/nmeth.3869) and [Deblur](10.1128/mSystems.00191-16).  ITSxpress processes QIIME artifacts of the type `SampleData[PairedEndSequencesWithQuality]` or `SampleData[SequencesWithQuality]`.

The plugin:
1. Merges reads (if paired-end) using [BBMerge](https://doi.org/10.1371/journal.pone.0185056)
2. Temporarily clusters highly similar sequences that are common in amplicon data using [VSEARCH](https://doi.org/10.7717/peerj.2584)
3. Identifies the ITS start and stop sites using  [Hmmsearch](https://doi.org/10.1371/journal.pcbi.1002195) on the representative sequences
4. Trims each original, merged sequence with quality scores, returning the merged or unmerged sequences with quality scores in a `.qza` file

ITSxpress speeds up the trimming of reads by a factor of 14-23 times on a 4-core computer by temporarily clustering highly similar sequences that are common in amplicon data and utilizing optimized parameters for Hmmsearch. For more information see [the paper](https://doi.org/10.12688/f1000research.15704.1).

ITSxpress is also available as a stand-alone software package from [Github](https://github.com/USDA-ARS-GBRU/itsxpress/), [PyPi](https://pypi.org/project/itsxpress/) and [Bioconda](https://bioconda.github.io/recipes/itsxpress/README.html).



## Installation

The instructions assume that you [installed QIIME 2 natively using Conda](https://docs.qiime2.org/2018.6/install/native/) and are using ITSxpress version 1.7.0.

Activate the QIIME 2 Conda environment.

```
source activate qiime2-2018.8
```

Install ITSxpress using Bioconda and Q2-itsxpress using pip. Be sure to install ITSxpress and Q2-itsxpress in the QIIME 2 environment, meaning you ran the **step above first**.

```
conda install -c bioconda itsxpress
pip install q2-itsxpress
```

In your QIIME2 environment, refresh the plugins.

```
qiime dev refresh-cache
```

Check to see if the ITSxpress plugin is installed. After running this command you should see a basic help menu.

```
qiime itsxpress
```

Download bbmap into the current directory (or Documents/anywhere) and follow this tutorial for your OS: https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/ 

Put bbmap folder into path, for example if you've downloaded bbmap into your current directory:
PATH="./bbmap/:${PATH}"
export PATH

If you receive an error like "No such file or directory: 'bbmerge.sh'" then your computer doesn't know the path to the bbmap folder and the previous step needs to be revisited. 

## Tutorial

> Note: this tutorial was updated for ITSxpress 1.7.0 on 08/13/2018.
> Recommendations about how to trim reads for use by Dada2 have changed.


This tutorial walks the user through the first portion of a typical ITS workflow:
1. Trimming the ITS region with ITSxpress
2. Calling sequence variants with Dada2 or Deblur
3. Training the QIIME 2 classifier
4. Classifying the sequences taxonomically

For this tutorial we will be starting with two paired-end samples than have already been demultiplexed into froward and reverse FASTQ files. A manifest file which lists the samples, files and read orientation is also used. *The example manifest uses the $PWD variable to complete the path for your computer. If you have issues you can replace it with the direct path.*

### Example data
We will be using data from two soil samples which have have their ITS1 region amplified with fungal primers. They have been subsampled to 10,000 read pairs for faster processing.

* [sample1_r1.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r1.fq.gz) and [sample1_r2.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r2.fq.gz)
* [sample2_r1.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r1.fq.gz) and [sample2_r2.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r2.fq.gz)
* A manifest file: [manifest.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/manifest.txt)
* A mapping file: [mapping.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/mapping.txt)

If you have the command line program `wget` you can download the data with these commands
```
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r1.fq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r2.fq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r1.fq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r2.fq.gz
wget https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/manifest.txt
wget https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/mapping.txt
```
### Import the sequence data

Make sure all the data files are in the same directory, then import the data into QIIME.

This step in the tutorial imports demultiplexed data into QIIME.

> NOTE: If you have multiplexed data in a format like `EMPPairedEndSequences` you will need to demultiplex it first using the `demux` plugin. For a paired-end example see example see [this tutorial.](https://docs.qiime2.org/2018.6/tutorials/atacama-soils/#atacama-demux)

```
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33\
  --input-path manifest.txt \
  --output-path sequences.qza
```
Run time: 4 seconds

* Output `sequences.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fsequences.qzv&type=html)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/sequences.qza)


We can see the quality of the data by running the summarize command.

```
qiime demux summarize \
  --i-data sequences.qza \
  --o-visualization sequences.qzv
```
Run time: 4 seconds

* Output `sequences.qzv` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fsequences.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/sequences.qzv)


### Trimming ITS samples with Q2-ITSxpress for Dada2

`ITSxpress trim-pair-output-unmerged` takes paired-end QIIME artifacts
`SampleData[PairedEndSequencesWithQuality]` for
trimming. It merges the reads, temporally clusters the reads, then looks for
the ends of the ITS region with Hmmsearch. HMM models are available for 18
different clades. `itsxpress trim-pair-output-unmerged` returns the unmerged, trimmed sequences.

```
qiime itsxpress trim-pair-output-unmerged\
  --i-per-sample-sequences sequences.qza \
  --p-region ITS1 \
  --p-taxa F \
  --o-trimmed trimmed.qza
```

```
qiime itsxpress trim-pair-output-unmerged\
  --i-per-sample-sequences sequences.qza \
  --p-region ITS1 \
  --p-taxa F \
  --p-cluster-id 1.0 \
  --p-threads 2 \
  --o-trimmed trimmed_exact.qza
  ```
Run time: 2 minutes 45 seconds

* Output `trimmed.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftrimmed.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/trimmed.qza)

### Use Dada2 to identify sequence variants

The trimmed sequences can be fed directly into Dada2 using the denoise-paired
command. Since BBmerge handled the merging and quality issues there is no need to trim or truncate the reads further. In this tutorial we have set a truncation length \ to 0 because the data quality was good. Be sure to examine the `sequences.qzv` file before deciding to hard trim your reads.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed.qza \
  --p-trunc-len-r 0 \
  --p-trunc-len-f 0 \
  --output-dir dada2out
```

```
qiime dada2 denoise-single \
  --i-demultiplexed-seqs trimmed.qza \
  --p-trunc-len 0 \
  --output-dir dada2wrongout
```
Run time: 1 minute

* Output:

  1. `dada2out/denoising_stats.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Fdenoising_stats.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/denoising_stats.qza)
  2. `dada2out/representative_sequences.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Frepresentative_sequences.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/representative_sequences.qza)
  3.  `dada2out/table.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Ftable.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/table.qza)

  ### Summarize the data for visual inspection:
  ```
  qiime feature-table summarize \
    --i-table dada2out/table.qza \
    --o-visualization tableviz.qzv
  ```
  Run time: 4 seconds

  * Output `tableviz.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftableviz.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/tableviz.qzv)


> Deblur is an alternative option for read correction. This tutorial uses Dada2 Because deblur requires uniform length reads, specified by the --p-trim-length flags, and ITS regions vary considerably in length. Tests across a range of trim lengths using Deblur yielded fewer sequence variants.


### Download reference data from UNITE for fungal classification

First download the newest [UNITE database for QIIME ](https://unite.ut.ee/repository.php) and unzip the file.

```
wget https://files.plutof.ut.ee/doi/0A/0B/0A0B25526F599E87A1E8D7C612D23AF7205F0239978CBD9C491767A0C1D237CC.zip
unzip 0A0B25526F599E87A1E8D7C612D23AF7205F0239978CBD9C491767A0C1D237CC.zip
```

### Import the latest UNITE data into QIIME 2:

Import the UNITE sequences for the smaller dataset selected with dynamic thresholds determined by fungal experts.

There has been discussion about whether trimming the database matters for classification. The QIIME team found that trimming the UNITE database [does not result in better classification](https://docs.qiime2.org/2018.6/tutorials/feature-classifier/) when untrimmed reads are used and recommended using the untrimmed developer database. Since we are using the trimmed ITS region, this tutorial recommends using the __trimmed__ database but this has not yet been systematically compared.


```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sh_refs_qiime_ver7_dynamic_01.12.2017.fasta \
  --output-path unite.qza
```
Run time 7 seconds

* Output `unite.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/unite.qza)

Import the associated UNITE taxonomy file.
```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path sh_taxonomy_qiime_ver7_dynamic_01.12.2017.txt \
  --output-path unite-taxonomy.qza
```
Run time 4 seconds

* Output `unite-taxonomy.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite-taxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/unite-taxonomy.qza)


### Train the QIIME classifier

QIIME provides its own naive Bayes classifier similar to [RDP](https://dx.doi.org/10.1128%2FAEM.00062-07) from the python package [SciKit Learn](http://scikit-learn.org/stable/modules/naive_bayes.html). Before using it the classifier must be trained using the data you just imported.  

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite.qza \
  --i-reference-taxonomy unite-taxonomy.qza \
  --o-classifier classifier.qza
```
Run time: 5 minutes

* Output `classifier.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fclassifier.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/classifier.qza)

### Classify the sequence variants

Once the classifier is trained sequences can be classified.

```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads dada2out/representative_sequences.qza \
  --o-classification taxonomy.qza
```
Run time: 1.5 minutes

* Output `taxonomy.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxonomy.qza)

### Summarize the results
Summarize the results for visualization in the QIIME 2 viewer.

```
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

Run time: 4 seconds

* Output `taxonomy.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxonomy.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxonomy.qzv)

### Create an interactive bar plot figure


```
qiime taxa barplot \
  --i-table dada2out/table.qza  \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file mapping.txt \
  --o-visualization taxa-bar-plots.qzv
```
Run time: 4 seconds

* Output `taxonomy.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxa-bar-plots.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxa-bar-plots.qzv)

This tutorial provides the basic process for analyzing ITS sequences. The data is now in a form where it can be analyzed further using many of the other methods provided by QIIME 2.

## Citation information for ITSxpress
* Rivers AR, Weber KC, Gardner TG et al. ITSxpress: Software to rapidly trim internally transcribed spacer sequences with quality scores for marker gene analysis [version 1; referees: awaiting peer review]. F1000Research 2018, 7:1418
[doi: 10.12688/f1000research.15704.1](https://doi.org/10.12688/f1000research.15704.1)

* ITSxpress software: [DOI:10.5281/zenodo.1304348](https://doi.org/10.5281/zenodo.1304348)

* ITSxpress QIIME 2 plugin: [DOI:10.5281/zenodo.1317578](https://doi.org/10.5281/zenodo.1317578)
