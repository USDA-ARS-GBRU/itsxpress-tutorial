# ITSxpress: a QIIME 2 plugin to trim ITS sequences
_Adam R. Rivers - USDA Agricultural Research Service_

_Sveinn V. Einarsson - Dep. Microbiology and Cell Science, U. Florida & USDA Agricultural Research Service_  

## Updates to ITSxpress ##
ITSxpress v2.0.0 is a major update to the ITSxpress package.
Release Highlights
- pip (PyPI) is no longer maintained for ITSxpress>=v2.0.0
    - Qiime2 plugin version of ITSxpress is now part of the standalone package of ITSxpress
    -   No longer need to install q2-itsxpress separately and will be installed if Qiime2 is already installed, otherwise only the standalone version will be installed
    - Seperate Qiime2 plugin version (PyPI) of ITSxpress is no longer maintained after q2-itsxpress v1.8.1
    - Package can be installed from Github, and Bioconda

- Removed BBmap dependency
    - BBmap scripts are no longer used in the pipeline, including:
        - reformat.sh (interleaved files no longer supported)
        - bbmerge.sh (merging of paired-end reads now done with Vsearch --fastq_mergepairs)
             - merging of paired-end reads is different between Vsearch and BBmerge, so results may differ
- Updated dereplication step for newer versions of Vsearch
    - Dereplication step now done using Vsearch --fastx_uniques (derep_fulllength command no longer supports fastq files)
- Pacbio sequences are now supported if fastq file scores are in Illumina format

Bug Fixes
- Fixed bug where the q2-itsxpress plugin was not handling single-end reads correctly, and was looking for a reverse read file
- Fixed a bug that could cause crashes when an intermediate file was empty

Previous release highlights since tutorial was written:
- 1.8.1 is the final version that uses BBmap scripts. This version is still available on the EOL-1.8.1 branch of the ITSxpress repository
- Fixed version of dependencies for version 1.8.1, to maintain compatibility with Qiime2 2022.8
- Updated pip install config file to pyproject.toml
- Updated readme usage section to reference compatible Qiime2 version
- Added read count output to log file
- Added support for primer sets in the reverse orientation

## Background

The internally transcribed spacer (ITS) region a widely used phylogenetic marker for fungi and other taxa. Previous work by [Nilsson et al. (2009)](https://doi.org/10.1111/j.1574-6968.2009.01618.x) showed that removing the conserved regions around the ITS results in more accurate taxonomic classification. An existing program, [ITSx](https://doi.org/10.1111/2041-210X.12073), can trim FASTA sequences by matching HMM profiles to the ends of the flanking conserved genes. ITSxpress is designed to extend this technique to trim the FASTQ files needed for the newer exact sequence variant methods used by in [QIIME 2](https://qiime2.org/): [Dada2](https://doi.org/10.1038/nmeth.3869) and [Deblur](10.1128/mSystems.00191-16).  ITSxpress processes QIIME artifacts of the type `SampleData[PairedEndSequencesWithQuality]` or `SampleData[SequencesWithQuality]`.

The plugin:
1. Merges reads (if paired-end) using [VSEARCH](https://doi.org/10.7717/peerj.2584)
2. Temporarily clusters highly similar sequences that are common in amplicon data also using [VSEARCH](https://doi.org/10.7717/peerj.2584)
3. Identifies the ITS start and stop sites using  [Hmmsearch](https://doi.org/10.1371/journal.pcbi.1002195) on the representative sequences
4. Trims each original, merged sequence with quality scores, returning the merged or unmerged sequences with quality scores in a `.qza` file

ITSxpress speeds up the trimming of reads by a factor of 14-23 times on a 4-core computer by temporarily clustering highly similar sequences that are common in amplicon data and utilizing optimized parameters for Hmmsearch. For more information see [the paper](#New link to V2 paper).

ITSxpress is installed as a standalone package with a Qiime2 plugin, IF Qiime2 is available. Package can be installed from [Github](https://github.com/USDA-ARS-GBRU/itsxpress/), and [Bioconda](https://bioconda.github.io/recipes/itsxpress/README.html).

PIP (PyPI) is no longer maintained for ITSxpress>=v2.0.0.


## Installation

The instructions assume that you [installed QIIME 2 natively using Mamba (or Conda)](https://docs.qiime2.org/2023.2/install/native/) and are using ITSxpress version 2.0.0.

Activate the QIIME 2 Conda environment.

```
mamba activate qiime2-2023.2
```

Install ITSxpress using Bioconda. Be sure to install ITSxpress in the QIIME 2 environment, meaning you ran the step above first. ITSxpress no longer has a seperate Qiime plugin, it is included with the standalone version IF you already have Qiime2 installed in your environment.

```
mamba install -c bioconda itsxpress
```

In your QIIME2 environment, refresh the plugins.

```
qiime dev refresh-cache
```

Check to see if the ITSxpress plugin is installed. After running this command you should see a basic help menu.

```
qiime itsxpress
```


> Note: this tutorial was updated for ITSxpress 2.0.0 on 03/06/2023.



This tutorial walks the user through the first portion of a typical ITS workflow:
1. Trimming the ITS region with ITSxpress
2. Calling sequence variants with Dada2 or Deblur
3. Training the QIIME 2 classifier
4. Classifying the sequences taxonomically

For this tutorial we will be starting with two paired-end samples than have already been demultiplexed into forward and reverse FASTQ files. A manifest file which lists the samples, files and read orientation is also used. *The example manifest uses the $PWD variable to complete the path for your computer. If you have issues you can replace it with the direct path.*

## Below the following tutorial you can find an experimental walkthrough of Pacbio ITS analysis. ##

### Example data
We will be using data from two soil samples which have have their ITS1 region amplified with fungal primers. They have been subsampled to 10,000 read pairs for faster processing.

* [sample1_r1.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r1.fastq.gz) and [sample1_r2.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r2.fastq.gz)
* [sample2_r1.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r1.fastq.gz) and [sample2_r2.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r2.fastq.gz)
* A manifest file: [manifest.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/manifest.txt)
* A mapping file: [mapping.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/mapping.txt)

If you have the command line program `wget` you can download the data with these commands
```
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r1.fastq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r2.fastq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r1.fastq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r2.fastq.gz
wget https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/manifest.txt
wget https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/mapping.txt
```
### Import the sequence data

Make sure all the data files are in the same directory, then import the data into QIIME.

This step in the tutorial imports demultiplexed data into QIIME.

> NOTE: If you have multiplexed data in a format like `EMPPairedEndSequences` you will need to demultiplex it first using the `demux` plugin. For a paired-end example see example see [this tutorial.](https://docs.qiime2.org/2023.5/tutorials/atacama-soils/#atacama-demux)

```
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --input-format PairedEndFastqManifestPhred33\
  --input-path manifest.txt \
  --output-path sequences.qza
```
Run time: 5 seconds

* Output `sequences.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fsequences.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/sequences.qza)


We can see the quality of the data by running the summarize command.

```
qiime demux summarize \
  --i-data sequences.qza \
  --o-visualization sequences.qzv
```
Run time: 9 seconds

* Output `sequences.qzv` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fsequences.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/sequences.qzv)


### Trimming ITS samples with ITSxpress for Dada2

`ITSxpress trim-pair-output-unmerged` takes paired-end QIIME artifacts
`SampleData[PairedEndSequencesWithQuality]` for
trimming. It merges the reads, temporally clusters the reads, then looks for
the ends of the ITS region with Hmmsearch. HMM models are available for 18
different clades. `itsxpress trim-pair-output-unmerged` returns the unmerged, trimmed sequences. `itsxpress trim-pair-output-merged` returns merged, trimmed sequences. You can adjust the --p-cluster-id value, which is the percent identity for clustering reads range [0.995-1.0], set to 1 for exact de-replication. Default is 1.

```
qiime itsxpress trim-pair-output-unmerged\
  --i-per-sample-sequences sequences.qza \
  --p-region ITS2 \
  --p-taxa F \
  --p-cluster-id 0.995 \
  --p-threads 16 \
  --o-trimmed trimmed.qza
```
Run time: 15 seconds

```
qiime itsxpress trim-pair-output-unmerged\
  --i-per-sample-sequences sequences.qza \
  --p-region ITS2 \
  --p-taxa F \
  --p-cluster-id 1.0 \
  --p-threads 16 \
  --o-trimmed trimmed_exact.qza
  ```
Run time: 23 seconds
* Output `trimmed.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftrimmed.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/trimmed.qza)

### Use Dada2 to identify sequence variants

The trimmed sequences can be fed directly into Dada2 using the denoise-paired
command. Since Vsearch handled the merging and quality issues there is no need to trim or truncate the reads further. In this tutorial we have set a truncation length \ to 0 because the data quality was good. Be sure to examine the `sequences.qzv` file before deciding to hard trim your reads.

```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs trimmed_exact.qza \
  --p-trunc-len-r 0 \
  --p-trunc-len-f 0 \
  --output-dir dada2out
```
Run time: 10 seconds

* Output:

  1. `dada2out/denoising_stats.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Fdenoising_stats.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/denoising_stats.qza)
  2. `dada2out/representative_sequences.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Frepresentative_sequences.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/representative_sequences.qza)
  3.  `dada2out/table.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Ftable.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/table.qza)


  ### Summarize the data for visual inspection:
  
  ```
  qiime feature-table summarize \
    --i-table dada2out/table.qza \
    --o-visualization tableviz.qzv
  ```
  
  Run time: 5 seconds

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

There has been discussion about whether trimming the database matters for classification. The QIIME team found that trimming the UNITE database [does not result in better classification](https://docs.qiime2.org/2023.5/tutorials/feature-classifier/) when untrimmed reads are used and recommended using the untrimmed developer database. Since we are using the trimmed ITS region, this tutorial recommends using the __trimmed__ database but this has not yet been systematically compared.


```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sh_refs_qiime_ver9_dynamic_29.11.2022.fasta \
  --output-path unite.qza
```
Run time 10 seconds

* Output `unite.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/unite.qza)

Import the associated UNITE taxonomy file.
```
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path sh_taxonomy_qiime_ver9_dynamic_29.11.2022.txt \
--output-path unite-taxonomy.qza
```
Run time 5 seconds

* Output `unite-taxonomy.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite-taxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/unite-taxonomy.qza)


### Train the QIIME classifier

QIIME provides its own naive Bayes classifier similar to [RDP](https://dx.doi.org/10.1128%2FAEM.00062-07) from the python package [SciKit Learn](http://scikit-learn.org/stable/modules/naive_bayes.html). Before using it the classifier must be trained using the data you just imported.  

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite.qza \
  --i-reference-taxonomy unite-taxonomy.qza \
  --o-classifier classifier.qza
```
Run time: 19 minutes

* Output `classifier.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fclassifier.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/classifier.qza)

### Note: If you come accross an error saying "[Errno 28] No space left on device". You may need to empty your qiime2 temp directory. See discussion here: [Temp directory space issue](https://forum.qiime2.org/t/is-there-something-unusual-about-how-the-tmpdir-is-being-set-in-2023-2/26402/8)

### Classify the sequence variants

Once the classifier is trained sequences can be classified.

```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads dada2out/representative_sequences.qza \
  --o-classification taxonomy.qza
```
Run time: 49 seconds

* Output `taxonomy.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxonomy.qza)

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



### EXPERIMENTAL - Below is a tutorial for Pacbio data. However, due to the variable nature of PACBIO read scoring and file types this is not maintained. This may be revisted in the future. ###

You can find the paper this test sample is from here: [Runnel et al. 2022] (https://doi.org/10.1111/1755-0998.13663)

### Example Pacbio data
We will be using data from amplifying fungal soil samples which have have the entire ITS1/5.8S/ITS2 region amplified with fungal primers.

* [pacbio_samp1.fastq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/pacbio/pacbio_samp1.fastq.gz)
* [pacbio_samp1_reformat.fastq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/pacbio/pacbio_samp1_reformat.fastq.gz)
* A manifest file: [manifest.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/pacbio/manifest_pacbio.txt)
* A mapping file: [mapping.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/pacbio/mapping_pacbio.txt)

If you have the command line program `wget` you can download the data with these commands
```
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/pacbio/pacbio_samp1.fastq.gz
wget https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/pacbio/pacbio_samp1_reformat.fastq.gz
wget https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/pacbio/manifest_pacbio.txt
wget https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/pacbio/mapping_pacbio.txt
```
### Import the sequence data

Make sure all the data files are in the same directory, then import the data into QIIME.

This step in the tutorial imports demultiplexed data into QIIME.

> NOTE: If you have multiplexed data in a format like `EMPPairedEndSequences` you will need to demultiplex it first using the `demux` plugin. For a paired-end example see example see [this tutorial.](https://docs.qiime2.org/2023.5/tutorials/atacama-soils/#atacama-demux)

> NOTE: The fastq file from Pacbio has had the Pacbio scoring reformated to Illumina scoring convention, using bbmap reformat.sh:

```
reformat.sh in=./pacbio/pacbio_samp1.fastq.gz out=pacbio_samp1_reformat.fastq.gz mincalledquality=2 maxcalledquality=41 qin=33
```

```
qiime tools import \
--input-path ./pacbio/manifest_pacbio.txt \
--input-format SingleEndFastqManifestPhred33 \
--type SampleData[SequencesWithQuality] \
--output-path ./pacbio/sequences_pacbio.qza
```
Run time: 5 seconds

* Output `sequences.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fsequences_pacbio.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/sequences_pacbio.qza)


We can see the quality of the data by running the summarize command.

```
qiime demux summarize   \
--i-data ./pacbio/sequences_pacbio.qza   \
--o-visualization ./pacbio/sequences_pacbio.qzv
```
Run time: 9 seconds

* Output `sequences.qzv` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fsequences_pacbio.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/sequences_pacbio.qzv)


### Trimming ITS samples with ITSxpress for Dada2

Since Pacbio long range sequencing allows for amplification of the entire ITS1/5.8S/ITS2 region, you can change the --p-region command below between ITS1, ITS2, and ALL to trim to the corresponding region and see the difference in taxonomic output.


```
qiime itsxpress trim-single \
--i-per-sample-sequences ./pacbio/sequences_pacbio.qza \
--p-region ITS2 \
--p-taxa F \
--p-cluster-id 1.0 \
--p-threads 16 \
--o-trimmed ./pacbio/trimmed_pacbio.qza
```
Run time: 15 seconds


* Output `trimmed.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Ftrimmed_pacbio.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/trimmed_pacbio.qza)

### Use Dada2 to identify sequence variants

The trimmed sequences can be fed directly into Dada2 using the denoise-single
command. Be sure to examine the `sequences_pacbio.qzv` file before deciding to hard trim your reads.

```
qiime dada2 denoise-single \
--i-demultiplexed-seqs ./pacbio/trimmed_pacbio.qza \
--p-trunc-len 0 \
--p-n-threads 16 \
--p-max-ee 20 \
--output-dir ./pacbio/dada2out
```
Run time: 10 seconds

* Output:

  1. `dada2out/denoising_stats.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fdada2out%2Fdenoising_stats.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/dada2out/denoising_stats.qza)
  2. `dada2out/representative_sequences.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fdada2out%2Frepresentative_sequences.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/dada2out/representative_sequences.qza)
  3.  `dada2out/table.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fdada2out%2Ftable.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/dada2out/table.qza)


  ### Summarize the data for visual inspection:
  
  ```
  qiime feature-table summarize \
    --i-table ./pacbio/dada2out/table.qza \
    --o-visualization ./pacbio/tableviz.qzv
  ```
  
  Run time: 5 seconds

  * Output `tableviz.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Ftableviz.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/tableviz.qzv)


> Deblur is an alternative option for read correction. This tutorial uses Dada2 Because deblur requires uniform length reads, specified by the --p-trim-length flags, and ITS regions vary considerably in length. Tests across a range of trim lengths using Deblur yielded fewer sequence variants.

### The following section is exactly the same as the above tutorial using paired-end example data. Feel free to use classifier developed above.

### Download reference data from UNITE for fungal classification

First download the newest [UNITE database for QIIME ](https://unite.ut.ee/repository.php) and unzip the file.

```
wget https://files.plutof.ut.ee/doi/0A/0B/0A0B25526F599E87A1E8D7C612D23AF7205F0239978CBD9C491767A0C1D237CC.zip
unzip 0A0B25526F599E87A1E8D7C612D23AF7205F0239978CBD9C491767A0C1D237CC.zip
```

### Import the latest UNITE data into QIIME 2:

Import the UNITE sequences for the smaller dataset selected with dynamic thresholds determined by fungal experts.

There has been discussion about whether trimming the database matters for classification. The QIIME team found that trimming the UNITE database [does not result in better classification](https://docs.qiime2.org/2023.5/tutorials/feature-classifier/) when untrimmed reads are used and recommended using the untrimmed developer database. Since we are using the trimmed ITS region, this tutorial recommends using the __trimmed__ database but this has not yet been systematically compared.


```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sh_refs_qiime_ver9_dynamic_29.11.2022.fasta \
  --output-path ./pacbio/unite.qza
```
Run time 10 seconds

* Output `unite.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/unite.qza)

Import the associated UNITE taxonomy file.
```
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path ./pacbio/sh_taxonomy_qiime_ver9_dynamic_29.11.2022.txt \
--output-path ./pacbio/unite-taxonomy.qza
```
Run time 5 seconds

* Output `unite-taxonomy.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite-taxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/unite-taxonomy.qza)


### Train the QIIME classifier

QIIME provides its own naive Bayes classifier similar to [RDP](https://dx.doi.org/10.1128%2FAEM.00062-07) from the python package [SciKit Learn](http://scikit-learn.org/stable/modules/naive_bayes.html). Before using it the classifier must be trained using the data you just imported.  

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads ./pacbio/unite.qza \
  --i-reference-taxonomy ./pacbio/unite-taxonomy.qza \
  --o-classifier ./pacbio/classifier.qza
```
Run time: 19 minutes

* Output `classifier.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fclassifier.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/classifier.qza)

### Note: If you come accross an error saying "[Errno 28] No space left on device". You may need to empty your qiime2 temp directory. See discussion here: [Temp directory space issue](https://forum.qiime2.org/t/is-there-something-unusual-about-how-the-tmpdir-is-being-set-in-2023-2/26402/8)

### Classify the sequence variants

Once the classifier is trained sequences can be classified.

```
qiime feature-classifier classify-sklearn \
  --i-classifier ./pacbio/classifier.qza \
  --i-reads ./pacbio/dada2out/representative_sequences.qza \
  --o-classification ./pacbio/taxonomy.qza
```
Run time: 49 seconds

* Output `taxonomy.qza` [View](https://view.qiime2.org/visualization/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Ftaxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/taxonomy.qza)

### Summarize the results
Summarize the results for visualization in the QIIME 2 viewer.

```
qiime metadata tabulate \
  --m-input-file ./pacbio/taxonomy.qza \
  --o-visualization ./pacbio/taxonomy.qzv
```

Run time: 4 seconds

* Output `taxonomy.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Ftaxonomy.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/taxonomy.qzv)

### Create an interactive bar plot figure


```
qiime taxa barplot \
  --i-table ./pacbio/dada2out/table.qza  \
  --i-taxonomy ./pacbio/taxonomy.qza \
  --m-metadata-file ./pacbio/mapping_pacbio.txt \
  --o-visualization ./pacbio/taxa-bar-plots.qzv
```
Run time: 4 seconds

* Output `taxa-bar-plots.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Ftaxa-bar-plots.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/taxa-bar-plots.qzv)


### The benefit of long read Pacbio analysis is the ability to see differences in taxonomy when trimming to different regions. The example above trimmed to the ITS1 region in the ITSxpress trim-single command, but you can do the same with the ITS2 and ALL regions and compare the taxa-bar-plots:

### ITS2
* Output `taxa-bar-plots.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fpacbio%2FITS2%2Ftaxa-bar-plots_ITS2.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/pacbio/ITS2/taxa-bar-plots_ITS2.qzv)

### ALL

* Output `taxa-bar-plots.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fpacbio%2Fpacbio%2FALL%2Ftaxa-bar-plots.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/pacbio/pacbio/ALL/taxa-bar-plots.qzv)

You'll see that including more conserved regions with the ALL trim command in ITSxpress, some of the finer scale variability is lost. This tutorial provides the basic process for analyzing ITS sequences. The data is now in a form where it can be analyzed further using many of the other methods provided by QIIME 2.

## Citation information for ITSxpress
* Rivers AR, Weber KC, Gardner TG et al. ITSxpress: Software to rapidly trim internally transcribed spacer sequences with quality scores for marker gene analysis [version 1; referees: awaiting peer review]. F1000Research 2018, 7:1418
[doi: 10.12688/f1000research.15704.1](https://doi.org/10.12688/f1000research.15704.1)

* ITSxpress software: [DOI:10.5281/zenodo.1304348](https://doi.org/10.5281/zenodo.1304348)
