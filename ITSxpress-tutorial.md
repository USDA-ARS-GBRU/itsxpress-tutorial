# Q2-ITSxpress a QIIME 2 plugin to trim ITS sequences

## Background

The internally transcribed spacer (ITS) region a widely used
phylogenetic marker for fungi and other taxa. Previous work by [Nilsson et al.
(2009)](https://doi.org/10.1111/j.1574-6968.2009.01618.x) showed that removing
the conserved regions around the ITS results in more accurate taxonomic
classification. An existing program,
[ITSx](https://doi.org/10.1111/2041-210X.12073), can trim  FASTA sequences by
matching HMM profiles to the ends of the flanking conserved genes. I Designed
ITSxpress to extend this technique to trim the FASTQ files needed for the newer
exact sequence variant methods used by in [QIIME 2](https://qiime2.org/):
[Dada2](https://doi.org/10.1038/nmeth.3869) and
[Deblur](10.1128/mSystems.00191-16). ITSxpress processes QIIME artifacts of the
type "paired-end sequences with quality" or "sequences with quality".

The plugin: 1. Merges reads (if paired-ended) using
[BBMerge](https://doi.org/10.1371/journal.pone.0185056) 2. Temporarily clusters
highly similar sequences that are common in amplicon data using
[VSEARCH](https://doi.org/10.7717/peerj.2584) 3. Identifies the ITS start and
stop sites using  [Hmmsearch](https://doi.org/10.1371/journal.pcbi.1002195) on
the representative sequences 4. Trims each original, merged sequence with
quality scores, returning the merged sequences with quality scores in a `.qza`
file

ITSxpress speeds up the trimming of reads by a factor of 14-23 times on a 4-core
computer by temporarily clustering highly similar sequences that are common in
amplicon data and utilizing optimized parameters for Hmmsearch.

ITSxpress is also available as a stand-alone software package from
[Github](https://github.com/USDA-ARS-GBRU/itsxpress/),
[PyPi](https://pypi.org/project/itsxpress/) and
[Bioconda](https://bioconda.github.io/recipes/itsxpress/README.html).



## Installation

The instructions assume that you [installed QIIME natively using the
Conda](https://docs.qiime2.org/2018.6/install/native/)

Activate the QIIME 2 Conda environment.

```
source activate qiime2-2018.6
```

Install Q2_itsxpress using BioConda. Be sure to install Q2_itsxpres in the QIIME
2 environment, meaning you ran the step above first.

```
conda config --add channels bioconda
conda install q2-itsxpress
```

In your QIIME2 environment, refresh the plugins.

```
qiime dev refresh-cache
```

Check to see if the ITSxpress plugin is installed. After running this command
you should see a basic help menu.

```
qiime itsxpress
```

## Tutorial
This tutorial walks the user through a typical workflow of importing sequences
and environmental data. This involves:
1. trimming the ITS region with ITSxpress
2. calling sequence variants with Dada2
3. training the Qiime Classifier
4. classifying the sequences Taxonomically

### Example data
We will be using two example soil samples from the Fungal ITS1 region. They have been
subsampled 10,000 read pairs for faster processing.

https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r2.fq.gz
* [sample1_r1.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r1.fq.gz) and [sample1_r2.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample1_r2.fq.gz)
* [sample2_r1.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r1.fq.gz) and [sample2_r2.fq.gz](https://github.com/USDA-ARS-GBRU/itsxpress-tutorial/raw/master/data/sample2_r2.fq.gz)
* A manifest file: [manifest.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/manifest.txt)
* A mapping file: [mapping.txt](https://raw.githubusercontent.com/USDA-ARS-GBRU/itsxpress-tutorial/master/data/mapping.txt)

### Workflow


For this tutorial we will be starting with two paired-end samples than have
already been demultiplexed into froward and reverse FASTQ files. A manifest
file which lists the samples, files and read orientation is also used. If you
have copied the examples files to your computer you will need to change the path
in the manifest to the complete path to your data.

#### Importing data

**Depending on where you downloaded the sample files you may need to edit the
manifest file to include the absolute path to your FASTQ files.**  Once you have
edited the files you can import the data into the `sequences.qza` file like this:

```
qiime tools import \
  --type SampleData[PairedEndSequencesWithQuality] \
  --source-format PairedEndFastqManifestPhred33\
  --input-path manifest.txt \
  --output-path sequences.qza
```
Run time: 5 seconds

* Output `sequences.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fsequences.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/sequences.qza)

#### Trimming ITS samples with Q2-ITSXpress

ITSxpress takes paired-end or single-end QIIME artifacts for trimming. It merges
reads (if paired-end), temporally clusters the reads, then looks for the ends of
the ITS region with  Hmmsearch. HMM models are available for 18 different
clades. For a full listing run `qiime itsxpress trim-pair --help`. ITSxpress
only returns merged sequences.

```
qiime itsxpress trim-pair\
  --i-per-sample-sequences sequences.qza \
  --p-region ITS1 \
  --p-taxa F \
  --p-threads 2 \
  --o-trimmed trimmed.qza
```
Run time: 1 minute

* Output `trimmed.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftrimmed.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/trimmed.qza)

#### Use Dada2 to identify sequence variants

The merged sequences can be fed directly into Dada2 using the denoise-single
option, even if the reads were paired ended to begin with.  Since BBmerge handled
the merging and quality issues there is no need to trim
or truncate the reads further.

```
time qiime dada2 denoise-single \
  --i-demultiplexed-seqs trimmed.qza \
  --p-trunc-len 0 \
  --p-n-threads 2 \
   --output-dir dada2out
```
Run time: 35 seconds

* Output:
  1.  `dada2out/table.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Ftable.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/table.qza)
  2. `dada2out/denoising_stats.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Fdenoising_stats.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/denoising_stats.qza)
  3. `dada2out/representative_sequences.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fdada2out%2Frepresentative_sequences.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/dada2out/representative_sequences.qza)

#### Summarize the data for visual inspection:
```
qiime feature-table summarize \
  --i-table dada2out/table.qza \
  --o-visualization tableviz.qzv
```
Run time: 4 seconds

* Output `tableviz.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftableviz.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/tableviz.qzv)

### Assigning fungal taxomomy

#### Download reference data from UNTIE

First download the newest [UNITE database for QIIME ](https://unite.ut.ee/repository.php) and unzip the file.

```
wget https://files.plutof.ut.ee/doi/0A/0B/0A0B25526F599E87A1E8D7C612D23AF7205F0239978CBD9C491767A0C1D237CC.zip
unzip 0A0B25526F599E87A1E8D7C612D23AF7205F0239978CBD9C491767A0C1D237CC.zip
```

#### Import the latest Unite data into UNITE:

Import the UNITE sequences for the smaller dataset selected with dynamic
thresholds determined by fungal experts.

There has been discussion about whether
trimming the database (as opposed to  the reads) matters for classification. The
QIIME team found that trimming the UNITE database [does not result in better
classification](https://docs.qiime2.org/2018.6/tutorials/feature-classifier/)
when untrimmed reads are used and recommended using the untrimmed developer
database. Since we are using the trimmed ITS region this tutorial recommends
using the _trimmed_ database but this has not yet been systematically compared.


```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path sh_refs_qiime_ver7_dynamic_01.12.2017.fasta \
  --output-path unite.qza
```
Run time 8 seconds

* Output `unite.qza` [View]https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/unite.qza)

Import the associated UNITE taxonomy file.
```
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --source-format HeaderlessTSVTaxonomyFormat \
  --input-path sh_taxonomy_qiime_ver7_dynamic_01.12.2017.txt \
  --output-path unite-taxonomy.qza
```
Run time 4 seconds

* Output `unite-taxonomy.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Funite-taxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/unite-taxonomy.qza)

#### Train the QIIME classifier

QIIME provides its own Naive Bayes classifier similar to
[RDP](https://dx.doi.org/10.1128%2FAEM.00062-07) from a python package called
[SciKit Learn](http://scikit-learn.org/stable/modules/naive_bayes.html). Before
using it  the classifier must be trained using the data you just imported.  

```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads unite.qza \
  --i-reference-taxonomy unite-taxonomy.qza \
  --o-classifier classifier.qza
```
Run time: 5 minutes

* Output `classifier.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Fclassifier.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/classifier.qza)

#### Classify the sequence variants

Once the classifier is trained sequences can be classified.

```
qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads dada2out/representative_sequences.qza \
  --o-classification taxonomy.qza
```
Run time: 1.5 minutes

* Output `taxonomy.qza` [View](https://view.qiime2.org/peek/?src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxonomy.qza)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxonomy.qza)

#### Summarize the results
Summarize the results for visualization in the QIIME viewer

```
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
```

Run time: 4 seconds

* Output `taxonomy.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxonomy.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxonomy.qzv)

#### Create barplot figures


```
qiime taxa barplot \
  --i-table dada2out/table.qza  \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file mapping.txt \
  --o-visualization taxa-bar-plots.qzv
```
Run time: 4 seconds

* Output `taxonomy.qzv` [View](https://view.qiime2.org/visualization/?type=html&src=https%3A%2F%2Fusda-ars-gbru.github.io%2Fitsxpress-tutorial%2Fdata%2Ftaxa-bar-plots.qzv)  \| [Download](https://usda-ars-gbru.github.io/itsxpress-tutorial/data/taxa-bar-plots.qzv)



## Citation information for ITSxpress
* Rivers, A. R., Weber K. C., Gardner T. G., Liu, S., Armstrong, S. D. 2018.
ITSxpress: Software to rapidly trim internally transcribed spacer sequences
with quality scores for marker gene analysis. F1000 Research. In press.

* ITSxpress software:
[DOI:10.5281/zenodo.1304348](https://doi.org/10.5281/zenodo.1304348)

* ITSxpress QIIME 2 plugin:
[DOI:10.5281/zenodo.1317578](https://doi.org/10.5281/zenodo.1317578)
