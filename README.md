# Combinatorial Dual Indexes Metabarcoding Tutorial
Complete tutorial for demultiplexing and analysing paired-end reads with combinatorial dual indexes.

### Requirements
- [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)
- [QIIME 2](https://docs.qiime2.org/2022.2/install/)

---

# Demultiplexing with Cutadapt

### adapters.fasta preparation
Create 2 FASTA files that list the adapters, one for the forward reads, the other for the reverse reads.  
The FASTA files must follow the format below:
```
>Adapter1 name
Adapter1_sequence
>Adapter2 name
Adapter2_sequence
```
For our particular type of adapter (non-internal 5'), an `X` must precede each adapter sequence. For other types of adapters, check the [Cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types).  
Furthermore, each adapter has a fixed number of N bases (from 2 to 4) ligated to the 5' end, so the final result will look something like this:
```
>adapter01
XNNacactgtg
>adapter02
XNNNctctgaga
>adapter03
XNNNNcgatcgat
```

## Let's Start Demultiplexing
Because forward and reverse reads are mixed together in R1 and R2, we'll run Cutadapt twice, once in one direction, once in the other.

### First run
```
cutadapt \
  -e 0 \ # maximum error rate
  -O 8 \ # minimum overlap (length of the adapter)
  --no-indels \ # no insertions or deletions
  --action=retain \ # the read is trimmed, but the adapter sequence itself is not removed
  -g file:barcode_F.fasta \ # FASTA with list of 5' fwd adapters to search
  -G file:barcode_R.fasta \ # FASTA with list of 5' rev adapters to search
  -o {name1}-{name2}.run1.1.fastq.gz \ # output R1
  -p {name1}-{name2}.run1.2.fastq.gz \ # output R2
  R1.fastq.gz \ # input R1
  R2.fastq.gz \ # input R2
  -j 0 # number of cores, -j 0 to automatically detect the number of cores
```
`barcode_F.fasta` and `barcode_R.fasta` are the FASTA files created in the previous section, `R1.fastq.gz` and `R2.fastq.gz` are our input data to be demultiplexed.  
`{name1}` will be replaced with the name of the best-matching R1 adapter and `{name2}` will be replaced with the name of the best-matching R2 adapter.
If there was no match of an R1 adapter, `{name1}` is set to `unknown`, and if there is no match of an R2 adapter, `{name2}` is set to `unknown`.

### Merge all the `unknown` R1 into a new file and do the same for R2
```
cat *unknown*.1.fastq.gz > unknown.1.fastq.gz
cat *unknown*.2.fastq.gz > unknown.2.fastq.gz
```
### Remove the original `unknown` reads
```
rm *unknown*.run1.*.fastq.gz
```

### Second run
```
cutadapt \
  -e 0 \
  -O 8 \
  --no-indels \
  --action=retain \
  -g file:barcode_F.fasta \
  -G file:barcode_R.fasta \
  -o {name1}-{name2}.run2.1.fastq.gz \ # run2 instead of run1
  -p {name1}-{name2}.run2.2.fastq.gz \ # run2 instead of run1
  unknown.2.fastq.gz \ # input R2 to search for forward reads
  unknown.1.fastq.gz \ # input R1 to search for reverse reads
  -j 0
```
### Remove `unknown` reads
```
rm *unknown*.fastq.gz
```

## Merge the output files of the 2 runs for each combination of adapters
Create a txt file (e.g. `combo`) that lists all the combinations and run the following command:
```
while read combo; do
  cat "$combo".run*.1.fastq.gz > "$combo".1.fastq.gz
  cat "$combo".run*.2.fastq.gz > "$combo".2.fastq.gz
done < combo
```
### Remove the `run` files
```
rm *run*.fastq.gz
```

## Rename files
Cutadapt does not have built-in functionality to rename multiple files at once, but you can use an external tool such as `mmv`.  
Create a txt file (e.g. `patterns`), `fwdindex1-revindex1` etc. are the names of the adapter combinations and `sampleA` etc. are the corrisponding sample names:
```
fwdindex1-revindex1.[12].fastq.gz sampleA.#1.fastq.gz
fwdindex1-revindex2.[12].fastq.gz sampleB.#1.fastq.gz
fwdindex1-revindex3.[12].fastq.gz sampleC.#1.fastq.gz
```
### Rename all files at once with
```
mmv < patterns
```

---

# Importing in QIIME2
### Get absolute path of forward/reverse reads
```
ls -1 "$PWD/"*.1.fastq.gz > R1.txt
ls -1 "$PWD/"*.2.fastq.gz > R2.txt
```
### Get identifiers (sample names), assuming 5 is the number of letters
```
find *.1.fastq.gz -printf "%.5f\n" > ids.txt
```
Edit `ids.txt` accordingly if the sample names have different lengths.
### Make manifest file
```
paste ids.txt R1.txt R2.txt > man.tsv
echo sample-id$'\t'forward-absolute-filepath$'\t'reverse-absolute-filepath > manifest.tsv
cat man.tsv >> manifest.tsv
```
Make sure the sample names are aligned correctly with their absolute file paths.

### Import fastq
```
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path manifest.tsv \
  --output-path seqs.qza \
  --input-format PairedEndFastqManifestPhred33V2
```
### Summary of sequence quality
```
qiime demux summarize \
  --i-data seqs.qza \
  --o-visualization seqs.qzv
```

# Denoising with DADA2 in QIIME2
### Paired-end reads denoising
```
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs seqs.qza \
  --p-trim-left-f 30 \ # base pairs to trim on left side of forward reads (adapter + primer length)
  --p-trim-left-r 32 \ # base pairs to trim on left side of reverse reads (adapter + primer length)
  --p-trunc-len-f 150 \ # length at which to truncate forward reads
  --p-trunc-len-r 150 \ # length at which to truncate reverse reads
  --p-max-ee-f 2.0 \ # reads with higher expected errors will be discarded
  --p-max-ee-r 2.0 \
  --p-chimera-method 'consensus' \ # 'consensus' (in samples individually), 'pooled' (in pooled reads), 'none'
  --p-n-reads-learn 1000000 \ # number of reads used for training error model
  --p-n-threads 0 \ # number of threads to use for multithreaded processing, if 0 all available core will be used
  --o-table table.qza \ # ASV table
  --o-representative-sequences rep-seqs.qza \ # sequences of the ASVs
  --o-denoising-stats denoising-stats.qza # summary of reads lost in each step of the DADA2 algorithm
```

### Visualizations of the outputs of DADA2
```
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```
Prepare [metadata.tsv](https://docs.qiime2.org/2022.2/tutorials/metadata/)

<details><summary><strong>In case you have a NovaSeq dataset</strong></summary>
<p>

# Import in QIIME2 DADA2 outputs with NovaSeq correction (R)
### Import ASVs sequences
```
qiime tools import \
  --input-path rep-seqs.fasta \
  --output-path rep-seqs.qza \
  --type 'FeatureData[Sequence]'
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv
```
### Import ASVs table
```
biom convert -i ashtab.tsv -o ashtab.biom --table-type="OTU table" --to-hdf5
```
```
qiime tools import \
  --input-path ashtab.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format \
  --output-path table.qza
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file metadata.tsv
```

</p>
</details>

# Taxonomic Classification
### Train classifier
```
qiime tools import \
  --type 'FeatureData[Sequence]' \
  --input-path 12S-seqs-derep-uniq.fasta \
  --output-path 12S-seqs-derep-uniq.qza
qiime tools import \
  --type 'FeatureData[Taxonomy]' \
  --input-format HeaderlessTSVTaxonomyFormat \
  --input-path 12S-tax-derep-uniq.txt \
  --output-path 12S-tax-derep-uniq.qza
```
```
qiime feature-classifier extract-reads \
  --i-sequences 12S-seqs-derep-uniq.qza \
  --p-f-primer AAACTCGTGCCAGCCACC \
  --p-r-primer GGGTATCTAATCCCAGTTTG \
  --p-trunc-len 400 \
  --p-min-length 100 \
  --p-max-length 700 \
  --p-n-jobs 8 \
  --o-reads 12-derep-uniq_EXTRACTED.qza \
  --verbose
qiime feature-table tabulate-seqs \
  --i-data 12-derep-uniq_EXTRACTED.qza \
  --o-visualization 12-derep-uniq_EXTRACTED.qzv
```
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 12-derep-uniq_EXTRACTED.qza \
  --i-reference-taxonomy 12S-tax-derep-uniq.qza \
  --p-verbose \
  --o-classifier NBclassifier.qza
```

### First classification
```
qiime feature-classifier classify-sklearn \
  --i-classifier NBclassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence=0 \
  --p-n-jobs 8 \
  --o-classification NB-unweighted-classification.qza
```

### Generate weights for classifier
```
qiime clawback generate-class-weights \
  --i-reference-taxonomy 12S-tax-derep-uniq.qza \
  --i-reference-sequences 12S-derep-uniq_EXTRACTED.qza \
  --i-samples table.qza \
  --i-taxonomy-classification NB-unweighted-classification.qza \
  --o-class-weight NB-weights.qza \
  --verbose
```

### Re-train classifier with weights
```
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads 12S-derep-uniq_EXTRACTED.qza \
  --i-reference-taxonomy 12S-tax-derep-uniq.qza \
  --i-class-weight NB-weights.qza \
  --o-classifier NB-WeightedClassifier.qza
```

### Re-classify and get final taxonomy
```
qiime feature-classifier classify-sklearn \
  --i-classifier NB-WeightedClassifier.qza \
  --i-reads rep-seqs.qza \
  --p-confidence 0.95 \
  --p-n-jobs 8 \
  --o-classification taxonomy-weighted_95.qza
qiime metadata tabulate \
  --m-input-file taxonomy-weighted_95.qza \
  --o-visualization taxonomy-weighted_95.qzv
```

### Barplot
```
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy-weighted_95.qza \
  --m-metadata-file metadata.tsv \
  --o-visualization Barplot.qzv
```
