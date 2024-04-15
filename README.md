# Seq2Vec - DNA sequence vectorization

This tool is intended to be used for data generation in Bioinformatics Machine Learning related tasks. You can use Seq2Vec to convert FASTA or FASTQ data sets into k-mer frequency vectors. We use memory mapped files to write faster and a multi-worker pipeline to vectorise the sequences.

### Citation for algorithms

```bibtex
@article{10.1093/bioinformatics/btaa441,
    author = {Wickramarachchi, Anuradha and Mallawaarachchi, Vijini and Rajan, Vaibhav and Lin, Yu},
    title = "{MetaBCC-LR: metagenomics binning by coverage and composition for long reads}",
    journal = {Bioinformatics},
    volume = {36},
    number = {Supplement_1},
    pages = {i3-i11},
    year = {2020},
    month = {07},
    abstract = "{Metagenomics studies have provided key insights into the composition and structure of microbial communities found in different environments. Among the techniques used to analyse metagenomic data, binning is considered a crucial step to characterize the different species of micro-organisms present. The use of short-read data in most binning tools poses several limitations, such as insufficient species-specific signal, and the emergence of long-read sequencing technologies offers us opportunities to surmount them. However, most current metagenomic binning tools have been developed for short reads. The few tools that can process long reads either do not scale with increasing input size or require a database with reference genomes that are often unknown. In this article, we present MetaBCC-LR, a scalable reference-free binning method which clusters long reads directly based on their k-mer coverage histograms and oligonucleotide composition.We evaluate MetaBCC-LR on multiple simulated and real metagenomic long-read datasets with varying coverages and error rates. Our experiments demonstrate that MetaBCC-LR substantially outperforms state-of-the-art reference-free binning tools, achieving ∼13\\% improvement in F1-score and ∼30\\% improvement in ARI compared to the best previous tools. Moreover, we show that using MetaBCC-LR before long-read assembly helps to enhance the assembly quality while significantly reducing the assembly cost in terms of time and memory usage. The efficiency and accuracy of MetaBCC-LR pave the way for more effective long-read-based metagenomics analyses to support a wide range of applications.The source code is freely available at: https://github.com/anuradhawick/MetaBCC-LR.Supplementary data are available at Bioinformatics online.}",
    issn = {1367-4803},
    doi = {10.1093/bioinformatics/btaa441},
    url = {https://doi.org/10.1093/bioinformatics/btaa441},
    eprint = {https://academic.oup.com/bioinformatics/article-pdf/36/Supplement\_1/i3/33488763/btaa441.pdf},
}
```

### Citation for this software

[![DOI](https://zenodo.org/badge/362989776.svg)](https://zenodo.org/badge/latestdoi/362989776)

```bibtex
@software{anuradha_wickramarachchi_2021_5515743,
  author       = {Anuradha Wickramarachchi},
  title        = {anuradhawick/seq2vec: release v1.0},
  month        = sep,
  year         = 2021,
  publisher    = {Zenodo},
  version      = {v1.0},
  doi          = {10.5281/zenodo.5515743},
  url          = {https://doi.org/10.5281/zenodo.5515743}
}
```


## Downloading and Compiling

First download the repository. The script named `build.sh` has all the required steps automated for easy compilation.

```
git clone https://github.com/anuradhawick/seq2vec.git
cd seq2vec
./build.sh
```

## Usage
Binary will be available at build/seq2vec. Help is available with `-h` command;

```
Seq2Vec fast sequence vectorization:
  -h [ --help ]              show help message
  -f [ --file ] arg          input file path
  -o [ --output ] arg        output vectors path
  -x [ --preset ] arg (=csv) output type, should be one of csv, tsv, or json
  -k [ --k-size ] arg (=3)   set k-mer size
  -t [ --threads ] arg (=8)  set thread count
```

## Output

A text file with the output will be generated at the output provided as the `-o` argument.

## Notes

* The default k-value is 3 and usually keep it under 8.
<!-- * The generated output directory will have several `*.txt` files containing the normalized vectors. Each line starts with sequence id (index starts at 1). You can process this output as you like. We provide the helper script `toH5.py` to sort-concatenate these vectors and to create an `H5` files (for ML tasks). Usage is as follows;

```
usage: toH5.py [-h] --seq2vec-outdir SEQ2VEC_OUTDIR --destination-file
               DESTINATION_FILE

This script of Seq2Vec helps you to convert the raw output to H5. Quite
helpful in machine learning work.

optional arguments:
  -h, --help            show this help message and exit
  --seq2vec-outdir SEQ2VEC_OUTDIR, -s2v SEQ2VEC_OUTDIR
                        Output directory of seq2vec containing all the *.txt
                        files.
  --destination-file DESTINATION_FILE, -h5 DESTINATION_FILE
                        Name of the destination h5 file.
```

You can find the vectors inside the `h5` file `vectors` dataset.

* You can also use the `gathered-sorted.txt` inside the **seq2vec** output folder generated by `toH5.py`. Note that each line starts with sequence id (index starts at 1). Dont forget to drop that column (use use it as pandas index).

* In linux use `cut -d' ' -f2- <SEQ2VEC OUTDIR>/gathered-sorted.txt > vectors.txt` to obtain vectors without seq ids. This can later be loaded to numpy as `np.loadtxt("vectors.txt")`. -->


***Have a good one!* 😃**
