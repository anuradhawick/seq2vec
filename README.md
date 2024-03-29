# Seq2Vec - DNA sequence vectorization

This tool is intended to be used for data generation in Bioinformatics Machine Learning related tasks. You can use Seq2Vec to convert FASTA or FASTQ data sets into k-mer frequency vectors. We use memory mapped files to write faster and a multi-worker pipeline to vectorise the sequences.

## Citing

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
