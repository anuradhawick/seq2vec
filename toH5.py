import h5py
from glob import glob
from tqdm import tqdm
import numpy as np
import argparse

def main(tmp_path, tmp_sort_path, h5_filename):
    ds_size = 0
    d_shape = 0

    for file in glob(tmp_path):
        for line in open(file):
            if d_shape == 0:
                d_shape = len(line.strip().split()) - 1
            ds_size += 1

    file_heads = [open(file) for file in glob(tmp_path)]
    sorted_output = open(tmp_sort_path, "w+")

    line_heads = [file.readline() for file in  file_heads]

    for x in tqdm(range(ds_size), total=ds_size, desc="Sorting vectors"):
        line_idx = np.array([int(line.strip().split()[0]) if line else float('inf') for line in line_heads])
        min_idx = np.argmin(line_idx)

        sorted_output.write(line_heads[min_idx])
        line_heads = [line if n!=min_idx else file_heads[n].readline() for n, line in enumerate(line_heads)]


    [file.close() for file in file_heads]


    h5 = h5py.File(h5_filename, 'w')
    ds = h5.create_dataset("vectors", (ds_size, d_shape), dtype='f')

    off_set = 0
    vecs = []

    for line in tqdm(open(tmp_sort_path), total=ds_size, desc="Writing to h5"):
        d = list(map(float, line.strip().split()))
        idx = str(d[0])
        data = np.array(d[1:])
        vecs.append(data)

        # when large enough
        if len(vecs) * d_shape > 1073741824:
            vecs = np.array(vecs)
            ds[off_set:off_set+len(vecs)] = vecs
            off_set += len(vecs)

            vecs = []
            
    vecs = np.array(vecs)
    ds[off_set:off_set+len(vecs)] = vecs
    off_set += len(vecs)

    vecs = []

    h5.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""This script of Seq2Vec helps you to convert the raw output to H5. \
        Quite helpful in machine learning work.""")

    parser.add_argument('--seq2vec-outdir', '-s2v',
                            help="Output directory of seq2vec containing all the *.txt files.",
                            type=str,
                            required=True)

    parser.add_argument('--destination-file', '-h5',
                            help="Name of the destination h5 file.",
                            type=str,
                            required=True)

    args = parser.parse_args()

    tmp_path = args.seq2vec_outdir + "/*.txt"
    tmp_sort_path = args.seq2vec_outdir + "/gathered-sorted.txt"
    h5_filename = args.destination_file

    if ".h5" not in h5_filename:
        h5_filename += ".h5"

    main(tmp_path, tmp_sort_path, h5_filename)