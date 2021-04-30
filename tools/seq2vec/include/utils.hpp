#pragma once
#include <gatb/gatb_core.hpp>
#include <iostream>
#include <string>

namespace s2v
{
    namespace utils
    {
        struct ProgressFunctor : public IteratorListener
        {
            u_int64_t done = 0;

            void inc(u_int64_t current)
            {
                done += current;
                cout << "Processed reads... " << done << "      \r" << flush;
            }

            void finish()
            {
                cout << "Processed reads... " << done << "      \n"
                     << flush;
            }
        };

        inline void clean_output_dir(string output_dir)
        {
            IFileSystem &fs = System::file();

            if (fs.doesExistDirectory(output_dir))
            {
                if (fs.doesExistDirectory(output_dir + "/tmp"))
                {
                    for (string s : fs.listdir(output_dir + "/tmp"))
                    {
                        if (strcmp(s.c_str(), ".") == 0 || strcmp(s.c_str(), "..") == 0)
                        {
                            continue;
                        }
                        fs.remove(output_dir + "/tmp/" + s);
                    }
                }
                for (string s : fs.listdir(output_dir))
                {
                    if (strcmp(s.c_str(), ".") == 0 || strcmp(s.c_str(), "..") == 0)
                    {
                        continue;
                    }
                    fs.remove(output_dir + "/" + s);
                }
                fs.rmdir(output_dir);
            }

            fs.mkdir(output_dir, 4095);
            fs.mkdir(output_dir + "/tmp", 4095);
        }
    };

    namespace kmerutils
    {
        inline pair<vector<u_int64_t>, u_int64_t> compute_kmer_inds(Kmer<>::ModelCanonical &kmermodel)
        {
            u_int64_t ind = 0;
            u_int64_t kmer_counts_length = 0;
            LargeInt<1> kmer_fc;
            LargeInt<1> kmer_rc;
            map<u_int64_t, u_int64_t> kmer_inds;
            vector<u_int64_t> kmer_inds_index;

            for (u_int64_t kmer = 0; kmer <= kmermodel.getKmerMax().getVal(); kmer++)
            {

                kmer_fc.setVal(kmer);
                kmer_rc = kmermodel.reverse(kmer_fc);

                if (kmer_inds.find(kmer_rc.getVal()) != kmer_inds.end())
                {
                    kmer_inds[kmer] = kmer_inds[kmer_rc.getVal()];
                }
                else
                {
                    kmer_inds[kmer] = ind;
                    kmer_counts_length++;
                    ind += 1;
                }
            }

            kmer_inds_index.resize(kmermodel.getKmerMax().getVal() + 1, 0);

            for (auto its = kmer_inds.begin(); its != kmer_inds.end(); its++)
            {
                kmer_inds_index[its->first] = its->second;
            }

            return make_pair(kmer_inds_index, kmer_counts_length);
        }
    };
};
