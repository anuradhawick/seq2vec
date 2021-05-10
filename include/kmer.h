#pragma once
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

using namespace std;

class KmerCounter
{
private:
    u_int64_t kmer_size = 0;
    u_int64_t kmer_counts_length = 0;
    vector<u_int64_t> kmer_inds_index;

    KmerCounter();

    u_int64_t rev_comp(u_int64_t x)
    {
        u_int64_t res = x;

        res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
        res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
        res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
        res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
        res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
        res = res ^ 0xAAAAAAAAAAAAAAAA;

        return (res >> (2 * (32 - kmer_size)));
    }

    void compute_kmer_inds()
    {
        u_int64_t ind = 0, kmer_fc, kmer_rc;
        map<u_int64_t, u_int64_t> kmer_inds;

        for (u_int64_t kmer = 0; kmer < (u_int64_t)pow(4, kmer_size); kmer++)
        {

            kmer_fc = kmer;
            kmer_rc = rev_comp(kmer_fc);

            if (kmer_inds.find(kmer_rc) != kmer_inds.end())
            {
                kmer_inds[kmer] = kmer_inds[kmer_rc];
            }
            else
            {
                kmer_inds[kmer] = ind;
                kmer_counts_length++;
                ind += 1;
            }
        }

        kmer_inds_index.resize((u_int64_t)pow(4, kmer_size), 0);

        for (auto its = kmer_inds.begin(); its != kmer_inds.end(); its++)
        {
            kmer_inds_index[its->first] = its->second;
        }
    }

public:
    KmerCounter(uint64_t kmer_size)
    {
        this->kmer_size = kmer_size;
        compute_kmer_inds();
    }

    vector<double> count_kmers(string seq)
    {
        vector<double> profile(kmer_counts_length, 0);
        double total = 0;
        long len = 0;
        u_int64_t val = 0;

        for (int i = 0; i < (int)seq.length(); i++)
        {
            if (!(seq[i] == 'A' || seq[i] == 'C' || seq[i] == 'G' || seq[i] == 'T'))
            {
                len = 0;
                val = 0;
                continue;
            }
            val = (val << 2);
            val = val & ((u_int64_t)pow(4, kmer_size) - 1);
            val += (seq[i] >> 1 & 3);
            len++;

            if (len == kmer_size)
            {
                // use val as the kmer for counting
                len--;
                profile[kmer_inds_index[val]]++;
                total++;
            }
        }

        for (size_t i = 0; i < profile.size(); i++)
        {
            profile[i] /= max(1.0, total);
        }

        return profile;
    }
};
