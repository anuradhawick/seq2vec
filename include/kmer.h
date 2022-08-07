#pragma once
#include <iostream>
#include <vector>
#include <valarray>
#include <map>
#include <cmath>

using namespace std;

class KmerCounter
{
private:
    u_int64_t kmer_size = 0;
    vector<u_int64_t> kmer_inds_index;
    // A=00, C=01, T=10, G=11
    u_int64_t reverse_complements[4] = { 0b10, 0b11, 0b00, 0b01 };

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
    u_int64_t kmer_counts_length = 0;

    KmerCounter(uint64_t kmer_size)
    {
        this->kmer_size = kmer_size;
        compute_kmer_inds();
    }

    valarray<double> count_kmers(string seq)
    {
        valarray<double> profile((double)0, kmer_counts_length);
        // val: forward and rval: reverse complements
        // len==k if the k-mer is correct with [acgtACGT]{k}
        u_int64_t val = 0, rval = 0, len = 0;
        char seqChar;

        for (size_t i = 0; i < seq.length(); i++)
        {
            seqChar = toupper(seq[i]);

            if (!(seqChar == 'A' || seqChar == 'C' || seqChar == 'G' || seqChar == 'T'))
            {
                len = 0;
                val = 0;
                rval = 0;
                continue;
            }
            
            val = (val << 2);
            val = val & ((u_int64_t)pow(4, kmer_size) - 1);
            val += (seqChar >> 1 & 3);
            rval += (reverse_complements[(seqChar >> 1 & 3)]<<(kmer_size*2));
            rval >>=2;
            rval = rval & ((u_int64_t)pow(4, kmer_size) - 1);
            len++;


            if (len == kmer_size)
            {
                // use val as the kmer for counting
                len--;
                profile[kmer_inds_index[val]]++;
            }
        }

        profile  = profile / max(1.0, profile.sum());

        return profile;
    }
};
