#pragma once
#include <string>
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

using namespace std;
using namespace boost;

class Seq
{
public:
    size_t seq_id;
    string seq_header;
    string seq_string;
};

class SeqReader
{
private:
    gzFile fp;
    kseq_t *ks;
    int ret;
    size_t seq_count = 0;
    size_t seq_id = 0;

public:
    SeqReader(string path)
    {
        fp = gzopen(path.c_str(), "r");
        ks = kseq_init(fp);
    }

    ~SeqReader()
    {
        kseq_destroy(ks);
        gzclose(fp);
    }

    size_t get_seq_count()
    {
        gzrewind(fp);
        kseq_rewind(ks);

        while((ret = kseq_read(ks)) >= 0)
        {
            seq_count++;
        }

        gzrewind(fp);
        kseq_rewind(ks);

        return seq_count;
    }

    bool get_seq(Seq &seq)
    {
        if ((ret = kseq_read(ks)) >= 0)
        {
            // using raw pointers can have unwated side effects
            // better to copy than debug
            seq.seq_string = string(ks->seq.s);
            seq.seq_header = string(ks->name.s);
            seq.seq_id = seq_id;
            seq_id++;
            
            return true;
        }
        return false;
    }
};