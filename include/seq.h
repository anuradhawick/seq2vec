#pragma once
#include <string>
#include <fstream>
#include <boost/iostreams/device/mapped_file.hpp>

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
    string path;
    ifstream file_handler;
    string seq_id = "";
    string seq_data = "";
    string line;
    int line_no = -1, seq_count = 0, seq_count_2 = 0;
    int mode;
    const char *fptr = NULL;   // pointer to file start
    const char *lptr = NULL;   // pointer to file end
    const char *tsptr = NULL;  // temporary pointer start
    const char *teptr = NULL;  // temporary pointer end
    const char *teptr2 = NULL; // temporary pointer end
    size_t fsize;
    iostreams::mapped_file mmap;
    bool _same_headers;

public:
    SeqReader(string path, bool same_headers = true)
    {
        path = path;
        file_handler.open(path, ios::in);
        mmap.open(path, iostreams::mapped_file::readonly);
        teptr2 = teptr = tsptr = fptr = mmap.const_data();
        lptr = fptr + mmap.size();
        fsize = mmap.size();
        assert(fsize > 0);

        if (static_cast<const char>(fptr[0]) == '>')
        {
            mode = 1;
        }
        else
        {
            mode = 2;
        }
        _same_headers = same_headers;
    }

    ~SeqReader()
    {
        file_handler.close();
        mmap.close();
    }

    size_t get_seq_count()
    {
        if (mode == 2)
        {
            while (teptr2 && teptr2 != lptr)
            {
                if ((teptr2 = static_cast<const char *>(memchr(teptr2, '\n', lptr - teptr2))))
                {
                    seq_count_2++;
                    teptr2++;
                }
            }
            return (seq_count_2 + 1) / 4;
        }
        else
        {
            while (teptr2 && teptr2 != lptr)
            {
                if ((teptr2 = static_cast<const char *>(memchr(teptr2, '>', lptr - teptr2))))
                {
                    seq_count_2++;
                    teptr2++;
                }
            }
            return seq_count_2;
        }
    }

    bool get_seq(Seq &seq)
    {

        if (mode == 1) // reading fasta
        {
            if (tsptr == lptr)
            {
                return false;
            }

            tsptr = static_cast<const char *>(memchr(tsptr, '>', lptr - tsptr));
            teptr = static_cast<const char *>(memchr(tsptr, '\n', lptr - tsptr));

            ++seq_count;

            if (_same_headers)
                seq.seq_header = string(tsptr, teptr);
            else
                seq.seq_header = "seq_" + to_string(seq_count);

            seq.seq_id = seq_count;
            teptr++;
            tsptr = teptr;

            seq.seq_string = "";

            // gather a sequence till next >

            while (tsptr != lptr)
            {
                // read a complete line
                teptr = static_cast<const char *>(memchr(tsptr, '\n', lptr - tsptr));

                // check if it is last line (with or w/o newline at end)
                seq.seq_string += string(tsptr, teptr != NULL ? teptr : lptr);

                if (teptr != NULL)
                {
                    // if not last line increment pointers
                    teptr++;
                    tsptr = teptr;
                }
                else
                {
                    // if last line set start to end
                    tsptr = lptr;
                    return true;
                }
                // if next line start with > sequence over
                if (static_cast<const char>(tsptr[0]) == '>')
                {
                    return true;
                }
            }
            // return true; handle end of file in next call to function
            return true;
        }
        else
        {
            if (tsptr == lptr)
            {
                return false;
            }

            teptr = static_cast<const char *>(memchr(tsptr, '\n', lptr - tsptr));

            ++seq_count;
            if (_same_headers)
                seq.seq_header = string(tsptr, teptr);
            else
                seq.seq_header = "seq_" + to_string(seq_count);
            
            seq.seq_id = seq_count;

            teptr++;
            tsptr = teptr;
            teptr = static_cast<const char *>(memchr(tsptr, '\n', lptr - tsptr));

            seq.seq_string = string(tsptr, teptr);

            teptr++;
            tsptr = teptr;
            teptr = static_cast<const char *>(memchr(tsptr, '\n', lptr - tsptr));

            teptr++;
            tsptr = teptr;
            teptr = static_cast<const char *>(memchr(tsptr, '\n', lptr - tsptr));

            tsptr = teptr;

            if (teptr != NULL)
            {
                teptr++;
                tsptr = teptr;
            }
            else
            {
                tsptr = lptr;
            }

            return true;
        }
    }
};
