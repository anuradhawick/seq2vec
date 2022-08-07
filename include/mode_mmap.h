#include <iostream>
#include <fstream>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <zlib.h>
#include <omp.h>
#include <thread>
#include <valarray>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/asio.hpp>

#include "./seq.h"
#include "./kmer.h"
#include "./progress.h"

using namespace std;

namespace po = boost::program_options;
namespace bio = boost::iostreams;
namespace basio = boost::asio;

namespace mmapkmers 
{
    void run(string &input, string &output, int &ksize, int &threads, char sep)
    {
        SeqReader reader(input);
        KmerCounter kc(ksize);

        cout << "Counting sequences" << endl;
        size_t total_reads = reader.get_seq_count();
        cout << total_reads <<  " sequences found" << endl;
        size_t per_line_size = kc.kmer_counts_length * (8 + 1);
        size_t estimated_file_size = total_reads * per_line_size; // sep + newline (9 ASCII chars per value)
        
        bio::mapped_file_params params;
        params.path = output;
        params.new_file_size = estimated_file_size;
        params.flags = bio::mapped_file::mapmode::readwrite;
        bio::mapped_file_sink mmout(params);

        asio::thread_pool pool(threads);
        ProgressDisplay pd(total_reads);
        mutex reader_mux;

        for (int _ = 0; _ < threads * 5; _++)
        {
            asio::post(pool, [&]() {
                bool has_read = true;
                Seq seq;
                auto sptr = mmout.begin();

                while (true)
                {
                    {
                        unique_lock<mutex> lock(reader_mux);
                        has_read = reader.get_seq(seq);
                        pd++;
                    }

                    if (has_read)
                    {
                        // process the read
                        valarray<double> dvec = kc.count_kmers(seq.seq_string);
                        ostringstream outss;
                        outss.precision(6);
                        outss << fixed;

                        for (size_t j = 0; j < dvec.size(); j++)
                        {
                            outss << dvec[j];
                            if (j < dvec.size() - 1)
                            {
                                outss << sep;
                            }
                        }
                        outss << '\n';

                        memcpy(sptr + seq.seq_id * per_line_size, outss.str().c_str(), outss.str().size());
                    }
                    else
                    {
                        break;
                    }
                }
            });
        }
        pool.join();
        mmout.close();
        pd.end();
    }
}
