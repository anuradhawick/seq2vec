#include <iostream>
#include <fstream>
#include <mutex>
#include <zlib.h>
#include <omp.h>
#include <thread>
#include <valarray>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/asio.hpp>

#include "./seq.h"
#include "./kmer.h"
#include "./progress.h"

using namespace std;

namespace po = boost::program_options;
namespace bio = boost::iostreams;
namespace basio = boost::asio;

namespace batchkmers 
{
    void run(string &input, string &output, int &ksize, int &threads, char sep)
    {
        SeqReader reader(input);
        KmerCounter kc(ksize);
        asio::thread_pool pool(threads);
        ProgressDisplay pd(0);
        mutex reader_mux;

        for (int _ = 0; _ < threads * 5; _++)
        {
            asio::post(pool, [&]() {
                bool has_read = true;
                Seq seq;

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
                        // TODO JSON implementation
                        // ostringstream outss;

                        // for (size_t j = 0; j < dvec.size(); j++)
                        // {
                        //     outss << dvec[j];
                        //     if (j < dvec.size() - 1)
                        //     {
                        //         outss << sep;
                        //     }
                        // }
                        // outss << '\n';
                    }
                    else
                    {
                        break;
                    }
                }
            });
        }
        pool.join();
        pd.end();
    }
}
