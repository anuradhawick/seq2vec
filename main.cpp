#include <iostream>
#include <fstream>
#include <mutex>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio.hpp>

#include "./include/seq.h"
#include "./include/kmer.h"

using namespace std;
namespace po = boost::program_options;
mutex mux;

void run(string input, string outdir, int ksize, int threads, bool use_mm)
{
    Reader *reader;
    if (use_mm)
    {
        reader = new SeqReaderMM(input);
    }
    else
    {
        reader = new SeqReader(input);
    }
    KmerCounter kc(ksize);
    asio::thread_pool pool(threads);

    for (size_t i = 0; i < threads; i++)
    {
        asio::post(pool, [&reader, &kc, &outdir]() {
            ofstream output(outdir + "/" + lexical_cast<string>(this_thread::get_id()) + ".txt", ios::out);
            bool found = true;
            string outstr = "";
            Seq seq;

            while (found)
            {
                {
                    unique_lock<mutex> lock(mux);
                    found = reader->get_seq(seq);
                }

                if (found)
                {
                    outstr = to_string(seq.seq_id);

                    for (double d : kc.count_kmers(seq.seq_string))
                    {
                        outstr += " " + to_string(d);
                    }
                    output << outstr << endl;
                    outstr = "";
                }
            }

            output.close();
        });
    }

    pool.join();
    delete reader;
}

int main(int ac, char **av)
{
    int ksize, threads;
    string input, outdir;
    bool use_mm = false;
    po::options_description desc("Seq2Vec fast sequence vectorization");
    desc.add_options()("help,h", "show help message");
    desc.add_options()("file,f", po::value<string>(&input)->required(), "input file path");
    desc.add_options()("outdir,o", po::value<string>(&outdir)->required(), "output directory");
    desc.add_options()("k-size,k", po::value<int>(&ksize)->default_value(3), "set k-mer size");
    desc.add_options()("threads,t", po::value<int>(&threads)->default_value(8), "set thread count");
    desc.add_options()("memory-mapped,m", "use memory mapped input (best for fastq)");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);

    if (vm.count("help") || ac == 1)
    {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("memory-mapped"))
    {
        use_mm = true;
    }

    po::notify(vm);

    run(input, outdir, ksize, threads, use_mm);

    return 0;
}
