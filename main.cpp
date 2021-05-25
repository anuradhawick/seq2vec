#include <iostream>
#include <fstream>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <iostream>
#include <zlib.h>
#include <omp.h>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/asio/thread_pool.hpp>
#include <boost/asio.hpp>

#include "./include/seq.h"
#include "./include/kmer.h"

using namespace std;

namespace po = boost::program_options;
mutex mux;
queue<string> reads_queue;
condition_variable condition;
volatile bool terminate_threads;

void off_load_process(string &output, KmerCounter &kc, int &threads)
{
    string seq;
    vector<string> batch;
    ofstream outfs(output, ios::out);

    while (true)
    {
        {
            unique_lock<mutex> lock(mux);

            while (reads_queue.size() > 0)
            {
                seq = reads_queue.front();
                batch.push_back(seq);
                reads_queue.pop();

                if (batch.size() == 10000)
                {
                    break;
                }
            }
        }

        condition.notify_all();

        if (batch.size() > 0)
        {
            vector<vector<double>> results(batch.size());
            ostringstream outss;

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
            for (size_t i = 0; i < batch.size(); i++)
            {
                results[i] = kc.count_kmers(seq);
            }

            for (auto dvec : results)
            {
                copy(dvec.begin(), dvec.end() - 1, ostream_iterator<int>(outss, "\t"));
                outss << dvec.back();
                outss << endl;
            }
            
            outfs << outss.str();
            outss.clear();

            batch.clear();
            results.clear();
        }

        {
            unique_lock<mutex> lock(mux);
            if (terminate_threads && reads_queue.size() == 0)
            {
                break;
            }
        }
    }
}

void io_thread(Reader &reader)
{
    Seq seq;
    int count = 0;

    while (reader.get_seq(seq))
    {
        {
            unique_lock<mutex> lock(mux);
            condition.wait(lock, [] { return reads_queue.size() < 50000; });
            reads_queue.push(seq.seq_string);
        }
        count++;

        cout << "Loaded Reads " << count << "       \r" << flush;
    }

    cout << endl;

    terminate_threads = true;
}

void run(string &input, string &output, int &ksize, int &threads, bool use_mm)
{
    Reader *reader;
    if (use_mm)
    {
        reader = new SeqReaderMM(input);
    }
    else
    {
        reader = new SeqReaderKS(input);
    }
    KmerCounter kc(ksize);
    asio::thread_pool pool(threads);

    thread iot(io_thread, std::ref(*reader));
    thread kct(off_load_process, std::ref(output), std::ref(kc), std::ref(threads));

    iot.join();
    kct.join();

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
    desc.add_options()("output,o", po::value<string>(&outdir)->required(), "output vectors path");
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
