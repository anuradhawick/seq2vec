#include <iostream>
#include <fstream>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <iostream>
#include <zlib.h>
#include <omp.h>
#include <thread>
#include <valarray>
#include <vector>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/asio.hpp>

#include "./include/seq.h"
#include "./include/kmer.h"

using namespace std;

namespace po = boost::program_options;
namespace bio = boost::iostreams;
namespace basio = boost::asio;

mutex mux;
queue<string> reads_queue;
condition_variable condition;
volatile bool terminate_threads;

class ProgressDisplay
{
private:
    int total = 0;
    int progress = 0;
public:
    ProgressDisplay(int total)
    {
        this->total = total;
    }

    void operator++(int)
    {
        this->operator++();
    }

    void operator++()
    {
        this->progress++;
        this->print();
    }

    void end()
    {
        this->progress=this->total;
        cout << "Completed " << fixed << setprecision(2) << 100.00 << "%       " << endl << flush;
    }

    void print()
    {
        float percentage = 100.0 * static_cast<float>(this->progress)/static_cast<float>(this->total);
        cout << "Completed " << fixed << setprecision(2) << percentage << "%             \r" << flush;
    }
};

void off_load_process(string &output, KmerCounter &kc, int &threads)
{
    string seq;
    vector<string> batch;
    ofstream outfs(output, ios::out);
    size_t batch_size;

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
            batch_size = batch.size();
            vector<valarray<double>> results(batch_size);
            ostringstream outss;
            outss.precision(6);
            outss << fixed;

#pragma omp parallel for num_threads(threads) schedule(dynamic, 1)
            for (size_t i = 0; i < batch_size; i++)
            {
                results[i] = kc.count_kmers(batch[i]);
            }

            for (auto dvec : results)
            {
                for (size_t j = 0; j < dvec.size(); j++)
                {
                    outss << dvec[j];
                    if (j < dvec.size() - 1)
                    {
                        outss << ' ';
                    }
                }
                outss << '\n';
            }

            outfs << outss.str();
            outss.clear();

            results.clear();
            batch.clear();
        }

        {
            unique_lock<mutex> lock(mux);
            if (terminate_threads && reads_queue.size() == 0)
            {
                break;
            }
        }
    }

    outfs.close();
}

void io_thread(SeqReader &reader)
{
    Seq seq;
    int count = 0;

    while (reader.get_seq(seq))
    {
        {
            unique_lock<mutex> lock(mux);
            condition.wait(lock, [] { return reads_queue.size() < 10000; });
            reads_queue.push(seq.seq_string);
        }
        count++;

        cout << "Loaded Reads " << count << "     \r" << flush;
    }

    cout << endl;

    terminate_threads = true;
}

void run(string &input, string &output, int &ksize, int &threads)
{
    SeqReader reader(input);
    KmerCounter kc(ksize);

    cout << "Counting sequences" << endl;
    size_t total_reads = reader.get_seq_count();
    cout << total_reads <<  " sequences found" << endl;
    size_t per_line_size = kc.kmer_counts_length * (8 + 1);
    size_t estimated_file_size = total_reads * per_line_size; // space after each value + newline (9 ASCII chars per value)
    
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
        asio::post(pool, [&reader_mux, &estimated_file_size, &reader, &pd, &mmout, &kc, &per_line_size]() {
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
                            outss << ' ';
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

int main(int ac, char **av)
{
    int ksize, threads;
    string input, output;

    po::options_description desc("Seq2Vec fast sequence vectorization");

    desc.add_options()("help,h", "show help message");
    desc.add_options()("file,f", po::value<string>(&input)->required(), "input file path");
    desc.add_options()("output,o", po::value<string>(&output)->required(), "output vectors path");
    desc.add_options()("k-size,k", po::value<int>(&ksize)->default_value(3), "set k-mer size");
    desc.add_options()("threads,t", po::value<int>(&threads)->default_value(8), "set thread count");

    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);

    if (vm.count("help") || ac == 1)
    {
        cout << desc << "\n";
        return 1;
    }

    po::notify(vm);

    cout << "Starting Seq2Vec sequence vectorization" << endl;
    run(input, output, ksize, threads);

    return 0;
}
