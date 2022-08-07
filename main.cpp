#include <iostream>
#include <thread>
#include <valarray>
#include <vector>
#include <iomanip>

#include <boost/program_options.hpp>

#include "./include/mode_mmap.h"

using namespace std;

int main(int ac, char **av)
{
    int ksize, threads;
    string input, output, type;

    po::options_description desc("Seq2Vec fast sequence vectorization");

    desc.add_options()("help,h", "show help message");
    desc.add_options()("file,f", po::value<string>(&input)->required(), "input file path");
    desc.add_options()("output,o", po::value<string>(&output)->required(), "output vectors path");
    desc.add_options()("preset,x", po::value<string>(&type)->default_value("csv"), "output type, should be one of csv, tsv, or json");
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

    if (type == "csv")
    {
        cout << "Starting Seq2Vec sequence vectorization: TSV output" << endl;
        mmapkmers::run(input, output, ksize, threads, ',');
    } 
    else if (type == "tsv")
    {
        cout << "Starting Seq2Vec sequence vectorization: TSV output" << endl;
        mmapkmers::run(input, output, ksize, threads, '\t');
    }

    return 0;
}
