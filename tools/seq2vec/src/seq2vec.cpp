#include <seq2vec.hpp>

using namespace std;

/********************************************************************************/

// We define some constant strings for names of command line parameters
// Use defaults

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
seq2vec::seq2vec() : Tool("seq2vec")
{
    // We add some custom arguments for command line interface
    getParser()->push_front(new OptionOneParam(STR_URI_FILE, "Path to reads", true));
    getParser()->push_front(new OptionOneParam(STR_URI_OUTPUT_DIR, "Output directory path", true));
    getParser()->push_front(new OptionOneParam(STR_KMER_SIZE, "Size of kmer for vectorization", false, "3", true));
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void seq2vec::execute()
{
    // We can do here anything we want.
    // For further information about the Tool class, please have a look
    // on the ToyTool snippet  (http://gatb-core.gforge.inria.fr/snippets_tools.html)

    // Program starts here-onwards
    string output_dir = getInput()->getStr(STR_URI_OUTPUT_DIR);
    u_int64_t nb_cores = getInput()->getInt(STR_NB_CORES);

    s2v::utils::clean_output_dir(output_dir);
    Kmer<>::ModelCanonical kmermodel(getInput()->getInt(STR_KMER_SIZE));

    cout << "Recoding all kmers      \r" << flush;
    pair<vector<u_int64_t>, u_int64_t> kmer_info = s2v::kmerutils::compute_kmer_inds(kmermodel);
    cout << "Recoding all kmers done   " << kmer_info.second << endl;

    vector<u_int64_t> kmer_inds_index = kmer_info.first;

    IBank *bank = Bank::open(getInput()->getStr(STR_URI_FILE));
    unordered_map<IThread::Id, ofstream> output_files;
    ISynchronizer *sync_map = System::thread().newSynchronizer();

    LOCAL(bank);
    LOCAL(sync_map);

    Dispatcher kmer_counter(nb_cores, 500);
    Iterator<Sequence> *bankit = bank->iterator();
    LOCAL(bankit);

    SubjectIterator<Sequence> itNotif(bankit, 1000);
    itNotif.addObserver(new s2v::utils::ProgressFunctor());

    cout << "Using Threads   " << nb_cores << endl;
    
    kmer_counter.iterate(itNotif, [&](Sequence &seq) {
        vector<double> profile(kmer_info.second, 0);
        Kmer<>::ModelCanonical::Iterator itkmer(kmermodel);
        itkmer.setData(seq.getData());
        double total_kmers = 0;
        IThread::Id tid = System::thread().getThreadSelf();
        size_t sid = seq.getIndex();
        string out_string = to_string(sid);

        for (itkmer.first(); !itkmer.isDone(); itkmer.next())
        {
            profile[kmer_inds_index[itkmer->value().getVal()]]++;
            total_kmers++;
        }

        for (size_t i = 0; i < profile.size(); i++)
        {
            profile[i] /= max(1.0, total_kmers);
            out_string += " " + to_string(profile[i]);
        }

        out_string += "\n";

        if (output_files.find(tid) == output_files.end())
        {
            sync_map->lock();
            output_files[tid] = ofstream(output_dir + "/tmp/" + to_string(tid) + ".txt", ios::out);
            sync_map->unlock();
        }

        output_files[tid] << out_string;
        out_string = "";
    });

    for (unordered_map<IThread::Id, ofstream>::iterator it = output_files.begin(); it != output_files.end(); ++it)
    {
        it->second.close();
    }

    
}
