#include "index.hpp"

#include "index_options.hpp"
#include "program_info.hpp"

#include <sys/utsname.h>
#include <string>
#include <thread>

using namespace std;

extern "C"
size_t getMemorySizeBytes();

void index_program_info(FILE* out)
{
    struct utsname _os_info_buf;
    struct utsname* os_info = nullptr;
    if (uname(&_os_info_buf) == 0) os_info = &_os_info_buf;

    size_t sys_mem_bytes = getMemorySizeBytes();
    string sys_mem = NStr::UInt8ToString_DataSize(sys_mem_bytes);

    int cpu_threads = thread::hardware_concurrency();

    fprintf(out, "\n");
    fprintf(out, "PROGRAM:\n");
    fprintf(out, "  Name:           %s\n", "falign2-index");
    fprintf(out, "  Version:        %s\n", HBN_PACKAGE_VERSION);
    fprintf(out, "  Description:    Genomic database index construction toolkit for falign2\n");
    fprintf(out, "  Contact:        chenying2016@gmail.com\n");

    fprintf(out, "\n");
    fprintf(out, "SYSTEM:\n");
    if (os_info) {
    fprintf(out, "  Computer:       %s\n", os_info->nodename);
    fprintf(out, "  Name:           %s\n", os_info->sysname);
    fprintf(out, "  Release:        %s\n", os_info->release);
    fprintf(out, "  Version:        %s\n", os_info->version);
    fprintf(out, "  Machine:        %s\n", os_info->machine);
    }
    fprintf(out, "  CPU threads:    %d\n", cpu_threads);
    fprintf(out, "  RAM:            %s\n", sys_mem.c_str());
    fprintf(out, "\n");
}

int build_index(int argc, char* argv[])
{
    if (argc < 2 || strcmp(argv[1], "index")) return false;
    IndexOptions opts;
    if (!opts.parse(argc, argv)) {
        opts.simple_usage(argc, argv);
        exit (1);
    }
    HbnProgramInfo program("falign-index", index_program_info);
    opts.dump(stderr);
    hbn_dfopen(out, opts.index_path.c_str(), "wb");

    HbnUnpackedDatabase updb(opts.reference_path);
    size_t num_kmer = 0;
    size_t num_hash = 0;
    {
        vector<pair<u64, u64>> hash_to_offset_list;
        build_lookup_table_mt(&updb, opts.kmer_size, opts.kmer_window,
            opts.kmer_rep_frac, opts.num_threads, hash_to_offset_list, num_kmer, out);
        pair<u64, u64>* a = hash_to_offset_list.data();
        num_hash = hash_to_offset_list.size();
        hbn_fwrite(a, sizeof(pair<u64, u64>), num_hash, out);
    }
    updb.dump(out);
    hbn_fwrite(&opts.kmer_size, sizeof(int), 1, out);
    hbn_fwrite(&num_kmer, sizeof(size_t), 1, out);
    hbn_fwrite(&num_hash, sizeof(size_t), 1, out);
    hbn_fclose(out);

    return 0;
}

void HbnIndex::load(const char* path) 
{
    std::string index_path = make_reference_index_path(path, nullptr);
    hbn_dfopen(in, index_path.c_str(), "rb");
    long int EndSize = sizeof(int) + sizeof(size_t) + sizeof(size_t);
    fseek(in, -EndSize, SEEK_END);
    hbn_fread(&M_kmer_size, sizeof(int), 1, in);
    hbn_fread(&M_kmer_list_size, sizeof(size_t), 1, in);
    size_t num_hash = 0;
    hbn_fread(&num_hash, sizeof(size_t), 1, in);
    fseek(in, 0, SEEK_SET);
    M_kmer_list = new KmerInfo[M_kmer_list_size];
    hbn_fread(M_kmer_list, sizeof(KmerInfo), M_kmer_list_size, in);
    M_hash_table = new HbnHashTable(M_kmer_size);
    pair<u64, u64> ip;
    for (size_t i = 0; i < num_hash; ++i) {
        hbn_fread(&ip, sizeof(pair<u64, u64>), 1, in);
        M_hash_table->add_one_hash_and_value(ip.first, ip.second);
    }
    M_updb = new HbnUnpackedDatabase();
    M_updb->load(in);

    hbn_fclose(in);
}
