#ifndef __INDEX_OPTIONS_HPP
#define __INDEX_OPTIONS_HPP

#include <string>

struct IndexOptions
{
    int kmer_size { 15 };
    int kmer_window { 3 };
    double kmer_rep_frac { 0.0003 };
    int num_threads { 8 };
    const char* reference_path { nullptr };
    std::string index_path;

    const char* _kmer_size { "-k" };
    const char* _kmer_window { "-w" };
    const char* _kmer_rep_frac { "-f" };
    const char* _num_threads { "-t" };
    const char* _index_path { "-o" };

    void simple_usage(int argc, char* argv[]);

    void full_usage(int argc, char* argv[]);

    bool parse(int argc, char* argv[]);

    void precheck_args(int argc, char* argv[]);

    void dump(FILE* out);
};

std::string make_reference_index_path(const char* reference_path, const char* user_specified_index_path);

#endif // __INDEX_OPTIONS_HPP
