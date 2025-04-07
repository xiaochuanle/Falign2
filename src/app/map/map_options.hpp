#ifndef __MAP_OPTIONS_HPP
#define __MAP_OPTIONS_HPP

#include <string>
#include <vector>

enum EOutputFmt {
    eOutputFmt_SAM = 0,
    eOutputFmt_FragSAM,
    eOutputFmt_BAM,
    eOutputFmt_FragBAM,
    eOutputFmt_PAF,
    eOutputFmt_Invalid
} ;

const char* output_format_name(EOutputFmt fmt);

EOutputFmt name_to_output_format(const char* name);

struct HbnMapOptions {

    double          ddf { 0.25 };
    int             kmer_dist { 400 };
    int             gap { 30 };
    int             chain_score { 2 };

    const char*     _ddf { "-f" };
    const char*     _kmer_dist { "-d" };
    const char*     _gap { "-g" };
    const char*     _chain_score { "-c" };

    /// enzyme
    const char*     enzyme { nullptr };
    size_t          ei_query_bases { 2000000000 };
    int             ei_num_hits { 5 };
    int             ei_frag_size { 200 };
    int             ei_end_dist { 200 };
    int             ei_flanking_bases { 30 };

    const char*     _enzyme { "-e" };
    const char*     _ei_query_bases { "-q" };
    const char*     _ei_num_hits { "-m" };
    const char*     _ei_frag_size { "-s" };
    const char*     _ei_end_dist { "-l" };
    const char*     _ei_flanking_bases { "-b" };

    /// Miscellaneous options

    int             num_threads { 1 };
    EOutputFmt      outfmt { eOutputFmt_PAF };

    const char*     reference { nullptr };
    const char*     output { "-" };

    const char*     _num_threads { "-t" };
    const char*     _outfmt { "-u" };
    const char*     _output { "-o" };

    /////////////////////////

    void precheck_args(int argc, char* argv[]);

    void simple_usage(int argc, char* argv[]);

    void full_usage(int argc, char* argv[]);

    bool parse(int argc, char* argv[], std::vector<std::string>& query_list);
};

void
s_dump_enzyme_names_and_seqs(FILE* out);

#endif // __MAP_OPTIONS_HPP
