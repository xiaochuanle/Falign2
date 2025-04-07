#include "../../corelib/arg_parse.hpp"
#include "../../corelib/fasta.hpp"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <limits>

using namespace std;

static int g_extracted_seqs = numeric_limits<int>::max();
static size_t g_extracted_bases = numeric_limits<size_t>::max();

static int
s_parse_arg(int argc, char* argv[])
{
    int i = 2;
    while (i < argc) {
        bool r = (argv[i][0] != '-') || (strlen(argv[i]) == 1);
        if (r) break;
        if (parse_data_size_arg_value(argc, argv, i, "-s", g_extracted_bases)) continue;
        if (parse_int_arg_value(argc, argv, i, "-n", g_extracted_seqs)) continue;
        fprintf(stderr, "ERROR: Unrecognised option '%s\n", argv[i]);
        return -1;
    }
    return i;
}

static void
s_dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "  %s %s [OPTIONS] fq-1 [... fq-n]\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS:\n");
    fprintf(stderr, "  -n <Integer>\n");
    fprintf(stderr, "    Number of sequences extracted\n");
    fprintf(stderr, "  -s <DataSize>\n");
    fprintf(stderr, "    Number of bases extracted\n");
}

int extract_fastq_main(int argc, char* argv[])
{
    int fq_idx = s_parse_arg(argc, argv);
    if (fq_idx == -1) {
        s_dump_usage(argc, argv);
        return 1;
    }
    int n_seq = 0;
    size_t n_base = 0;
    for (; fq_idx < argc; ++fq_idx) {
        HBN_LOG("Extract sequences from %s", argv[fq_idx]);
        HbnFastaReader in(argv[fq_idx]);
        bool is_done = false;
        while (in.ReadOneSeq() != -1) {
            bool has_qual = !in.qual().empty();

            if (has_qual) fprintf(stdout, "@"); else fprintf(stdout, ">");
            fprintf(stdout, "%s", in.name().c_str());
            if (!in.comment().empty()) fprintf(stdout, " %s", in.comment().c_str());
            fprintf(stdout, "\n");

            fprintf(stdout, "%s\n", in.sequence().c_str());

            if (has_qual) fprintf(stdout, "+\n");

            if (has_qual) fprintf(stdout, "%s\n", in.qual().c_str());
        
            ++n_seq;
            n_base += in.sequence().size();
            if (n_seq >= g_extracted_seqs || n_base >= g_extracted_bases) {
                is_done = true;
                break;
            }
        }
        if (is_done) break;
    }

    string size = NStr::UInt8ToString_DataSize(n_base);
    HBN_LOG("Extract %d sequences (%s)", n_seq, size.c_str());

    return 0;
}