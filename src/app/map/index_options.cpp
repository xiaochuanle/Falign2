#include "index_options.hpp"

#include "../../corelib/arg_parse.hpp"
#include "program_info.hpp"

#include <string>

using namespace std;

std::string make_reference_index_path(const char* reference_path, const char* user_specified_index_path)
{
    if (user_specified_index_path) return string(user_specified_index_path);
    if (strcmp(reference_path, "-") == 0) return string("falign.fi");
    string path = reference_path;
    path += ".fi";
    return path;
}

void IndexOptions::simple_usage(int argc, char* argv[])
{
    FILE* out = stderr;
    fprintf(out, "\n");
    fprintf(out, "USAGE\n");
    fprintf(out, "  %s %s [OPTIONS] ref.fa\n", argv[0], argv[1]);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION\n");
    fprintf(out, "  Genomic database index construction toolkit for falign\n");

    fprintf(out, "\n");
    fprintf(out, "VERSION\n");
    const char* version = HBN_PACKAGE_VERSION;
    fprintf(out, "  %s\n", version);

    fprintf(out, "\n");
    fprintf(out, "Type option '-help' to see details of optional arguments\n");    
}

void IndexOptions::full_usage(int argc, char* argv[])
{
    string size;
    FILE* out = stderr;
    fprintf(out, "\n");
    fprintf(out, "USAGE\n");
    fprintf(out, "  %s %s [OPTIONS] reference\n", argv[0], argv[1]);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION\n");
    fprintf(out, "  Genomic database index construction toolkit for falign\n");

    fprintf(out, "\n");
    fprintf(out, "VERSION\n");
    const char* version = HBN_PACKAGE_VERSION;
    fprintf(out, "  %s\n", version);

    fprintf(out, "\n");
    fprintf(out, "OPTIONAL ARGUMENTS\n");

    fprintf(out, "  %s\n", "-h");
    fprintf(out, "    Print USAGE and DESCRIPTION; ignore all other parameters\n");
    fprintf(out, "  %s\n", "-help");
    fprintf(out, "    Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\n");

    fprintf(out, "  %s <Integer>\n", _kmer_size);
    fprintf(out, "    Kmer size (length of best perfect match)\n");
    fprintf(out, "    Default = '%d'\n", kmer_size);
    fprintf(out, "  %s <Integer>\n", _kmer_window);
    fprintf(out, "    Kmer sampling window size in reference sequences\n");
    fprintf(out, "    Default = '%d'\n", kmer_window);
    fprintf(out, "  %s <Real>\n", _kmer_rep_frac);
    fprintf(out, "    Filter out this fraction of most frequently occurring kmers on repeat reference regions\n");
    fprintf(out, "    Default = '%g'\n", kmer_rep_frac);
    fprintf(out, "  %s <Path>\n", _index_path);
    fprintf(out, "    Output genomic index to this path\n");

    fprintf(out, "  %s <Integer>\n", _num_threads);
    fprintf(out, "    Number of threads (CPUs) to use in the search\n");
    fprintf(out, "    Default = '%d'\n", num_threads);
}

void IndexOptions::precheck_args(int argc, char* argv[])
{
    for (int i = 2; i < argc; ++i) {
        bool r = (argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (strcmp(argv[i] + 1, "h") == 0) {
            simple_usage(argc, argv);
            exit(0);
        }
        if (strcmp(argv[i] + 1, "help") == 0) {
            full_usage(argc, argv);
            exit(0);
        }
    }
}

bool IndexOptions::parse(int argc, char* argv[])
{
    precheck_args(argc, argv);
    const char* c_index_path = nullptr;
    int i = 2;
    while (i < argc) {
        bool r = (argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (parse_int_arg_value(argc, argv, i, _kmer_size, kmer_size)) continue;
        if (parse_int_arg_value(argc, argv, i, _kmer_window, kmer_window)) continue;
        if (parse_real_arg_value(argc, argv, i, _kmer_rep_frac, kmer_rep_frac)) continue;
        if (parse_int_arg_value(argc, argv, i, _num_threads, num_threads)) continue;
        if (parse_string_arg_value(argc, argv, i, _index_path, c_index_path)) {
            index_path = c_index_path;
            continue;
        }

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    reference_path = argv[i];

    index_path = make_reference_index_path(reference_path, c_index_path);

    return true;
}

void IndexOptions::dump(FILE* out)
{
    fprintf(out, "\n");
    fprintf(out, "\n");
    fprintf(out, "====> Genomic index construction params:\n");
    fprintf(out, "kmer: %d\n", kmer_size);
    fprintf(out, "kmer-window: %d\n", kmer_window);
    fprintf(out, "repeat-kmer-frac: %g\n", kmer_rep_frac);
    fprintf(out, "CPU-threads: %d\n", num_threads);
    fprintf(out, "Reference: %s\n", reference_path);
    fprintf(out, "index: %s\n", index_path.c_str());
    fprintf(out, "\n\n");
}