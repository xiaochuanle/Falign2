#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <iostream>

#include "../../corelib/fasta.hpp"
#include "../../corelib/seq_name2id_map.hpp"
#include "../../corelib/split_string_by_char.hpp"
#include "../../htslib/sam.h"

using namespace std;

static void
s_parse_one_gt_frag_map(const char* s, const int sl, int& sfrom, int& sto)
{
    int i = 0, j = 0;
    while (j < sl && s[j] != '-') ++j;
    hbn_assert(j < sl);
    hbn_assert(s[j] == '-');

    i = j + 1;
    j = i;
    sfrom = 0;
    hbn_assert(i < sl);
    while (j < sl) {
        if (s[j] == '-') break;
        sfrom = sfrom * 10 + s[j] - '0';
        ++j;
    }

    hbn_assert(j < sl);
    hbn_assert(s[j] == '-');
    i = j + 1;
    j = i;
    sto = 0;
    while (j < sl) {
        if (s[j] == '-') break;
        sto = sto * 10 + s[j] - '0';
        ++j;
    }
}

static void
s_load_sim_read_info(const char* fq_path, 
    size_t& num_reads,
    size_t& num_frags,
    size_t& num_bases,
    size_t& num_contacts)
{
    ifstream in(fq_path);
    string line;
    int read_id = 0;
    while (getline(in, line)) {
        if (line[0] != '@') continue;
        const size_t N = line.size();
        size_t i = 1, j = i;

        while (j < N) {
            if (isspace(line[j])) break;
            ++j;
        }
        hbn_assert(j < N);
        hbn_assert(isspace(line[j]));

        i = j + 1;
        j = i;
        hbn_assert(i < N);
        while (j < N) {
            if (isspace(line[j])) break;
            ++j;
        }
        hbn_assert(j < N);
        hbn_assert(isspace(line[j]));

        i = j + 1;
        hbn_assert(i < N);
        int frag = 0;
        while (i < N) {
            j = i + 1;
            while (j < N && line[j] != ';') ++j;
            int sfrom, sto;
            s_parse_one_gt_frag_map(line.c_str() + i, j - i, sfrom, sto);
            ++num_frags;
            num_bases += (sto - sfrom);
            ++frag;
            i = j + 1;
        }
        num_contacts += frag * (frag-1)/2;
        ++num_reads;
    }
}

int sim_fastq_stats(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s sim-fq\n", argv[0], argv[1]);
        return 1;
    }
    const char* sim_fq_path = argv[2];

    size_t num_reads = 0;
    size_t num_frags = 0;
    size_t num_bases = 0;
    size_t num_contacts = 0;
    s_load_sim_read_info(sim_fq_path, num_reads, num_frags, num_bases, num_contacts);
    fprintf(stderr, "Reads: %zu\n", num_reads);
    fprintf(stderr, "Frags: %zu\n", num_frags);
    fprintf(stderr, "Bases: %zu\n", num_bases);
    fprintf(stderr, "Contacts: %zu\n", num_contacts);
    fprintf(stderr, "Avg-read-length: %d\n", num_bases / num_reads);
    fprintf(stderr, "Avg-frag-length: %d\n", num_bases / num_frags);

    return 0;
}
