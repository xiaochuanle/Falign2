#include "simulate_pore_c_read.hpp"

#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <random>
#include <sstream>
#include <vector>

using namespace std;

const size_t kTargetReads = 20000;

static const int kShortRangeContactWeight = 4;
static const int kLongRangeContactWeight = 2;
static const int kInterChromContactWeight = 1;

struct ReadOrderPerc
{
    int min_order;
    int max_order;
    int perc;
};

#if 0
static const int kMaxReadOrder = 15;
static const int kReadOrderPercListSize = 4;
static const ReadOrderPerc kReadOrderPercList[] = {
    { 1, 1,  5  },
    { 2, 2,  15  },
    { 3, 5,  45  }, //45
    { 6, 15,  35  } //35
};
#else
static const int kMaxReadOrder = 15;
static const int kReadOrderPercListSize = 5;
static const ReadOrderPerc kReadOrderPercList[] = {
    { 1, 1,  8  },
    { 2, 2,  8 },
    { 3, 5,  24 },
    { 6, 10, 45 },
    { 11, 15, 15 }
};
#endif

static int gen_read_order(pair<int, double>* a, int c)
{
    double f = gen_one_prob();
    for (int i = 0; i < c - 1; ++i) {
        if (f <= a[i].second) return a[i].first;
        f -= a[i].second;
    }
    return a[c-1].first;
}

static void
s_parse_read_order_frac(const ReadOrderPerc* a, int c, vector<pair<int, double>>& order_frac_list)
{
    int sum = 0;
    for (int i = 0; i < c; ++i) sum += a[i].perc;
    hbn_assert(sum == 100);

    for (int i = 0; i < c; ++i) {
        int m = a[i].min_order;
        int M = a[i].max_order;
        hbn_assert(M >= m);
        double f = 1.0 * a[i].perc / (M - m + 1) / 100.0;
        for (int s = m; s <= M; ++s) order_frac_list.emplace_back(s, f);
    }

    fprintf(stderr, "Read-order dist:\n");
    for (auto& p : order_frac_list) fprintf(stderr, "%d\t%g\n", p.first, p.second);
}

bool s_is_overlapped_frag_pair(const SimFrag& x, const SimFrag& y)
{
    if (x.chr_id != y.chr_id) return false;
    if (x.from <= y.from && x.to > y.from) return true;
    if (y.from <= x.from && y.to > x.from) return true;
    return false;
}

static void
s_sample_inter_chrom_frag(SimFragLibrary& library, int main_chr_id, int num_frag, vector<SimFrag>& frags)
{
    int num_chr = library.chr_list.size();
    if (num_chr < 2) return;
    int added = 0;
    while (added < num_frag) {
        int ci = gen_one_positive_int() % num_chr;
        ChrFragInfo& chr = library.chr_list[ci];
        if (chr.chr_id == main_chr_id) continue;
        int ki = gen_one_positive_int() % chr.sfcc;
        SimFragChunk& chunk = chr.sfca[ki];
        int fi = gen_one_positive_int() % chunk.sfc;
        SimFrag& frag = chunk.sfa[fi];
        bool is_added = false;
        for (auto& s : frags) {
            if (s_is_overlapped_frag_pair(frag, s)) {
                is_added = true;
                break;
            }
        }
        if (is_added) continue;
        frags.push_back(frag);
        ++added;
    }
}

static void
s_sample_long_range_frag(ChrFragInfo& chr, int main_chunk_id, int num_frag, vector<SimFrag>& frags)
{
    if (chr.sfcc < 2) return;
    int num_chunks = chr.sfcc;
    int added = 0;
    while (added < num_frag) {
        int ki = gen_one_positive_int() % num_chunks;
        if (ki == main_chunk_id) continue;
        SimFragChunk& chunk = chr.sfca[ki];
        int fi = gen_one_positive_int() % chunk.sfc;
        SimFrag& frag = chunk.sfa[fi];
        bool is_added = false;
        for (auto& s : frags) {
            if (s_is_overlapped_frag_pair(s, frag)) {
                is_added = true;
                break;
            }
        }
        if (is_added) continue;
        frags.push_back(frag);
        ++added;
    }
}

static void
s_sample_short_range_frag(SimFragChunk& chunk, int main_frag_id, int num_frag, vector<SimFrag>& frags)
{
    int added = 0;
    while (added < num_frag) {
        int fi = gen_one_positive_int() % chunk.sfc;
        if (fi == main_frag_id) continue;
        SimFrag& frag = chunk.sfa[fi];
        bool is_added = false;
        for (auto& s : frags) {
            if (s_is_overlapped_frag_pair(frag, s)) {
                is_added = true;
                break;
            }
        }
        if (is_added) continue;
        frags.push_back(frag);
        ++added;
    }
}

static void
s_gen_fasta_header(HbnUnpackedDatabase* updb, int main_chr_id, int read_id, SimFrag* a, int c, ostringstream& hdr)
{
    hdr.str("");
    hdr << ">p0_" << updb->target_name(main_chr_id) << '_' << read_id << ' ';
    int strand_flag[] = { 0, 16 };
    for (int i = 0; i < c; ++i) {
        hdr << updb->target_name(a[i].chr_id) << '-' << a[i].from << '-' << a[i].to << '-'
            << strand_flag[a[i].strand] << '-' << '0';
        if (i < c - 1) hdr << ';';
    }
}

static char
s_complement_base(char c)
{
    switch (c) {
        case 'a':
        case 'A':
            return 'T';
        case 'c':
        case 'C':
            return 'G';
        case 'g':
        case 'G':
            return 'C';
        case 't':
        case 'T':
            return 'A';
        default:
            HBN_ERR("Illegal DNA base %c", c);
    }
}

static void
s_gen_fasta_sequence(HbnUnpackedDatabase* updb, const int enzyme_size, SimFrag* frags, int n, string& seq)
{
    seq.clear();
    for (int i = 0; i < n; ++i) {
        SimFrag& frag = frags[i];
        const char* chr = updb->target_sequence(frag.chr_id, FWD);
        if (frag.strand == FWD) {
            for (int i = frag.from; i < frag.to - enzyme_size; ++i) {
                seq += toupper(chr[i]);
            }
        } else {
		for (int i = frag.to; i > frag.from + enzyme_size; --i) {
			seq += s_complement_base(chr[i-1]);
		}
        }
    }
}

static void
s_describe_length_dist(int* a, int ac)
{
    if (!ac) return;
    sort(a, a + ac);
    size_t sum = 0;
    for (int i = 0; i < ac; ++i) sum += a[i];
    fprintf(stderr, "min: %d\n", a[0]);
    fprintf(stderr, "max: %d\n", a[ac-1]);
    fprintf(stderr, "avg: %d\n", sum/ac);
    fprintf(stderr, "median: %d\n", a[ac/2]);
    fprintf(stderr, "N25: %d\n", a[(int)(ac*0.25)]);
    fprintf(stderr, "N50: %d\n", a[(int)(ac*.5)]);
    fprintf(stderr, "N75: %d\n", a[(int)(ac*0.75)]);
}

struct SimReadInfo
{
    int chr_id;
    int read_id;
    size_t frag_offset;
    int frag_cnt;
    int length;
    bool has_short_frag;
    bool is_valid;
};

void gen_simulated_pore_c_reads(SimFragLibrary& frag_library, HbnUnpackedDatabase* updb, 
    const char* enzyme_seq, const char* output)
{
    vector<pair<int, double>> order_frac_list;
    s_parse_read_order_frac(kReadOrderPercList, kReadOrderPercListSize, order_frac_list);

    vector<SimFrag> srf, lrf, icf, mf, frags, read_frag_list;
    vector<SimReadInfo> read_list;
    int cnt = 0;
    int order_stats[kMaxReadOrder+1];
    fill(order_stats, order_stats + kMaxReadOrder + 1, 0);
    int strand_stats[] = { 0, 0 };
    RestrictEnzyme enzyme;
    RestrictEnzyme_Init(enzyme_seq, &enzyme);
    random_device rd;
    mt19937 gen(rd());
    hbn_dfopen(out, output, "w");
    ostringstream hdr;
    string seq;
    size_t num_reads = 0, num_bases = 0;
    vector<int> read_length_list, frag_length_list;

    int num_chr = frag_library.chr_list.size();
    for (int i = 0; i < num_chr; ++i) {
        ChrFragInfo& chr = frag_library.chr_list[i];
        int read_id = 0;
        for (int j = 0; j < chr.sfcc; ++j) {
            SimFragChunk& chunk = chr.sfca[j];
            for (int k = 0; k < chunk.sfc; ++k) {
                int order = gen_read_order(order_frac_list.data(), order_frac_list.size());
                if (order == 1) {
                    frags.clear();
                    frags.push_back(chunk.sfa[k]);
                } else {
                    srf.clear(); s_sample_short_range_frag(chunk, k, order * kShortRangeContactWeight, srf);
                    lrf.clear(); s_sample_long_range_frag(chr, j, order * kLongRangeContactWeight, lrf);
                    icf.clear(); s_sample_inter_chrom_frag(frag_library, i, order * kInterChromContactWeight, icf);
                    mf.clear();
                    mf.insert(mf.end(), srf.begin(), srf.end());
                    mf.insert(mf.end(), lrf.begin(), lrf.end());
                    mf.insert(mf.end(), icf.begin(), icf.end());

                    shuffle(mf.begin(), mf.end(), gen);
                    shuffle(mf.begin(), mf.end(), gen);
                    shuffle(mf.begin(), mf.end(), gen);
                    shuffle(mf.begin(), mf.end(), gen);
                    shuffle(mf.begin(), mf.end(), gen);

                    int n = mf.size();
                    frags.clear();
                    for (int s = 0; s < order - 1 && s < n; ++s) frags.push_back(mf[s]);
                    frags.push_back(chunk.sfa[k]);

                    sort(frags.begin(), frags.end(), [](const SimFrag& x, const SimFrag& y) {
                        return (x.chr_id < y.chr_id) || (x.chr_id == y.chr_id && x.from < y.from);
                    });
                }

                SimReadInfo read;
                read.chr_id = i;
                read.read_id = read_id++;
                read.frag_offset = read_frag_list.size();
                read.frag_cnt = frags.size();
                read.has_short_frag = false;
                read.length = 0;
                read.is_valid = false;
                for (auto& f : frags) {
                    read.length += f.to - f.from;
                    if (f.to - f.from < 200) read.has_short_frag = true;
                }
                read_list.push_back(read);
                read_frag_list.insert(read_frag_list.end(), frags.begin(), frags.end());
            }
        }
    }

    shuffle(read_list.begin(), read_list.end(), gen);
    shuffle(read_list.begin(), read_list.end(), gen);
    shuffle(read_list.begin(), read_list.end(), gen);
    shuffle(read_list.begin(), read_list.end(), gen);
    shuffle(read_list.begin(), read_list.end(), gen);
    int n = 0;
    for (auto& r : read_list) {
        if (r.has_short_frag == false) continue;
        ++n;
        r.is_valid = true;
        if (n == kTargetReads) break;
    }
    if (n < kTargetReads) {
        for (auto& r : read_list) {
	    if (r.is_valid) continue;
            ++n;
            r.is_valid = true;
            if (n == kTargetReads) break;
        }
    }
    sort(read_list.begin(), read_list.end(), [](const SimReadInfo& x, const SimReadInfo& y) {
        return (x.chr_id < y.chr_id) || (x.chr_id == y.chr_id && x.read_id < y.read_id);
    });
    for (auto& r : read_list) {
        if (!r.is_valid) continue;
        SimFrag* sfa = read_frag_list.data() + r.frag_offset;
        int sfc = r.frag_cnt;
        for (int i = 0; i < sfc; ++i) {
            frag_length_list.push_back(sfa[i].to - sfa[i].from);
            ++strand_stats[sfa[i].strand];
        }
        ++order_stats[sfc];
        s_gen_fasta_header(updb, r.chr_id, r.read_id, sfa, sfc, hdr);
        s_gen_fasta_sequence(updb, enzyme.enzyme_size, sfa, sfc, seq);
        fprintf(out, "%s\n", hdr.str().c_str());
        fprintf(out, "%s\n", seq.c_str());
        ++num_reads;
        num_bases += seq.size();
        read_length_list.push_back(seq.size());        
    }
    hbn_fclose(out);
    fprintf(stderr, "====> Read order-stats:\n");
    for (int i = 1; i <= kMaxReadOrder; ++i) fprintf(stderr, "%d\t%d\n", i, order_stats[i]);
    fprintf(stderr, "====> Strand-stats:\n");
    for (int i = 0; i < 2; ++i) fprintf(stderr, "%d\t%d\n", i, strand_stats[i]);
    fprintf(stderr, "====> Read length stats:\n");
    s_describe_length_dist(read_length_list.data(), read_length_list.size());
    fprintf(stderr, "====> Frag length stats:\n");
    s_describe_length_dist(frag_length_list.data(), frag_length_list.size());
    string size = NStr::UInt8ToString_DataSize(num_bases);
    HBN_LOG("Extract %zu reads (%s)", num_reads, size.c_str());
}
