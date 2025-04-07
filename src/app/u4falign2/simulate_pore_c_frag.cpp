#include "simulate_pore_c_frag.hpp"

#include "../map/restrict_enzyme_loci_list.hpp"

#include <algorithm>
#include <random>

using namespace std;

random_device G_rd;
mt19937 G_gen(G_rd());
uniform_real_distribution<> G_real_gen(0.0, 1.0);
uniform_int_distribution<> G_int_gen(0);

static const int kFragOrderList[] = { 1, 2, 3 };
static const double kFragOrderDist[] = { 0.8, 0.1, 0.1 };
//static const double kFragOrderDist[] = { 0.7, 0.15, 0.15 };
static const int kFragOrderListSize = 3;
static const int kMaxFragOrder = 3;

static const int kMinFragSize = 100;
static const int kMaxFragSize = 5000;

double gen_one_prob()
{
    return G_real_gen(G_gen);
}

int gen_one_positive_int()
{
    return G_int_gen(G_gen);
}

static int gen_frag_order()
{
    double f = gen_one_prob();
    for (int i = 0; i < kFragOrderListSize - 1; ++i) {
        if (f <= kFragOrderDist[i]) return kFragOrderList[i];
        f -= kFragOrderDist[i];
    }
    return kFragOrderList[kFragOrderListSize-1];
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

static void
s_gen_frag_for_one_chr(const int chr_id, const char* chr_name, const char* chr, const int chr_size, 
    RestrictEnzyme* enzyme, vector<SimFrag>& frag_list)
{
    vector<int> enzyme_pos_list;
    extract_enzyme_loci_list_for_one_seq(chr, chr_size, enzyme, enzyme_pos_list);
    int* ea = enzyme_pos_list.data() + 1;
    int ec = enzyme_pos_list.size();
    if (ec < 1000) return;
    ec -= 3;
    hbn_assert(ea[0] > 0);
    hbn_assert(ea[ec-1] + enzyme->enzyme_size <= chr_size);
    fprintf(stderr, "number of enzyme pos: %d\n", ec);

    SimFrag frag;
    int order_stats[kMaxFragOrder+1];
    fill(order_stats, order_stats + kMaxFragOrder + 1, 0);
    int strand_cnts[] = { 0, 0 };
    vector<int> length_list;
    for (int i = 0; i < ec - 1; ++i) {
        int order = gen_frag_order();
        int j = min(i+order, ec-1);
        int from = ea[i];
        int to = ea[j] + enzyme->enzyme_size;
        order = j - i;

        if (to - from < kMinFragSize) continue;
        if (to - from > kMaxFragSize) continue;
        bool has_amb_base = false;
        for (int s = from; s < to; ++s) {
            int c = chr[s];
            if (nst_nt4_table[c] > 3) {
                has_amb_base = true;
                break;
            }
        }
        if (has_amb_base) continue;
	int m = 0;
	for (int s = from; s < to; ++s) {
		if (chr[s] >= 'a' && chr[s] <= 'z') ++m;
	}
	if (m > (to - from) * 0.4) {
		fprintf(stderr, "repeat frag %s %d-%d\n", chr_name, from, to);
		for (int s = from; s < to; ++s) fprintf(stderr, "%c", chr[s]); 
		fprintf(stderr, "\n");
		has_amb_base = true;
	}
	if (has_amb_base) continue;

        ++order_stats[order];
        frag.chr_id = chr_id;
        frag.chunk_id = -1;
        frag.from = from;
        frag.to = to;
        frag.strand = gen_one_positive_int()%2;
        ++strand_cnts[frag.strand];
        frag_list.push_back(frag);
        length_list.push_back(to - from);
    }

    fprintf(stderr, "Frag-order stats for %d:%s:%d\n", chr_id, chr_name, chr_size);
    for (int i = 1; i <= kMaxFragOrder; ++i) fprintf(stderr, "%d\t%d\n", i, order_stats[i]);
    fprintf(stderr, "Strand-stats:\n");
    for (int i = 0; i < 2; ++i) fprintf(stderr, "%d\t%d\n", i, strand_cnts[i]);
    s_describe_length_dist(length_list.data(), length_list.size());
}

static bool
s_is_unlocalized_or_unplaced_or_alternate_chromosome(const char* name)
{
    const int N = 3;
    const char* kTargetList[] = { "UN", "RANDOM", "ALT" };
    string sname = name;
    for (auto& c : sname) c = toupper(c);
    for (int i = 0; i < N; ++i) if (strstr(sname.c_str(), kTargetList[i])) return true;
    return false;
}

SimFragLibrary::SimFragLibrary(HbnUnpackedDatabase* updb, const char* enzyme_seq)
{
    RestrictEnzyme enzyme;
    RestrictEnzyme_Init(enzyme_seq, &enzyme);

    const int num_chr = updb->num_targets();
    for (int i = 0; i < num_chr; ++i) {
        const char* name = updb->target_name(i);
        if (s_is_unlocalized_or_unplaced_or_alternate_chromosome(name)) continue;
        const int chr_size = updb->target_size(i);
        const char* chr = updb->target_sequence(i, FWD);
        s_gen_frag_for_one_chr(i, name, chr, chr_size, &enzyme, frag_list);
    }

    SimFrag* sfa = frag_list.data();
    size_t sfc = frag_list.size();
    size_t sfi = 0;
    while (sfi < sfc) {
        size_t sfj = sfi + 1;
        while (sfj < sfc && sfa[sfi].chr_id == sfa[sfj].chr_id) ++sfj;

        ChrFragInfo chr;
        chr.chr_id = sfa[sfi].chr_id;
        chr.sfa = sfa + sfi;
        chr.sfc = sfj - sfi;
        chr.__frag_offset = sfi;
        chr.__chunk_offset = chunk_list.size();
        chr.sfca = nullptr;
        chr.sfcc = 0;

        int chunk_size = (chr.sfc + 9) / 10;
        int chunk_id = 0;
        size_t i = sfi;
        while (i < sfj) {
            size_t j = min(i+ chunk_size, sfj);
            for (size_t k = i; k < j; ++k) sfa[k].chunk_id = chunk_id;

            SimFragChunk chunk;
            chunk.chr_id = sfa[sfi].chr_id;
            chunk.chunk_id = chunk_id;
            chunk.__frag_offset = i;
            chunk.sfa = sfa + i;
            chunk.sfc = j - i;
            chunk_list.push_back(chunk);
            ++chr.sfcc;
            ++chunk_id;
            i = j;
        }

        chr_list.push_back(chr);
        sfi = sfj;
    }
    for (auto& chr : chr_list) chr.sfca = chunk_list.data() + chr.__chunk_offset;

#if 0
    for (auto& chr : chr_list) {
        fprintf(stderr, "chr-id: %d, frag-offset: %zu, frag-cnt: %d, chunk-offset: %zu, chunk-cnt: %d\n",
            chr.chr_id, chr.__frag_offset, chr.sfc, chr.__chunk_offset, chr.sfcc);
        for (int i = 0; i < chr.sfcc; ++i) {
            SimFragChunk* chunk = chr.sfca + i;
            fprintf(stderr, "\tchunk-id: %d, frag-offset: %zu, frag-cnt: %d\n", i, chunk->__frag_offset, chunk->sfc);
        }
    }
#endif
}
