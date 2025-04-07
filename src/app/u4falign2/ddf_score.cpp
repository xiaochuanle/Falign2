#include "../../corelib/frag_id.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/unpacked_seqdb.hpp"
#include "../../htslib/sam.h"

using namespace std;

static constexpr int kMinFragSize = 0;
static constexpr int kMaxFragSize = 2000;
static constexpr int kFragInc = 100;

static bam1_t*
s_read_next_bam(const char* sam_path, samFile* fp, sam_hdr_t* hdr)
{
    bam1_t* bam = bam_init1();
    int r = sam_read1(fp, hdr, bam);
    if (r == -1) {
        bam_destroy1(bam);
        return nullptr;
    }
    if (r < 0) {
        bam_destroy1(bam);
        HBN_ERR("ERROR at reading BAM record from %s", sam_path);
    }
    return bam;
}

int ddf_score_stats(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s bam\n", argv[0], argv[1]);
        return 1;
    }
    const char* bam_path = argv[2];

    int N = (kMaxFragSize - kMinFragSize) / kFragInc;
    pair<size_t, size_t>* stats = new pair<size_t, size_t>[N];
    for (int i = 0; i < N; ++i) stats[i].first = stats[i].second = 0;

    samFile* in = sam_open(bam_path, "rb");
    if (!in) HBN_ERR("Could not open file %s for reading", bam_path);
    hts_set_threads(in, 8);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", bam_path);
    bam1_t* bam = bam_init1();
    char tag[2] = { 'd', 's' };
    while (1) {
        int r = sam_read1(in, hdr, bam);
        if (r == -1) break;
        if (r < 0) break;
        int L = bam->core.l_qseq;
	if (L >= kMaxFragSize) continue;
	if (bam->core.qual < 5) continue;
        uint8_t* stag = bam_aux_get(bam, tag);
        hbn_assert(stag);
        int sc = bam_aux2i(stag);
        int p = L / kFragInc;
	if (p >= N) continue;
        hbn_assert(p < N);
        ++stats[p].first;
        stats[p].second += sc;
    }
    sam_hdr_destroy(hdr);
    sam_close(in);

    int L = 50;
    for (int i = 0; i < N; ++i) {
        double sc = 1.0 * stats[i].second / stats[i].first;
        fprintf(stdout, "%d\t%g\n", L, sc);
	L += kFragInc;
    }

    delete[] stats;
    return 0;
}
