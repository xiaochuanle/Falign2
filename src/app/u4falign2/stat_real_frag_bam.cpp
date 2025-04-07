#include "../../corelib/frag_id.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/unpacked_seqdb.hpp"
#include "../../htslib/sam.h"
#include "../map/restrict_enzyme_loci_list.hpp"

#include <cctype>

using namespace std;

static size_t g_num_reads = 0;
static size_t g_num_frags = 0;
static size_t g_num_contacts = 0;
static size_t g_enzyme_frag_ends = 0;

static const int kMinMapQ = 5;
static const int kConfidentEndDist = 20;

size_t low_pi = 0;

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

static void
s_stat_one_bam_list(const int enzyme_size, sam_hdr_t* hdr, vector<bam1_t*>& bam_list)
{
    ++g_num_reads;
    int N = bam_list.size();
    for (int i = 0; i < N; ++i) {
        bam1_t* bi = bam_list[i];
        int qi = bi->core.qual;
        if (qi >= kMinMapQ) ++g_num_frags;
        for (int j = i + 1; j < N; ++j) {
            bam1_t* bj = bam_list[j];
            if (bi->core.tid != bj->core.tid) continue;
            int qj = bj->core.qual;
            if (qj >= kMinMapQ) ++g_num_contacts;
        }
    }

    char tag[2];
    uint8_t* stag;
    for (int i = 0; i < N; ++i) {
        bam1_t* bam = bam_list[i];
        if (bam->core.qual < kMinMapQ) continue;

        tag[0] = 'q'; tag[1] = 's';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int qb = bam_aux2i(stag);

        tag[0] = 'q'; tag[1] = 'S';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int qS = bam_aux2i(stag);

        tag[0] = 'q'; tag[1] = 'e';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int qe = bam_aux2i(stag);

        tag[0] = 'q'; tag[1] = 'E';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int qE = bam_aux2i(stag);

        tag[0] = 'q'; tag[1] = 'l';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int ql = bam_aux2i(stag);

        tag[0] = 's'; tag[1] = 's';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int sb = bam_aux2i(stag);

        tag[0] = 'v'; tag[1] = 'S';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int vS = bam_aux2i(stag);

        tag[0] = 's'; tag[1] = 'e';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int se = bam_aux2i(stag);

        tag[0] = 'v'; tag[1] = 'E';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int vE = bam_aux2i(stag);

        tag[0] = 's'; tag[1] = 'l';
        stag = bam_aux_get(bam, tag);
        if (!stag) HBN_ERR("Could not extract tag %c%c in BAM record", tag[0], tag[1]);
        int sl = bam_aux2i(stag);

        bool r = (qb <= enzyme_size) || (sb <= enzyme_size) || (abs(qb - qS) <= enzyme_size) || (abs(sb - vS) <= enzyme_size);
        if (r) ++g_enzyme_frag_ends;

        r = (abs(ql - qe) <= enzyme_size) || (abs(sl - se <= enzyme_size)) || (abs(qe - qE) <= enzyme_size) || (abs(se - vE) <= enzyme_size);
        if (r) ++g_enzyme_frag_ends;
    }
}

int frag_and_contact_stats_main(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s %s enzyme bam\n", argv[0], argv[1]);
        exit(1);
    }
    const char* enzyme_seq = argv[2];
    const char* bam_path = argv[3];

    RestrictEnzyme enzyme;
    RestrictEnzyme_Init(enzyme_seq, &enzyme);
    samFile* in = sam_open(bam_path, "rb");
    if (!in) HBN_ERR("Could not open file %s for reading", bam_path);
    hts_set_threads(in, 8);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", bam_path);

    vector<bam1_t*> bam_list;
    bam1_t* bam = s_read_next_bam(bam_path, in, hdr);
    while (bam) {
        bam_list.clear();
        int rid1, fid1, sid1, soff1;
        int rid2, fid2, sid2, soff2;
        const char* qname = bam_get_qname(bam);
        extract_frag_id_from_name(qname, bam->core.l_qname, &rid1, &fid1, &sid1, &soff1);
        bam_list.push_back(bam);
        while ((bam = s_read_next_bam(bam_path, in, hdr))) {
            const char* next_qname = bam_get_qname(bam);
            extract_frag_id_from_name(next_qname, bam->core.l_qname, &rid2, &fid2, &sid2, &soff2);
            if (rid1 != rid2) break;
            bam_list.push_back(bam);
        }
        s_stat_one_bam_list(enzyme.enzyme_size, hdr, bam_list);
        for (bam1_t* b : bam_list) bam_destroy1(b);
        bam_list.clear();
    }
    sam_hdr_destroy(hdr);
    sam_close(in);

    fprintf(stderr, "Reads: %zu\n", g_num_reads);
    fprintf(stderr, "Frags: %zu\n", g_num_frags);
    fprintf(stderr, "Contacts: %zu\n", g_num_contacts);
    //double p1 = 1.0 * g_confident_frag_ends / (2 * g_num_frags);
    //fprintf(stderr, "Confident frag ends: %g\n", p1);
    double p2 = 1.0 * g_enzyme_frag_ends / (2 * g_num_frags);
    fprintf(stderr, "Enzyme frag ends: %g\n", p2);

    return 0;
}
