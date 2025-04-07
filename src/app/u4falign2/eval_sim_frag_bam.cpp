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

struct FragMapInfo
{
    int read_id;
    int fwd_qb, fwd_qe;
    int ctg_id;
    int from;
    int to;
    int mapQ;
    bool is_found;
    bool is_correct;

    const char* chr_name;
    const char* qname;
};

std::ostream& operator<<(std::ostream& os, const FragMapInfo& fmi)
{
    return os << fmi.qname << '\t' << fmi.chr_name << '\t' << fmi.from << '\t' << fmi.to;
}

static bool 
s_parse_one_frag_bam(sam_hdr_t* hdr, bam1_t* bam, SeqName2IdMap& read_name2id, FragMapInfo& fmi)
{
    static int cnt = 0;
    fmi.ctg_id = bam->core.tid;
    fmi.from = bam->core.pos;
    fmi.to = fmi.from;

    if (bam->core.n_cigar == 0) return false;
    uint32_t* cigar = bam_get_cigar(bam);
    for (int i = 0; i < bam->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int num = bam_cigar_oplen(cigar[i]);
        if (op == 'M' || op == '=' || op == 'X' || op == 'D') fmi.to += num;
    }

    const char* qname = bam_get_qname(bam);
    int qnl = strlen(qname);
    while (qnl && qname[qnl-1] != '_') --qnl;
    --qnl;
    hbn_assert(qnl > 0 && qname[qnl] == '_');
    NStr::CTempString cqname(qname, qnl);
    fmi.read_id = read_name2id.GetIdFromNameSafe(cqname);
    fmi.mapQ = bam->core.qual;

    fmi.qname = read_name2id.GetSeqName(fmi.read_id);
    fmi.chr_name = sam_hdr_tid2name(hdr, fmi.ctg_id);

    char tag[2];
    /// qdir
    tag[0] = 'q'; tag[1] = 'd';
    uint8_t* stag = bam_aux_get(bam, tag);
    if (!stag) HBN_ERR(" Could not find tag %c%c in BAM record", tag[0], tag[1]);
    int qdir = bam_aux2i(stag);
    hbn_assert(qdir == FWD || qdir == REV);

    /// qb
    tag[0] = 'q'; tag[1] = 's';
    stag = bam_aux_get(bam, tag);
    if (!stag) HBN_ERR(" Could not find tag %c%c in BAM record", tag[0], tag[1]);
    int qb = bam_aux2i(stag);

    /// qe
    tag[0] = 'q'; tag[1] = 'e';
    stag = bam_aux_get(bam, tag);
    if (!stag) HBN_ERR(" Could not find tag %c%c in BAM record", tag[0], tag[1]);
    int qe = bam_aux2i(stag);

    /// ql
    tag[0] = 'q'; tag[1] = 'l';
    stag = bam_aux_get(bam, tag);
    if (!stag) HBN_ERR(" Could not find tag %c%c in BAM record", tag[0], tag[1]);
    int ql = bam_aux2i(stag);

    if (qdir == FWD) {
        fmi.fwd_qb = qb;
        fmi.fwd_qe = qe;
    } else {
        fmi.fwd_qb = ql - qe;
        fmi.fwd_qe = ql - qb;
    }

    //if (cnt < 10) fprintf(stderr, "read-id = %d, chr-id = %d, from = %d, to = %d, mapQ = %d\n",
     //   fmi.read_id, fmi.ctg_id, fmi.from, fmi.to, fmi.mapQ);
    ++cnt;

    return true;
}

static void
s_parse_one_gt_frag_map(const char* s, const int sl, sam_hdr_t* hdr, FragMapInfo& fmi)
{
    int i = 0, j = 0;
    while (j < sl && s[j] != '-') ++j;
    hbn_assert(j < sl);
    hbn_assert(s[j] == '-');
    string chr_name(s+i, s+j);
    fmi.ctg_id = sam_hdr_name2tid(hdr, chr_name.c_str());
    if (fmi.ctg_id == -1) HBN_ERR("Unknown chromosome %s\n", chr_name.c_str());
    if (fmi.ctg_id == -2) HBN_ERR("Could not parse BAM header");

    i = j + 1;
    j = i;
    fmi.from = 0;
    hbn_assert(i < sl);
    while (j < sl) {
        if (s[j] == '-') break;
        fmi.from = fmi.from * 10 + s[j] - '0';
        ++j;
    }

    hbn_assert(j < sl);
    hbn_assert(s[j] == '-');
    i = j + 1;
    j = i;
    fmi.to = 0;
    while (j < sl) {
        if (s[j] == '-') break;
        fmi.to = fmi.to * 10 + s[j] - '0';
        ++j;
    }
}

static void
s_load_gt_read_maps(const char* fq_path, 
    sam_hdr_t* hdr,
    vector<char>& read_name_list,
    SeqName2IdMap& read_name2id, 
    vector<pair<size_t, int>>& read_map_info,
    vector<FragMapInfo>& frag_maps)
{
    ifstream in(fq_path);
    string line;
    FragMapInfo fmi; fmi.is_found = false;
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
        read_name_list.insert(read_name_list.end(), line.c_str() + i, line.c_str() + j);
        read_name_list.push_back('\0');

        i = j + 1;
        j = i;
        hbn_assert(i < N);
        while (j < N) {
            if (isspace(line[j])) break;
            ++j;
        }
        hbn_assert(j < N);
        hbn_assert(isspace(line[j]));

        pair<size_t, int> ip(frag_maps.size(), 0);
        fmi.read_id = read_id++;

        i = j + 1;
        hbn_assert(i < N);
        while (i < N) {
            j = i + 1;
            while (j < N && line[j] != ';') ++j;
            s_parse_one_gt_frag_map(line.c_str() + i, j - i, hdr, fmi);
            ++ip.second;
            frag_maps.push_back(fmi);
            i = j + 1;
        }
        read_map_info.push_back(ip);
    }

    size_t N = read_name_list.size();
    size_t i = 0;
    int add = 0;
    while (i < N) {
        size_t j = i + 1;
        while (read_name_list[j]) ++j;
        read_name2id.add_one_name(read_name_list.data() + i, j - i);
        ++add;
        i = j + 1;
    }
    HBN_LOG("add %d read names", add);
    read_name2id.build_name2id_map();
}

static void
s_load_cmp_read_maps(samFile* in, sam_hdr_t* hdr, SeqName2IdMap& read_name2id, 
    vector<FragMapInfo>& cmp_list)
{
    bam1_t* bam = bam_init1();
    while (1) {
        int r = sam_read1(in, hdr, bam);
        if (r == -1) break;
        if (r < 0) HBN_ERR("Could not parse BAM record");
        FragMapInfo cmp;
        s_parse_one_frag_bam(hdr, bam, read_name2id, cmp);
        cmp_list.push_back(cmp);
    }
    bam_destroy1(bam);
}

static void
s_compute_metrics(vector<FragMapInfo>& gt_frag_maps, vector<pair<size_t, int>>& read_map_info, 
    vector<FragMapInfo>& cmp_frag_maps, const int min_mapQ)
{
    const double kMinOvlpFrac = 0.8;
    const int kConfidentEndDist = 20, kCorrectEndDist = 4;
    static int cnt = 0;
    int n_cmp = 0, n_tp = 0, n_fp = 0, n_found = 0;
    int n_gt_end = 0, n_confident_end = 0, n_tp_end = 0;
    for (auto& cfm : cmp_frag_maps) {
        cfm.is_correct = false;
        if (cfm.mapQ < min_mapQ) continue;
        ++n_cmp;
        FragMapInfo* Gfma = gt_frag_maps.data() + read_map_info[cfm.read_id].first;
        int Gfmc = read_map_info[cfm.read_id].second;
        hbn_assert(Gfmc > 0);
        bool is_tp = false;
        FragMapInfo* TP = nullptr;
        for (int i = 0; i < Gfmc; ++i) {
            FragMapInfo* F = Gfma + i;
            hbn_assert(F->read_id == cfm.read_id);
            if (F->ctg_id != cfm.ctg_id) continue;
            int ovlp = 0;
            int sb = max(F->from, cfm.from);
            int se = min(F->to, cfm.to);
            if (se > sb) ovlp = se - sb;
            if (!ovlp) continue;
            F->is_found = true;
            ++n_found;
            if (ovlp >= (cfm.to - cfm.from) * kMinOvlpFrac && ovlp >= (F->to - F->from) * kMinOvlpFrac) {
                is_tp = true;
                TP = F;
                break;
            }
        }
        //if (!is_tp) cerr << cfm << '\n';
        if (is_tp) {
            ++n_tp;
            cfm.is_correct = true;
            n_gt_end += 2;
            if (abs(cfm.from - TP->from) <= kCorrectEndDist) ++n_tp_end;
            if (abs(cfm.to - TP->to) <= kCorrectEndDist) ++n_tp_end;
            if (abs(cfm.from - TP->from) <= kConfidentEndDist) ++n_confident_end;
            if (abs(cfm.to - TP->to) <= kConfidentEndDist) ++n_confident_end;
        } else {
            ++n_fp;
        }
    }
    //fprintf(stderr, "**** Metrics for mapQ>=%d:\n", min_mapQ);
    fprintf(stderr, "**** Mapping position metrics:\n");
    int n_gt = gt_frag_maps.size();
    HBN_LOG("n-gt = %d, n-cmp = %d, n-found = %d, n-tp = %d, n-fp = %d", n_gt, n_cmp, n_found, n_tp, n_fp);
    double acc = 1.0 * n_tp / n_cmp;
    double recall = 1.0 * n_tp / n_gt;
    HBN_LOG("acc: %g", acc);
    HBN_LOG("recall: %g", recall);
    //fprintf(stderr, "**** Confident fragment boundary metrics:\n");
    //recall = 1.0 * n_confident_end / n_gt_end;
    //HBN_LOG("n-gt: %d, n-tp: %d", n_gt_end, n_confident_end);
    //HBN_LOG("acc: %g", recall);
    fprintf(stderr, "**** Correct fragment boundary metrics:\n");
    recall = 1.0 * n_tp_end / n_gt_end;
    HBN_LOG("n-gt: %d, n-tp: %d", n_gt_end, n_tp_end);
    HBN_LOG("acc: %g", recall);
}

static void
s_compute_contact_metrics(vector<FragMapInfo>& gt_frag_maps, vector<FragMapInfo>& cmp_frag_maps, int min_mapQ)
{
    FragMapInfo* a = cmp_frag_maps.data();
    size_t c = cmp_frag_maps.size();
    int n_cmp = 0, n_tp = 0, n_fp = 0, n_found = 0;
    size_t i = 0;
    while (i < c) {
        size_t j = i + 1;
        while (j < c && a[i].read_id == a[j].read_id) ++j;

        for (size_t x = i; x < j; ++x) {
            if (a[x].mapQ < min_mapQ) continue;
            for (size_t y = x + 1; y < j; ++y) {
                if (a[y].mapQ < min_mapQ) continue;
                ++n_cmp;
                if (a[x].is_correct && a[y].is_correct) {
                    ++n_tp;
                } else {
                    ++n_fp;
                }
            }
        }

        i = j;
    }

    a = gt_frag_maps.data();
    c = gt_frag_maps.size();
    i = 0;
    size_t n_gt = 0;
    while (i < c) {
        size_t j = i + 1;
        while (j < c && a[i].read_id == a[j].read_id) ++j;

        int n = j - i;
        n_gt += n * (n-1) / 2;

        i = j;
    }

    fprintf(stderr, "**** Pairwise contact metrics:\n");
    HBN_LOG("n-gt = %d, n-cmp = %d, n-tp = %d, n-fp = %d", n_gt, n_cmp, n_found, n_tp, n_fp);
    double acc = 1.0 * n_tp / n_cmp;
    double recall = 1.0 * n_tp / n_gt;
    HBN_LOG("acc: %g", acc);
    HBN_LOG("recall: %g", recall);
}

int eval_sim_frag_bam(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s reads frag-bam\n", argv[0], argv[1]);
        return 1;
    }
    const char* read_path = argv[2];
    const char* bam_path = argv[3];

    samFile* in = sam_open(bam_path, "rb");
    sam_hdr_t* hdr = sam_hdr_read(in);

    vector<char> read_name_list;
    SeqName2IdMap read_name2id;
    vector<pair<size_t, int>> read_map_info;
    vector<FragMapInfo> gt_frag_maps;
    s_load_gt_read_maps(read_path, hdr, read_name_list, read_name2id, read_map_info, gt_frag_maps);
    vector<FragMapInfo> cmp_frag_maps;
    s_load_cmp_read_maps(in, hdr, read_name2id, cmp_frag_maps);

    s_compute_metrics(gt_frag_maps, read_map_info, cmp_frag_maps, 5);

    s_compute_contact_metrics(gt_frag_maps, cmp_frag_maps, 5);

    sam_hdr_destroy(hdr);
    sam_close(in);

    return 0;
}