#include "gapped_align.hpp"

#include "../../corelib/pdqsort.h"
#include "../../sw/hbn_traceback_aux.h"

#include <vector>

using namespace std;

static const int kMaxW = 50;
static const double kEpsilon = 0.25;
static const int kEnzymeNbhd = 50;

static const int kMinCov = 3;

static size_t gExtendedHits = 0;
void dump_extended_hits()
{
    fprintf(stderr, "Extended hits: %zu\n", gExtendedHits);
}

static BOOL 
chain_seed_list_is_contained_in_align(PoreCInitHit* align_list, const int align_cnt,
    HbnChainInfo* chains, HbnChainInfo* chain)
{

	for (int i = 0; i < align_cnt; ++i) {
		PoreCInitHit* ai = align_list + i;
		if (ai->sid != chain->sid || ai->sdir != chain->sdir) continue;
		int qb = max(ai->qoff, chain->qb);
		int qe = min(ai->qend, chain->qe);
		if (qb > qe) continue;
		int sb = max(ai->soff, chain->soff);
		int se = min(ai->send, chain->send);
		if (sb > se) continue;
		return TRUE;
	}

    const int E = 0;
    int min_cov = kMinCov;
    for (int i = 0; i < align_cnt; ++i) {
        HbnChainInfo* c = chains + align_list[i].parent;
	if (c->qb <= chain->qb + E && c->qe + E >= chain->qe) {
	    if (c->sc > chain->sc * 2) return TRUE;
            if (align_list[i].id >= min_cov) return TRUE;
            ++align_list[i].id;
            return FALSE;
        }
    }

    return FALSE;
}

int s_extend_one_mg_list_a(PoreCQuery* query,
    int sid, 
    int sdir,
    const char* subject,
    const int subject_size,
    const pu64_t* a, int n_a,
    HbnTracebackData* tbck_data,
    PoreCInitHit& align, bool xxx)
{
    if (!tbck_data->map_one_chain(query->fwd_qs, query->size, subject, subject_size, a, n_a, xxx)) return 0;

    align.qoff = tbck_data->qoff;
    align.qend = tbck_data->qend;
    align.sid = sid;
    align.sdir = sdir;
    align.soff = tbck_data->soff;
    align.send = tbck_data->send;
    align.map_score = tbck_data->score;
    align.pi = tbck_data->pi;

   return 1;
}

static BOOL subject_subseq_cov_is_full(u8* cov_stats, int soff, int send, const int kMinCov, const int kMaxCov)
{
    for (int i = soff; i < send; ++i) if (cov_stats[i] < kMinCov) return FALSE;

    int n = 0;
    for (int i = soff; i < send; ++i) 
        if (cov_stats[i] >= kMaxCov) ++n;
    if (send - soff > n) return FALSE;
    return TRUE; 
}

static void
s_invalidate_hits(const u8* cov_stats, const int min_cov, HbnChainInfo* hita, int hitc, int subject_size, int hit_i)
{
    int from = hita[hit_i].qb, to = hita[hit_i].qe;
    while (from && cov_stats[from-1] >= min_cov) --from;
    while (to < subject_size && cov_stats[to] >= min_cov) ++to;

    for (int i = hit_i + 1; i < hitc; ++i) {
        if (hita[i].sc < 0) continue;
        if (hita[i].qb >= from && hita[i].qe <= to) hita[i].sc = -hita[i].sc;
    }
}

void extend_ddf_chain_list(HbnTracebackData* tbck_data,
    const HbnUnpackedDatabase* updb,
    PoreCQuery* query,
    HbnChainInfo* ac_a,
    int ac_n,
    pu64_t* ac_km,
    std::vector<u8>& cov_list,
    std::vector<PoreCInitHit>& align_list,
    std::string& align_strings)
{
    align_list.clear();
    align_strings.clear();
    cov_list.clear();
    cov_list.resize(query->size); fill(cov_list.begin(), cov_list.end(), 0);

	int added = 0;
    PoreCInitHit align;
    for (int i = 0; i < ac_n; ++i) {
        if (ac_a[i].sc <= 0) continue;
        pu64_t* F = ac_km + ac_a[i].offset;
        int n_mg = ac_a[i].sc;
        int qb = ac_a[i].qb;
        int qe = ac_a[i].qe;
        int sb = ac_a[i].soff;
        int se = ac_a[i].send;
        const char* subject = updb->target_sequence(ac_a[i].sid, ac_a[i].sdir);
        const int subject_size = updb->target_size(ac_a[i].sid);

        if (chain_seed_list_is_contained_in_align(align_list.data(), align_list.size(), ac_a, ac_a+i)) continue;
        if (subject_subseq_cov_is_full(cov_list.data(), qb, qe, kMinCov, kMinCov)) {
            s_invalidate_hits(cov_list.data(), kMinCov, ac_a, ac_n, query->size, i);
            continue;
        }
        //++gExtendedHits;
        for (int s = qb; s < qe; ++s) cov_list[s] = min<u8>(100, cov_list[s]+1);

        //cerr << "extend hit " << i << '\t' << lc[i] << '\n';

        int r = s_extend_one_mg_list_a(query, ac_a[i].sid, ac_a[i].sdir, subject, subject_size, F, n_mg, tbck_data, align, true);
        if (!r) continue;
        tbck_data->cigar_to_align_string(query->fwd_rqs, subject);
        align.qas_offset = align_strings.size();
        align_strings.insert(align_strings.end(), tbck_data->qas, tbck_data->qae);
        align.sas_offset = align_strings.size();
        align_strings.insert(align_strings.end(), tbck_data->sas, tbck_data->sae);
        align.as_size = tbck_data->qae - tbck_data->qas;
        align.ddf_score = n_mg;
        align.parent = i;
        align.id = 0;
        align.ssize = subject_size;
        align_list.push_back(align);
	    if (++added == 1000) break;
    }

    PoreCInitHit* fasa = align_list.data();
    int fasc = align_list.size();
    pdqsort(fasa, fasa + fasc, [](const PoreCInitHit& a, const PoreCInitHit& b) { 
        return (a.sid < b.sid) || (a.sid == b.sid && a.sdir < b.sdir) || (a.sid == b.sid && a.sdir == b.sdir && a.soff < b.soff);
    });
    int score_list[fasc]; for (int i = 0; i < fasc; ++i) score_list[i] = fasa[i].map_score;
    int pred_list[fasc]; fill(pred_list, pred_list + fasc, -1);
    int succ_list[fasc]; fill(succ_list, succ_list + fasc, 0);
    int avail_list[fasc]; fill(avail_list, avail_list + fasc, 1);

    for (int i = 1; i < fasc; ++i) {
        PoreCInitHit* fi = fasa + i;
        int max_score = fasa[i].map_score;
        int max_j = -1;
        for (int j = i - 1; j >= 0; --j) {
            PoreCInitHit* fj = fasa + j;
            if (fi->sid != fj->sid) break;
            if (fi->sdir != fj->sdir) break;
            if (fj->qend > fi->qoff || fj->send > fi->soff) continue;
	        if (fi->qoff - fj->qend > 500 && fi->soff - fj->send > 500) continue;
            int q_d = fi->qend - fj->qoff;
            int s_d = fi->send - fj->soff;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > 0.25) continue;
            int score = score_list[j] + fasa[i].map_score;
            if (score > max_score) {
                max_score = score;
                max_j = j;
            }
        }
        score_list[i] = max_score;
        pred_list[i] = max_j;
    }

    for (int i = 0; i < fasc; ++i) {
        if (pred_list[i] >= 0) succ_list[pred_list[i]] = 1;
    }
    vector<pair<int, int>> idx_and_score_list;
    for (int i = 0; i < fasc; ++i) {
        if (succ_list[i] == 0 && score_list[i] > 0) {
            idx_and_score_list.emplace_back(i, score_list[i]);
        }
    }
    pdqsort(idx_and_score_list.begin(), 
         idx_and_score_list.end(),
         [](const pair<int, int>& a, const pair<int, int>& b)->bool {
             return (a.second > b.second) || (a.second == b.second && a.first < b.first);
         });

    vector<PoreCInitHit> colin_align_list;
    vector<pu64_t> seeds;
    for (auto& ias : idx_and_score_list) {
        if (!avail_list[ias.first]) continue;
        int p = ias.first;
        colin_align_list.clear();
        while (p >= 0) {
            avail_list[p] = 0;
            colin_align_list.push_back(align_list[p]);
            p = pred_list[p];
        }
        reverse(colin_align_list.begin(), colin_align_list.end());
        if (colin_align_list.size() < 2) continue;

        seeds.clear();
        int last_qoff = -1, last_soff = -1;
        for (auto& ca : colin_align_list) {
            HbnChainInfo* caa = ac_a + ca.parent;
            pu64_t* ka = ac_km + caa->offset;
            int kc = caa->sc;
            for (int s = 0; s < kc; ++s) {
                int qoff = (int32_t)ka[s].y + 1 - (ka[s].y>>32&0xff);
                int soff = (int32_t)ka[s].x + 1 - (ka[s].y>>32&0xff);
                if (qoff < last_qoff || soff < last_soff) continue;
                if (qoff == last_qoff && soff == last_soff) continue;
                seeds.push_back(ka[s]);
                last_qoff = qoff;
                last_soff = soff;
            }
        }

        int sid = colin_align_list.front().sid;
        int sdir = colin_align_list.front().sdir;
        const char* subject = updb->target_sequence(sid, sdir);
        const int subject_size = updb->target_size(sid);
        int r = s_extend_one_mg_list_a(query, sid, sdir, subject, subject_size, seeds.data(), seeds.size(), tbck_data, align, true);
        if (!r) continue;
        tbck_data->cigar_to_align_string(query->fwd_rqs, subject);
        align.qas_offset = align_strings.size();
        align_strings.insert(align_strings.end(), tbck_data->qas, tbck_data->qae);
        align.sas_offset = align_strings.size();
        align_strings.insert(align_strings.end(), tbck_data->sas, tbck_data->sae);
        align.as_size = tbck_data->qae - tbck_data->qas;
        align.ddf_score = 0;
        for (auto& a : colin_align_list) align.ddf_score += a.ddf_score;
        align_list.push_back(align);
    }
    return;

    fasa = align_list.data();
    fasc = align_list.size();
    pdqsort(fasa, fasa + fasc, [](const PoreCInitHit& a, const PoreCInitHit& b) { return a.map_score > b.map_score; });
    for (int i = 0; i < fasc; ++i) {
        PoreCInitHit* fi = fasa + i;
        if (fi->sid == -1) continue;
        for (int j = i + 1; j < fasc; ++j) {
            PoreCInitHit* fj = fasa + j;
            if (fj->sid == -1 || fi->sid != fj->sid || fi->sdir != fj->sdir) continue;
            const int E = 10;
            int r = (fj->qoff + E >= fi->qoff) && (fj->qend <= fi->qend + E) 
                    && 
                    (fj->soff + E >= fi->soff) && (fj->send <= fi->send + E);
            if (r) fj->sid = -1;
        }
    }
    int n = 0;
    for (int i = 0; i < fasc; ++i) if (fasa[i].sid >= 0) fasa[n++] = fasa[i];
    fasc = n;
    align_list.resize(n);
}

BOOL
refine_pi_for_one_pca(PrelimSearch* prelim,
    HbnTracebackData* tbck_data,
    const u8* read,
    const char* raw_read,
    const char* subject,
    PoreCAlign* pca)
{
    //cerr << "refined " << *pca << '\n';
    HbnChainInfo hit;
    int km_off;
    pu64_t* km_a;
    int km_c;
    if (!prelim->R_find_init_hit(pca->qoff, pca->qend, subject, pca->soff, pca->send, hit, km_off, km_a, km_c)) return FALSE;
    //cerr << hit << '\n';

    if (!tbck_data->refined_map_one_chain_edlib(read, pca->qoff, pca->qend, subject, pca->soff, pca->send, km_a, km_c)) return FALSE;
    
    tbck_data->cigar_to_align_string(raw_read, subject);
    //tbck_data->dump(stderr);
    int as_size = tbck_data->qae - tbck_data->qas;
    validate_aligned_string_2(HBN_LOG_ARGS_DEFAULT, 0, read, pca->qoff, pca->qend, tbck_data->qas,
        0, subject, pca->soff, pca->send, tbck_data->sas, as_size, TRUE);
    pca->map_score = tbck_data->score;
    pca->pi = tbck_data->pi;
    pca->is_modified = false;
    pca->km_off = km_off;
    pca->km_c = km_c;
    return TRUE;
}

struct enzyme_offset {
    int qoff, soff;
    int enzyme_qoff, enzyme_soff;

    enzyme_offset() {
        qoff = -1;
        soff = -1;
        enzyme_qoff = -1;
        enzyme_soff = -1;
    }

    void dump(FILE* out, const char* sep = NULL) {
        fprintf(out, "[%d (%d), %d (%d)]", qoff, enzyme_qoff, soff, enzyme_soff);
        if (sep) fprintf(out, "%s", sep);
    }
};

static void
s_add_nearby_enzyme_pos(const int* reloci_array,
    const int intv_cnt,
    const int intv_i,
    const int offset,
    const int min_enzyme_offset,
    const int max_enzyme_offset,
    const int max_dist,
    vector<int>& offset_list)
{
    offset_list.clear();
    hbn_assert(offset >= reloci_array[intv_i] && offset <= reloci_array[intv_i + 1]);
    hbn_assert(offset >= min_enzyme_offset && offset <= max_enzyme_offset);

    int i = intv_i;
    while (i >= 0) {
        if (offset - reloci_array[i] > max_dist) break;
        if (reloci_array[i] < min_enzyme_offset) break;
        offset_list.push_back(reloci_array[i]);
        --i;
    }

    i = intv_i;
    while (i < intv_cnt) {
        if (reloci_array[i+1] - offset > max_dist) break;
        if (reloci_array[i+1] > max_enzyme_offset) break;
        offset_list.push_back(reloci_array[i+1]);
        ++i;
    }
    pdqsort(offset_list.begin(), offset_list.end());
}

static bool 
s_resolve_unfixed_start_offset(HbnTracebackData* A_tbck_data,
    int qoff,
    int qend,
    int soff,
    int send,
    const char* qas,
    const char* sas,
    const int as_size,
    const u8* query,
    const char* subject,
    enzyme_offset* eop)
{
    vector<u8>& qfrag = A_tbck_data->qfrag;
    vector<u8>& sfrag = A_tbck_data->sfrag;
    hbn_assert(eop->qoff == -1 || eop->soff == -1);
    if (eop->qoff == -1) {
        hbn_assert(eop->soff < send);
        if (eop->soff >= soff) {
            int qi = qoff;
            int si = soff;
            int ai = 0;
            while (si < eop->soff && ai < as_size) {
                if (qas[ai] != GAP_CHAR) ++qi;
                if (sas[ai] != GAP_CHAR) ++si;
                ++ai;
            }
            if (si != eop->soff) return false;
            hbn_assert(qi < qend);
            eop->qoff = qi;
            return true;
        }
        hbn_assert(eop->soff < soff);
        int sf = eop->soff;
        int st = soff + ALIGN_END_MATCH;
        hbn_assert(sf < st);
        int qt = qoff + ALIGN_END_MATCH;
        int qf = qt - (st - sf) - 50; qf = max(0, qf);
        hbn_assert(sf >= 0 && sf < st && st <= send);
        hbn_assert(qf >= 0 && qf < qt && qt <= qend);
        qfrag.clear();
        for (int i = qt - 1; i >= qf; --i) {
            qfrag.push_back(query[i]);
        }
        sfrag.clear();
        for (int i = st - 1; i >= sf; --i) {
            int c = subject[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        int qe = 0, se = 0;
        //HBN_LOG("1 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        run_shw(A_tbck_data->small_edlib, A_tbck_data->edlib, qfrag.data(), qfrag.size(),
            sfrag.data(), sfrag.size(), &qe, &se, nullptr, nullptr);
        if (sfrag.size() == se) {
            eop->qoff = qoff + ALIGN_END_MATCH - qe;
            hbn_assert(eop->qoff >= 0);
            return true;
        }
    }
    if (eop->soff == -1) {
        hbn_assert(eop->qoff < qend);
        if (eop->qoff >= qoff) {
            int qi = qoff;
            int si = soff;
            int ai = 0;
            while (qi < eop->qoff && ai < as_size) {
                if (qas[ai] != GAP_CHAR) ++qi;
                if (sas[ai] != GAP_CHAR) ++si;
                ++ai;
            }
            //fprintf(stderr, "====> qi = %d, si = %d, qoff = %d, soff = %d, eqoff = %d, as_size = %d\n", 
            //    qi, si, qoff, soff, eop->qoff, as_size);
            if (qi != eop->qoff) return false;
            hbn_assert(si < send, "si = %d, soff = %d, send = %d, qi = %d, qoff = %d, qend = %d", si, soff, send, qi, qoff, qend);
            eop->soff = si;
            return true;
        }

        int qf = eop->qoff;
        int qt = qoff + ALIGN_END_MATCH;
        hbn_assert(qf < qt);
        int st = soff + ALIGN_END_MATCH;
        int sf = st - (qt - qf) - 50; sf = max(0, sf);
        hbn_assert(sf >= 0 && sf < st && st <= send);
        hbn_assert(qf >= 0 && qf < qt && qt <= qend);
        qfrag.clear();
        for (int i = qt - 1; i >= qf; --i) qfrag.push_back(query[i]);
        sfrag.clear();
        for (int i = st - 1; i >= sf; --i) {
            int c = subject[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        int qe = 0, se = 0;
        //HBN_LOG("2 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        run_shw(A_tbck_data->small_edlib, A_tbck_data->edlib, qfrag.data(), qfrag.size(),
            sfrag.data(), sfrag.size(), &qe, &se, nullptr, nullptr);
        if (qfrag.size() == qe) {
            eop->soff = soff + ALIGN_END_MATCH - se;
            return true; 
        }     
    }
    return false;
}

static void
s_add_start_offset_pairs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    const int* q_reloci_array,
    const int q_reloci_cnt,
    const u8* query,
    const char* subject,
    const int sid,
    const int sdir,
    const int qoff,
    const int qend,
    const int soff,
    const int send,
    const char* qas,
    const char* sas,
    const int as_size,
    vector<enzyme_offset>& leop_list)
{
    const int vid = sid * 2 + sdir;
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[vid].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[vid].enzyme_loci_cnt;
    int l_q_intv_c = 0;
    int l_q_intv_i = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, qoff, &l_q_intv_c);  
    int l_s_intv_c = 0;
    int l_s_intv_i = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, soff, &l_s_intv_c);

    vector<int> ql_list;
    vector<int> sl_list;
    s_add_nearby_enzyme_pos(q_reloci_array, l_q_intv_c, l_q_intv_i, qoff, 0, qend - 1, kEnzymeNbhd, ql_list);
    s_add_nearby_enzyme_pos(s_reloci_array, l_s_intv_c, l_s_intv_i, soff, 0, send - 1, kEnzymeNbhd, sl_list);
    for (auto qi : ql_list) hbn_assert(qi < qend);
    for (auto si : sl_list) hbn_assert(si < send);

#if 0
    fprintf(stderr, "====> qoff = %d, soff = %d\n", qoff, soff);
    fprintf(stderr, "qoff:");for (auto qi : ql_list) fprintf(stderr, "\t%d", qi); fprintf(stderr, "\n");
    fprintf(stderr, "qoff:");for (auto si : sl_list) fprintf(stderr, "\t%d", si); fprintf(stderr, "\n");
#endif

    for (auto& qi : ql_list) {
        int q_d = qend - abs(qi);
        hbn_assert(q_d > 0);
        for (auto& si : sl_list) {
            int s_d = send - abs(si);
            hbn_assert(s_d > 0);
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            enzyme_offset eop;
            eop.qoff = abs(qi);
            eop.soff = abs(si);
            eop.enzyme_qoff = abs(qi);
            eop.enzyme_soff = abs(si);
            leop_list.push_back(eop);   
            qi = -abs(qi);
            si = -abs(si);        
        }
    }
    for (auto& qi : ql_list) {
        if (qi < 0) continue;
        enzyme_offset eop;
        eop.qoff = qi;
        eop.enzyme_qoff = qi;
        eop.soff = -1;
        eop.enzyme_soff = ((soff - s_reloci_array[l_s_intv_i]) <= (s_reloci_array[l_s_intv_i+1] - soff)) ? s_reloci_array[l_s_intv_i] : s_reloci_array[l_s_intv_i+1];
        //fprintf(stderr, "fix soff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_start_offset(tbck_data, qoff, qend, soff, send, qas, sas, as_size, query, subject, &eop)) {
            //fprintf(stderr, "fixed soff for\t"); eop.dump(stderr, "\n");
            if (eop.soff != eop.enzyme_soff) {
                int xxc = 0;
                int xxi = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, eop.soff, &xxc);
                int ld = eop.soff - s_reloci_array[xxi];
                hbn_assert(ld >= 0);
                int rd = s_reloci_array[xxi+1] - eop.soff;
                hbn_assert(rd >= 0);
                eop.enzyme_soff = (ld <= rd) ? s_reloci_array[xxi] : s_reloci_array[xxi+1];
            }
            leop_list.push_back(eop);
        }        
    }
    for (auto& si : sl_list) {
        if (si < 0) continue;
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qoff - q_reloci_array[l_q_intv_i]) <= (q_reloci_array[l_q_intv_i+1] - qoff)) ? q_reloci_array[l_q_intv_i] : q_reloci_array[l_q_intv_i+1];
        eop.soff = si;
        eop.enzyme_soff = si;
        //fprintf(stderr, "fix qoff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_start_offset(tbck_data, qoff, qend, soff, send, qas, sas, as_size, query, subject, &eop)) {
            //fprintf(stderr, "fixed qoff for\t"); eop.dump(stderr, "\n");
            if (eop.qoff != eop.enzyme_qoff) {
                int xxc = 0;
                int xxi = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, eop.qoff, &xxc);
                int ld = eop.qoff -  q_reloci_array[xxi];
                hbn_assert(ld >= 0);
                int rd = q_reloci_array[xxi+1] - eop.qoff;
                hbn_assert(rd >= 0);
                eop.enzyme_qoff = (ld <= rd) ? q_reloci_array[xxi] : q_reloci_array[xxi+1];
            }
            leop_list.push_back(eop);
        }          
    }

    if (leop_list.empty()) {
        enzyme_offset eop;
        eop.qoff = qoff;
        eop.enzyme_qoff = ((qoff - q_reloci_array[l_q_intv_i]) <= (q_reloci_array[l_q_intv_i+1] - qoff)) ? q_reloci_array[l_q_intv_i] : q_reloci_array[l_q_intv_i+1];
        eop.soff = soff;
        eop.enzyme_soff = ((soff - s_reloci_array[l_s_intv_i]) <= (s_reloci_array[l_s_intv_i+1] - soff)) ? s_reloci_array[l_s_intv_i] : s_reloci_array[l_s_intv_i+1];
        leop_list.push_back(eop);
    }

    int n = leop_list.size();
    for (int i = 0; i < n; ++i) {
        if (leop_list[i].qoff == -1) continue;
        int iq = leop_list[i].qoff;
        int is = leop_list[i].soff;
        for (int j = i + 1; j < n; ++j) {
            if (leop_list[j].qoff == -1) continue;
            int jq = leop_list[j].qoff;
            int js = leop_list[j].soff;
            if (jq == iq && js == is) leop_list[j].qoff = -1;
        }
    }
    int m = 0;
    for (int i = 0; i < n; ++i) if (leop_list[i].qoff>=0) leop_list[m++] = leop_list[i];
    leop_list.resize(m);
}

static bool 
s_resolve_unfixed_end_offset(HbnTracebackData* A_tbck_data,
    const u8* query,
    const int query_size,
    const char* subject,
    const int subject_size,
    int qoff,
    int qend,
    int soff,
    int send,
    const char* qas,
    const char* sas,
    const int as_size,
    enzyme_offset* eop)
{   
    vector<u8>& sfrag = A_tbck_data->sfrag;
    hbn_assert(eop->qoff == -1 || eop->soff == -1);
    if (eop->qoff == -1) {
        hbn_assert(eop->soff > soff);
        if (eop->soff <= send) {
            int qi = qend;
            int si = send;
            int ai = as_size;
            while (ai) {
                if (sas[ai-1] != GAP_CHAR && si == eop->soff) break;
                if (qas[ai-1] != GAP_CHAR) --qi;
                if (sas[ai-1] != GAP_CHAR) --si;
                --ai;
            }
            if (si != eop->soff) return false;
            hbn_assert(qi <= qend);
            eop->qoff = qi;
            return true;
        }
        int sf = send - ALIGN_END_MATCH;
        int st = eop->soff;
        hbn_assert(sf < st);
        int qf = qend - ALIGN_END_MATCH;
        int qt = qf + (st - sf) + 50; qt = min(qt, query_size);
        hbn_assert(sf >= soff && sf < st && st <= subject_size);
        hbn_assert(qf >= qoff && qf < qt && qt <= query_size);
        int qe = 0, se = 0;
        sfrag.clear();
        for (int i = sf; i < st; ++i) {
            int c = subject[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        //HBN_LOG("3 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        run_shw(A_tbck_data->small_edlib, A_tbck_data->edlib, query+qf, qt - qf, sfrag.data(), st - sf, &qe, &se, nullptr, nullptr);
        //if (se != (st - sf)) return 0;
        if (se == (st - sf)) {
            eop->qoff = qend - ALIGN_END_MATCH + qe;
            return true;
        }
    }
    if (eop->soff == -1) {
        hbn_assert(eop->qoff > qoff);
        if (eop->qoff <= qend) {
            int qi = qend;
            int si = send;
            int ai = as_size;
            while (ai) {
                if (qas[ai-1] != GAP_CHAR && qi == eop->qoff) break;
                if (qas[ai-1] != GAP_CHAR) --qi;
                if (sas[ai-1] != GAP_CHAR) --si;
                --ai;
            }
            if (qi != eop->qoff) return false;
            hbn_assert(si <= send);
            eop->soff = si;
            return true;            
        }
        int qf = qend - ALIGN_END_MATCH;
        int qt = eop->qoff;
        hbn_assert(qf < qt);
        int sf = send - ALIGN_END_MATCH;
        int st = sf + (qt - qf) + 50; st = min(st, subject_size);
        hbn_assert(sf >= soff && sf < st && st <= subject_size);
        hbn_assert(qf >= qoff && qf < qt && qt <= query_size);
        sfrag.clear();
        for (int i = sf; i < st; ++i) {
            int c = subject[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        int qe = 0, se = 0;
        //HBN_LOG("4 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        run_shw(A_tbck_data->small_edlib, A_tbck_data->edlib, query+qf, qt - qf, sfrag.data(), st - sf, &qe, &se, nullptr, nullptr);
        if (qe == (qt - qf)) {
            eop->soff = send - ALIGN_END_MATCH + se;
            return true;
        }
    }
    return false;
}

static void
s_add_end_offset_pairs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    const int* q_reloci_array,
    const int q_reloci_cnt,
    const u8* query,
    const int query_size,
    const char* subject,
    const int subject_size,
    const int sid,
    const int sdir,
    const int qoff,
    const int qend,
    const int soff,
    const int send,
    const char* qas,
    const char* sas,
    const int as_size,
    vector<enzyme_offset>& reop_list)
{
    const int vid = sid * 2 + sdir;
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[vid].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[vid].enzyme_loci_cnt;

    int r_q_intv_c = 0;
    int r_q_intv_i = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, qend, &r_q_intv_c);  
    int r_s_intv_c = 0;
    int r_s_intv_i = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, send, &r_s_intv_c);
    vector<int> qr_list;
    vector<int> sr_list;
    s_add_nearby_enzyme_pos(q_reloci_array, r_q_intv_c, r_q_intv_i, qend, qoff + 1, query_size, kEnzymeNbhd, qr_list);
    s_add_nearby_enzyme_pos(s_reloci_array, r_s_intv_c, r_s_intv_i, send, soff + 1, subject_size, kEnzymeNbhd, sr_list);

    for (auto qi : qr_list) hbn_assert(qi > qoff);
    for (auto si : sr_list) hbn_assert(si > soff);

#if 0
    fprintf(stderr, "====> qend = %d, send = %d\n", qend, send);
    fprintf(stderr, "qend: "); for (auto qi : qr_list) fprintf(stderr, "%d\t", qi); fprintf(stderr, "\n");
    fprintf(stderr, "send: "); for (auto si : sr_list) fprintf(stderr, "%d\t", si); fprintf(stderr, "\n");
#endif

    for (auto& qi : qr_list) {
        int q_d = abs(qi) - qoff;
        //HBN_LOG("qi = %d, q_d = %d, b_qoff = %d", qi, q_d, b_qoff);
        hbn_assert(q_d > 0);
        for (auto& si : sr_list) {
            int s_d = abs(si) - soff;
            //HBN_LOG("\tsi = %d, s_d = %d, b_soff = %d", si, s_d, b_soff);
            hbn_assert(s_d > 0);
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            enzyme_offset eop;
            eop.qoff = abs(qi);
            eop.enzyme_qoff = abs(qi);
            eop.soff = abs(si);
            eop.enzyme_soff = abs(si);
            reop_list.push_back(eop);
            qi = -abs(qi);
            si = -abs(si);
        }
    } 
    for (auto& qi : qr_list) {
        if (qi < 0) continue;
        enzyme_offset eop;
        eop.qoff = qi;
        eop.enzyme_qoff = qi;
        eop.soff = -1;
        eop.enzyme_soff = ((send - s_reloci_array[r_s_intv_i]) < (s_reloci_array[r_s_intv_i+1] - send)) ? s_reloci_array[r_s_intv_i] : s_reloci_array[r_s_intv_i+1];
        //fprintf(stderr, "resolve soff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, qoff, qend,
            soff, send, qas, sas, as_size, &eop)) {
            //fprintf(stderr, "resolved soff for\t"); eop.dump(stderr, "\n");
            if (eop.soff != eop.enzyme_soff) {
                int xxc = 0;
                int xxi = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, eop.soff, &xxc);
                int ld = eop.soff - s_reloci_array[xxi];
                hbn_assert(ld >= 0);
                int rd = s_reloci_array[xxi+1] - eop.soff;
                hbn_assert(rd >= 0);
                eop.enzyme_soff = (ld <= rd) ? s_reloci_array[xxi] : s_reloci_array[xxi+1];
            }
            reop_list.push_back(eop);
        }        
    }
    for (auto& si : sr_list) {
        if (si < 0) continue;
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qend - q_reloci_array[r_q_intv_i]) < (q_reloci_array[r_q_intv_i+1] - qend)) ? q_reloci_array[r_q_intv_i] : q_reloci_array[r_q_intv_i+1];
        eop.soff = si;
        eop.enzyme_soff = si;
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, qoff, qend,
            soff, send, qas, sas, as_size, &eop)) {
            if (eop.qoff != eop.enzyme_qoff) {
                int xxc = 0;
                int xxi = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, eop.qoff, &xxc);
                int ld = eop.qoff -  q_reloci_array[xxi];
                hbn_assert(ld >= 0);
                int rd = q_reloci_array[xxi+1] - eop.qoff;
                hbn_assert(rd >= 0);
                eop.enzyme_qoff = (ld <= rd) ? q_reloci_array[xxi] : q_reloci_array[xxi+1];
            }
            reop_list.push_back(eop);
        }        
    }

    if (reop_list.empty()) {
        enzyme_offset eop;
        eop.qoff = qend;
        eop.enzyme_qoff = ((qend - q_reloci_array[r_q_intv_i]) < (q_reloci_array[r_q_intv_i+1] - qend)) ? q_reloci_array[r_q_intv_i] : q_reloci_array[r_q_intv_i+1];
        eop.soff = send;
        eop.enzyme_soff = ((send - s_reloci_array[r_s_intv_i]) < (s_reloci_array[r_s_intv_i+1] - send)) ? s_reloci_array[r_s_intv_i] : s_reloci_array[r_s_intv_i+1];
        reop_list.push_back(eop);
    }

    int n = reop_list.size();
    for (int i = 0; i < n; ++i) {
        if (reop_list[i].qoff == -1) continue;
        int iq = reop_list[i].qoff;
        int is = reop_list[i].soff;
        for (int j = i + 1; j < n; ++j) {
            if (reop_list[j].qoff == -1) continue;
            int jq = reop_list[j].qoff;
            int js = reop_list[j].soff;
            if (jq == iq && js == is) reop_list[j].qoff = -1;
        }
    }
    int m = 0;
    for (int i = 0; i < n; ++i) if (reop_list[i].qoff>=0) reop_list[m++] = reop_list[i];
    reop_list.resize(m);
}

static void
s_fix_enzyme_align_ends(const HbnUnpackedDatabase* updb,
    PoreCQuery* query,
    PrelimSearch* prelim,
    HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    PoreCInitHit& align,
    const int align_idx,
    const int ddf_score,
    const int mapQ,
    const char* qas,
    const char* sas,
    const int as_size,
    std::vector<PoreCAlign>& pca_list)
{
    vector<enzyme_offset> leop_list;
    const char* subject = updb->target_sequence(align.sid, align.sdir);
    const int subject_size = updb->target_size(align.sid);
    s_add_start_offset_pairs(tbck_data, reloci_list, query->fwd_enzyme_pos, query->fwd_enzyme_pos_cnt, query->fwd_qs, 
        subject, align.sid, align.sdir, align.qoff, align.qend, align.soff, align.send, qas, sas, as_size, leop_list);

    vector<enzyme_offset> reop_list;
    s_add_end_offset_pairs(tbck_data, reloci_list, query->fwd_enzyme_pos, query->fwd_enzyme_pos_cnt, 
        query->fwd_qs, query->size, subject, subject_size,
        align.sid, align.sdir, align.qoff, align.qend, align.soff, align.send, qas, sas, as_size, reop_list);

    PoreCAlign pca;
    pca.qid = query->id;
    pca.qsize = query->size;
    pca.sid = align.sid;
    pca.sdir = align.sdir;
    pca.ssize = subject_size;
    pca.raw_aln_idx = align_idx;
    pca.ddf_score = ddf_score;
    pca.map_q = mapQ;
    for (auto& sx : leop_list) {
        pca.qoff = sx.qoff;
        pca.enzyme_qoff = sx.enzyme_qoff;
        pca.soff = sx.soff;
        pca.enzyme_soff = sx.enzyme_soff;
        for (auto& ex : reop_list) {
            pca.qend = ex.qoff;
            pca.enzyme_qend = ex.enzyme_qoff;
            pca.send = ex.soff;
            pca.enzyme_send = ex.enzyme_soff;
            if (pca.qend - pca.qoff < 20) continue;
            if (!refine_pi_for_one_pca(prelim, tbck_data, query->fwd_qs, query->fwd_rqs, subject, &pca)) continue;
            //if (!refine_pi_for_one_pca(tbck_data, query->fwd_qs, query->fwd_rqs, subject, &pca)) continue;
            if (!set_pca_chain_offset(&pca, reloci_list->enzyme.enzyme_size)) continue;

            pca_list.push_back(pca);
        }
    }
}

void fix_enzyme_align_ends(HbnUnpackedDatabase* updb,
    RestrictEnzymeLociList* reloci_list,
    PoreCQuery* query,
    std::vector<PoreCInitHit>& align_list,
    std::string& align_strings,
    PrelimSearch* prelim,
    HbnTracebackData* tbck_data,
    std::vector<PoreCAlign>& pca_list)
{
    prelim->R_init(query->fwd_qs, query->size);
    int na = align_list.size();
    for (int i = 0; i < na; ++i) {
        //cerr << "fix enzyme offsets for " << i << '\t' << align_list[i] << '\n';
        PoreCInitHit& align = align_list[i];
	    //cerr << i << '\t' << align << '\n';
        //if (align.id != align.parent) continue;
        const char* qas = align_strings.data() + align.qas_offset;
        const char* sas = align_strings.data() + align.sas_offset;
        const int as_size = align.as_size;
        s_fix_enzyme_align_ends(updb, query, prelim, tbck_data, reloci_list, align, i, align.ddf_score,
            align.mapQ, qas, sas, as_size, pca_list);
    }
}

////////////////////

void PoreCGapAlign::align_ddf_chains(PrelimSearch* prelim,
    PoreCQuery* query,
    HbnUnpackedDatabase* updb,
    RestrictEnzymeLociList* enzyme,
    HbnChainInfo* ac_a,
    int ac_n,
    pu64_t* ac_km)
{
    M_align_list.clear();
    M_align_strings.clear();
    extend_ddf_chain_list(M_tbck_data, updb, query, ac_a, ac_n, ac_km, M_cov_list, M_align_list, M_align_strings);
    set_mapQ_for_init_hits(M_align_list.data(), M_align_list.size());

    M_pca_list.clear();
    fix_enzyme_align_ends(updb, enzyme, query, M_align_list, M_align_strings, prelim, M_tbck_data, M_pca_list);
}
