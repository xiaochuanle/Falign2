#include "trim_overlap_subseq.hpp"

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "../../corelib/pdqsort.h"
#include "../../sw/hbn_traceback_aux.h"

using namespace std;

static const int kMinFragSize = 50;

void
get_pca_align_string(HbnTracebackData* tbck_data,
    PrelimSearch* prelim_search,
    RestrictEnzymeLociList* reloci_list,
    HbnUnpackedDatabase* updb,
    PoreCQuery* query,
    PoreCAlign* pca,
    TrimPcaList* tpca_list)
{
    HbnChainInfo hit;
    int km_off;
    pu64_t* km_a;
    int km_c;
    const char* subject = updb->target_sequence(pca->sid, pca->sdir);
    if (!pca->is_modified) {
        km_a = prelim_search->R_get_chain_seeds(pca->km_off);
        km_c = pca->km_c;
    } else {
        if (!prelim_search->R_find_init_hit(pca->qoff, pca->qend, subject, pca->soff, pca->send, hit, km_off, km_a, km_c)) return;
    }

    if (!tbck_data->refined_map_one_chain(query->fwd_qs, pca->qoff, pca->qend, subject, pca->soff, pca->send, km_a, km_c)) return;
    
    tbck_data->cigar_to_align_string(query->fwd_rqs, subject);
    //tbck_data->dump(stderr);
    int as_size = tbck_data->qae - tbck_data->qas;
    validate_aligned_string_2(HBN_LOG_ARGS_DEFAULT, 0, query->fwd_qs, pca->qoff, pca->qend, tbck_data->qas,
        0, subject, pca->soff, pca->send, tbck_data->sas, as_size, TRUE);
    pca->map_score = tbck_data->score;
    pca->pi = tbck_data->pi;
    tpca_list->add(pca, tbck_data->qas, tbck_data->sas, as_size);
}

static void
s_trim_qend(RestrictEnzymeLociList* gfa_enzyme, TrimPcaList* tpca_list, TrimPca* tpca, int trim_qend_to)
{
    hbn_assert(tpca->pca.qend > trim_qend_to);
    if (tpca->pca.qoff >= trim_qend_to) {
        tpca->is_valid = 0;
        return;
    }
    int qend = tpca->pca.qend;
    int vtrim = 0;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_size = tpca->as_size;
    int ai = as_size;
    while (ai) {
        if (qas[ai-1] != GAP_CHAR && qend == trim_qend_to) break;
        if (qas[ai-1] != GAP_CHAR) --qend;
        if (sas[ai-1] != GAP_CHAR) ++vtrim;
        --ai;
    }
    hbn_assert(qend == trim_qend_to);
    tpca->pca.qend = qend;
    tpca->pca.send -= vtrim;
    tpca->as_size = ai;

    int vid = tpca->pca.sid * 2 + tpca->pca.sdir;
    const int* reloci_array = gfa_enzyme->reloci_array + gfa_enzyme->seq_reloci_info_array[vid].enzyme_loci_offset;
    const int reloci_cnt = gfa_enzyme->seq_reloci_info_array[vid].enzyme_loci_cnt;
    tpca->pca.enzyme_send = get_enzyme_pos(reloci_array, reloci_cnt, tpca->pca.send);
}

static void
s_trim_qoff(RestrictEnzymeLociList* gfa_enzyme, TrimPcaList* tpca_list, TrimPca* tpca, int trim_qoff_to)
{
    hbn_assert(tpca->pca.qoff < trim_qoff_to);
    if (tpca->pca.qend <= trim_qoff_to) {
        tpca->is_valid = 0;
        return;
    }

    int qoff = tpca->pca.qoff;
    int vtrim = 0;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int ai = 0;   
    while (ai < tpca->as_size) {
        if (qas[ai] != GAP_CHAR) ++qoff;
        if (sas[ai] != GAP_CHAR) ++vtrim;
        ++ai;
        if (qoff == trim_qoff_to) break;
    }    
    hbn_assert(qoff == trim_qoff_to);
    hbn_assert(ai <= tpca->as_size);
    tpca->pca.qoff = qoff;
    tpca->qas_offset += ai;
    tpca->sas_offset += ai;
    tpca->as_size -= ai;

    tpca->pca.qoff = qoff;
    tpca->pca.soff += vtrim;
    int vid = tpca->pca.sid * 2 + tpca->pca.sdir;
    const int* reloci_array = gfa_enzyme->reloci_array + gfa_enzyme->seq_reloci_info_array[vid].enzyme_loci_offset;
    const int reloci_cnt = gfa_enzyme->seq_reloci_info_array[vid].enzyme_loci_cnt;
    tpca->pca.enzyme_soff = get_enzyme_pos(reloci_array, reloci_cnt, tpca->pca.soff);
}

static void
s_trim_send(TrimPcaList* tpca_list, TrimPca* tpca, int trim_send_to)
{
    hbn_assert(tpca->pca.send > trim_send_to);
    if (tpca->pca.soff >= trim_send_to) {
        tpca->is_valid = 0;
        return;
    }
    int qend = tpca->pca.qend;
    int send = tpca->pca.send;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_i = tpca->as_size;
    while (as_i) {
        if (sas[as_i-1] != GAP_CHAR && send == trim_send_to) break;
        if (qas[as_i-1] != GAP_CHAR) --qend;
        if (sas[as_i-1] != GAP_CHAR) --send;
        --as_i;
    }
    hbn_assert(send == trim_send_to, "qend = %d, send = %d, ai = %d, as = %d", qend, send, as_i, tpca->as_size);
    tpca->pca.qend = qend;
    tpca->pca.send = send;
    tpca->as_size = as_i;
}

static void
s_trim_soff(TrimPcaList* tpca_list, TrimPca* tpca, int trim_soff_to)
{
    hbn_assert(tpca->pca.soff < trim_soff_to);
    if (tpca->pca.send <= trim_soff_to) {
        tpca->is_valid = 0;
        return;
    }
    int qoff = tpca->pca.qoff;
    int soff = tpca->pca.soff;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_i = 0;
    while (as_i < tpca->as_size) {
        if (qas[as_i] != GAP_CHAR) ++qoff;
        if (sas[as_i] != GAP_CHAR) ++soff;
        ++as_i;
        if (soff == trim_soff_to) break;
    }    
    if (as_i >= tpca->as_size || soff != trim_soff_to) return;
    tpca->pca.qoff = qoff;
    tpca->pca.soff = soff;
    tpca->qas_offset += as_i;
    tpca->sas_offset += as_i;
    tpca->as_size -= as_i;
}

static void
s_trim_query_overlap_subseq(RestrictEnzymeLociList* gfa_enzyme, TrimPcaList* tpca_list)
{
    int n_tpca = tpca_list->tpca_list.size();
    for (int i = 0; i < n_tpca - 1; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        TrimPca* pj = &tpca_list->tpca_list[i+1];
        if (pi->is_valid == 0 || pj->is_valid == 0) continue;
        if (pi->pca.qend <= pj->pca.qoff) continue;
	    if ((*pi) < (*pj)) {
		    pi->is_valid = 0;
		    break;
	    }
	    if ((*pj) < (*pi)) {
		    pj->is_valid = 0;
		    continue;
        }

        //cerr << "trim qoff for" << '\n' << pi->pca << '\n' << pj->pca << '\n';

        int i_q_e = (pi->pca.qend == pi->pca.enzyme_qend);
        int i_s_e = (pi->pca.send == pi->pca.enzyme_send);
        int j_q_s = (pj->pca.qoff == pj->pca.enzyme_qoff);
        int j_s_s = (pj->pca.soff == pj->pca.enzyme_soff);

        if (i_q_e && i_s_e && j_q_s && j_s_s) {
            if (pi->pca.qend - pi->pca.qoff > pj->pca.qend - pj->pca.qoff) {
                s_trim_qoff(gfa_enzyme, tpca_list, pj, pi->pca.qend);
            } else {
                s_trim_qend(gfa_enzyme, tpca_list, pi, pj->pca.qoff);
            }
            continue;
        }

        if (i_q_e && i_s_e) {
            s_trim_qoff(gfa_enzyme, tpca_list, pj, pi->pca.qend);
            continue;
        }

        if (j_q_s && j_s_s) {
            s_trim_qend(gfa_enzyme, tpca_list, pi, pj->pca.qoff);
            continue;
        }

        if ((i_q_e || i_s_e) && (j_q_s == 0 && j_s_s == 0)) {
            s_trim_qoff(gfa_enzyme, tpca_list, pj, pi->pca.qend);
            continue;
        }

        if ((i_q_e == 0 && i_s_e == 0) && (j_q_s || j_s_s)) {
            s_trim_qend(gfa_enzyme, tpca_list, pi, pj->pca.qoff);
            continue;
        }

        if (pi->pca.qend - pi->pca.qoff > pj->pca.qend - pj->pca.qoff) {
            s_trim_qoff(gfa_enzyme, tpca_list, pj, pi->pca.qend);
        } else {
            s_trim_qend(gfa_enzyme, tpca_list, pi, pj->pca.qoff);
        }
    }

    for (int i = 0; i < n_tpca; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        if (!pi->is_valid) continue;
        hbn_assert(pi->pca.qoff < pi->pca.qend);
        if (pi->pca.qend - pi->pca.qoff < kMinFragSize) pi->is_valid = 0;
    }
}

static void
s_trim_subject_overlap_subseq(TrimPcaList* tpca_list)
{
    const int n_tpca = tpca_list->tpca_list.size();
    for (int i = 0; i < n_tpca; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        if (!pi->is_valid) continue;
        for (int j = i + 1; j < n_tpca; ++j) {
            TrimPca* pj = &tpca_list->tpca_list[j];
            if (!pj->is_valid) continue;

            TrimPca* pl = nullptr;
            TrimPca* pr = nullptr;
            if (pi->pca.sid == pj->pca.sid && pi->pca.sdir == pj->pca.sdir 
                && pi->pca.soff < pj->pca.soff && pi->pca.send > pj->pca.soff) {
                pl = pi;
                pr = pj;
            }   
            if (pj->pca.sid == pi->pca.sid && pj->pca.sdir == pi->pca.sdir 
                && pj->pca.soff < pi->pca.soff && pj->pca.send > pi->pca.soff) {
                pl = pj;
                pr = pi;
            }
            if (pl == nullptr || pr == nullptr) continue;

            //cerr << "trim subject offset for" << '\n' << pl->pca << '\n' << pr->pca << '\n';

            int l_q_e = (pl->pca.qend == pl->pca.enzyme_qend);
            int l_s_e = (pl->pca.send == pl->pca.enzyme_send);
            int r_q_s = (pr->pca.qoff == pr->pca.enzyme_qoff);
            int r_s_s = (pr->pca.soff == pr->pca.enzyme_soff);

            if (l_q_e && l_s_e && r_q_s && r_s_s) {
                if (pl->pca.send - pl->pca.soff < pr->pca.send - pr->pca.soff) {
                    s_trim_send(tpca_list, pl, pr->pca.soff);
                } else {
                    s_trim_soff(tpca_list, pr, pl->pca.send);
                }
                continue;
            }

            if (l_q_e && l_s_e) {
                s_trim_soff(tpca_list, pr, pl->pca.send);
                continue;
            }

            if (r_q_s && r_s_s) {
                s_trim_send(tpca_list, pl, pr->pca.soff);
                continue;
            }

            if ((l_q_e || l_s_e) && (r_q_s == 0 && r_s_s)) {
                s_trim_soff(tpca_list, pr, pl->pca.send);
                continue;
            }

            if ((l_q_e == 0 && l_s_e == 0) && (r_q_s || r_s_s)) {
                s_trim_send(tpca_list, pl, pr->pca.soff);
                continue;
            }

            if (pl->pca.send - pl->pca.soff > pr->pca.send - pr->pca.soff) {
                s_trim_soff(tpca_list, pr, pl->pca.send);
            } else {
                s_trim_send(tpca_list, pl, pr->pca.soff);
            }
        }
    }    

    for (int i = 0; i < n_tpca; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        if (!pi->is_valid) continue;
        if (pi->pca.qend - pi->pca.qoff < kMinFragSize) pi->is_valid = 0;
    }
}

bool s_is_complete_chain(const int* vdfa, const int vdfc, const PoreCAlign* pca_a, const int pca_c);

void
set_trim_pca_list_for_fwd_query(TrimPcaList* src, TrimPcaList* dst, const int enzyme_size)
{
    dst->align_strings.clear();
    dst->tpca_list.clear();
    TrimPca dpca;
    for (auto& tpca : src->tpca_list) {
        if (!tpca.is_valid) continue;
        {
            dpca = tpca;
            const char* qas = src->align_strings.c_str() + tpca.qas_offset;
            const char* sas = src->align_strings.c_str() + tpca.sas_offset;
            const int as_size = tpca.as_size;
            dpca.qas_offset = dst->align_strings.size();
            dst->align_strings.append(qas, as_size);
            dpca.sas_offset = dst->align_strings.size();
            dst->align_strings.append(sas, as_size);
            dpca.as_size = as_size;
            dst->tpca_list.push_back(dpca);
            continue;
        }
    }
}

void
set_trim_pca_list_for_fwd_subject(TrimPcaList* src, TrimPcaList* dst, const int enzyme_size)
{
    dst->align_strings.clear();
    dst->tpca_list.clear();
    TrimPca dpca;
    for (auto& tpca : src->tpca_list) {
        if (!tpca.is_valid) {
            //HBN_LOG("invalid pca");
            //dump_pca(fprintf, stderr, tpca.pca, -1);
            continue;
        }
        {
            dpca = tpca;
            const char* qas = src->align_strings.c_str() + tpca.qas_offset;
            const char* sas = src->align_strings.c_str() + tpca.sas_offset;
            const int as_size = tpca.as_size;
            dpca.qas_offset = dst->align_strings.size();
            dst->align_strings.append(qas, as_size);
            dpca.sas_offset = dst->align_strings.size();
            dst->align_strings.append(sas, as_size);
            dpca.as_size = as_size;
            dst->tpca_list.push_back(dpca);
            continue;
        }
    }
}

static void
s_normalize_pca_offsets(TrimPcaList* tpca_list,
    PoreCQuery* query,
    HbnUnpackedDatabase* updb)
{
    for (auto& tpca : tpca_list->tpca_list) {
        if (tpca.pca.sdir == FWD) continue;
        PoreCAlign& pca = tpca.pca;
        //cerr << "normalize\t" << pca << '\n';

        {
            const char* qas = tpca_list->align_strings.c_str() + tpca.qas_offset;
            const char* sas = tpca_list->align_strings.c_str() + tpca.sas_offset;
            const int as_size = tpca.as_size;
            const char* ss = updb->target_sequence(pca.sid, pca.sdir);
            //dump_align_string(qas, sas, as_size, stderr);
            validate_aligned_string_2(__FILE__, __FUNCTION__, __LINE__,
                pca.qid, query->fwd_qs, pca.qoff, pca.qend, qas,
                pca.sid, ss, pca.soff, pca.send, sas, as_size, 1);
        }
        
        int x = pca.qsize - pca.qend;
        int y = pca.qsize - pca.qoff;
        pca.qoff = x;
        pca.qend = y;

        x = pca.qsize - pca.enzyme_qend;
        y = pca.qsize - pca.enzyme_qoff;
        pca.enzyme_qoff = x;
        pca.enzyme_qend = y;

        x = pca.ssize - pca.send;
        y = pca.ssize - pca.soff;
        pca.soff = x;
        pca.send = y;

        x = pca.ssize - pca.enzyme_soff;
        y = pca.ssize - pca.enzyme_send;
        pca.enzyme_soff = x;
        pca.enzyme_send = y;

        size_t qa = tpca.qas_offset;
        size_t sa = tpca.sas_offset;
        const int as_size = tpca.as_size;
        reverse(&tpca_list->align_strings[qa], &tpca_list->align_strings[qa+as_size]);
        for (int i = 0; i < as_size; ++i) {
            char& c = tpca_list->align_strings[qa+i];
            if (c != GAP_CHAR) c = complement_base_table[(int)c];
        }
        reverse(&tpca_list->align_strings[sa], &tpca_list->align_strings[sa+as_size]);
        for (int i = 0; i < as_size; ++i) {
            char& c = tpca_list->align_strings[sa+i];
            if (c != GAP_CHAR) c = complement_base_table[(int)c];
        }

        const u8* qs = (pca.sdir) ? query->rev_qs : query->fwd_qs;
        const char* ss = updb->target_sequence(pca.sid, FWD);
        const char* qas = tpca_list->align_strings.c_str() + qa;
        const char* sas = tpca_list->align_strings.c_str() + sa;
        //dump_align_string(qas, sas, as_size, stderr);
        validate_aligned_string_2(__FILE__, __FUNCTION__, __LINE__,
            query->id, qs, pca.qoff, pca.qend, qas,
            pca.sid, ss, pca.soff, pca.send, sas, as_size, 1);
    }
}

void
trim_overlap_subseqs(HbnTracebackData* tbck_data,
    PrelimSearch* prelim_search,
    HbnUnpackedDatabase* updb,
    RestrictEnzymeLociList* reloci_list,
    PoreCQuery* query,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    PoreCAlign* pca_a,
    int pca_c,
    const EChainType chain_type,
    TrimPcaList& trim_pca_list)
{
    if (1)
    {
        if (1) {
            for (int i = 0; i < pca_c; ++i) {
                if (pca_a[i].qend - pca_a[i].qoff >= 200) continue;
                const int E = 20;
                int ld = abs(pca_a[i].soff - pca_a[i].enzyme_soff) <= E;
                int rd = abs(pca_a[i].send - pca_a[i].enzyme_send) <= E;
		        int lsd = abs(pca_a[i].qoff - pca_a[i].enzyme_qoff) <= E;
		        int rsd = abs(pca_a[i].qend - pca_a[i].enzyme_qend) <= E;
                if (ld || rd || lsd || rsd) continue;
                //fprintf(stderr, "garbage pca\n");
                //dump_chain_pca(fprintf, stderr, pca_a[i], i);
                pca_a[i].qid = -1;
            }
        } else if (1) {
		    for (int i = 0; i < pca_c; ++i) {
                	    if (pca_a[i].qend - pca_a[i].qoff >= 100) continue;
			    int lsd = abs(pca_a[i].soff - pca_a[i].enzyme_soff) <= 20;	
			    int rsd = abs(pca_a[i].send - pca_a[i].enzyme_send) <= 20;
			    int lqd = 0;//abs(pca_a[i].qoff - pca_a[i].enzyme_qoff) <= 20;	
			    int rqd = 0;//abs(pca_a[i].qend - pca_a[i].enzyme_qend) <= 20;
			    if (lsd == 0 && rsd == 0 && lqd == 0 && rqd == 0) {
				    //fprintf(stderr, "garbage complte pca\n");
				    //dump_chain_pca(fprintf, stderr, pca_a[i], i);
				    pca_a[i].qid = -1;
			    }
		    }
	    }
        int n = 0;
        for (int i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
        pca_c = n;
    }

    TrimPcaList all_tpca_list;
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        get_pca_align_string(tbck_data, prelim_search, reloci_list, updb, query, pca, &all_tpca_list);
    }

    TrimPcaList qry_tpca_list;
    set_trim_pca_list_for_fwd_query(&all_tpca_list, &qry_tpca_list, reloci_list->enzyme.enzyme_size);
    pdqsort(qry_tpca_list.tpca_list.begin(),
        qry_tpca_list.tpca_list.end(),
        [](const TrimPca& a, const TrimPca& b)->bool { return a.pca.qoff < b.pca.qoff; });
    s_trim_query_overlap_subseq(reloci_list, &qry_tpca_list);

    TrimPcaList sbj_tpca_list;
    set_trim_pca_list_for_fwd_subject(&qry_tpca_list, &sbj_tpca_list, reloci_list->enzyme.enzyme_size);
    int n_pca = sbj_tpca_list.tpca_list.size();
    s_trim_subject_overlap_subseq(&sbj_tpca_list);
    n_pca = sbj_tpca_list.tpca_list.size();

    trim_pca_list.clear();
    trim_pca_list.align_strings = sbj_tpca_list.align_strings;
    for (int i = 0; i < n_pca; ++i) {
        TrimPca* tpca = &sbj_tpca_list.tpca_list[i];
        if (!tpca->is_valid) continue;
        trim_pca_list.tpca_list.push_back(*tpca);
    }

    s_normalize_pca_offsets(&trim_pca_list, query, updb);
}
