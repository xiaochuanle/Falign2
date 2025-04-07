#include "smooth_pca_list.hpp"

#include "../../corelib/pdqsort.h"

#include <algorithm>
#include <iostream>
#include <cmath>

using namespace std;

static void
s_merge_two_pca(PoreCAlign* L, PoreCAlign* R, PoreCAlign* M)
{
    hbn_assert(L->sid == R->sid && L->sdir == R->sdir);

    M->qid = L->qid;
    M->qoff = L->qoff;
    M->enzyme_qoff = L->enzyme_qoff;
    M->qend = R->qend;
    M->enzyme_qend = R->enzyme_qend;
    M->qsize = L->qsize;

    M->sid = L->sid;
    M->sdir = L->sdir;
    M->soff = L->soff;
    M->enzyme_soff = L->enzyme_soff;
    M->send = R->send;
    M->enzyme_send = R->enzyme_send;
    M->ssize = L->ssize;

    M->ddf_score = L->ddf_score + R->ddf_score;
    M->map_score = L->map_score + R->map_score;
    M->pi = (L->pi + R->pi) / 2;

    M->is_modified = true;
    M->km_off = 0;
    M->km_c = 0;
}

static bool 
s_gf_p2_left_extend(PoreCAlign* L,
    PoreCAlign* R,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    if (L->qend != L->enzyme_qend) return false;
    int lqe = L->qend;
    int rqe = R->qend;

    int qf = lqe;
    int qt = rqe;
    int sf = R->enzyme_soff;
    int st = R->send;

    //HBN_LOG("le qf = %d, qt = %d, sf = %d, st = %d\n", qf, qt, sf, st);
    int q_d = qt - qf;
    int s_d = st - sf;
    if (s_d <= 0) return false;
    double ddf = fabs(1.0 - 1.0 * q_d / s_d);
    if (ddf > 0.2) return false;

    const u8* read = fwd_read + qf;
    int rl = qt - qf;
    const char* sbj = updb->target_sequence(R->sid, R->sdir) + sf;
    int sl = st - sf;
    vector<u8>& sfrag = tbck_data->sfrag;
    sfrag.clear();
    for (int i = 0; i < sl; ++i) {
        int c = sbj[i];
        c = nst_nt4_table[c];
        if (c > 3) c = 0;
        sfrag.push_back(c);
    }
    run_nw(read, rl, sfrag.data(), sl, tbck_data->small_edlib, tbck_data->edlib, 
        tbck_data->ext_qabuf, tbck_data->ext_sabuf);

    const char* qas = tbck_data->ext_qabuf.c_str(); 
    const char* sas = tbck_data->ext_sabuf.c_str(); 
    int as_size = tbck_data->ext_qabuf.size();
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    //fprintf(stderr, "pi = %g\n", ident);
    if (ident < R->pi - 5.0 && ident < 75.0) return false;

    R->qoff = lqe;
    R->enzyme_qoff = lqe;
    R->soff = sf;
    set_pca_chain_offset(R, enzyme_size);
    R->is_modified = true;
    R->km_off = 0;
    R->km_c = 0;

    return true;
}

static bool 
s_gf_p2_right_extend(PoreCAlign* L,
    PoreCAlign* R,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    if (R->qoff != R->enzyme_qoff) return false;
    int lqb = L->qoff;
    int rqb = R->qoff;

    int qf = lqb;
    int qt = rqb;
    int sf = L->soff;
    int st = L->enzyme_send;
    //HBN_LOG("re qf = %d, qt = %d, sf = %d, st = %d\n", qf, qt, sf, st);
    int q_d = qt - qf;
    int s_d = st - sf;
    if (s_d <= 0) return false;
    double ddf = fabs(1.0 - 1.0 * q_d / s_d);
    //HBN_LOG("q_d = %d, s_d = %d, ddf = %g", q_d, s_d, ddf);
    if (ddf > 0.2) return false;

    const u8* read = fwd_read + qf;
    int rl = qt - qf;
    const char* sbj = updb->target_sequence(L->sid, L->sdir) + sf;
    int sl = st - sf;
    vector<u8>& sfrag = tbck_data->sfrag;
    sfrag.clear();
    for (int i = 0; i < sl; ++i) {
        int c = sbj[i];
        c = nst_nt4_table[c];
        if (c > 3) c = 0;
        sfrag.push_back(c);
    }
    run_nw(read, rl, sfrag.data(), sl, tbck_data->small_edlib, tbck_data->edlib, 
        tbck_data->ext_qabuf, tbck_data->ext_sabuf);

    const char* qas = tbck_data->ext_qabuf.c_str(); 
    const char* sas = tbck_data->ext_sabuf.c_str(); 
    int as_size = tbck_data->ext_qabuf.size();
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    //fprintf(stderr, "pi: %g ---> %f\n", L->pi, ident);
    if (ident < L->pi - 5.0 && ident < 75.0) return false;

    L->qend = qt;
    L->enzyme_qend = qt;
    L->send = st;
    set_pca_chain_offset(L, enzyme_size);
    L->is_modified = true;
    L->km_off = 0;
    L->km_c = 0;

    return true;
}

static bool 
s_gf_p2_left_most_extend(PoreCAlign* pca,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    vector<u8>& qfrag = tbck_data->qfrag;
    vector<u8>& sfrag = tbck_data->sfrag;
    const u8* read = fwd_read;
    const char* subject = updb->target_sequence(pca->sid, pca->sdir);

    {
        int qblk = min(pca->qoff, pca->soff + 50);
        int sblk = min(pca->qoff + 50, pca->soff);
        qfrag.clear(); for (int i = 1; i <= qblk; ++i) qfrag.push_back(read[pca->qoff-i]);
        sfrag.clear(); 
        for (int i = 1; i <= sblk; ++i) {
            int c = subject[pca->soff - i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        //fprintf(stderr, "lme 1 qf = %d, sf = %d, qblk = %d, sblk = %d\n", pca->qoff, pca->soff, qblk, sblk);
        //if (qblk > 5000 || sblk > 5000) fprintf(stderr, "lme 1 qf = %d, sf = %d, qblk = %d, sblk = %d\n", pca->qoff, pca->soff, qblk, sblk);
        if (qblk > 300 || sblk > 300) return false;
        int qfae = 0, sfae = 0;
        run_shw(tbck_data->small_edlib, tbck_data->edlib, qfrag.data(), qblk, sfrag.data(), sblk, &qfae, &sfae, 0, 0);
        //fprintf(stderr, "**** qfae = %d, sfae = %d\n", qfae, sfae);
        //HBN_LOG("LM fix "); cerr << *pca << '\n';
        pca->qoff -= qfae;
        pca->soff -= sfae;
        if (pca->qoff == 0) pca->enzyme_qoff = pca->qoff;
        //fprintf(stderr, "into\t"); cerr << *pca << '\n';
        pca->is_modified = true;
        pca->km_off = 0;
        pca->km_c = 0;
	    return true;
    } 
}

static bool 
s_gf_p2_right_most_extend(PoreCAlign* pca,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    vector<u8>& sfrag = tbck_data->sfrag;
    const u8* read = fwd_read;
    const char* subject = updb->target_sequence(pca->sid, pca->sdir);

    {
        int qt = pca->qend;
        int st = pca->send;
        int dq = pca->qsize - qt;
        int ds = pca->ssize - st;
        int qblk = min(ds + 50, dq);
        int sblk = min(dq + 50, ds);
	    if (qblk == 0 || sblk == 0) return false;
        //HBN_LOG("rme1 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, pca->qsize, pca->rv_vsize);
        //if (qblk > 5000 || sblk > 5000) HBN_LOG("rme1 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, pca->qsize, pca->ssize);
        if (qblk > 300 || sblk > 300) return false;
	    if (qblk == 0 || sblk == 0) return false;
        //HBN_LOG("RM fix"); cerr << *pca << '\n';
        sfrag.clear();
        for (int i = 0; i < sblk; ++i) {
            int c = subject[st+i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        int qfae = 0, sfae = 0;
        run_shw(tbck_data->small_edlib, tbck_data->edlib, read+qt, qblk, sfrag.data(), sblk, &qfae, &sfae, 0, 0);
        //HBN_LOG("qfae = %d, sfae = %d, qt = %d, st = %d", qfae, sfae, qt, st);
        qt += qfae;
        st += sfae;
        hbn_assert(qt <= pca->qsize);
        hbn_assert(st <= pca->ssize);
        pca->qend = qt;
        pca->qend = qt;
        pca->send = st;
        if (pca->qend == pca->qsize) pca->enzyme_qend = pca->qend;
        //fprintf(stderr, "into\t"); cerr << *pca << '\n';
        pca->is_modified = true;
        pca->km_off = 0;
        pca->km_c = 0;
	    return true;
    }
}

static bool 
s_gf_p2_left_extend_gap(PoreCAlign* L,
    PoreCAlign* R,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    int sf = R->enzyme_soff;
    int st = R->soff;
    int qt = R->qoff;
    int qf = qt - (st - sf + 50);
    qf = max(0, qf);

    //fprintf(stderr, "[%s] [%d, %d] x [%d, %d]\n", __FUNCTION__, qf, qt, sf, st);

    vector<u8>& qfrag = tbck_data->qfrag;
    qfrag.clear();
    for (int i = qt - 1; i >= qf; --i) {
        qfrag.push_back(fwd_read[i]);
        int qc = fwd_read[i];
        qc = DECODE_RESIDUE(qc);
        //fprintf(stderr, "%c", qc);
    }
    //fprintf(stderr, "\n");

    vector<u8>& sfrag = tbck_data->sfrag;
    sfrag.clear();
    const char* sbj = updb->target_sequence(R->sid, R->sdir);
    for (int i = st - 1; i >= sf; --i) {
        int c = sbj[i];
        //fprintf(stderr, "%c", c);
        c = nst_nt4_table[c];
        if (c > 3) c = 0;
        sfrag.push_back(c);
    }
    //fprintf(stderr, "\n");

    int qfae = 0, sfae = 0;
    int ql = qfrag.size(), sl = sfrag.size();
    run_shw(tbck_data->small_edlib, tbck_data->edlib, qfrag.data(), ql, sfrag.data(), sl, &qfae, &sfae, 
        &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
    //fprintf(stderr, "qfae = %d, sfae = %d, ql = %d, sl = %d\n", qfae, sfae, ql, sl);

    const char* qas = tbck_data->ext_qabuf.c_str(); 
    const char* sas = tbck_data->ext_sabuf.c_str(); 
    int as_size = tbck_data->ext_qabuf.size();
    //dump_align_string(qas, sas, as_size, stderr);
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    //fprintf(stderr, "pi = %g\n", ident);
    if (ident < 60.0) return false;

    R->qoff -= qfae;
    R->soff -= sfae;
    set_pca_chain_offset(R, enzyme_size);
    R->is_modified = true;
    R->km_off = 0;
    R->km_c = 0;

    return true;
}

static bool 
s_gf_p2_right_extend_gap(PoreCAlign* L,
    PoreCAlign* R,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int read_size,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    int sf = L->send;
    int st = L->enzyme_send;
    int qf = L->qend;
    int qt = qf + (st - sf) + 50;
    qt = min(qt, read_size);

    //fprintf(stderr, "[%s] [%d, %d] x [%d, %d]\n", __FUNCTION__, qf, qt, sf, st);

    vector<u8>& qfrag = tbck_data->qfrag;
    qfrag.clear();
    for (int i = qf; i < qt; ++i) {
        qfrag.push_back(fwd_read[i]);
        int qc = fwd_read[i];
        qc = DECODE_RESIDUE(qc);
        //fprintf(stderr, "%c", qc);
    }
    //fprintf(stderr, "\n");

    vector<u8>& sfrag = tbck_data->sfrag;
    sfrag.clear();
    const char* sbj = updb->target_sequence(L->sid, L->sdir);
    for (int i = sf; i < st; ++i) {
        int c = sbj[i];
        //fprintf(stderr, "%c", c);
        c = nst_nt4_table[c];
        if (c > 3) c = 0;
        sfrag.push_back(c);
    }
    //fprintf(stderr, "\n");

    int qfae = 0, sfae = 0;
    int ql = qfrag.size(), sl = sfrag.size();
    run_shw(tbck_data->small_edlib, tbck_data->edlib, qfrag.data(), ql, sfrag.data(), sl, &qfae, &sfae, 
        &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
    //fprintf(stderr, "qfae = %d, sfae = %d, ql = %d, sl = %d\n", qfae, sfae, ql, sl);

    const char* qas = tbck_data->ext_qabuf.c_str(); 
    const char* sas = tbck_data->ext_sabuf.c_str(); 
    int as_size = tbck_data->ext_qabuf.size();
    //dump_align_string(qas, sas, as_size, stderr);
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    //fprintf(stderr, "pi = %g\n", ident);
    if (ident < 60.0) return false;

    L->qend += qfae;
    L->send += sfae;
    set_pca_chain_offset(L, enzyme_size);
    L->is_modified = true;
    L->km_off = 0;
    L->km_c = 0;

    return true;
}

static bool 
s_gf_p2_left_extend_gap_2(PoreCAlign* L,
    PoreCAlign* R,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    int qf = L->qend;
    int qt = R->qoff;
    int ql = qt - qf;
    int st = R->soff;
    int sf = max(0, st - ql - 30);
    int sl = st - sf;
    if (ql == 0 || sl == 0) return false;

    //cerr << *L << '\n' << *R << '\n';
    //fprintf(stderr, "[%s] [%d, %d] x [%d, %d]\n", __FUNCTION__, qf, qt, sf, st);

    vector<u8>& qfrag = tbck_data->qfrag;
    qfrag.clear();
    for (int i = qt - 1; i >= qf; --i) {
        qfrag.push_back(fwd_read[i]);
        int qc = fwd_read[i];
        qc = DECODE_RESIDUE(qc);
    }

    vector<u8>& sfrag = tbck_data->sfrag;
    sfrag.clear();
    const char* sbj = updb->target_sequence(R->sid, R->sdir);
    for (int i = st - 1; i >= sf; --i) {
        int c = sbj[i];
        c = nst_nt4_table[c];
        if (c > 3) c = 0;
        sfrag.push_back(c);
    }

    int qfae = 0, sfae = 0;
    run_shw(tbck_data->small_edlib, tbck_data->edlib, qfrag.data(), ql, sfrag.data(), sl, &qfae, &sfae, 
        &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
    //fprintf(stderr, "qfae = %d, sfae = %d, ql = %d, sl = %d\n", qfae, sfae, ql, sl);

    const char* qas = tbck_data->ext_qabuf.c_str(); 
    const char* sas = tbck_data->ext_sabuf.c_str(); 
    int as_size = tbck_data->ext_qabuf.size();
    while (as_size) {
        int qc = qas[as_size-1];
        int sc = sas[as_size-1];
        if (qc == sc) break;
        if (qc != GAP_CHAR) --qfae;
        if (sc != GAP_CHAR) --sfae;
        --as_size;
    }
    if (as_size == 0) return false;
    if (ql - qfae > 10) return false;
    //dump_align_string(qas, sas, as_size, stderr);
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    //fprintf(stderr, "pi = %g\n", ident);
    if (ident < 60.0) return false;

    R->qoff -= qfae;
    R->soff -= sfae;
    set_pca_chain_offset(R, enzyme_size);
    R->is_modified = true;
    R->km_off = 0;
    R->km_c = 0;

    return true;
}

static bool 
s_gf_p2_right_extend_gap_2(PoreCAlign* L,
    PoreCAlign* R,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int read_size,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    int qf = L->qend;
    int qt = R->qoff;
    int ql = qt - qf;
    int sf = L->send;
    int st = min(L->ssize, sf + ql + 30);
    int sl = st - sf;
    if (ql == 0 || sl == 0) return false;

    //cerr << *L << '\n' << *R << '\n';
    //fprintf(stderr, "[%s] [%d, %d] x [%d, %d]\n", __FUNCTION__, qf, qt, sf, st);

    vector<u8>& qfrag = tbck_data->qfrag;
    qfrag.clear();
    for (int i = qf; i < qt; ++i) {
        qfrag.push_back(fwd_read[i]);
        int qc = fwd_read[i];
        qc = DECODE_RESIDUE(qc);
    }

    vector<u8>& sfrag = tbck_data->sfrag;
    sfrag.clear();
    const char* sbj = updb->target_sequence(L->sid, L->sdir);
    for (int i = sf; i < st; ++i) {
        int c = sbj[i];
        c = nst_nt4_table[c];
        if (c > 3) c = 0;
        sfrag.push_back(c);
    }

    int qfae = 0, sfae = 0;
    run_shw(tbck_data->small_edlib, tbck_data->edlib, qfrag.data(), ql, sfrag.data(), sl, &qfae, &sfae, 
        &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
    //fprintf(stderr, "qfae = %d, sfae = %d, ql = %d, sl = %d\n", qfae, sfae, ql, sl);

    const char* qas = tbck_data->ext_qabuf.c_str(); 
    const char* sas = tbck_data->ext_sabuf.c_str(); 
    int as_size = tbck_data->ext_qabuf.size();
    while (as_size) {
        int qc = qas[as_size-1];
        int sc = sas[as_size-1];
        if (qc == sc) break;
        if (qc != GAP_CHAR) --qfae;
        if (sc != GAP_CHAR) --sfae;
        --as_size;
    }
    if (as_size == 0) return false;
    if (ql - qfae > 10) return false;
    //dump_align_string(qas, sas, as_size, stderr);
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    //fprintf(stderr, "pi = %g\n", ident);
    if (ident < 60.0) return false;

    L->qend += qfae;
    L->send += sfae;
    set_pca_chain_offset(L, enzyme_size);
    L->is_modified = true;
    L->km_off = 0;
    L->km_c = 0;

    return true;
}


void
smooth_pca_list_pass2(std::vector<PoreCAlign>& chain,
    EChainType& chain_type,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int* vdfa,
    const int vdfc,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    PoreCAlign* pcaa = chain.data();
    int pcac = chain.size();

    {
        PoreCAlign* first = pcaa;
        bool r = (first->enzyme_qoff == 0 && first->qoff > 0)
                 ||
                 (first->qoff > 0 && first->qoff <= 200);
        if (r) s_gf_p2_left_most_extend(pcaa, updb, fwd_read, enzyme_size, tbck_data);
    }

    for (int i = 0; i < pcac - 1; ++i) {
        PoreCAlign* pi = pcaa + i;
        PoreCAlign* pj = pcaa + i + 1;
        if (pi->qend >= pj->qoff) continue;

        if (s_gf_p2_left_extend(pi, pj, updb, fwd_read, enzyme_size, tbck_data)) continue;
        if (s_gf_p2_right_extend(pi, pj, updb, fwd_read, enzyme_size, tbck_data)) continue;
        bool r = (pi->enzyme_send > pi->send && pi->enzyme_send - pi->send <= 200) && (pi->qend != pi->enzyme_qend);
        if (r) {
            s_gf_p2_right_extend_gap(pi, pj, updb, fwd_read, pi->qsize, enzyme_size, tbck_data);
        }
        r = (pj->qoff != pj->enzyme_qoff) && (pj->soff > pj->enzyme_soff && pj->soff - pj->enzyme_soff <= 200);
        if (r) {
            s_gf_p2_left_extend_gap(pi, pj, updb, fwd_read, enzyme_size, tbck_data);
        }        
    }

    if (1)
    for (int i = 0; i < pcac - 1; ++i) {
        PoreCAlign* pi = pcaa + i;
        PoreCAlign* pj = pcaa + i + 1;
        //cerr << i << '\n' << *pi << '\n' << *pj << '\n';
        if (pi->qend >= pj->qoff || pj->qoff - pi->qend < 20 || pj->qoff - pi->qend > 200) continue;  

        if (s_gf_p2_left_extend_gap_2(pi, pj, updb, fwd_read, enzyme_size, tbck_data)) continue;
        if (s_gf_p2_right_extend_gap_2(pi, pj, updb, fwd_read, pi->qsize, enzyme_size, tbck_data)) continue;
    }

    {
        PoreCAlign* last = pcaa + pcac - 1;
        bool r = (last->enzyme_qend == last->qsize && last->qend < last->qsize)
                 ||
                 (last->qend < last->qsize && last->qsize - last->qend <= 200);
        if (r) s_gf_p2_right_most_extend(pcaa + pcac - 1, updb, fwd_read, enzyme_size, tbck_data);
    }

    pdqsort(chain.begin(), chain.end(), [](const PoreCAlign& x, const PoreCAlign& y) { return x.qoff < y.qoff; });
    chain_type = pca_chain_type(vdfa, vdfc, chain.data(), chain.size());
}
/////////////////

void
s_merge_adjacent_pca(std::vector<PoreCAlign>& pca_list,
    HbnUnpackedDatabase* updb,
    const u8* fwd_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    PoreCAlign* pca_a = pca_list.data();
    int pca_c = pca_list.size();
    int i = 0;
    while (i < pca_c) {
        PoreCAlign* pi = pca_a + i;
        if (pi->qid == -1) { ++i; continue; }
        int j = i + 1;
        while (j < pca_c) {
            PoreCAlign* pj = pca_a + j;
            if (pj->qid == -1) { ++j; continue; }
            if (pi->sid != pj->sid || pi->sdir != pj->sdir) break;

            int iqb = pi->qoff;
            int isb = pi->soff;
            int jqe = pj->qend;
            int jse = pj->send;

            int q_d = jqe - iqb;
            int s_d = jse - isb;
            if (s_d <= 0) break;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > 0.2) break;
            const char* sbj = updb->target_sequence(pi->sid, pi->sdir);
            vector<u8>& sfrag = tbck_data->sfrag;
            sfrag.clear();
            for (int q = 0; q < s_d; ++q) {
                int c = sbj[isb + q];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* read = fwd_read;
            run_nw(read+iqb, q_d, sfrag.data(), s_d, tbck_data->small_edlib, tbck_data->edlib, 
                tbck_data->ext_qabuf, tbck_data->ext_sabuf);
            const char* qas = tbck_data->ext_qabuf.c_str(); 
            const char* sas = tbck_data->ext_sabuf.c_str(); 
            int as_size = tbck_data->ext_qabuf.size();
            double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
	        double avg_pi = (pi->pi + pj->pi) / 2;
	        if (ident < avg_pi - 4.0 && ident < 75.0) { ++j; continue; }
            PoreCAlign pca = *pi;
            s_merge_two_pca(pi, pj, &pca);
            set_pca_chain_offset(&pca, enzyme_size);
            //HBN_LOG("merge, pi = %g, avg_pi = %g", ident, avg_pi); cerr << i << '\t' << *pi << '\n' << j << '\t' << *pj << '\n';
            //fprintf(stderr, "into"); dump_chain_pca(fprintf, stderr, pca, -1);
            *pi = pca;
            pj->qid = -1;
            ++j;
        }
        i = j;
    }
    int n = 0;
    for (i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
    pca_list.resize(n);
}

static void
x_update_chain(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    PoreCAlign* pa = chain.data();
    int pc = chain.size();
    for (auto& p : pca_list) {
        for (int i = 0; i < pc - 1; ++i) {
            if (pa[i].qid == -1 || pa[i+1].qid == -1) continue;
            if (p.qoff > pa[i].qoff) break;
            if (p.qoff == pa[i].qoff && p.qend == pa[i+1].qend) {
                double avg_pi = (pa[i].pi + pa[i+1].pi) / 2;
                bool r1 = (p.pi > avg_pi - 4.0) || p.pi >= 80.0;
                //bool r2 = (p.pi > avg_pi - 4.0) && ((!pa[i].lqm) || (!pa[i].rqm) || (!pa[i+1].lqm) || (!pa[i+1].rqm));
                if (r1) {
                    //HBN_LOG("%s replace", read_name); dump_chain_pca(fprintf, stderr, pa[i], i); dump_chain_pca(fprintf, stderr, pa[i+1], i+1);
                    //fprintf(stderr, "by\n"); dump_chain_pca(fprintf, stderr, p, -1); 
                    pa[i].qid = -1;
                    pa[i+1] = p;
                    break;
                }
            }
        }
    }
    int n = 0;
    for (int i = 0; i < pc; ++i) if (pa[i].qid >= 0) pa[n++] = pa[i];
    chain.resize(n);
}

static void
x_update_chain_p2(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    PoreCAlign* pa = chain.data();
    int pc = chain.size();
    for (auto& pca : pca_list) {
        int f = -1, t = -1;
        for (int i = 0; i < pc; ++i) {
            if (pa[i].qoff == pca.qoff) {
                f =i;
                break;
            }
        }
        if (f == -1) continue;
        for (int i = 0; i < pc; ++i) {
            if (pa[i].qend == pca.qend) {
                t = i;
                break;
            }
        }
        if (t == -1) continue;
        if (f >= t) continue;

        double avg_pi = 0;
        for (int i = f; i <= t; ++i) avg_pi += pa[i].pi;
        avg_pi /= (t - f + 1);
        if (pca.pi < avg_pi - 4.0 && pca.pi < 80.0) continue;
        pa[f] = pca;
        for (int i = f + 1; i <= t; ++i) pa[i].qid = -1;
        int n = 0;
        for (int i = 0; i < pc; ++i) if (pa[i].qid >= 0) pa[n++] = pa[i];
        pc = n;
    }
    chain.resize(pc);
}

static void
x_update_chain_p3(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    PoreCAlign* pa = chain.data();
    int pc = chain.size();
    for (auto& pca : pca_list) {
        for (int i = 0; i < pc; ++i) {
	        if (pca.qoff == pa[i].qoff && pca.qend == pa[i].qend && pca.pi > pa[i].pi) {
                pa[i] = pca;
                break;
            }
        }
    }
}

void
update_chain(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    pdqsort(chain.begin(), chain.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    x_update_chain(chain, pca_list, read_name);
    pdqsort(chain.begin(), chain.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    x_update_chain_p2(chain, pca_list, read_name);
    x_update_chain_p3(chain, pca_list, read_name);
}

void
smooth_pca_list(std::vector<PoreCAlign>& pca_list,
    EChainType& chain_type,
    std::vector<PoreCAlign>& all_pca_list,
    HbnUnpackedDatabase* updb,
    PoreCQuery* query,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    pdqsort(pca_list.begin(), pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    s_merge_adjacent_pca(pca_list, updb, query->fwd_qs, enzyme_size, tbck_data);
    pdqsort(pca_list.begin(), pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    update_chain(pca_list, all_pca_list, query->name);
    smooth_pca_list_pass2(pca_list, chain_type, updb, query->fwd_qs,
        query->fwd_enzyme_pos, query->fwd_enzyme_pos_cnt, enzyme_size, tbck_data);
}
