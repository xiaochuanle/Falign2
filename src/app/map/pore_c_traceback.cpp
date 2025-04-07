#include "pore_c_traceback.hpp"

#include "../../sw/hbn_traceback_aux.h"

using namespace std;

void
run_nw(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    string& qaln,
    string& saln)
{
    if (query_length < SMALL_EDLIB_MAX_SEQ_SIZE && subject_length < SMALL_EDLIB_MAX_SEQ_SIZE) {
        small_edlib_nw(query, query_length, subject, subject_length, small_edlib, qaln, saln);
    } else {
        edlib_nw(edlib, query, query_length, subject, subject_length, qaln, saln);
    }
}

void
run_nw_cigar(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    vector<u64>& cigar)
{
    if (query_length < SMALL_EDLIB_MAX_SEQ_SIZE && subject_length < SMALL_EDLIB_MAX_SEQ_SIZE) {
        small_edlib_nw_cigar(query, query_length, subject, subject_length, small_edlib, cigar);
    } else {
        edlib_nw_cigar(edlib, query, query_length, subject, subject_length, cigar);
    }
}

int run_nw_dist(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    int tolerence,
    int* score)
{
    if (query_length < SMALL_EDLIB_MAX_SEQ_SIZE && subject_length < SMALL_EDLIB_MAX_SEQ_SIZE) {
        return small_edlib_nw_dist(query, query_length, subject, subject_length, small_edlib, tolerence, score);
    } else {
        return edlib_nw_dist(edlib, query, query_length, subject, subject_length, tolerence, score);
    }
}

int run_shw(small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    int* _qend,
    int* _send,
    std::string* qaln,
    std::string* saln)
{
    if (query_length > SMALL_EDLIB_MAX_SEQ_SIZE || subject_length > SMALL_EDLIB_MAX_SEQ_SIZE) {
        if (query_length <= subject_length) {
             return edlib_shw(edlib, query, query_length, subject, subject_length, _qend, _send, qaln, saln);
        } else {
             return edlib_shw(edlib, subject, subject_length, query, query_length, _send, _qend, saln, qaln);
        }
    } else {
        if (query_length <= subject_length) {
            return small_edlib_shw(small_edlib, query, query_length, subject, subject_length, _qend, _send, qaln, saln);
        } else {
            return small_edlib_shw(small_edlib, subject, subject_length, query, query_length, _send, _qend, saln, qaln);   
        }
    }
}

void
right_extend_align(small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    const u8* query,
    int* qend,
    const int query_length,
    const char* subject,
    vector<u8>& sfrag,
    int* send,
    const int subject_length,
    string& qabuf,
    string& sabuf,
    string& qas,
    string& sas)
{
    qas.clear();
    sas.clear();

    int qt = *qend;
    int st = *send;
    bool last_block = false;
    bool done = false;
    while (1) {
        int qblk, sblk;
        if (query_length - qt <= ALIGN_MAX_EXT_SEG || subject_length - st <= ALIGN_MAX_EXT_SEG) {
            qblk = (subject_length - st) + 30;;
            qblk = min<int>(qblk, query_length - qt);
            sblk = (query_length - qt) + 30;
            sblk = min<int>(sblk, subject_length - st);
            last_block = true;
        } else {
            qblk = ALIGN_EXT_SEG;
            sblk = ALIGN_EXT_SEG;
            last_block = false;
        }
        if (qblk == 0 || sblk == 0) break;
        //HBN_LOG("1 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, query_length, subject_length);
        int qfae = 0, sfae = 0;
        qabuf.clear();
        sabuf.clear();
        sfrag.clear();
        for (int i = 0; i < sblk; ++i) {
            int c = subject[st+i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
#if 1
        run_shw(small_edlib, edlib, query + qt, qblk, sfrag.data(), sblk, &qfae, &sfae, &qabuf, &sabuf);
	if (qblk - qfae > 10 || sblk - sfae > 10) done = true;
#else
        run_nw(query + qt, qblk, sfrag.data(), sblk, small_edlib, edlib, qabuf, sabuf);
        qfae = qblk;
        sfae = sblk;
#endif
        int acnt = 0, qcnt = 0, scnt = 0;
        int as_size = qabuf.size(), k = as_size - 1, m = 0;
        while (k >= 0) {
            char qc = qabuf[k];
            char sc = sabuf[k];
            if (qc != GAP_CHAR) ++qcnt;
            if (sc != GAP_CHAR) ++scnt;
            m = (qc == sc) ? (m+1) : 0;
            ++acnt;
            if (m == ALIGN_END_MATCH) break;
            --k;
        }
        //HBN_LOG("qblk = %d, qfae = %d, sblk = %d, sfae = %d, m = %d, k = %d, qcnt = %d, scnt = %d, acnt = %d", qblk, qfae, sblk, sfae, m, k, qcnt, scnt, acnt);
        if (m != ALIGN_END_MATCH || k < 1) {
            as_size = 0;
            for (int i = 0; i < qblk && i < sblk; ++i, ++qt, ++st) {
                int qc = query[qt];
                qc = DECODE_RESIDUE(qc);
                int sc = sfrag[i];
                sc = DECODE_RESIDUE(sc);
                if (qc != sc) break;
                qas += qc;
                sas += sc;
                ++as_size;
            }
            done = true;
        } else {
            as_size -= acnt;
            qt += (qfae - qcnt);
            st += (sfae - scnt);
            if (done) {
                as_size += ALIGN_END_MATCH;
                qt += ALIGN_END_MATCH;
                st += ALIGN_END_MATCH;
            }
            for (int i = 0; i < as_size; ++i) {
                qas += qabuf[i];
                sas += sabuf[i];
            }
        }
        //HBN_LOG("qt = %d, st = %d", qt, st);
        if (done) break;
    }

    *qend = qt;
    *send = st;
}

void
left_extend_align(small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    const u8* query,
    vector<u8>& qfrag,
    int* qoff,
    const char* subject,
    vector<u8>& sfrag,
    int* soff,
    string& qabuf,
    string& sabuf,
    string& qas,
    string& sas)
{
    qas.clear();
    sas.clear();

    int qf = *qoff;
    int sf = *soff;
    bool last_block = false;
    bool done = false;
    while (1) {
        int qblk, sblk;
        if (qf <= ALIGN_MAX_EXT_SEG || sf <= ALIGN_MAX_EXT_SEG) {
            qblk = min<int>(sf + 30, qf);
            sblk = min<int>(qf + 30, sf);
            last_block = true;
        } else {
            qblk = ALIGN_EXT_SEG;
            sblk = ALIGN_EXT_SEG;
            last_block = false;
        }
        if (qblk == 0 || sblk == 0) break;
        qfrag.clear(); 
        for (int i = 1; i <= qblk; ++i) {
            qfrag.push_back(query[qf-i]);
        }
        sfrag.clear(); 
        for (int i = 1; i <= sblk; ++i) {
            int c = subject[sf-i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        //fprintf(stderr, "1 qf = %d, sf = %d, qblk = %d, sblk = %d\n", qf, sf, qblk, sblk);

        int qfae = 0, sfae = 0;
        qabuf.clear();
        sabuf.clear();
        run_shw(small_edlib, edlib, qfrag.data(), qblk, sfrag.data(), sblk, &qfae, &sfae, &qabuf, &sabuf);
	    if (qblk - qfae > 10 || sblk - sfae > 10) done = true;
        int acnt = 0, qcnt = 0, scnt = 0;
        int as_size = qabuf.size(), k = as_size - 1, m = 0;
        while (k >= 0) {
            char qc = qabuf[k];
            char sc = sabuf[k];
            if (qc != GAP_CHAR) ++qcnt;
            if (sc != GAP_CHAR) ++scnt;
            m = (qc == sc) ? (m+1) : 0;
            ++acnt;
            if (m == ALIGN_END_MATCH) break;
            --k;
        }
        //HBN_LOG("qfae = %d, safe = %d, m = %d, k = %d, qcnt = %d, scnt = %d, acnt = %d", qfae, sfae, m, k, qcnt, scnt, acnt);
        if (m != ALIGN_END_MATCH || k < 1) {
            for (int i = 0; i < qblk && i < sblk; ++i) {
                int qc = qfrag[i];
                qc = DECODE_RESIDUE(qc);
                int sc = sfrag[i];
                sc = DECODE_RESIDUE(sc);
                if (qc != sc) break;
                qas += qc;
                sas += sc;
                --qf; --sf;
            }
            done = true;
        } else {
            as_size -= acnt;
            qf -= (qfae - qcnt);
            sf -= (sfae - scnt);
            if (done) {
                as_size += ALIGN_END_MATCH;
                qf -= ALIGN_END_MATCH;
                sf -= ALIGN_END_MATCH;
            }
            for (int i = 0; i < as_size; ++i) {
                qas += qabuf[i];
                sas += sabuf[i];
            }
        }
        //HBN_LOG("sf = %d, qf = %d", qf, sf);
        if (done) break;
    }

    *qoff = qf;
    *soff = sf;
}

BOOL
truncate_align_bad_ends(const char* qaln,
    const char* saln,
    const int aln_size,
    int* qoff,
    int* qend,
    int* soff,
    int* send,
    const char** qas_,
    const char** qae_,
    const char** sas_,
    const char** sae_)
{
    const char* qas = qaln;
    const char* qae = qaln + aln_size;
    const char* sas = saln;
    const char* sae = saln + aln_size;
    int qcnt = 0, scnt = 0;
    const char* qa = qas;
    const char* sa = sas;
    int m = 0;

    while (qa < qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == ALIGN_END_MATCH) break;
    }
    if (m < ALIGN_END_MATCH || qa == qae) return 0;
    *qas_ = qa - ALIGN_END_MATCH;
    *sas_ = sa - ALIGN_END_MATCH;
    *qoff += qcnt - ALIGN_END_MATCH;
    *soff += scnt - ALIGN_END_MATCH;  

    qcnt = 0;
    scnt = 0;
    qa = qae;
    sa = sae;
    m = 0;
    while (qa > qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == ALIGN_END_MATCH) break;
    }
    hbn_assert(m == ALIGN_END_MATCH, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == qas) return 0;

    *qae_ = qa + ALIGN_END_MATCH;
    *sae_ = sa + ALIGN_END_MATCH;
    *qend -= qcnt - ALIGN_END_MATCH;
    *send -= scnt - ALIGN_END_MATCH; 
    return TRUE;
}

static void
s_append_cigar1(u64 op, u64 len, vector<u64>& v)
{
    if ((!v.empty()) && (v.back()&0xf) == op) {
        v.back() += (len<<4);
    } else {
        v.push_back((len<<4)|op);
    }
}

static void
s_append_cigar(int n_cigar, u32* cigar, vector<u64>& v)
{
    if (n_cigar == 0) return;
    s_append_cigar1(cigar[0]&0xf, cigar[0]>>4, v);
    for (int i = 1; i < n_cigar; ++i) v.push_back(cigar[i]);
}

static void
s_append_alignment_to_cigar(const char* qas, const char* sas, const int as_size, int qb, int qe, int sb, int se, vector<u64>& v)
{
    int qi = qb, si = sb;
    for (int i = 0; i < as_size; ++i) {
        if (qas[i] != GAP_CHAR) ++qi;
        if (sas[i] != GAP_CHAR) ++si;
    }
    if (qi != qe || si != se) dump_align_string(qas, sas, as_size, stderr);
    hbn_assert(qi == qe, "qi = %d, qb = %d, qe = %d, si = %d, sb = %d, se = %d", qi, qb, qe, si, sb, se);
    hbn_assert(si == se);
    qi = qb;
    si = sb;
    int i = 0;
    while (i < as_size) {
        u64 op = 0, len = 0;
        int j = i;
        if (qas[i] == sas[i]) {
            op = 7;
            while (j < as_size && qas[j] == sas[j]) ++j;
            qi += (j - i);
            si += (j - i);
        } else if (qas[i] == GAP_CHAR) {
            op = 2;
            while (j < as_size && qas[j] == GAP_CHAR) ++j;
            si += (j - i);
        } else if (sas[i] == GAP_CHAR) {
            op = 1;
            while (j < as_size && sas[j] == GAP_CHAR) ++j;
            qi += (j - i);
        } else {
            op = 8;
            while (j < as_size) {
                if (qas[j] == GAP_CHAR || sas[j] == GAP_CHAR || qas[j] == sas[j]) break;
                ++j;
            }
            qi += (j - i);
            si += (j - i);
        }
        len = j - i;
        s_append_cigar1(op, len, v);
        i = j;
    }
    hbn_assert(qi == qe, "qi = %d, qb = %d, qe = %d, si = %d, sb = %d, se = %d", qi, qb, qe, si, sb, se);
    hbn_assert(si == se);
}

static bool
s_truncate_cigar_bad_ends(vector<u64>& cigar, int& qb, int& qe, int& sb, int& se)
{
    int nc = cigar.size();
    if (!nc) return false;
    int fi = 0;
    while (fi < nc) {
        int op = cigar[fi]&0xf;
        int len = cigar[fi]>>4;
        if (op == 7 && len >= ALIGN_END_MATCH) break;
        if (op == 7 || op == 8) {
            qb += len;
            sb += len;
        } else if (op == 1) {
            qb += len;
        } else if (op == 2) {
            sb += len;
        }
        ++fi;
    }
    if (fi >= nc) return false;

    int ri = nc;
    while (ri) {
        int op = cigar[ri-1]&0xf;
        int len = cigar[ri-1]>>4;
        if (op == 7 && len >= ALIGN_END_MATCH) break;
        if (op == 7 || op == 8) {
            qe -= len;
            se -= len;
        } else if (op == 1) {
            qe -= len;
        } else if (op == 2) {
            se -= len;
        }
        --ri;
    }
    hbn_assert(ri > fi);

    int n = 0;
    if (fi > 0) {
        for (int i = fi; i < ri; ++i) cigar[n++] = cigar[i];
    } else {
        n = ri;
    }
    cigar.resize(n);
    return true;
}

static void
s_validate_cigar_and_seq_len(const int n_cigar, const u64* cigar, 
    int qb, int qe, int sb, int se)
{
    int qi = qb, si = sb;
    for (int i = 0; i < n_cigar; ++i) {
        int op = cigar[i]&0xf;
        int len = cigar[i]>>4;
        if (op == 7) {
            qi += len;
            si += len;
        } else if (op == 8) {
            qi += len;
            si += len;
        } else if (op == 1) {
            qi += len;
        } else if (op == 2) {
            si += len;
        }
    }
    hbn_assert(qi == qe, "qi = %d, qb = %d, qe = %d, si = %d, sb = %d, se = %d", qi, qb, qe, si, sb, se);
    hbn_assert(si == se);
}

static void
s_compute_alignment_score_from_cigar(const int n_cigar, const u64* cigar, 
    int qb, int qe, int sb, int se,
    int* score, int* match, float* pi, float* epi)
{
    int qi = qb, si = sb;
    int mat = 0, eff_len = 0, as_size = 0;
    int sc = 0;
    int a, b, go, ge;
    extract_sw_scoring_params(&a, &b, &go, &ge, nullptr, nullptr);
    for (int i = 0; i < n_cigar; ++i) {
        int op = cigar[i]&0xf;
        int len = cigar[i]>>4;
        as_size += len;
        if (op == 7) {
            sc += a * len;
            qi += len;
            si += len;
            mat += len;
            eff_len += len;
        } else if (op == 8) {
            sc -= b * len;
            qi += len;
            si += len;
            eff_len += len;
        } else if (op == 1) {
            sc -= ge * len;
            qi += len;
            if (len <= ALIGN_END_MATCH) eff_len += len;
        } else if (op == 2) {
            sc -= ge * len;
            si += len;
            if (len <= ALIGN_END_MATCH) eff_len += len;
        }
    }
    hbn_assert(qi == qe, "qi = %d, qb = %d, qe = %d, si = %d, sb = %d, se = %d", qi, qb, qe, si, sb, se);
    hbn_assert(si == se);
    if (score) *score = sc;
    if (match) *match = mat;
    if (pi) *pi = 100.0 * mat / as_size;
    if (epi) *epi = 100.0 * mat / eff_len;
}

void HbnTracebackData::cigar_to_align_string(const char* FQ, const char* S)
{
    vqas.clear();
    vsas.clear();
    int qi = qoff;
    int si = soff;
    for (auto c : cigar) {
        int op = c&0xf;
        int len = c>>4;
        if (op == 7) {
            for (int s = 0; s < len; ++s, ++qi, ++si) {
                vqas.push_back(FQ[qi]);
                vsas.push_back(S[si]);
            }
        } else if (op == 8) {
            for (int s = 0; s < len; ++s, ++qi, ++si) {
                vqas.push_back(FQ[qi]);
                vsas.push_back(S[si]);
            }
        } else if (op == 1) {
            for (int s = 0; s < len; ++s, ++qi) {
                vqas.push_back(FQ[qi]);
                vsas.push_back(GAP_CHAR);
            }
        } else if (op == 2) {
            for (int s = 0; s < len; ++s, ++si) {
                vqas.push_back(GAP_CHAR);
                vsas.push_back(S[si]);
            }
        }
    }
    hbn_assert(qi == qend);
    hbn_assert(si == send);
    qas = vqas.data();
    qae = qas + vqas.size();
    sas = vsas.data();
    sae = sas + vsas.size();
}

bool HbnTracebackData::map_one_chain(const u8* EQ, int qs, const char* S, int ss, const pu64_t* a, int n_a, bool xxx)
{
    qas = qae = nullptr;
    sas = sae = nullptr;
    cigar.clear();
    
    int kmer = a->y>>32&0xff;
    int qb = (int32_t)a->y + 1 - kmer;
    int sb = (int32_t)a->x + 1 - kmer;
    qoff = qb;
    soff = sb;
    if (xxx)
    if (qb && sb) {
        left_extend_align(small_edlib, edlib, EQ, qfrag, &qoff, S, sfrag, &soff, ext_qabuf, ext_sabuf, qabuf, sabuf);
        reverse(qabuf.begin(), qabuf.end());
        reverse(sabuf.begin(), sabuf.end());
        s_append_alignment_to_cigar(qabuf.c_str(), sabuf.c_str(), qabuf.size(), qoff, qb, soff, sb, cigar);
        s_validate_cigar_and_seq_len(cigar.size(), cigar.data(), qoff, qb, soff, sb);
    }

    s_append_cigar1(7, kmer, cigar);
    int j0 = 0;
    for (int j = 1; j < n_a; ++j) {
        const pu64_t* p = a + j;
        const pu64_t* q = a + j0;
        int l_seq = (int32_t)p->x - (int32_t)q->x;
        const char* vss = S + (int32_t)q->x + 1;
        int qlen = (int32_t)p->y - (int32_t)q->y;
        const u8* qss = EQ + (int32_t)q->y + 1;
        if (l_seq == 0) s_append_cigar1(1, qlen, cigar);
        else if (qlen == 0) s_append_cigar1(2, l_seq, cigar);
        else if (l_seq == qlen && qlen <= (q->y>>32&0xff)) s_append_cigar1(7, qlen, cigar);
        else {
            sfrag.clear();
            for (int s = 0; s < l_seq; ++s) {
                int c = vss[s];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            run_nw_cigar(qss, qlen, sfrag.data(), l_seq, small_edlib, edlib, cigar);
        }
        j0 = j;
    }
    int qe = (int32_t)a[n_a-1].y + 1;
    int se = (int32_t)a[n_a-1].x + 1;
    qend = qe;
    send = se;
    {
        int m = 0;
        while (qend < qs && send < ss) {
            int qc = EQ[qend];
            qc = DECODE_RESIDUE(qc);
            int sc = S[send];
            if (qc != sc) break;
            ++m;
            ++qend;
            ++send;
        }
        if (m) s_append_cigar1(7, m, cigar);
    }
    if (!s_truncate_cigar_bad_ends(cigar, qoff, qend, soff, send)) return false;
    qe = qend;
    se = send;

    if (xxx)
    if (qend < qs && send < ss) {
        right_extend_align(small_edlib, edlib, EQ, &qend, qs, S, sfrag, &send, ss, ext_qabuf, ext_sabuf, qabuf, sabuf);
        s_append_alignment_to_cigar(qabuf.c_str(), sabuf.c_str(), qabuf.size(), qe, qend, se, send, cigar);
    }

    s_compute_alignment_score_from_cigar(cigar.size(), cigar.data(), qoff, qend, soff, send, &score, &match, &pi, &epi);
    s_validate_cigar_and_seq_len(cigar.size(), cigar.data(), qoff, qend, soff, send);
    return true;
}

bool HbnTracebackData::map_one_chain(const u8* EQ, int qb, int qe, int qs, const char* S, int sb, int se, int ss)
{
    qas = qae = nullptr;
    sas = sae = nullptr;
    cigar.clear();

    while (qb && sb) {
        int qc = DECODE_RESIDUE(EQ[qb-1]);
        int sc = toupper(S[sb-1]);
        if (qc != sc) break;
        --qb;
        --sb;
    }
    while (qe < qs && se < ss) {
        int qc = DECODE_RESIDUE(EQ[qe]);
        int sc = toupper(S[se]);
        if (qc != sc) break;
        ++qe;
        ++se;
    }

        sfrag.clear();
        for (int i = sb; i < se; ++i) {
            int c = S[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            sfrag.push_back(c);
        }
        const u8* qss = EQ + qb;
        int ql = qe - qb;
        const u8* vss = sfrag.data();
        int sl = sfrag.size();
        run_nw(qss, ql, vss, sl, small_edlib, edlib, ext_qabuf, ext_sabuf);
        int as_size = ext_qabuf.size();
        const char* xqas = ext_qabuf.c_str();
        const char* xqae = xqas + as_size;
        const char* xsas = ext_sabuf.c_str();
        const char* xsae = xsas + as_size;
        if (!truncate_align_bad_ends(xqas, xsas, as_size, &qb, &qe, &sb, &se, &xqas, &xqae, &xsas, &xsae)) return false;
        if (qe - qb < 20 || se - sb < 20) return false;
        as_size = xqae - xqas;
    validate_aligned_string_2(__FILE__, __FUNCTION__, __LINE__, 
        0, EQ, qb, qe, xqas,
        0, S, sb, se, xsas, as_size, TRUE);
        qabuf.assign(xqas, xqae);
        sabuf.assign(xsas, xsae);

    if (qb && sb) {
        left_extend_align(small_edlib, edlib, EQ, qfrag, &qb, S, sfrag, &sb, edlib->sqaln, edlib->staln, ext_qabuf, ext_sabuf);
        reverse(ext_qabuf.begin(), ext_qabuf.end());
        reverse(ext_sabuf.begin(), ext_sabuf.end());
        ext_qabuf += qabuf;
        ext_sabuf += sabuf;
        qabuf = ext_qabuf;
        sabuf = ext_sabuf;

        validate_aligned_string_2(__FILE__, __FUNCTION__, __LINE__, 
        0, EQ, qb, qe, qabuf.c_str(),
        0, S, sb, se, sabuf.c_str(), qabuf.size(), TRUE);
    } 


    if (qe < qs && se < ss) {
        right_extend_align(small_edlib, edlib, EQ, &qe, qs, S, sfrag, &se, ss, edlib->sqaln, edlib->staln, ext_qabuf, ext_sabuf);
        xqas = ext_qabuf.c_str();
        xsas = ext_sabuf.c_str();
        as_size = ext_qabuf.size();
        //dump_align_string(xqas, xsas, as_size, stderr);
        qabuf += ext_qabuf;
        sabuf += ext_sabuf;
    }
    
    validate_aligned_string_2(__FILE__, __FUNCTION__, __LINE__, 
        0, EQ, qb, qe, qabuf.c_str(), 
        0, S, sb, se, sabuf.c_str(), qabuf.size(), TRUE);

    qoff = qb;
    qend = qe;
    soff = sb;
    send = se;
    s_append_alignment_to_cigar(qabuf.c_str(), sabuf.c_str(), qabuf.size(), qb, qe, sb, se, cigar);

    s_compute_alignment_score_from_cigar(cigar.size(), cigar.data(), qoff, qend, soff, send, &score, &match, &pi, &epi);
    s_validate_cigar_and_seq_len(cigar.size(), cigar.data(), qoff, qend, soff, send);
    return true;
}

////////////////////

static void
s_appedn_ksw_cigar(Ksw2Data* ksw, const u8* fwd_qs, int qb, int qe, const char* S, int sb, int se, vector<u64>& cigar)
{
    int qi = qb, si = sb;
    for (int i = 0; i < ksw->ez.n_cigar; ++i) {
        int len = ksw->ez.cigar[i]>>4;
        int op = ksw->ez.cigar[i]&0xf;
        if (op == 1) {
            s_append_cigar1(op, len, cigar);
            qi += len;
        } else if (op == 2) {
            s_append_cigar1(op, len, cigar);
            si += len;
        } else {
            hbn_assert(op == 0);
            int s = 0;
            while (s < len) {
                int qc = DECODE_RESIDUE(fwd_qs[qi+s]);
                int sc = toupper(S[si+s]);
                if (qc == sc) {
                    int t = s + 1;
                    while (t < len) {
                        qc = DECODE_RESIDUE(fwd_qs[qi+t]);
                        sc = toupper(S[si+t]);
                        if (qc != sc) break;
                        ++t;
                    }
                    s_append_cigar1(7, t - s, cigar);
                    s = t;
                } else {
                    int t = s + 1;
                    while (t < len) {
                        qc = DECODE_RESIDUE(fwd_qs[qi+t]);
                        sc = toupper(S[si+t]);
                        if (qc == sc) break;
                        ++t;
                    }
                    s_append_cigar1(8, t - s, cigar);
                    s = t;
                }
            }
            qi += len;
            si += len;
        }
    }
    hbn_assert(qi == qe);
    hbn_assert(si == se);
}

bool HbnTracebackData::refined_map_one_chain(const u8* fwd_qs, int iqb, int iqe, const char* S, int isb, int ise, const pu64_t* a, int n_a)
{
    cigar.clear();
    qas = qae = sas = sae = nullptr;

    qoff = iqb;
    soff = isb;
    int qb = (int32_t)a[0].y + 1 - (a[0].y>>32&0xff);
    int sb = (int32_t)a[0].x + 1 - (a[0].y>>32&0xff);
    {
        bool need_nw = (qb >= iqb && sb >= isb) && (qb > iqb || sb > isb);
        if (need_nw && iqb == qb) {
            s_append_cigar1(2, sb - isb, cigar);
        } else if (need_nw && isb == sb) {
            s_append_cigar1(1, qb - iqb, cigar);
        } else if (need_nw) {
            //fprintf(stderr, "num-mm: %d\n", n_a);
            //fprintf(stderr, "init [%d, %d] x [%d, %d]\n", iqb, qb, isb, sb);
            const u8* qss = fwd_qs + iqb;
            int ql = qb - iqb;
            sfrag.clear();
            for (int i = isb; i < sb; ++i) {
                int c = S[i];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* sss = sfrag.data();
            int sl = sb - isb;
            ksw2_nw(ksw, qss, ql, sss, sl);
            s_appedn_ksw_cigar(ksw, qss, 0, ql, S, isb, sb, cigar);
        } else {
            qoff = qb;
            soff = sb;
        }
    }

    s_append_cigar1(7, a[0].y>>32&0xff, cigar);
    int j0 = 0;
    for (int j = 1; j < n_a; ++j) {
        const pu64_t* p = a + j;
        const pu64_t* q = a + j0;
        const char* vs = S + (int32_t)q->x + 1;
        int sl = (int32_t)p->x - (int32_t)q->x;
        const u8* qss = fwd_qs + (int32_t)q->y + 1;
        int ql = (int32_t)p->y - (int32_t)q->y;
        if (sl == 0) s_append_cigar1(1, ql, cigar);
        else if (ql == 0) s_append_cigar1(2, sl, cigar);
        else if (ql == sl && ql <= (q->y>>32&0xff)) s_append_cigar1(7, ql, cigar);
        else {
            sfrag.clear();
            for (int s = 0; s < sl; ++s) {
                int c = vs[s];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* sss = sfrag.data();
            ksw2_nw(ksw, qss, ql, sss, sl);
            s_appedn_ksw_cigar(ksw, qss, 0, ql, vs, 0, sl, cigar);
        }
        j0 = j;
    }
    qend = iqe;
    send = ise;
    int qe = (int32_t)a[n_a-1].y + 1;
    int se = (int32_t)a[n_a-1].x + 1;
    {
        bool need_nw = (iqe >= qe && ise >= se) && (iqe > qe || ise > se);
        if (need_nw && qe == iqe) {
            s_append_cigar1(2, ise - se, cigar);
        } else if (need_nw && se == ise) {
            s_append_cigar1(1, iqe - qe, cigar);
        } else if (need_nw) {
            const u8* qss = fwd_qs + qe;
            int ql = iqe - qe;
            sfrag.clear();
            for (int i = se; i < ise; ++i) {
                int c = S[i];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* sss = sfrag.data();
            int sl = ise - se;
            ksw2_nw(ksw, qss, ql, sss, sl);
            s_appedn_ksw_cigar(ksw, qss, 0, ql, S, se, ise, cigar);
        } else {
            qend = qe;
            send = se;
        }
    }

    s_compute_alignment_score_from_cigar(cigar.size(), cigar.data(), qoff, qend, soff, send, &score, &match, &pi, &epi);
    s_validate_cigar_and_seq_len(cigar.size(), cigar.data(), qoff, qend, soff, send);
    return true;
}

bool HbnTracebackData::refined_map_one_chain_edlib(const u8* fwd_qs, int iqb, int iqe, 
    const char* S, int isb, int ise, const pu64_t* a, int n_a)
{
    cigar.clear();
    qas = qae = sas = sae = nullptr;

    qoff = iqb;
    soff = isb;
    int qb = (int32_t)a[0].y + 1 - (a[0].y>>32&0xff);
    int sb = (int32_t)a[0].x + 1 - (a[0].y>>32&0xff);
    {
        bool need_nw = (qb >= iqb && sb >= isb) && (qb > iqb || sb > isb);
        if (need_nw && iqb == qb) {
            s_append_cigar1(2, sb - isb, cigar);
        } else if (need_nw && isb == sb) {
            s_append_cigar1(1, qb - iqb, cigar);
        } else if (need_nw) {
            //fprintf(stderr, "num-mm: %d\n", n_a);
            //fprintf(stderr, "init [%d, %d] x [%d, %d]\n", iqb, qb, isb, sb);
            const u8* qss = fwd_qs + iqb;
            int ql = qb - iqb;
            sfrag.clear();
            for (int i = isb; i < sb; ++i) {
                int c = S[i];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* sss = sfrag.data();
            int sl = sb - isb;
            run_nw_cigar(qss, ql, sss, sl, small_edlib, edlib, cigar);
        } else {
            qoff = qb;
            soff = sb;
        }
    }

    s_append_cigar1(7, a[0].y>>32&0xff, cigar);
    int j0 = 0;
    for (int j = 1; j < n_a; ++j) {
        const pu64_t* p = a + j;
        const pu64_t* q = a + j0;
        const char* vs = S + (int32_t)q->x + 1;
        int sl = (int32_t)p->x - (int32_t)q->x;
        const u8* qss = fwd_qs + (int32_t)q->y + 1;
        int ql = (int32_t)p->y - (int32_t)q->y;
        if (sl == 0) s_append_cigar1(1, ql, cigar);
        else if (ql == 0) s_append_cigar1(2, sl, cigar);
        else if (ql == sl && ql <= (q->y>>32&0xff)) s_append_cigar1(7, ql, cigar);
        else {
            sfrag.clear();
            for (int s = 0; s < sl; ++s) {
                int c = vs[s];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* sss = sfrag.data();
            run_nw_cigar(qss, ql, sss, sl, small_edlib, edlib, cigar);
        }
        j0 = j;
    }
    qend = iqe;
    send = ise;
    int qe = (int32_t)a[n_a-1].y + 1;
    int se = (int32_t)a[n_a-1].x + 1;
    {
        bool need_nw = (iqe >= qe && ise >= se) && (iqe > qe || ise > se);
        if (need_nw && qe == iqe) {
            s_append_cigar1(2, ise - se, cigar);
        } else if (need_nw && se == ise) {
            s_append_cigar1(1, iqe - qe, cigar);
        } else if (need_nw) {
            const u8* qss = fwd_qs + qe;
            int ql = iqe - qe;
            sfrag.clear();
            for (int i = se; i < ise; ++i) {
                int c = S[i];
                c = nst_nt4_table[c];
                if (c > 3) c = 0;
                sfrag.push_back(c);
            }
            const u8* sss = sfrag.data();
            int sl = ise - se;
            run_nw_cigar(qss, ql, sss, sl, small_edlib, edlib, cigar);
        } else {
            qend = qe;
            send = se;
        }
    }

    s_compute_alignment_score_from_cigar(cigar.size(), cigar.data(), qoff, qend, soff, send, &score, &match, &pi, &epi);
    s_validate_cigar_and_seq_len(cigar.size(), cigar.data(), qoff, qend, soff, send);
    return true;
}
