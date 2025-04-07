#include "prelim_search.hpp"

using namespace std;

static void
x_extract_kmer_matches(PoreCQuery* query, HbnIndex* index, vector<pu64_t>& km_list)
{
    km_list.clear();

    pu64_t H;
    KmerInfo* ol;
    u64 oll;
    const int kmer_size = index->kmer_size();
    const u64 kHashMask = (U64_ONE << (kmer_size * 2)) - 1;
    const u64 kRevShift = 2 * (kmer_size - 1);

    int pos = 0;
    u64 fhash = 0, rhash = 0;
    for (int i = 0; i < query->size; ++i) {
        fhash = ((fhash << 2) | query->fwd_qs[i]) & kHashMask;
        rhash = (rhash >> 2) | (3ULL^query->fwd_qs[i]) << kRevShift;
        if (i + 1 - pos == kmer_size) {
            if (fhash < rhash) {
                index->extract_offset_list(fhash, &ol, &oll);
                for (u64 s = 0; s < oll; ++s) {
                    int soff = ol[s].seq_offset/2;
                    int sdir = ol[s].seq_offset&1;
                    int qoff = pos;
                    if (sdir == REV) qoff = query->size - 1 - (pos + kmer_size - 1);
                    int strand = (sdir == FWD) ? FWD : REV;
                    set_kmer_match(qoff + kmer_size - 1, ol[s].seq_id, strand, soff + kmer_size - 1, kmer_size, &H);
                    km_list.push_back(H);
                }
            } else if (fhash > rhash) {
                index->extract_offset_list(rhash, &ol, &oll);
                for (u64 s = 0; s < oll; ++s) {
                    int soff = ol[s].seq_offset/2;
                    int sdir = ol[s].seq_offset&1;
                    int qoff = pos;
                    if (sdir == FWD) qoff = query->size - 1 - (pos + kmer_size - 1);
                    int strand = (sdir == REV)? FWD : REV;
                    set_kmer_match(qoff + kmer_size - 1, ol[s].seq_id, strand, soff + kmer_size - 1, kmer_size, &H);
                    km_list.push_back(H);
                }
            }
            ++pos;
        }
    }

    pu64_t* A = km_list.data();
    int NA = km_list.size();
    sort_pu64_x(A, A + NA);
    {
        int i = 0;
        while (i < NA) {
            int j = i + 1;
            while (j < NA && A[i].x == A[j].x) ++j;
            if (j - i > 10) for (int k = i; k < j; ++k) A[k].x = U64_MAX;
            i = j;
        }
        int n = 0;
        for (i = 0; i < NA; ++i) if (A[i].x != U64_MAX) A[n++] = A[i];
        NA = n;
        km_list.resize(NA);
    }    
}

void x_ddf_scoring(pu64_t* A, int NA, int query_size, const HbnUnpackedDatabase* updb,
    const int kMaxKmerDist,
    const int kMaxGapSize,
    const double kMaxDDF,
    const int kMinDDFScore,
    vector<int>& D_pred_list,
    vector<int>& D_score_list,
    vector<pair<int, int>>& D_head_list,
    vector<u8>& D_avail_list,
    vector<int>& D_seed_idx_list,
    vector<HbnChainInfo>& M_chain_list,
    vector<pu64_t>& M_chain_seed_list)
{
    M_chain_list.clear();
    M_chain_seed_list.clear();
    if (NA == 0) return;

    D_pred_list.clear();
    D_pred_list.resize(NA);
    D_score_list.clear();
    D_score_list.resize(NA);
    int* pred_list = D_pred_list.data();
    int* score_list = D_score_list.data();

    int st = 0;
    for (int i = 0; i < NA; ++i) {
        int max_j = -1;
        int max_f = 1;
        while (st < i && (A[i].x>>32 != A[st].x>>32 || A[i].x > A[st].x + kMaxKmerDist)) ++st;
        int iqb = (int32_t)A[i].y;
        int isb = (int32_t)A[i].x;
	    int j = i - 1;
        int i_span = A[i].y>>32&0xff;
        for (; j >= st; --j) {
            int jqb = (int32_t)A[j].y;
            int jsb = (int32_t)A[j].x;
            hbn_assert(jsb <= isb);
            if (jsb == isb || jqb >= iqb || iqb - jqb > kMaxKmerDist) continue;
            int q_d = iqb - jqb;
            int s_d = isb - jsb;
            int d_d = abs(q_d - s_d); // distance difference
	        if (jqb + i_span > iqb && d_d == 0) {
		        max_f = score_list[j] + 1;
		        max_j = j;
		        break;
	        }
            bool r = (jqb + i_span > iqb || jsb + i_span > isb) && d_d;
            if (r) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kMaxDDF) continue;
            if (d_d > kMaxGapSize) continue;
            int sc = score_list[j] + 1;
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
            } 
        }
        score_list[i] = max_f;
        pred_list[i] = max_j;
    }

    D_avail_list.clear();
    D_avail_list.resize(NA); 
    u8* avail_list = D_avail_list.data();
    fill(avail_list, avail_list + NA, 0);
    for (int i = 0; i < NA; ++i) {
        if (pred_list[i] >= 0) avail_list[pred_list[i]] = 1;
    }
    D_head_list.clear();
    for (int i = 0; i < NA; ++i) {
        if (avail_list[i] == 0 && score_list[i] >= kMinDDFScore) {
            D_head_list.emplace_back(i, score_list[i]);
        }
    }
    pdqsort(D_head_list.begin(), D_head_list.end(),
         [](const pair<int, int>& a, const pair<int, int>& b)->bool {
             return (a.second > b.second) || (a.second == b.second && a.first < b.first);
         });

    const int num_vtx = updb->num_targets() * 2;
    int added = 0;
    for (auto& ias : D_head_list) {
        if (avail_list[ias.first] == 2) continue;
        int p = ias.first;
        D_seed_idx_list.clear();
        while (p >= 0) {
            D_seed_idx_list.push_back(p);
            p = pred_list[p];
        }
        reverse(D_seed_idx_list.begin(), D_seed_idx_list.end());

        const pu64_t* F = A + D_seed_idx_list.front();
        const pu64_t* T = A + D_seed_idx_list.back();
        int f_span = F->y>>32&0xff;
        int qb = (int32_t)F->y + 1 - f_span;
        int qe = (int32_t)T->y + 1;
        int vid = F->x>>32;
        hbn_assert(vid >= 0 && vid < num_vtx);
        int vb = (int32_t)F->x + 1 - f_span;
        int ve = (int32_t)T->x + 1;
        if (1) {
            HbnChainInfo c;
            c.qb = qb;
            c.qe = qe;
            c.sid = vid/2;
            c.sdir = vid&1;
            c.soff = vb;
            c.send = ve;
            c.offset = M_chain_seed_list.size();
            c.sc = D_seed_idx_list.size();
            for (auto i : D_seed_idx_list) M_chain_seed_list.push_back(A[i]);
            M_chain_list.push_back(c);
        }

        for (auto i : D_seed_idx_list) avail_list[i] = 2;
	    if (++added == 10000) break;
    }
    {
        int nc = M_chain_list.size();
        for (int i = 0; i < nc; ++i) {
            HbnChainInfo& hci = M_chain_list[i];
            if (hci.sdir == FWD) continue;
            int x = query_size - hci.qe;
            int y = query_size - hci.qb;
            hci.qb = x;
            hci.qe = y;
            int ssize = updb->target_size(hci.sid);
            x = ssize - hci.send;
            y = ssize - hci.soff;
            hci.soff = x;
            hci.send = y;

            pu64_t* SA = M_chain_seed_list.data() + hci.offset;
            for (int k = 0; k < hci.sc; ++k) {
                int kmer_size = SA[k].y>>32&0xff;
                int qoff = (int32_t)SA[k].y + 1 - kmer_size;
                int soff = (int32_t)SA[k].x + 1 - kmer_size;
                qoff = query_size - 1 - qoff;
                soff = ssize - 1 - soff;
                SA[k].x = ((u64)hci.sid<<33) | ((u64)hci.sdir<<32) | ((u64)soff);
                SA[k].y = ((u64)kmer_size<<32) | ((u64)qoff);
            }
            reverse(SA, SA + hci.sc);

            //cerr << i << '\t' << D_chain_list[i] << '\n';
        }
    }
}

void PrelimSearch::find_init_hits(PoreCQuery* query, HbnIndex* index, HbnChainInfo*& chains, int& num_chains, pu64_t*& seeds)
{
    chains = nullptr;
    num_chains = 0;

    x_extract_kmer_matches(query, index, D_seed_list);
    if (D_seed_list.empty()) return;

    pu64_t* A = D_seed_list.data();
    size_t N = D_seed_list.size();
    x_ddf_scoring(A, N, query->size, index->updb(),
        M_kmer_dist,
        M_gap_size,
        M_ddf,
        M_ddf_score,
        D_pred_list,
        D_score_list,
        D_head_list,
        D_avail_list,
        D_seed_idx_list,
        M_chain_list,
        M_chain_seed_list);
    if (M_chain_list.empty()) return;

    chains = M_chain_list.data();
    num_chains = M_chain_list.size();
    seeds = M_chain_seed_list.data();
}

/////////////////////////

void PrelimSearch::R_init(const u8* fwd_query, const int query_size)
{
    R_hash_table.clear();
    R_kmer_list.clear();

    int pos = 0;
    u32 hash = 0;
    pu32_t qk;
    for (int i = 0; i < query_size; ++i) {
        u32 c = fwd_query[i];
        hash = (hash << 2 | c) & R_hash_mask;
        if (i + 1 - pos == R_kmer_size) {
            qk.x = hash;
            qk.y = pos + R_kmer_size - 1;
	        //if ((pos%4)==0)
            R_kmer_list.push_back(qk);
            ++pos;
        }
    }

    pu32_t* a = R_kmer_list.data();
    size_t c = R_kmer_list.size();
    size_t i = 0;
    while (i < c) {
        size_t j = i + 1;
        while (j < c && a[i].x == a[j].x) ++j;
        //if (j - i > 10) { i = j; continue; }
        u64 cnt = j - i;
        u64 u = lktbl_pack_hash_value(i, cnt);
        R_hash_table.insert(pair<u64, u64>(a[i].x, u));
        i = j;
    }

    M_chain_seed_list.clear();
}

static void x_refined_ddf_scoring2(pu64_t* A, int NA,
    const int kMaxKmerDist,
    const int kMaxGapSize,
    const double kMaxDDF,
    const int kMinDDFScore,
    vector<int>& D_pred_list,
    vector<int>& D_score_list,
    vector<u8>& D_avail_list,
    HbnChainInfo& hit,
    vector<pu64_t>& chain)
{
    chain.clear();
    if (NA == 0) return;

    D_pred_list.clear();
    D_pred_list.resize(NA);
    D_score_list.clear();
    D_score_list.resize(NA);
    int* pred_list = D_pred_list.data();
    int* score_list = D_score_list.data();

    int st = 0;
    for (int i = 0; i < NA; ++i) {
        int max_j = -1;
        int max_f = 1;
        while (st < i && (A[i].x>>32 != A[st].x>>32 || A[i].x > A[st].x + kMaxKmerDist)) ++st;
        int iqb = (int32_t)A[i].y;
        int isb = (int32_t)A[i].x;
	    int j = i - 1;
        int i_span = A[i].y>>32&0xff;
        for (; j >= st; --j) {
            int jqb = (int32_t)A[j].y;
            int jsb = (int32_t)A[j].x;
            hbn_assert(jsb <= isb);
            if (jsb == isb || jqb >= iqb || iqb - jqb > kMaxKmerDist) continue;
            int q_d = iqb - jqb;
            int s_d = isb - jsb;
            int d_d = abs(q_d - s_d); // distance difference
	        if (jqb + i_span > iqb && d_d == 0) {
		        max_f = score_list[j] + 1;
		        max_j = j;
		        break;
	        }
            bool r = (jqb + i_span > iqb || jsb + i_span > isb) && d_d;
            if (r) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kMaxDDF) continue;
            if (d_d > kMaxGapSize) continue;
            int sc = score_list[j] + 1;
            if (sc > max_f) {
                max_f = sc;
                max_j = j;
            } 
        }
        score_list[i] = max_f;
        pred_list[i] = max_j;
    }

    D_avail_list.clear();
    D_avail_list.resize(NA); 
    u8* avail_list = D_avail_list.data();
    fill(avail_list, avail_list + NA, 0);
    for (int i = 0; i < NA; ++i) {
        if (pred_list[i] >= 0) avail_list[pred_list[i]] = 1;
    }
    int max_sc = 0;
    int max_i = -1;
    for (int i = 0; i < NA; ++i) {
        if (avail_list[i] == 0 && score_list[i] >= kMinDDFScore) {
            if (score_list[i] > max_sc) {
                max_sc = score_list[i];
                max_i = i;
            }
        }
    }
    if (!max_sc) return;

    while (max_i >= 0) {
        chain.push_back(A[max_i]);
        max_i = pred_list[max_i];
    }
    reverse(chain.begin(), chain.end());

    const pu64_t* F = &chain.front();
    const pu64_t* T = &chain.back();
    int f_span = F->y>>32&0xff;
    int qb = (int32_t)F->y + 1 - f_span;
    int qe = (int32_t)T->y + 1;
    int vb = (int32_t)F->x + 1 - f_span;
    int ve = (int32_t)T->x + 1;
    hit.qb = qb;
    hit.qe = qe;
    hit.sid = 0;
    hit.sdir = FWD;
    hit.soff = vb;
    hit.send = ve;
    hit.sc = chain.size();
}

bool PrelimSearch::R_find_init_hit(int qb, int qe, const char* subject, int sb, int se, HbnChainInfo& chain, int& km_off, pu64_t*& km_a, int& km_c)
{
    const int min_qoff = qb + R_kmer_size - 1;
    D_seed_list.clear();

    int pos = sb;
    u64 hash = 0;
    pu64_t M;
    for (int i = sb; i < se; ++i) {
        u64 c = subject[i];
        c = nst_nt4_table[c];
        if (c > 3) { hash = 0; pos = i + 1; continue; }
        hash = ((hash << 2) | c) & R_hash_mask;
        if (i + 1 - pos == R_kmer_size) {
            auto iter = R_hash_table.find(hash);
            if (iter != R_hash_table.end()) {
                u64 offset = lktbl_extract_offset(iter->second);
                u64 cnt = lktbl_extract_cnt(iter->second);
                for (u64 s = 0; s < cnt; ++s) {
                    pu32_t& qk = R_kmer_list[offset+s];
                    if (qk.y < min_qoff || qk.y >= qe) continue;
                    M.x = pos + R_kmer_size - 1;
                    M.y = ((u64)R_kmer_size<<32) | qk.y;
                    D_seed_list.push_back(M);
                }
            }
            ++pos;
        }
    }
    if (D_seed_list.empty()) return false;

    pu64_t* A = D_seed_list.data();
    size_t N = D_seed_list.size();
    sort_pu64_x(A, A + N);
    {
        int i = 0;
        while (i < N) {
            int j = i + 1;
            while (j < N && A[i].x == A[j].x) ++j;
            if (j - i > 20) for (int k = i; k < j; ++k) A[k].x = U64_MAX;
            i = j;
        }
        int n = 0;
        for (i = 0; i < N; ++i) if (A[i].x != U64_MAX) A[n++] = A[i];
        N = n;
    }

    x_refined_ddf_scoring2(A, N,
        M_kmer_dist,
        M_gap_size,
        M_ddf,
        M_ddf_score,
        D_pred_list,
        D_score_list,
        D_avail_list,
        chain,
        R_chain);
    if (R_chain.empty()) return false;

    km_off = M_chain_seed_list.size();
    km_c = R_chain.size();
    M_chain_seed_list.insert(M_chain_seed_list.end(), R_chain.begin(), R_chain.end());
    km_a = M_chain_seed_list.data() + km_off;
    return true;
}