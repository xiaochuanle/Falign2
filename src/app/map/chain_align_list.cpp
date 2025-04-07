#include "chain_align_list.hpp"

using namespace std;

#include <algorithm>
#include <array>
#include <vector>

using namespace std;

#define PC_MAX_ALIGN_OVLP 20
#define kMaxMissedCovBase 200

const char* get_chain_type_name(EChainType type)
{
    static const char* name[] = { "perfect", "complete", "complete", "incomplete" };
    return name[type];
}

static BOOL 
s_is_perfect_chain_offset(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    if (al->is_perfect == false || ar->is_perfect == false) return FALSE;
    if (al->qend != ar->qoff) return FALSE;
    return TRUE;
}

static BOOL
s_is_validate_soff_relation_general(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    bool r = (al->sid == ar->sid && al->sdir == ar->sdir);
    if (!r) return TRUE;

    if (al->soff <= ar->soff) {
        if (al->send > ar->send) {
            if (al->send - ar->soff > max_align_ovlp) {
                return FALSE;
            }
        }
    } else {
        if (ar->send > al->soff) {
            if (ar->send - al->soff > max_align_ovlp) {
                return FALSE;
            }
        }
    }

    return TRUE;
}

static BOOL
s_is_pseudo_perfect_chain_offset(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    hbn_assert(al->qoff <= ar->qoff);
    if (ar->qend - al->qend <= 5 || ar->qoff - al->qoff <= 5) return FALSE;
    if (al->is_pseudo_perfect == false || ar->is_pseudo_perfect == false) return FALSE;

    if (al->rcp && ar->lcp && al->qend != ar->qoff) return FALSE;

    if (al->rcp && al->qend > ar->qoff) return FALSE;

    if (ar->lcp && al->qend > ar->qoff) return FALSE;

    if (al->qend > ar->qoff && al->qend - ar->qoff > max_align_ovlp) {
        return FALSE;
    }

    return s_is_validate_soff_relation_general(al, ar, max_align_ovlp);
}

BOOL 
is_validate_chain_offset_general(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    hbn_assert(al->qoff <= ar->qoff);
    if (ar->qend - al->qend <= 5 || ar->qoff - al->qoff <= 5) return FALSE;
    
    if (al->rcp && ar->lcp) {
        if (al->qend > ar->qoff) return FALSE;
        if (al->sid == ar->sid && al->sdir == ar->sdir) {
            if (al->soff <= ar->soff && al->send > ar->soff) return FALSE;
            if (ar->soff <= al->soff && ar->send > al->soff) return FALSE;
        }
        return TRUE;
    }
    
    int l_ext = ar->qoff - al->qoff;
    int r_ext = ar->qend - al->qend;
    if (l_ext < 50 && r_ext < 50 &&  al->qend > ar->qoff && al->qend - ar->qoff > max_align_ovlp) {
        return FALSE;
    }

    //return s_is_validate_soff_relation_general(al, ar, max_align_ovlp);
    return TRUE;
}

static int
s_pca_list_relation(vector<int>& a_idx_list, vector<int>& b_idx_list)
{
    int* l_a = NULL;
    int l_c = 0;
    int* s_a = NULL;
    int s_c = 0;
    bool reverse_relation = false;
    if (a_idx_list.size() >= b_idx_list.size()) {
        l_a = a_idx_list.data();
        l_c = a_idx_list.size();
        s_a = b_idx_list.data();
        s_c = b_idx_list.size();
    } else {
        l_a = b_idx_list.data();
        l_c = b_idx_list.size();
        s_a = a_idx_list.data();
        s_c = a_idx_list.size();
        reverse_relation = true;
    }

    int l_i = 0, s_i = 0, n_eq = 0;
    while (l_i < l_c && s_i < s_c) {
        while (l_i < l_c && l_a[l_i] < s_a[s_i]) ++l_i;
        if (l_i >= l_c) break;
        while (s_i < s_c && s_a[s_i] < l_a[l_i]) ++s_i;
        if (s_i >= s_c) break;
        if (l_a[l_i] == s_a[s_i]) {
            ++n_eq;
            ++l_i;
            ++s_i;
        }
    }

    if (n_eq != s_c) return 0;
    if (reverse_relation) return -1;
    return 1;
}

bool s_is_complete_chain(const int* vdfa, const int vdfc, const PoreCAlign* pca_a, const int pca_c)
{
    if (vdfc <= 2) return true;

    if (pca_a[0].qoff > 50) return false;
    /// left gap
    {
        int l_gap = 0;
        for (int i = 1; i < vdfc; ++i) {
            if (vdfa[i] > pca_a[0].qoff) break;
            l_gap = vdfa[i];
        }
        //HBN_LOG("l_gap = %d", l_gap);
        if (l_gap > 50) return false;
    }

    if (pca_a[pca_c-1].qsize - pca_a[pca_c-1].qend > 50) return false;
    /// right gap
    {
        int qsize = pca_a[0].qsize;
        int r_gap = 0;
        for (int i = vdfc - 3; i >= 0; --i) {
            if (vdfa[i] < pca_a[pca_c-1].qend) break;
            r_gap = qsize - vdfa[i];
        }
        //HBN_LOG("r_gap = %d", r_gap);
        if (r_gap > 50) return false;
    }

    int gap_size = pca_a[0].qoff;
    for (int i = 0; i < pca_c - 1; ++i) {
        int gap = pca_a[i+1].qoff - pca_a[i].qend;
        if (gap > 50) return false;
        if (gap > 0) gap_size += gap;
    }
    gap_size += (pca_a[pca_c-1].qsize - pca_a[pca_c-1].qend);

    return (gap_size < kMaxMissedCovBase);
}

bool s_is_perfect_chain(PoreCAlign* pca_a, int pca_c)
{
    for (int i = 0; i < pca_c - 1; ++i) {
        if (!s_is_perfect_chain_offset(pca_a + i, pca_a + i + 1, PC_MAX_ALIGN_OVLP)) return FALSE;
    }
    return TRUE;
}

bool s_is_pseudo_perfect_chain(PoreCAlign* pca_a, int pca_c)
{
    for (int i = 0; i < pca_c - 1; ++i) {
        if (!s_is_pseudo_perfect_chain_offset(pca_a + i, pca_a + i + 1, PC_MAX_ALIGN_OVLP)) return FALSE;
    }
    return TRUE;
}

EChainType pca_chain_type(const int* vdfa, const int vdfc, PoreCAlign* pca_a, int pca_c)
{
    if (s_is_perfect_chain(pca_a, pca_c)) return ePerfectChain;
    if (s_is_pseudo_perfect_chain(pca_a, pca_c)) return ePseudoPerfectChain;
    if (s_is_complete_chain(vdfa, vdfc, pca_a, pca_c)) return eCompleteChain;
    return eMaxCovChain;
}

static inline int calc_gap_penalty(int gap_size)
{
    if (gap_size == 0) return 0;
#if 1
    int go1 = 4;
    int ge1 = 2;
    int go2 = 24;
    int ge2 = 1; 
#else
	int go1 = 5;
	int ge1 = 4;
	int go2 = 56;
	int ge2 = 1;
#endif
    return go1 + ge1 * gap_size;
    //return min(go1 + ge1 * gap_size, go2 + ge2 * gap_size);
}

static bool 
s_is_valid_edge(PoreCAlign* pi, PoreCAlign* pj)
{
    hbn_assert(pi->qoff <= pj->qoff);
    if (pj->qoff - pi->qoff < 50) return false;
    if (pj->qend - pi->qend < 50) return false;
    return true;
}

bool shortest_path_max_cov(PoreCAlign* pca_array, int pca_count, PoreCAlignChainData& data)
{
    vector<PoreCAlign>& chain = data.chain;
	chain.clear();
	if (pca_count == 0) return false;
	if (pca_count == 1) {
		chain.push_back(pca_array[0]);
		return true;
	}
    const int kMaxDist = 100000000;
    int num_nodes = pca_count + 2;
    data.dist_list.resize(num_nodes);
    int* dist = data.dist_list.data();
    fill(dist, dist + num_nodes, kMaxDist);
    dist[0] = 0;
    //const int min_aln_score = -100;//10;

    vector<GraphEdge>& edges = data.edge_list;
    edges.clear();
    for (int i = 0; i < pca_count; ++i) {
        //if (pca_array[i].map_score < min_aln_score) continue;
        GraphEdge e;
        e.start = 0;
        e.end = i+1;
        e.weight = calc_gap_penalty(pca_array[i].qoff) - pca_array[i].map_score;
        edges.push_back(e);
        //fprintf(stderr, "++++ add start edge, weight = %d\n", e.weight); dump_chain_pca(fprintf, stderr, pca_array[i], i);
    }
    for (int i = 0; i < pca_count; ++i) {
        //if (pca_array[i].map_score < min_aln_score) continue;
        GraphEdge e;
        e.start = i+1;
        e.end = pca_count + 1;
        e.weight = calc_gap_penalty(pca_array[i].qsize - pca_array[i].qend);
        edges.push_back(e);
        //fprintf(stderr, "---- add end edge, weight = %d\n", e.weight); dump_chain_pca(fprintf, stderr, pca_array[i], i);
    }
    for (int i = 0; i < pca_count; ++i) {
        PoreCAlign* pi = pca_array + i;
        //if (pca_array[i].map_score < min_aln_score) continue;
        for (int j = i + 1; j < pca_count; ++j) {
            //if (pca_array[j].map_score < min_aln_score) continue;
            PoreCAlign* pj = pca_array + j;
            hbn_assert(pi->qoff <= pj->qoff);
            if (!s_is_valid_edge(pi, pj)) continue;
            int weight = 0;
            weight = calc_gap_penalty(abs(pj->qoff - pi->qend)) - pj->map_score;
            GraphEdge e;
            e.start = i + 1;
            e.end = j + 1;
            e.weight = weight;
            edges.push_back(e);
            //fprintf(stderr, "Add edge, weight = %d\n", weight); dump_chain_pca(fprintf, stderr, *pi, i); dump_chain_pca(fprintf, stderr, *pj, j);
        }
    }

    const int n_edge = edges.size();
    data.pred_list.resize(num_nodes);
    int* pred_list = data.pred_list.data();
    fill(pred_list, pred_list + num_nodes, kMaxDist);
    pred_list[0] = 0;

    for (int i = 0; i < num_nodes - 1; ++i) {
        int has_update = 0;
        for (int j = 0; j < n_edge; ++j) {
            if (dist[edges[j].end] > dist[edges[j].start] + edges[j].weight) {
                dist[edges[j].end] = dist[edges[j].start] + edges[j].weight;
                has_update = 1;
                pred_list[edges[j].end] = edges[j].start;
            }
        }
        if (has_update == 0) break;
    }
    if (dist[num_nodes-1] == kMaxDist) return 0;

    int p = pred_list[pca_count+1];
    while (p != 0 && p != kMaxDist) {
        chain.push_back(pca_array[p-1]);
        p = pred_list[p];
    }

    reverse(chain.begin(), chain.end());
    //fprintf(stderr, "***************** find perfect path:\n");
    //for (auto& pca : chain) dump_chain_pca(fprintf, stderr, pca, -1);
    return true;
}

bool select_pca_chain( PoreCAlign* pca_array, int pca_count, PoreCAlignChainData& data)
{
	data.chain.clear();
    if (pca_count == 0) return 0;
    sort(pca_array, pca_array + pca_count, [](const PoreCAlign& a, const PoreCAlign& b)->bool { return a.qoff < b.qoff; });

#if 0
    HBN_LOG("chaining pca list");
    for (int i = 0; i < pca_count; ++i) {
	    cerr << i << '\t' << pca_array[i] << '\n';
    }
#endif

    return shortest_path_max_cov(pca_array, pca_count, data);
}
