#ifndef __PORE_C_AUX_HPP
#define __PORE_C_AUX_HPP

#include "../../corelib/pdqsort.h"
#include "lookup_table.hpp"
#include "restrict_enzyme_loci_list.hpp"

#include <sstream>

inline void set_kmer_match(u64 qoff, u64 sid, u64 sdir, u64 soff, u64 kmer, pu64_t* H)
{
    H->x = (sid << 33) | (sdir << 32) | soff;
    H->y = (kmer << 32) | qoff;
}

inline void sort_pu64_x(pu64_t* beg, pu64_t* end)
{
    pdqsort(beg, end, [](const pu64_t& a, const pu64_t& b) { return a.x < b.x; });
}

inline void sort_pu32_x(pu32_t* beg, pu32_t* end)
{
    pdqsort(beg, end, [](const pu32_t& a, const pu32_t& b) { return a.x < b.x; });
}

struct HbnChainInfo
{
    int qb, qe;
    int sid, sdir, soff, send;
    
    int offset;
    int sc;

    bool is_hom;
    bool is_valid;
    int cnt;
    int id;
    int parent;
    int mapQ;
};

void set_mapQ_for_ddf_chains(HbnChainInfo* a, int c);

inline std::ostream& operator<<(std::ostream& os, const HbnChainInfo& hit)
{
    os << '[' << hit.qb << ", " << hit.qe << ']'
       << " x "
       << '[' << hit.sid << ", " << hit.sdir << ", " << hit.soff << ", " << hit.send << ']'
       << ", " << hit.sc;
    return os;
}

struct PoreCInitHit
{
    int qoff, qend;

    int sid, sdir, soff, send, ssize;

    int ddf_score;
    int map_score;
    double pi;
    int mapQ;
    int ddf_mapQ;
    int map_mapQ;

    int id;
    int parent;
    bool is_hom;

    int qas_offset;
    int sas_offset;
    int as_size;
};

inline std::ostream& operator<<(std::ostream& os, const PoreCInitHit& hit)
{
    os << '[' << hit.qoff << ", " << hit.qend << ']'
       << " x "
       << '[' << hit.sid << ", " << hit.sdir << ", " << hit.soff << ", " << hit.send << ']'
       << ", " << hit.ddf_score << ", " << hit.map_score << ", " << hit.mapQ;
    return os;
}

void set_mapQ_for_init_hits(PoreCInitHit* a, int c);

struct PoreCAlign {
    int qid, qoff, qend, qsize;
    int enzyme_qoff, enzyme_qend;
    int qc;

    int sid, sdir, soff, send, ssize;
    int enzyme_soff, enzyme_send;

    bool lqm, rqm, lsm, rsm;
    bool lp, lpp, rp, rpp;
    bool lcp, lcpp, rcp, rcpp;
    bool is_perfect;
    bool is_pseudo_perfect;
    bool is_homologous;

    bool is_modified;
    int km_off;
    int km_c;

    bool is_hom;
    int map_q;
    double pi;
    int map_score;
    int ddf_score;
    int raw_aln_idx;
};

inline std::ostream& operator<<(std::ostream& os, const PoreCAlign& pca)
{
    os << '[' << pca.qoff << ", " << pca.qend << ']'
       << " x "
       << '[' << pca.sid << ", " << pca.sdir << ", " << pca.soff << ", " << pca.send << ']'
       << ", q-enzyme = [" << pca.enzyme_qoff << ", " << pca.enzyme_qend << ']'
       << ", v-enzyme = [" << pca.enzyme_soff << ", " << pca.enzyme_send << ']'
       << ", pi = " << pca.pi << ", score = " << pca.map_score << ", ddf-score = " << pca.ddf_score;
    return os;
}

inline BOOL set_pca_chain_offset(PoreCAlign* pca, const int enzyme_size)
{
    pca->lqm = (pca->qoff == 0) || (pca->qoff == pca->enzyme_qoff);
    pca->rqm = (pca->qend == pca->qsize) || (pca->qend == pca->enzyme_qend);
    pca->lsm = (pca->soff == 0) || (pca->soff == pca->enzyme_soff);
    pca->rsm = (pca->send == pca->ssize) || (pca->send == pca->enzyme_send);
    
    bool l_perfect = (pca->lqm && pca->lsm) || (pca->qoff == 0) || (pca->soff == 0);
    bool r_perfect = (pca->rqm && pca->rsm) || (pca->qend == pca->qsize) || (pca->send == pca->ssize);
    pca->is_perfect = l_perfect && r_perfect;
    pca->lp = l_perfect;
    pca->rp = r_perfect;

    bool l_pseudo_perfect = l_perfect || pca->lqm || pca->lsm;
    bool r_pseudo_perfect = r_perfect || pca->rqm || pca->rsm;
    pca->is_pseudo_perfect = l_pseudo_perfect && r_pseudo_perfect;
    pca->lpp = l_pseudo_perfect;
    pca->rpp = r_pseudo_perfect;

    pca->qc = pca->qend - pca->qoff;
    if (pca->qoff >= pca->qend) return FALSE;
    if (pca->soff >= pca->send || pca->enzyme_soff > pca->enzyme_send) return FALSE;
    return TRUE;
}

struct TrimPca {
    int is_valid;
    PoreCAlign pca;
    int qas_offset;
    int sas_offset;
    int as_size;
    int frag_id;
};

static inline bool operator < (const TrimPca& lhs, const TrimPca& rhs)
{
    if (lhs.pca.qoff >= rhs.pca.qoff && lhs.pca.qend <= rhs.pca.qend) return true;

    if ((lhs.pca.sid == rhs.pca.sid)
        &&
        (lhs.pca.sdir == rhs.pca.sdir)
        &&
        (lhs.pca.soff >= rhs.pca.soff) && (lhs.pca.send <= rhs.pca.send)) return true;

    return false;
}

struct TrimPcaList
{
    std::vector<TrimPca> tpca_list;
    std::string align_strings;

    void add(PoreCAlign* pca, const char* qas, const char* sas, const int as_size) 
    {
        TrimPca tpca;
        tpca.is_valid = true;
        tpca.pca = *pca;
        tpca.qas_offset = align_strings.size();
        align_strings.append(qas, as_size);
        tpca.sas_offset = align_strings.size();
        align_strings.append(sas, as_size);
        tpca.as_size = as_size;
        tpca_list.push_back(tpca);
    }

    void clear() {
        tpca_list.clear();
        align_strings.clear();
    }
};

#endif // __PORE_C_AUX_HPP
