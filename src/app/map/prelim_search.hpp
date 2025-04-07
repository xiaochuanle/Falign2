#ifndef __PRELIM_SEARCH_HPP
#define __PRELIM_SEARCH_HPP

#include "index.hpp"
#include "map_options.hpp"
#include "pore_c_aux.hpp"
#include "query_input.hpp"

#include <vector>

class PrelimSearch
{
public:
    using hash_t = large_hash_type;

    PrelimSearch(HbnMapOptions* options)
    {
        M_kmer_dist = options->kmer_dist;
        M_gap_size = options->gap;
        M_ddf = options->ddf;
        M_ddf_score = options->chain_score;

        R_kmer_size = 15;
        R_hash_mask = (U64_ONE << (2 * R_kmer_size)) - 1;
    }

    void find_init_hits(PoreCQuery* query, HbnIndex* index, HbnChainInfo*& chains, int& num_chains, pu64_t*& seeds);

    std::vector<u8>& cov_list() {
        return D_avail_list;
    }

public:
    void R_init(const u8* fwd_query, const int query_size);

    bool R_find_init_hit(int qb, int qe, const char* subject, int sb, int se, HbnChainInfo& chain, int& km_off, pu64_t*& km_a, int& km_c);

    pu64_t* R_get_chain_seeds(int km_off) 
    {
        return M_chain_seed_list.data() + km_off;
    }

private:
    std::vector<pu64_t>              D_seed_list;

    std::vector<int>                    D_pred_list;
    std::vector<int>                    D_score_list;
    std::vector<std::pair<int, int>>    D_head_list;
    std::vector<u8>                     D_avail_list;
    std::vector<int>                    D_seed_idx_list;
    
    std::vector<HbnChainInfo>           M_chain_list;
    std::vector<pu64_t>              M_chain_seed_list;

private:
    std::vector<pu32_t>                 R_kmer_list;
    std::vector<pu64_t>              R_chain;
    hash_t                              R_hash_table;

    int                                 R_kmer_size;
    u32                                 R_hash_mask;

private:
    int     M_kmer_dist;
    int     M_gap_size;
    double  M_ddf;
    int     M_ddf_score;
};

#endif // __PRELIM_SEARCH_HPP
