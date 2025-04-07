#include "align_one_read.hpp"

#include "smooth_pca_list.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

using namespace std;

bool select_pca_chain( PoreCAlign* pca_array, int pca_count, std::vector<PoreCAlign>& chain);

void
align_one_read(HbnIndex* index,
    PoreCQuery* query,
    HbnMapOptions* options,
    PrelimSearch* prelim_search,
    PoreCGapAlign* gapped_search,
    PoreCAlignChainData* pca_chain_data,
    TrimPcaList* trim_pca_list)
{
    trim_pca_list->clear();

    HbnChainInfo* ac_a = nullptr;
    int ac_n = 0;
    pu64_t* ac_km = nullptr;
    prelim_search->find_init_hits(query, index, ac_a, ac_n, ac_km);
    if (!ac_n) return;
    HbnUnpackedDatabase* updb = index->updb();
    gapped_search->align_ddf_chains(prelim_search, query, updb, index->enzyme_pos_list(), ac_a, ac_n, ac_km);

    EChainType chain_type;
    vector<PoreCAlign>& all_pca_list = gapped_search->get_pca_list();
    const int* vdfa = query->fwd_enzyme_pos;
    const int vdfc = query->fwd_enzyme_pos_cnt;
    PoreCAlign* pa = all_pca_list.data();
    int pc = all_pca_list.size();

    select_pca_chain(pa, pc, *pca_chain_data);
    vector<PoreCAlign>& chain = pca_chain_data->chain;
    if (chain.empty()) return;
    chain_type = pca_chain_type(vdfa, vdfc, chain.data(), chain.size());

    RestrictEnzymeLociList* updb_enzyme = index->enzyme_pos_list();

    smooth_pca_list(chain, chain_type, all_pca_list, updb, query, query->enzyme->enzyme_size, gapped_search->tbck_data());

    trim_overlap_subseqs(gapped_search->tbck_data(),
        prelim_search,
        updb,
        index->enzyme_pos_list(),
        query,
        all_pca_list.data(),
        all_pca_list.size(),
        chain.data(),
        chain.size(),
        chain_type,
        *trim_pca_list);

    for (auto& tpca : trim_pca_list->tpca_list) {
        auto& pca = tpca.pca;
        int* ea = (pca.sdir == FWD) ? query->fwd_enzyme_pos : query->rev_enzyme_pos;
        int ec = (pca.sdir == FWD) ? query->fwd_enzyme_pos_cnt : query->rev_enzyme_pos_cnt;
        pca.enzyme_qoff = get_enzyme_pos(ea, ec, pca.qoff);
        pca.enzyme_qend = get_enzyme_pos(ea, ec, pca.qend);
        ea = updb_enzyme->reloci_array + updb_enzyme->seq_reloci_info_array[pca.sid*2].enzyme_loci_offset;
        ec = updb_enzyme->seq_reloci_info_array[pca.sid*2].enzyme_loci_cnt;
        pca.enzyme_soff = get_enzyme_pos(ea, ec, pca.soff);
        pca.enzyme_send = get_enzyme_pos(ea, ec, pca.send);
    }
}
