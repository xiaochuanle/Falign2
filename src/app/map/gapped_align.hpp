#ifndef __GAPPED_ALIGN_HPP
#define __GAPPED_ALIGN_HPP

#include "index.hpp"
#include "pore_c_aux.hpp"
#include "pore_c_traceback.hpp"
#include "prelim_search.hpp"
#include "query_input.hpp"

void extend_ddf_chain_list(HbnTracebackData* tbck_data,
    const HbnUnpackedDatabase* updb,
    PoreCQuery* query,
    HbnChainInfo* ac_a,
    int ac_n,
    pu64_t* ac_km,
    std::vector<u8>& cov_list,
    std::vector<PoreCInitHit>& align_list,
    std::string& align_strings);

class PoreCGapAlign
{
public:
    PoreCGapAlign()
    {
        M_tbck_data = new HbnTracebackData();
    }

    ~PoreCGapAlign()
    {
        delete M_tbck_data;
    }

    void align_ddf_chains(PrelimSearch* prelim,
        PoreCQuery* query,
        HbnUnpackedDatabase* updb,
        RestrictEnzymeLociList* enzyme,
        HbnChainInfo* ac_a,
        int ac_n,
        pu64_t* ac_km);

    std::vector<PoreCAlign>& get_pca_list()
    {
        return M_pca_list;
    }

    std::vector<PoreCInitHit>& get_init_align_list()
    {
        return M_align_list;
    }

    HbnTracebackData* tbck_data()
    {
        return M_tbck_data;
    }

    u8* cov_list() { return M_cov_list.data(); }

private:
    HbnTracebackData*         M_tbck_data;
    std::vector<u8>             M_cov_list;

    std::vector<PoreCInitHit>   M_align_list;
    std::string                 M_align_strings;

    std::vector<PoreCAlign>     M_pca_list;
};

#endif // __GAPPED_ALIGN_HPP
