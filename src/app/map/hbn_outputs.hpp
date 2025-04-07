#ifndef __HBN_OUTPUTS_HPP
#define __HBN_OUTPUTS_HPP

#include "../../corelib/hbn_aux.h"
#include "../../corelib/unpacked_seqdb.hpp"
#include "../../htslib/sam.h"
#include "map_options.hpp"
#include "pore_c_aux.hpp"
#include "restrict_enzyme_loci_list.hpp"

#include <mutex>

#define SAM_VERSION "1.6"

class HbnOutputs
{
public:
    HbnOutputs(int argc, char* argv[], HbnUnpackedDatabase* reference, const char* output, const EOutputFmt outfmt) {
        M_outfmt = outfmt;
        M_out = nullptr;
        M_sam_hdr = nullptr;
        M_sam_out = nullptr;

        if (outfmt == eOutputFmt_PAF) {
            hbn_fopen(M_out, output, "w");
        } else if (outfmt == eOutputFmt_SAM || outfmt == eOutputFmt_FragSAM) {
            M_sam_out = sam_open(output, "w");
            hts_set_threads(M_sam_out, 8);
            x_init_sam_hdr(argc, argv, reference);
        } else if (outfmt == eOutputFmt_BAM || outfmt == eOutputFmt_FragBAM) {
            M_sam_out = sam_open(output, "wb");
            hts_set_threads(M_sam_out, 8);
            x_init_sam_hdr(argc, argv, reference);
        }
    }

    ~HbnOutputs() {
        if (M_out) hbn_fclose(M_out);
        if (M_sam_hdr) sam_hdr_destroy(M_sam_hdr);
        if (M_sam_out) sam_close(M_sam_out);
    }

    void dump(RestrictEnzymeLociList* reloci_list, 
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList& pca_list);

private:
    void x_init_sam_hdr(int argc, char* argv[], HbnUnpackedDatabase* reference);

private:
    EOutputFmt  M_outfmt;
    FILE*       M_out;
    sam_hdr_t*  M_sam_hdr;
    samFile*    M_sam_out;
    std::mutex  M_out_mutex;
};

#endif // __HBN_OUTPUTS_HPP
