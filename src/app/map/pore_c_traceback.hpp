#ifndef __PORE_C_TRACEBACK_HPP
#define __PORE_C_TRACEBACK_HPP

#include "../../sw/edlib_wrapper.hpp"
#include "../../sw/ksw2_wrapper.hpp"
#include "../../sw/small_edlib_align.hpp"
#include "lookup_table.hpp"
#include "pore_c_aux.hpp"

#define ALIGN_END_MATCH     8
#define ALIGN_EXT_SEG       100
#define ALIGN_MAX_EXT_SEG   150

int run_nw_dist(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    int tolerence,
    int* score);

void
run_nw(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    std::string& qaln,
    std::string& saln);

int run_shw(small_edlib_align_struct* small_edlib,
    EdlibAlignData* edlib,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    int* _qend,
    int* _send,
    std::string* qaln,
    std::string* saln);

struct HbnTracebackData
{
public:
    int qdir, qoff, qend, qsize;
    int soff, send, ssize;
    int score, match;
    float pi, epi;

    std::string qabuf;
    std::string sabuf;
    std::string ext_qabuf;
    std::string ext_sabuf;
    std::vector<u8> qfrag;
    std::vector<u8> sfrag;
    std::vector<u64> cigar;

    std::vector<char> vqas;
    std::vector<char> vsas;
    char* qas;
    char* qae;
    char* sas;
    char* sae;

    EdlibAlignData* edlib;
    Ksw2Data* ksw;
    small_edlib_align_struct* small_edlib;

public:
    large_hash_type                 M_lktbl;
    std::vector<KmerHashAndOffset>  M_query_khao_list;
    
    std::vector<pu64_t>            M_query_kmer_list;
    std::vector<pu64_t>            M_subject_kmer_list;
    
    std::vector<pu64_t>            M_seed_list;
    std::vector<pu64_t>            M_chain;

public:

    HbnTracebackData() {
        edlib = new EdlibAlignData();
        ksw = new Ksw2Data();
        small_edlib = small_edlib_align_struct_new();
    }

    ~HbnTracebackData() {
        small_edlib_align_struct_free(small_edlib);
        delete edlib;
        delete ksw;
    }

    bool map_one_chain(const u8* EQ, int qs, const char* S, int ss, const pu64_t* a, int n_a, bool xxx);

    bool map_one_chain(const u8* EQ, int qb, int qe, int qs, const char* S, int sb, int se, int ss);

    bool refined_map_one_chain(const u8* fwd_qs, int iqb, int iqe, const char* S, int isb, int ise, const pu64_t* a, int n_a);

    bool refined_map_one_chain_edlib(const u8* fwd_qs, int iqb, int iqe, 
        const char* S, int isb, int ise, const pu64_t* a, int n_a);

    void dump(FILE* out) {
        fprintf(out, "[%d, %d] x [%d, %d], score = %d, pi = %g, epi = %g\n", qoff, qend, soff, send, score, pi, epi);
    }

    void cigar_to_align_string(const char* FQ, const char* S);
};

#endif // __UL_TRACEBACK_HPP
