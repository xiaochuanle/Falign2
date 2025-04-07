#ifndef __QUERY_INPUT_HPP
#define __QUERY_INPUT_HPP

#include "../../corelib/fasta.hpp"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "restrict_enzyme_loci_list.hpp"

#include <mutex>

struct PoreCQuery
{
public:
    PoreCQuery(RestrictEnzyme* re): enzyme(re) {
        ks_initialize(&I_hdr);
        ks_initialize(&I_seq);
        ks_initialize(&I_qv);
    }

    ~PoreCQuery() {
        ks_free(&I_hdr);
        ks_free(&I_seq);
        ks_free(&I_qv);
    }

    void trim_end_spaces(kstring_t* line) {
        while (line->l && isspace(line->s[line->l-1])) --line->l;
    }

    void init()
    {
        trim_end_spaces(&I_hdr);
        trim_end_spaces(&I_seq);
        trim_end_spaces(&I_qv);

        M_name.clear();
        for (size_t i = 1; i < I_hdr.l; ++i) if (isspace(I_hdr.s[i])) break; else M_name += I_hdr.s[i];
        name = M_name.c_str();
        size = I_seq.l;

        size = I_seq.l;
        M_fwd_qs.clear();
        M_fwd_qs.resize(size);
        M_rev_rqs.clear();
        M_rev_rqs.resize(size);
        M_rev_qs.clear();
        M_rev_qs.resize(size);
        for (int fi = 0, ri = size - 1; fi < size; ++fi, --ri) {
            I_seq.s[fi] = toupper(I_seq.s[fi]);
            int fc = I_seq.s[fi];
            fc = nst_nt4_table[fc];
            if (fc > 3) fc = 0;
            M_fwd_qs[fi] = fc;

            int rc = 3 - fc;
            M_rev_qs[ri] = rc;
            rc = DECODE_RESIDUE(rc);
            M_rev_rqs[ri] = rc;
        }
        fwd_rqs = I_seq.s;
        fwd_qs = M_fwd_qs.data();
        rev_rqs = M_rev_rqs.data();
        rev_qs = M_rev_qs.data();

        fwd_qv = nullptr;
        rev_qv = nullptr;
        if (I_qv.l) {
            fwd_qv = I_qv.s;
	        for (int i = 0; i < size; ++i) I_qv.s[i] -= 33;
            M_rev_qv.assign(fwd_qv, fwd_qv + size);
            reverse(M_rev_qv.begin(), M_rev_qv.end());
            rev_qv = M_rev_qv.data();
        }

        M_fwd_enzyme_pos.clear();
        M_rev_enzyme_pos.clear();
        if (!enzyme) return;
        extract_enzyme_loci_list_for_one_seq(fwd_rqs, size, enzyme, M_fwd_enzyme_pos);
        fwd_enzyme_pos = M_fwd_enzyme_pos.data();
        fwd_enzyme_pos_cnt = M_fwd_enzyme_pos.size();

        extract_enzyme_loci_list_for_one_seq(rev_rqs, size, enzyme, M_rev_enzyme_pos);
        rev_enzyme_pos = M_rev_enzyme_pos.data();
        rev_enzyme_pos_cnt = M_rev_enzyme_pos.size();
    }

    RestrictEnzyme* enzyme;
    int id;
    const char* name;
    int size;
    const char* fwd_rqs;
    const char* rev_rqs;
    const u8* fwd_qs;
    const u8* rev_qs;
    const char* fwd_qv;
    const char* rev_qv;
    int* fwd_enzyme_pos;
    int fwd_enzyme_pos_cnt;
    int* rev_enzyme_pos;
    int rev_enzyme_pos_cnt;

private:
    std::vector<u8>     M_fwd_qs;
    std::vector<u8>     M_rev_qs;
    std::string         M_rev_rqs;
    std::string         M_rev_qv;
    std::vector<int>    M_fwd_enzyme_pos;
    std::vector<int>    M_rev_enzyme_pos;

public:
    std::string         M_name;

    kstring_t           I_hdr;
    kstring_t           I_seq;
    kstring_t           I_qv;
};

class QueryReader
{

public:
    QueryReader(std::vector<std::string>& query_fn_list)
    {
        M_query_fn_list = query_fn_list;
        M_num_fn = query_fn_list.size();
        M_fn_idx = 0;
        M_in = nullptr;

        M_mapped_bases = 0;
        M_mapped_queries = 0;

        B_num_queries = 0;
        B_query_id = 0;
    }

    QueryReader(const char* fn, size_t batch_bases)
    {
        M_query_fn_list.push_back(fn);
        M_num_fn = 1;
        M_fn_idx = 0;
        M_in = nullptr;

        M_mapped_bases = 0;
        M_mapped_queries = 0;

        B_num_queries = 0;
        B_query_id = 0;
    }

    ~QueryReader()
    {
        delete M_in;
    }

    bool get_next_query_file()
    {
        delete M_in;
        M_in = nullptr;
        if (M_fn_idx >= M_num_fn) return false;

        M_in = new HbnLineReader(M_query_fn_list[M_fn_idx].c_str());
        ++M_fn_idx;
        return true;
    }

    void dump_mapped_stats()
    {
        std::string size = NStr::UInt8ToString_DataSize(M_mapped_bases);
        HBN_LOG("%10zu queries (%s) mapped", M_mapped_queries, size.c_str());
    }

    bool get_next_batch_query(PoreCQuery* query) 
    {
        query->I_hdr.l = 0;
        query->I_seq.l = 0;
        query->I_qv.l = 0;
        std::lock_guard<std::mutex> __(M_query_mutex);

        if (B_query_id < B_num_queries) {
            BatchQueryInfo& bqi = B_query_list[B_query_id];
            const char* h = B_hdr_list.data() + bqi.hdr_offset;
            kputsn(h, strlen(h), &query->I_hdr); 
            const char* s = B_seq_list.data() + bqi.seq_offset;
            kputsn(s, bqi.seq_size, &query->I_seq);
            if (!B_qual_list.empty()) {
                s = B_qual_list.data() + bqi.seq_offset;
                kputsn(s, bqi.seq_size, &query->I_qv);
            }
            //fprintf(stderr, "%zu\t%zu\t%zu\n", query->M_hdr.l, query->M_fwd_raw_seq.l, query->M_fwd_qv.l);

            query->id = B_query_id;
            ++B_query_id;
            return true;
        }
        return false;  
    }

    bool get_next_query(PoreCQuery* query)
    {
        query->I_hdr.l = 0;
        query->I_seq.l = 0;
        query->I_qv.l = 0;
        std::lock_guard<std::mutex> __(M_query_mutex);

        if (B_query_id < B_num_queries) {
            BatchQueryInfo& bqi = B_query_list[B_query_id];
            const char* h = B_hdr_list.data() + bqi.hdr_offset;
            kputsn(h, strlen(h), &query->I_hdr); 
            const char* s = B_seq_list.data() + bqi.seq_offset;
            kputsn(s, bqi.seq_size, &query->I_seq);
            if (!B_qual_list.empty()) {
                s = B_qual_list.data() + bqi.seq_offset;
                kputsn(s, bqi.seq_size, &query->I_qv);
            }

            ++B_query_id;
            query->id = M_mapped_queries++;
            M_mapped_bases += bqi.seq_size;
            if ((M_mapped_queries % 100000) == 0) dump_mapped_stats();
            return true;
        }

        if (!M_in->ReadOneLine(&query->I_hdr)) return false;
        bool r = M_in->ReadOneLine(&query->I_seq);
        if (!r) HBN_ERR("Unexpected end of query file");
        hbn_assert(query->I_hdr.s[0] == '@' || query->I_hdr.s[0] == '>');
        if (query->I_hdr.s[0] == '@') {
            r = M_in->ReadOneLine(&query->I_qv);
            if (!r) HBN_ERR("Unexpected end of query file");
            hbn_assert(query->I_qv.s[0] == '+');
            r = M_in->ReadOneLine(&query->I_qv);
            if (!r) HBN_ERR("Unexpected end of query file");
        }
        query->id = M_mapped_queries++;
        M_mapped_bases += query->I_seq.l;
        if ((M_mapped_queries % 100000) == 0) dump_mapped_stats();
        return true;
    }

    void load_query_batch(const size_t num_bases)
    {
        B_query_list.clear();
        B_hdr_list.clear();
        B_seq_list.clear();
        B_qual_list.clear();

        kstring_t hdr; ks_initialize(&hdr);
        kstring_t seq; ks_initialize(&seq);
        kstring_t plus; ks_initialize(&plus);
        kstring_t qual; ks_initialize(&qual);

        B_num_queries = 0;
        size_t loaded_bases = 0;
        BatchQueryInfo bqi;
        while (1) {
            if (!M_in->ReadOneLine(&hdr)) break;
            bool r = M_in->ReadOneLine(&seq);
            if (!r) HBN_ERR("Unexpected end of query file");
            hbn_assert(hdr.s[0] == '@' || hdr.s[0] == '>');
            if (hdr.s[0] == '@') {
                r = M_in->ReadOneLine(&plus) && M_in->ReadOneLine(&qual);
                if (!r) HBN_ERR("Unexpected end of query file");
                hbn_assert(plus.s[0] == '+');
                hbn_assert(seq.l == qual.l);
            }

            bqi.hdr_offset = B_hdr_list.size();
            B_hdr_list.insert(B_hdr_list.end(), hdr.s, hdr.s + hdr.l);
            B_hdr_list.push_back('\0');
            bqi.seq_offset = B_seq_list.size();
            bqi.seq_size = seq.l;
            B_seq_list.insert(B_seq_list.end(), seq.s, seq.s + seq.l);
            if (qual.l) B_qual_list.insert(B_qual_list.end(), qual.s, qual.s + qual.l);
            B_query_list.push_back(bqi);

            ++B_num_queries;
            loaded_bases += seq.l;
            if (loaded_bases >= num_bases) break;
        }

        ks_free(&hdr);
        ks_free(&seq);
        ks_free(&plus);
        ks_free(&qual);

        std::string size = NStr::UInt8ToString_DataSize(loaded_bases);
        HBN_LOG("Load %d (%s) queries", B_num_queries, size.c_str());
    }

    void reset_batch_query_id()
    {
        B_query_id = 0;
    }

    const char* get_curr_query_fn() const {
        return M_fn_idx ? M_query_fn_list[M_fn_idx-1].c_str() : nullptr;
    }

private:

    struct BatchQueryInfo
    {
        size_t hdr_offset;
        size_t seq_offset;
        int seq_size;
    };

private:

    std::vector<std::string>    M_query_fn_list;
    size_t                      M_fn_idx;
    size_t                      M_num_fn;

    size_t                      M_mapped_queries;
    size_t                      M_mapped_bases;
    HbnLineReader*              M_in;
    std::mutex                  M_query_mutex;

    std::vector<BatchQueryInfo> B_query_list; 
    int                         B_num_queries;
    int                         B_query_id;
    std::mutex                  B_query_id_mutex;
    std::vector<char>           B_hdr_list;
    std::vector<char>           B_seq_list;
    std::vector<char>           B_qual_list;
};

#endif // __QUERY_INPUT_HPP
