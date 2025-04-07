#ifndef __INDEX_HPP
#define __INDEX_HPP

#include "../../corelib/unpacked_seqdb.hpp"
#include "index_options.hpp"
#include "lookup_table.hpp"
#include "restrict_enzyme_loci_list.hpp"

int build_index(int argc, char* argv[]);

class HbnIndex
{
public:
    HbnIndex()
    {
        M_updb = nullptr; 
        
        M_kmer_size = 0;
        M_hash_table = nullptr;
        M_kmer_list = nullptr;

        M_enzyme_pos_list = nullptr;
    }

    ~HbnIndex()
    {
        delete M_updb; 

        delete M_hash_table;
        delete[] M_kmer_list;

        if (M_enzyme_pos_list) RestrictEnzymeLociListFree(M_enzyme_pos_list);
    }

    HbnUnpackedDatabase* updb() {
        return M_updb;
    }

    inline int kmer_size() const {
        return M_kmer_size;
    }

    inline void 
    extract_offset_list(const u64 hash, KmerInfo** ol, u64* cnt) const {
        *ol = nullptr;
        *cnt = 0;
        u64 u = M_hash_table->extract_value(hash);
        if (!u) return;
        u64 offset = lktbl_extract_offset(u);
        *cnt = lktbl_extract_cnt(u);
        *ol = M_kmer_list + offset;
    }

    inline RestrictEnzymeLociList* enzyme_pos_list() {
        return M_enzyme_pos_list;
    }

    void load(const char* path);

    void init_enzyme(const char* enzyme, const int num_threads)
    {
        HBN_LOG("Build enzyme loci list");
        if (M_enzyme_pos_list) M_enzyme_pos_list = RestrictEnzymeLociListFree(M_enzyme_pos_list);
        M_enzyme_pos_list = RestrictEnzymeLociListNew(M_updb, enzyme, num_threads);
        HBN_LOG("Done");
    }

private:
    HbnUnpackedDatabase*    M_updb;
    
    int                     M_kmer_size;
    HbnHashTable*           M_hash_table;
    KmerInfo*               M_kmer_list;
    size_t                  M_kmer_list_size;

    RestrictEnzymeLociList* M_enzyme_pos_list;
};

#endif // __INDEX_HPP
