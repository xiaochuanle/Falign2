#ifndef __LOOKUP_TABLE_HPP
#define __LOOKUP_TABLE_HPP

#include "../../corelib/sparse_map.h"
#include "../../corelib/unpacked_seqdb.hpp"

#define LKTBL_OFFSET_BIT    48
#define LKTBL_CNT_BIT   16

#define HASH_BUCKETS 16
#define HASH_BUCKET_BITS 4 // 2^4 = 16

inline u64 lktbl_pack_hash_value(u64 offset, u64 cnt)
{
    return (offset << LKTBL_CNT_BIT) | cnt;
}

inline u64 lktbl_extract_offset(u64 u) 
{
    return u >> LKTBL_CNT_BIT;
}

inline u64 lktbl_extract_cnt(u64 u) 
{
    u64 mask = 1;
    mask = (mask << LKTBL_CNT_BIT) - 1;
    return u & mask;
}

struct BobHenkin_u64tou64
{
    u64 operator()(uint64_t key) const {
            key = (~key) + (key << 21); // key = (key << 21) - key - 1; 
            key = key ^ (key >> 24); 
            key = (key + (key << 3)) + (key << 8); // key * 265 
            key = key ^ (key >> 14); 
            key = (key + (key << 2)) + (key << 4); // key * 21 
            key = key ^ (key >> 28); 
            key = key + (key << 31); 
            return key;
    }
};

using large_hash_type = tsl::sparse_map<u64, u64, BobHenkin_u64tou64>;

struct KmerInfo
{
    int seq_id;
    int seq_offset;
};

struct KmerHashAndOffset
{
    u64 hash;
    KmerInfo ki;

    bool operator<(const KmerHashAndOffset& rhs) const {
        return hash < rhs.hash;
    }
};

class HbnHashTable
{
public:
    HbnHashTable(const int kmer_size) {
        hbn_assert(kmer_size > 0 && kmer_size <= kMaxKmerSize, "kmer = %d", kmer_size);
        M_max_hash_value = U64_ONE << (2 * kmer_size);
        if (kmer_size <= kLargeTableKmerCutoff) {
            HBN_LOG("Use small hash table for k = %d", kmer_size);
            M_small_hash = (u64*)calloc(M_max_hash_value, sizeof(u64));
            M_large_hash = nullptr;   
        } else {
            HBN_LOG("Use large hash table for k = %d", kmer_size);
            M_small_hash = nullptr;
            M_large_hash = new large_hash_type;
        }
    }

    ~HbnHashTable() {
        if (M_small_hash) free(M_small_hash);
        if (M_large_hash) delete M_large_hash;
    }

    u64 extract_value(const u64 hash) const {
        hbn_assert(hash < M_max_hash_value);
        if (M_small_hash) {
            return M_small_hash[hash];
        } else {
            auto pos = M_large_hash->find(hash);
            return (pos == M_large_hash->end()) ? 0 : pos->second;
        }
    }

    void add_one_hash_and_value(const u64 hash, const u64 value) {
        hbn_assert(hash < M_max_hash_value);
        if (M_small_hash) {
            M_small_hash[hash] = value;
        } else {
            M_large_hash->insert(std::pair<u64, u64>(hash, value));
        }
    }

    void reserve(const u64 num_hash) {
        if (M_large_hash) M_large_hash->reserve(num_hash);
    }

private:
    static const int kLargeTableKmerCutoff = 15;
    static const int kMaxKmerSize = 31;

    u64                 M_max_hash_value;
    u64*                M_small_hash;
    large_hash_type*    M_large_hash;
};

void build_lookup_table_mt(const HbnUnpackedDatabase* updb,
        const int kmer_size,
        const int kmer_window,
        const double repeat_kmer_frac,
        const int num_threads,
        std::vector<std::pair<u64, u64>>& hash_to_offset_list,
        size_t& num_kmer,
        FILE* out);

void dump_lookup_table(int kmer_size, std::pair<u64, u64>* hashs, size_t num_hash, KmerInfo* kmers, size_t num_kmers, FILE* out);

void load_lookup_table(HbnHashTable*& hash_tbale, KmerInfo*& kmer_list, int& kmer_size, FILE* in);

#endif // __LOOKUP_TABLE_HPP
