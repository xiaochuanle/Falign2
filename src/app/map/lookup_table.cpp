#include "lookup_table.hpp"

#include "../../corelib/parasort.h"
#include "../../corelib/pdqsort.h"

#include <mutex>

using namespace std;

struct KmerStats
{
    u64 total_kmers;
    u64 distinct_kmers;
    u64 removed_kmers;
    u64 removed_distinct_kmers;

    KmerStats() {
        total_kmers = 0;
        distinct_kmers = 0;
        removed_kmers = 0;
        removed_distinct_kmers = 0;
    }

    void add(u64 tk, u64 dk, u64 rk, u64 rdk) {
        total_kmers += tk;
        distinct_kmers += dk;
        removed_kmers += rk;
        removed_distinct_kmers += rdk;
    }

    void stats() {
        double p = total_kmers ? 100.0 * removed_kmers / total_kmers : 0.0;
        fprintf(stderr, "Total kmers: %zu, %zu (%g%%) are removed\n", total_kmers, removed_kmers, p);
        p = distinct_kmers ? 100.0 * removed_distinct_kmers / distinct_kmers : 0.0;
        fprintf(stderr, "Distinct kmers: %zu, %zu (%g%%) are removed\n", distinct_kmers, removed_distinct_kmers, p);
    }
};

class KmerExtractThreadData
{
public:
    KmerExtractThreadData(const HbnUnpackedDatabase* updb, int kmer_size, int kmer_window)
    {
        M_updb = updb;
        M_num_subjects = M_updb->num_targets();
        
        M_kmer_size = kmer_size;
        M_kmer_window = kmer_window;

        u64 max_hash_value = U64_ONE << (kmer_size << 1);
        u64 hash_per_bucket = max_hash_value / kNumBuckets;

        M_bucket_boundaries[0] = 0;
        for (int i = 1; i < kNumBuckets; ++i) M_bucket_boundaries[i] = M_bucket_boundaries[i-1] + hash_per_bucket;
        M_bucket_boundaries[kNumBuckets] = U64_MAX;

        M_bucket_id = -1;
    }

    void reset_bucket_id()
    {
        M_bucket_id = -1;
    }

    bool set_next_bucket()
    {
        ++M_bucket_id;
        if (M_bucket_id >= kNumBuckets) return false;

        M_sid = 0;
        M_min_hash = M_bucket_boundaries[M_bucket_id];
        M_max_hash = M_bucket_boundaries[M_bucket_id+1];

        M_kmers.clear();

        return true;
    }

    void extract_kmers(KmerHashAndOffset*& A, size_t& NA)
    {
        A = M_kmers.data();
        NA = M_kmers.size();
    }

    void append_kmers(KmerHashAndOffset* a, size_t ac)
    {
        std::lock_guard<std::mutex> __(M_khao_list_mutex);
        M_kmers.insert(M_kmers.end(), a, a + ac);
    }    

    int next_sid()
    {
        std::lock_guard<std::mutex> __(M_sid_mutex);
        int i = -1;
        if (M_sid < M_num_subjects) i = M_sid++;
        return i;
    }

    const HbnUnpackedDatabase* updb() const
    {
        return M_updb;
    }

    u64 min_hash()
    {
        return M_min_hash;
    }

    u64 max_hash()
    {
        return M_max_hash;
    }

    int kmer_size()
    {
        return M_kmer_size;
    }

    int kmer_window()
    {
        return M_kmer_window;
    }

private:
    static constexpr int kNumBuckets = 16;
    u64 M_bucket_boundaries[kNumBuckets+1];
    int M_bucket_id;

private:
    const HbnUnpackedDatabase*      M_updb;
    int                             M_num_subjects;
    int                             M_sid;
    std::mutex                      M_sid_mutex;

    u64                             M_min_hash;
    u64                             M_max_hash;

    int                             M_kmer_size;
    int                             M_kmer_window;

    std::mutex                      M_khao_list_mutex;
    std::vector<KmerHashAndOffset>             M_kmers;
};

static void*
s_kmer_extract_thread(void* params)
{
    KmerExtractThreadData* data = (KmerExtractThreadData*)(params);
    const HbnUnpackedDatabase* updb = data->updb();
    u64 min_hash = data->min_hash();
    u64 max_hash = data->max_hash();
    int kmer_size = data->kmer_size();
    int kmer_window = data->kmer_window();
    const u64 kHashMask = (U64_ONE << (2 * kmer_size)) - 1;
    const u64 kRevShift = 2 * (kmer_size - 1);
    vector<KmerHashAndOffset> khao_list;
    KmerHashAndOffset khao;

    int sid;
    while ((sid = data->next_sid()) != -1) {
        const char* chr = updb->target_sequence(sid, FWD);
        const int chr_size = updb->target_size(sid);
        //HBN_LOG("=========> %d:%d, %d, %d", sid, chr_size, kmer_size, kmer_window);
        int pos = 0;
        u64 fhash = 0, rhash = 0;
        for (int i = 0; i < chr_size; ++i) {
            u64 c = chr[i];
            c = nst_nt4_table[c];
            if (c > 3) { pos = i + 1; fhash = 0; rhash = 0; continue; }
            fhash = ((fhash << 2) | c) & kHashMask;
            rhash = (rhash >> 2) | (3ULL^c) << kRevShift;
            if (i + 1 - pos == kmer_size) {
                u64 hash = min(fhash, rhash);
                if (fhash != rhash && hash >= min_hash && hash < max_hash && (pos % kmer_window) == 0) {
                    khao.ki.seq_id = sid;
                    if (fhash < rhash) {
                        khao.ki.seq_offset = pos * 2;
                    } else {
                        khao.ki.seq_offset = pos * 2 + 1;
                    }
                    khao.hash = hash;
                    khao_list.push_back(khao);
                    if (khao_list.size() == 10000) {
                        data->append_kmers(khao_list.data(), khao_list.size());
                        khao_list.clear();
                    }
                }
                ++pos;
            }
        }
    }
    if (!khao_list.empty()) data->append_kmers(khao_list.data(), khao_list.size());

    return nullptr;
}

void extract_gfa_kmers(KmerExtractThreadData* data, const int num_threads)
{
    pthread_t jobs[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobs + i, nullptr, s_kmer_extract_thread, data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], nullptr);
    }
}

int
resolve_repeat_kmer_occ_cutoff(KmerHashAndOffset* A,
    u64 NA, 
    int num_threads,
    int kmer_size,
    double rep_frac, 
    int _max_kmer_occ)
{
    u64 dnk = 0;
    u64 i = 0;
    while (i < NA) {
        u64 j = i + 1;
        while (j < NA && A[i].hash == A[j].hash) ++j;
        ++dnk;
        i = j;
    }
    u64 cnt_list_size = dnk;
    int* cnt_list = new int[cnt_list_size];
    const u64 kKmerOccT = numeric_limits<int>::max() / 2;
    u64 cnt_list_idx = 0;
    i = 0;
    while (i < NA) {
        u64 j = i + 1;
        while (j < NA && A[i].hash == A[j].hash) ++j;
        u64 cnt = j - i;
        cnt = min<u64>(cnt, kKmerOccT);
        if (cnt) cnt_list[cnt_list_idx++] = cnt;
        i = j;
    }
    hbn_assert(cnt_list_idx == cnt_list_size);

    int* a = cnt_list;
    u64 c = cnt_list_size;
    pdqsort(a, a + c, greater<int>());
    const u64 kRepC = cnt_list_size * rep_frac;
    int cutoff = a[kRepC];
    HBN_LOG("Total kmers: %llu", NA);
    HBN_LOG("distinct kmers: %llu, rep_frac = %g, cutoff = %d", c, rep_frac, cutoff);
    cutoff = hbn_max(cutoff, _max_kmer_occ);
    cutoff = hbn_max(100, cutoff);
    fprintf(stderr, "max kmer occ: %llu, rep_c = %llu, cutoff = %d\n", a[0], kRepC, cutoff);

    delete[] cnt_list;
    return cutoff;
}

void build_lookup_table_mt(const HbnUnpackedDatabase* updb,
        const int kmer_size,
        const int kmer_window,
        const double repeat_kmer_frac,
        const int num_threads,
        std::vector<std::pair<u64, u64>>& hash_to_offset_list,
        size_t& num_kmer,
        FILE* out)
{
    KmerExtractThreadData data(updb, kmer_size, kmer_window);
    size_t ki = 0;
    KmerStats all_stats;
    num_kmer = 0;

    data.reset_bucket_id();
    while (data.set_next_bucket()) {
        extract_gfa_kmers(&data, num_threads);
        KmerHashAndOffset* A;
        size_t NA;
        data.extract_kmers(A, NA);
        if (NA == 0) continue;
        parasort(NA, A, num_threads);
        size_t cutoff = resolve_repeat_kmer_occ_cutoff(A, NA, num_threads, kmer_size, repeat_kmer_frac, 1);
        KmerStats stats;
        size_t i = 0;
        while (i < NA) {
            size_t j = i + 1;
            while (j < NA && A[i].hash == A[j].hash) ++j;
            size_t cnt = j - i;
            if (cnt > cutoff) { 
                stats.add(cnt, 1, cnt, 1);
                all_stats.add(cnt, 1, cnt, 1);
                i = j; 
                continue; 
            }
            stats.add(cnt, 1, 0, 0);
            all_stats.add(cnt, 1, 0, 0);

            u64 offset = ki;
            u64 u = lktbl_pack_hash_value(offset, cnt);
            hash_to_offset_list.emplace_back(A[i].hash, u);
            for (size_t k = i; k < i + cnt; ++k, ++ki) {
                hbn_fwrite(&A[k].ki, sizeof(KmerInfo), 1, out);
            }
            i = j;
        }
        stats.stats();
    } 
    all_stats.stats();
    num_kmer = ki;
}

void dump_lookup_table(int kmer_size, std::pair<u64, u64>* hashs, size_t num_hash, KmerInfo* kmers, size_t num_kmers, FILE* out)
{
    void* s = nullptr;
    size_t cnt = 0;
    size_t size = 0;

    s = (void*)(&kmer_size);
    cnt = 1;
    size = sizeof(int);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(&num_hash);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(hashs);
    cnt = num_hash;
    size = sizeof(hashs[0]);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(&num_kmers);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fwrite(s, size, cnt, out);

    s = (void*)(kmers);
    cnt = num_kmers;
    size = sizeof(KmerInfo);
    hbn_fwrite(s, size, cnt, out);
}

void load_lookup_table(HbnHashTable*& hash_table, KmerInfo*& kmer_list, int& kmer_size, FILE* in)
{
    void* s = nullptr;
    size_t cnt = 0;
    size_t size = 0;

    s = (void*)(&kmer_size);
    cnt = 1;
    size = sizeof(int);
    hbn_fread(s, size, cnt, in);

    size_t num_hash;
    s = (void*)(&num_hash);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fread(s, size, cnt, in);

    pair<u64, u64>* hashs = new pair<u64, u64>[num_hash];
    s = (void*)(hashs);
    cnt = num_hash;
    size = sizeof(hashs[0]);
    hbn_fread(s, size, cnt, in);

    hash_table = new HbnHashTable(kmer_size);
    hash_table->reserve(num_hash);
    for (size_t i = 0; i < num_hash; ++i) hash_table->add_one_hash_and_value(hashs[i].first, hashs[i].second);
    delete[] hashs;

    size_t num_kmers;
    s = (void*)(&num_kmers);
    cnt = 1;
    size = sizeof(size_t);
    hbn_fread(s, size, cnt, in);

    kmer_list = new KmerInfo[num_kmers];
    s = (void*)(kmer_list);
    cnt = num_kmers;
    size = sizeof(KmerInfo);
    hbn_fread(s, size, cnt, in); 
}
