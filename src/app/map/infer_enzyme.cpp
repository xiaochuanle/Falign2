#include "infer_enzyme.hpp"

#include "../../corelib/pdqsort.h"
#include "gapped_align.hpp"
#include "prelim_search.hpp"

#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <vector>

using namespace std;

#ifndef U32_MAX
constexpr const u32 U32_MAX = numeric_limits<u32>::max();
#endif

static const char* kEnzymeNameList[] = {
    "DpnII",
    "HindIII",
    "NcoI",
    "NlaIII"
};

static const char* kEnzymeList[] = {
    "^GATC",
    "A^AGCTT",
    "C^CATGG",
    "CATG^"
};
static constexpr const int kEnzymeListSize = 4;

class EnzymeList
{
public:
    EnzymeList(HbnMapOptions* options);
    ~EnzymeList();

    RestrictEnzyme** get_enzyme_list() {
        return M_enzyme_list.data();
    }

    void add_one_enzyme_count(u32 hash) {
        std::lock_guard<std::mutex> lg(M_mutex);
        auto pos = M_enzyme_stats.find(hash);
        hbn_assert(pos != M_enzyme_stats.end());
        ++pos->second.second;
    }

    void dump_stats();

    void find_enzyme(const u8* query,
        const int qoff, 
        const int qend,
        const int qsize,
        const char* subject,
        const int soff,
        const int send,
        const int ssize);

    std::string infer_enzyme() {
        dump_stats();
        
        pair<int, size_t> cnts[kEnzymeListSize];
        for (int i = 0; i < kEnzymeListSize; ++i) {
            u32 hash = M_enzyme_list[i]->enzyme_hash;
            auto pos = M_enzyme_stats.find(hash);
            hbn_assert(pos != M_enzyme_stats.end());
            cnts[i] = pair<int, size_t>(i, pos->second.second);
        }
        pdqsort(cnts, cnts + kEnzymeListSize, [](const pair<int, size_t>& x, const pair<int, size_t>& y) { return x.second > y.second; });
        int c0 = cnts[0].second;
        int c1 = cnts[1].second;
        hbn_assert(c0 >= c1);
        string enzyme;
        if (c0 > c1 * 2) enzyme = kEnzymeList[cnts[0].first];
        return enzyme; 
    }

private:
    HbnMapOptions*  M_options;
    std::mutex  M_mutex;
    std::vector<RestrictEnzyme*> M_enzyme_list;
    std::map<u32, std::pair<RestrictEnzyme*, size_t>>   M_enzyme_stats;
};

EnzymeList::EnzymeList(HbnMapOptions* options)
{
    M_options = options;
    for (int i = 0; i < kEnzymeListSize; ++i) {
        RestrictEnzyme* enzyme = new RestrictEnzyme();
        RestrictEnzyme_Init(kEnzymeList[i], enzyme);
        M_enzyme_list.push_back(enzyme);
        M_enzyme_stats[enzyme->enzyme_hash] = pair<RestrictEnzyme*, size_t>(enzyme, 0);
    }
}

EnzymeList::~EnzymeList()
{
    for (auto enzyme : M_enzyme_list) delete enzyme;
}

void EnzymeList::dump_stats()
{
    char buf[256];
    for (int i = 0; i < kEnzymeListSize; ++i) {
        u32 hash = M_enzyme_list[i]->enzyme_hash;
        auto pos = M_enzyme_stats.find(hash);
        hbn_assert(pos != M_enzyme_stats.end());
        int cnt = pos->second.second;
        snprintf(buf, 256, "%s(%s)", kEnzymeNameList[i], kEnzymeList[i]);
        print_fixed_width_string(stderr, buf, 20);
        fprintf(stderr, "%d\n", cnt);
    }
}

static bool 
s_enzyme_exists(RestrictEnzyme* enzyme,
    const int kEndDist,
    const int kEnzymeFlankingBases,
    const u8* query,
    int qoff,
    int qsize,
    const char* subject,
    int soff,
    int ssize)
{
    bool r = (qoff < kEndDist) || (qsize - qoff < kEndDist)
        || (soff < kEndDist) || (ssize - soff < kEndDist);
    if (r) return false;
    
    int qfrom = qoff - kEnzymeFlankingBases;
    int qto = qoff + kEnzymeFlankingBases;
    int p = qfrom;
    bool q_exists = false;
    u32 hash = 0;
    for (; p < qfrom + enzyme->enzyme_size - 1; ++p) hash = (hash << 2) | query[p];
    for (; p < qto; ++p) {
        hash = (hash << 2) | query[p];
        hash = hash & enzyme->enzyme_mask;
        if (hash == enzyme->enzyme_hash) {
            q_exists = true;
            break;
        }
    }
    if (!q_exists) return false;

    int sfrom = soff - kEnzymeFlankingBases;
    int sto = soff + kEnzymeFlankingBases;
    p = sfrom;
    bool s_exists = false;
    hash = 0;
    for (; p < sfrom + enzyme->enzyme_size - 1; ++p) {
        u64 c = subject[p];
        c = nst_nt4_table[c];
        if (c > 3) { hash = 0; continue; }
        hash = (hash << 2) | c;
    }
    for (; p < sto; ++p) {
        u64 c = subject[p];
        c = nst_nt4_table[c];
        if (c > 3) { hash = 0; continue; }
        hash = (hash << 2) | c;
        hash = hash & enzyme->enzyme_mask;
        if (hash == enzyme->enzyme_hash) {
            s_exists = true;
            break;
        }        
    }
    if (!s_exists) return false;

    return true;
}

void EnzymeList::find_enzyme(const u8* query,
    const int qoff, 
    const int qend,
    const int qsize,
    const char* subject,
    const int soff,
    const int send,
    const int ssize)
{
    RestrictEnzyme** enzymes = get_enzyme_list();
    u32 L_exist[kEnzymeListSize]; fill(L_exist, L_exist + kEnzymeListSize, U32_MAX);
    u32 R_exist[kEnzymeListSize]; fill(R_exist, R_exist + kEnzymeListSize, U32_MAX);
    for (int i = 0; i < kEnzymeListSize; ++i) {
	    if (s_enzyme_exists(enzymes[i], M_options->ei_end_dist, M_options->ei_flanking_bases, query, qoff, qsize, subject, soff, ssize)) L_exist[i] = enzymes[i]->enzyme_hash;
    }
    for (int i = 0; i < kEnzymeListSize; ++i) {
	    if (s_enzyme_exists(enzymes[i], M_options->ei_end_dist, M_options->ei_flanking_bases, query, qend, qsize, subject, send, ssize)) R_exist[i] = enzymes[i]->enzyme_hash;
    }

    for (int i = 0; i < kEnzymeListSize; ++i) {
	    if (L_exist[i] == U32_MAX || R_exist[i] == U32_MAX) continue;
	    add_one_enzyme_count(L_exist[i]);
    }
}

////////////////////////////

struct EnzymeInferenceThreadWorkData
{
    int thread_id;
    HbnIndex* index;
    QueryReader* queries;
    HbnMapOptions* opts;
    EnzymeList* enzyme_list;
};

EnzymeInferenceThreadWorkData* 
EnzymeInferenceThreadWorkDataNew(int thread_id,
    HbnIndex* index,
    QueryReader* queries,
    HbnMapOptions* options,
    EnzymeList* enzyme_list)
{
    EnzymeInferenceThreadWorkData* data = new EnzymeInferenceThreadWorkData();
    data->thread_id = thread_id;
    data->index = index;
    data->queries = queries;
    data->opts = options;
    data->enzyme_list = enzyme_list;

    return data;
}

EnzymeInferenceThreadWorkData*
EnzymeInferenceThreadWorkDataFree(EnzymeInferenceThreadWorkData* data)
{
    delete data;
    return NULL;
}

#if 0
static void*
s_enzyme_inference_thread(void* params)
{
    EnzymeInferenceThreadWorkData* data = (EnzymeInferenceThreadWorkData*)(params);
    HbnUnpackedDatabase* updb = data->index->updb();
    PrelimSearch prelim_search(data->opts);
    HbnTracebackData tbck_data;
    PoreCQuery query(nullptr);
    vector<PoreCInitHit> align_list;
    string align_strings;

    while (data->queries->get_next_batch_query(&query)) {
        query.init();
        HbnChainInfo* chains = nullptr;
        int num_chains = 0;
        pu64_t* seeds;
        prelim_search.find_init_hits(&query, data->index, chains, num_chains, seeds);
        for (int i = 0; i < num_chains && i < data->opts->ei_num_hits; ++i) {
            const char* subject = updb->target_sequence(chains[i].sid, chains[i].sdir);
            const int subject_size = updb->target_size(chains[i].sid);
            bool r = tbck_data.map_one_chain(query.fwd_qs, chains[i].qb, chains[i].qe, query.size,
                subject, chains[i].soff, chains[i].send, subject_size);
            if (!r) continue;
	    if (tbck_data.qend - tbck_data.qoff < data->opts->ei_frag_size) continue;
            data->enzyme_list->find_enzyme(query.fwd_qs, tbck_data.qoff, tbck_data.qend, query.size, subject, tbck_data.soff, tbck_data.send, subject_size);
        }
    }
    return NULL;
}
#else
static void*
s_enzyme_inference_thread(void* params)
{
    EnzymeInferenceThreadWorkData* data = (EnzymeInferenceThreadWorkData*)(params);
    HbnUnpackedDatabase* updb = data->index->updb();
    PrelimSearch prelim_search(data->opts);
    HbnTracebackData tbck_data;
    PoreCQuery query(nullptr);
    vector<PoreCInitHit> align_list;
    string align_strings;
    vector<u8> cov_list;

    while (data->queries->get_next_batch_query(&query)) {
        query.init();
        HbnChainInfo* chains = nullptr;
        int num_chains = 0;
        pu64_t* seeds;
        prelim_search.find_init_hits(&query, data->index, chains, num_chains, seeds);
        if (!num_chains) continue;
        extend_ddf_chain_list(&tbck_data, updb, &query, chains, num_chains, seeds, cov_list, align_list, align_strings);
        if (align_list.empty()) continue;
        set_mapQ_for_init_hits(align_list.data(), align_list.size());
        for (auto& a : align_list) {
            if (a.mapQ < 10) continue;
            if (a.qend - a.qoff < data->opts->ei_frag_size) continue;
            const char* subject = updb->target_sequence(a.sid, a.sdir);
            const int subject_size = updb->target_size(a.sid);
            data->enzyme_list->find_enzyme(query.fwd_qs, a.qoff, a.qend, query.size, subject, a.soff, a.send, subject_size);
        }
    }
    return NULL;
}
#endif

///////////////////////////

std::string infer_enzyme_mt(HbnMapOptions* options, HbnIndex* index, QueryReader* queries)
{
    queries->load_query_batch(options->ei_query_bases);
    queries->reset_batch_query_id();
    const int num_threads = options->num_threads;
    EnzymeList enzyme_list(options);
    EnzymeInferenceThreadWorkData* data_array[num_threads];
    pthread_t jobids[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = EnzymeInferenceThreadWorkDataNew(i, index, queries, options, &enzyme_list);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobids + i, NULL, s_enzyme_inference_thread, data_array[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobids[i], NULL);
    }
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = EnzymeInferenceThreadWorkDataFree(data_array[i]);
    }
    queries->reset_batch_query_id();
    return enzyme_list.infer_enzyme();
}
