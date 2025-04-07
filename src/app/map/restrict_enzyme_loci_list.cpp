#include "restrict_enzyme_loci_list.hpp"

#include <mutex>
#include <vector>
#include <cstring>

using namespace std;

void
RestrictEnzyme_Init(const char* enzyme, RestrictEnzyme* re)
{
    int l = strlen(enzyme);
    //if (l > MAX_ENZYME_SIZE) HBN_ERR("Restrict enzyme is too long: %s", enzyme);

    re->enzyme_size = 0;
    re->break_loci = -1;
    for (int i = 0; i < l; ++i) {
        int c = enzyme[i];
        if (c == '^') {
            re->break_loci = i;
            continue;
        }
        c = nst_nt16_table[c];
        if (c > 3) {
            HBN_ERR("Illegal character '%c' in restrict enzyme '%s'", enzyme[i], enzyme);
        }
        re->enzyme[re->enzyme_size] = enzyme[i];
        re->encoded_enzyme[re->enzyme_size] = c;
        ++re->enzyme_size;
    }
    if (re->break_loci == -1) {
        HBN_ERR("'^' is missing from restrict enzyme '%s'", enzyme);
    }
    re->enzyme_mask = (U64_ONE << (2 * re->enzyme_size)) - 1;
    re->enzyme_hash = 0;
    for (int i = 0; i < re->enzyme_size; ++i) {
        re->enzyme_hash = (re->enzyme_hash << 2) | re->encoded_enzyme[i];
    }

    for (int i = 0; i < re->enzyme_size; ++i) {
        int fc = re->encoded_enzyme[i];
        int rc = re->encoded_enzyme[re->enzyme_size - 1 - i];
        if (fc + rc != 3) {
            HBN_ERR("Illegal restrict enzyme '%s'", enzyme);
        }
    }
}

void extract_enzyme_loci_list_for_one_seq(const char* subject,
    const int subject_size,
    RestrictEnzyme* re,
    std::vector<int>& reloci_list)
{
    const u64 kHashMask = re->enzyme_mask;
    const u64 kEnzymeHash = re->enzyme_hash;
    const int enzyme_size = re->enzyme_size;
    reloci_list.clear();
    reloci_list.push_back(0);

    u64 hash = 0;
    int pos = 1;
    for (int i = 1; i < subject_size; ++i) {
        u64 c = subject[i];
        c = nst_nt4_table[c];
        if (c > 3) { hash = 0; pos = i + 1; continue; }
        hash = ((hash << 2) | c) & kHashMask;
        if (i + 1 - pos == enzyme_size) {
            if (hash == kEnzymeHash) reloci_list.push_back(pos);
            ++pos;
        }
    }
    reloci_list.push_back(subject_size);
    reloci_list.push_back(subject_size + 1);
}

struct RestrictEnzymeListThreadData
{
public:
    RestrictEnzymeListThreadData(const HbnUnpackedDatabase* updb, RestrictEnzymeLociList* enzyme_pos_list)
    {
        M_updb = updb;
        M_vtx_id = 0;
        M_num_vtx = M_updb->num_targets() * 2;
        M_enzyme_pos_list = enzyme_pos_list;
    }

    int get_next_vtx_id()
    {
        lock_guard<mutex> __(M_vtx_id_mutex);
        return M_vtx_id++;
    }

public:
    const HbnUnpackedDatabase*    M_updb;
    int                     M_vtx_id;
    int                     M_num_vtx;
    mutex                   M_vtx_id_mutex;
    RestrictEnzymeLociList* M_enzyme_pos_list;
};

static void*
s_chr_enzyme_pos_cnt_thread(void* params)
{
    RestrictEnzymeListThreadData* data = (RestrictEnzymeListThreadData*)(params);
    vector<int> enzyme_pos_list;
    int vtx_id;
    while (1) {
        vtx_id = data->get_next_vtx_id();
        if (vtx_id >= data->M_num_vtx) break;
        int sid = vtx_id / 2;
        int sdir = vtx_id % 2;
        const int chr_size = data->M_updb->target_size(sid);

        enzyme_pos_list.clear();
        const char* chr_seq = data->M_updb->target_sequence(sid, sdir);
        extract_enzyme_loci_list_for_one_seq(chr_seq, chr_size, &data->M_enzyme_pos_list->enzyme, enzyme_pos_list);
        data->M_enzyme_pos_list->seq_reloci_info_array[vtx_id].enzyme_loci_cnt = enzyme_pos_list.size();
    }

    return nullptr;
}

static void
s_compute_chr_enzyme_pos_counts_mt(RestrictEnzymeListThreadData* data, const int num_threads)
{
    pthread_t jobs[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobs + i, nullptr, s_chr_enzyme_pos_cnt_thread, data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], nullptr);
    }
}

static void*
s_chr_enzyme_pos_fill_thread(void* params)
{
    RestrictEnzymeListThreadData* data = (RestrictEnzymeListThreadData*)(params);
    vector<int> enzyme_pos_list;
    int vtx_id;
    while (1) {
        vtx_id = data->get_next_vtx_id();
        if (vtx_id >= data->M_num_vtx) break;
        int sid = vtx_id / 2;
        int sdir = vtx_id % 2;
        const int chr_size = data->M_updb->target_size(sid);

        enzyme_pos_list.clear();
        const char* chr_seq = data->M_updb->target_sequence(sid, sdir);
        extract_enzyme_loci_list_for_one_seq(chr_seq, chr_size, &data->M_enzyme_pos_list->enzyme, enzyme_pos_list);
        int n = enzyme_pos_list.size();
        hbn_assert(n == data->M_enzyme_pos_list->seq_reloci_info_array[vtx_id].enzyme_loci_cnt);
        int* dst = data->M_enzyme_pos_list->reloci_array + data->M_enzyme_pos_list->seq_reloci_info_array[vtx_id].enzyme_loci_offset;
        copy(enzyme_pos_list.begin(), enzyme_pos_list.end(), dst);
    }

    return nullptr;
}

static void
s_fill_chr_enzyme_pos_list_mt(RestrictEnzymeListThreadData* data, const int num_threads)
{
    pthread_t jobs[num_threads];
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobs + i, nullptr, s_chr_enzyme_pos_fill_thread, data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], nullptr);
    }
}

RestrictEnzymeLociList*
RestrictEnzymeLociListNew(const HbnUnpackedDatabase* updb, const char* enzyme, const int num_threads)
{
    RestrictEnzymeLociList* list = (RestrictEnzymeLociList*)calloc(1, sizeof(RestrictEnzymeLociList));
    RestrictEnzyme_Init(enzyme, &list->enzyme);

    const int n_seq = updb->num_targets();
    list->seq_reloci_info_array = (SeqRestrictEnzymeLociInfo*)calloc(2 * n_seq, sizeof(SeqRestrictEnzymeLociInfo));

    RestrictEnzymeListThreadData data(updb, list);
    data.M_vtx_id = 0;
    s_compute_chr_enzyme_pos_counts_mt(&data, num_threads);
    
    data.M_vtx_id = 0;
    list->reloci_cnt = 0;
    for (int i = 0; i < 2 * n_seq; ++i) {
        list->seq_reloci_info_array[i].enzyme_loci_offset = list->reloci_cnt;
        list->reloci_cnt += list->seq_reloci_info_array[i].enzyme_loci_cnt;
        //HBN_LOG("%d\t%zu\t%d", i, list->seq_reloci_info_array[i].enzyme_loci_offset, list->seq_reloci_info_array[i].enzyme_loci_cnt);
    }
    list->reloci_array = (int*)calloc(list->reloci_cnt, sizeof(int));

    s_fill_chr_enzyme_pos_list_mt(&data, num_threads);

    return list;
}

RestrictEnzymeLociList*
RestrictEnzymeLociListFree(RestrictEnzymeLociList* list)
{
    free(list->seq_reloci_info_array);
    free(list->reloci_array);
    free(list);
    return NULL;
}

int 
offset_to_enzyme_intv_idx(const int* loci_array, const int loci_cnt, const int offset, int* intv_cnt)
{
    hbn_assert(loci_cnt >= 3);
    int left = 0, right = loci_cnt, mid = 0;
    while (left < right) {
        mid = (left + right) / 2;
        if (offset >= loci_array[mid]) {
            if (mid == loci_cnt - 1) break;
            if (offset < loci_array[mid+1]) break;
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    //fprintf(stderr, "mid = %d, cnt = %d\n", mid, loci_cnt);
    if (mid == loci_cnt - 1 || mid == loci_cnt - 2) {
        mid = loci_cnt - 3;
    }
    hbn_assert(mid >= 0);
    hbn_assert(offset >= loci_array[mid] && offset <= loci_array[mid+1], "mid = %d, cnt = %d, offset = %d, s1 = %d, s2 = %d", mid, loci_cnt, offset,
		    loci_array[mid], loci_array[mid+1]);

    if (intv_cnt) *intv_cnt = loci_cnt - 1;
    return mid;
}

int get_enzyme_pos(const int* loci_array, const int loci_cnt, const int offset)
{
    int x = offset_to_enzyme_intv_idx(loci_array, loci_cnt, offset, nullptr);
    int ld = offset - loci_array[x];
    int rd = loci_array[x+1] - offset;
    return (ld <= rd) ? loci_array[x] : loci_array[x+1];
}
