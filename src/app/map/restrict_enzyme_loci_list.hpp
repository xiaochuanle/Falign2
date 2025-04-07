#ifndef __RESTRICT_ENZYME_LOCI_LIST_HPP
#define __RESTRICT_ENZYME_LOCI_LIST_HPP

#include "../../corelib/hbn_aux.h"
#include "../../corelib/unpacked_seqdb.hpp"

#include <vector>

#define MAX_ENZYME_SIZE 6

#define LEFT_ENZYME_LOCI_MATCH 1
#define RIGHT_ENZYME_LOCI_MATCH 2

typedef struct {
    char enzyme[MAX_ENZYME_SIZE];
    u8 encoded_enzyme[MAX_ENZYME_SIZE];
    int enzyme_size;
    int break_loci;
    u64 enzyme_hash;
    u64 enzyme_mask;
} RestrictEnzyme;

void
RestrictEnzyme_Init(const char* enzyme, RestrictEnzyme* re);

typedef struct {
    size_t enzyme_loci_offset;
    int enzyme_loci_cnt;
} SeqRestrictEnzymeLociInfo;

typedef struct {
    RestrictEnzyme enzyme;
    SeqRestrictEnzymeLociInfo* seq_reloci_info_array;
    int* reloci_array;
    size_t reloci_cnt;
} RestrictEnzymeLociList;

RestrictEnzymeLociList*
RestrictEnzymeLociListFree(RestrictEnzymeLociList* list);

RestrictEnzymeLociList*
RestrictEnzymeLociListNew(const HbnUnpackedDatabase* updb, const char* enzyme, const int num_threads);

int 
offset_to_enzyme_intv_idx(const int* loci_array, const int loci_cnt, const int offset, int* intv_cnt);

int get_enzyme_pos(const int* loci_array, const int loci_cnt, const int offset);

void extract_enzyme_loci_list_for_one_seq(const char* subject,
    const int subject_size,
    RestrictEnzyme* re,
    std::vector<int>& reloci_list);

#endif // __RESTRICT_ENZYME_LOCI_LIST_HPP