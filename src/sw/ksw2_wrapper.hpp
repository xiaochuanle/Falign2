#ifndef KSW2_WRAPPER_H
#define KSW2_WRAPPER_H

#include "../corelib/hbn_aux.h"
#include "kalloc.h"
#include "ksw2.h"

#include <string>
#include <vector>

struct Ksw2Data {
    void* km;
    ksw_extz_t ez;

    int reward, penalty, ambi_penalty;
    int go, ge;
    int go1, ge1;
    int zdrop;
    int band_width;
    int end_bonus;
    int8_t mat[25];

    Ksw2Data();

    ~Ksw2Data() {
        if (ez.cigar) kfree(km, ez.cigar);
        if (km) km_destroy(km);
    }

    void free_cigar()
    {
        if (ez.cigar) kfree(km, ez.cigar);
        ez.n_cigar = 0;
        ez.cigar = 0;
    }
};

void extract_sw_scoring_params(int* match_reward, int* mismatch_penalty, int* gap_open, int* gap_extend, int* gap_open1, int* gap_extend1);

void ksw2_nw(Ksw2Data* data, const u8* seq1, const int sl1, const u8* seq2, const int sl2);

#endif // KSW2_WRAPPER_H
