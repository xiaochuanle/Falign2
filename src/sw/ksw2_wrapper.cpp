#include "ksw2_wrapper.hpp"

#include "ksw2.h"
#include "hbn_traceback_aux.h"

#include <cstring>

using namespace std;

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	sc_ambi = sc_ambi > 0? -sc_ambi : sc_ambi;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = sc_ambi;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = sc_ambi;
}

Ksw2Data::Ksw2Data() {
    //km = km_init();
    km = nullptr;
    memset(&ez, 0, sizeof(ksw_extz_t));

#if 0
    reward = 2;
    penalty = 4;
    go = 4;
    ge = 2;
    go1 = 24;
    ge1 = 1;
#else
    //-o 5 -O 56 -e 4 -E 1 -A 2 -B 5
    reward = 2;
    penalty = 5;
    go = 5;
    ge = 4;
    go1 = 56;
    ge1 = 1;
#endif 
    ambi_penalty = 1;
    zdrop = -1;
    end_bonus = -1;
    band_width = -1;
    ksw_gen_simple_mat(5, mat, reward, penalty, ambi_penalty);
}

void extract_sw_scoring_params(int* match_reward, int* mismatch_penalty, int* gap_open, int* gap_extend, int* gap_open1, int* gap_extend1)
{
#if 0
    if (match_reward) *match_reward = 2;
    if (mismatch_penalty) *mismatch_penalty = 4;
    if (gap_open) *gap_open = 4;
    if (gap_extend) *gap_extend = 2;
    if (gap_open1) *gap_open1 = 24;
    if (gap_extend1) *gap_extend1 = 1;
#else
    if (match_reward) *match_reward = 2;
    if (mismatch_penalty) *mismatch_penalty = 5;
    if (gap_open) *gap_open = 5;
    if (gap_extend) *gap_extend = 4;
    if (gap_open1) *gap_open1 = 56;
    if (gap_extend1) *gap_extend1 = 1;
#endif
}

static void
s_run_ksw2_nw(Ksw2Data* data, const u8* seq1, const int sl1, const u8* seq2, const int sl2,
    const int flag, const int zdrop, const int bw, ksw_extz_t* ez)
{
    if (data->go == data->go1 && data->ge == data->ge1) {
        ksw_extz2_sse(data->km, sl1, seq1, sl2, seq2, 5, data->mat,
            data->go, data->ge, bw, zdrop, data->end_bonus, flag, ez);
    } else {
        ksw_extd2_sse(data->km, sl1, seq1, sl2, seq2, 5, data->mat,
            data->go, data->ge, data->go1, data->ge1, bw, zdrop, data->end_bonus, flag, ez);
    }   
}

void ksw2_nw(Ksw2Data* data, const u8* seq1, const int sl1, const u8* seq2, const int sl2)
{
    data->free_cigar();
    ksw_extz_t* ez = &data->ez;
    memset(ez, 0, sizeof(ksw_extz_t));
    if (sl1 == 0 || sl2 == 0) return;
    int flag = 0;
    int zdrop = -1;
    int bw = 2 * abs(sl1 - sl2);
    bw = max(32, bw);

    s_run_ksw2_nw(data, seq1, sl1, seq2, sl2, flag, zdrop, bw, ez);
    if (ez->n_cigar == 0) {
        if (ez->cigar) kfree(data->km, ez->cigar);
        bw = max(sl1, sl2);
        bw = max(32, bw);
        s_run_ksw2_nw(data, seq1, sl1, seq2, sl2, flag, zdrop, bw, ez);
    }
    hbn_assert(ez->n_cigar > 0);
}
