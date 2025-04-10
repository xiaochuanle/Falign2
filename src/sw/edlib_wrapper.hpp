#ifndef __EDLIB_WRAPPER_HPP
#define __EDLIB_WRAPPER_HPP

#include "../corelib/hbn_aux.h"
#include "kalloc.h"

#include <string>
#include <vector>

struct EdlibAlignData {
    std::vector<char>       vqaln;
    std::vector<char>       vtaln;
    std::string             sqaln;
    std::string             staln;
    int                     tolerance;
    int                     do_traceback;
    int                     dist;
    int                     score;
    void*                   km;

    EdlibAlignData() {
        tolerance = -1;
        do_traceback = 0;
        dist = -1;
        score = -1;
        //km = km_init();
	km = nullptr;
    }

    ~EdlibAlignData() {
        //km_destroy(km);
    }
};

int
edlib_nw(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    std::string& qaln,
    std::string& taln);

int
edlib_nw_cigar(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    std::vector<u64>& cigar);

int
edlib_nw_dist(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    int tolerance,
    int* score);

int
edlib_shw(EdlibAlignData* data,
    const u8* query,
    int query_size,
    const u8* target,
    int target_size,
    int* qend,
    int* tend,
    std::string* qaln,
    std::string* taln);

void
edlib_extend(EdlibAlignData* data,
    const u8* query,
    const int query_size,
    const u8* target,
    const int target_size,
    const int block_size,
    const BOOL right_extend,
    std::vector<u8>& qfrag,
    std::vector<u8>& tfrag,
    int* qend,
    int* tend,
    std::string& qaln,
    std::string& taln);

#endif // __EDLIB_WRAPPER_H
