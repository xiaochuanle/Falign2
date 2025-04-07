#include "edlib_wrapper.hpp"

#include "edlib.h"
#include "hbn_traceback_aux.h"

using namespace std;

static void x_edlibAlignmentToCigar(const unsigned char* alignment, 
							 int alignmentLength, 
							 int tgtStart, int tgtEnd, 
							 int qryStart, int qryEnd, 
        vector<u64>& cigar) 
{
    int i = 0, qi = qryStart, si = tgtStart;
    while (i < alignmentLength) {
        unsigned char opuc = alignment[i];
        int j = i + 1;
        while (j < alignmentLength && alignment[j] == opuc) ++j;
        u64 Cop = 0;
        u64 Ccnt = j - i;
        if (opuc == EDLIB_EDOP_MATCH) {
            Cop = 7;
            qi += Ccnt;
            si += Ccnt;
        } else if (opuc == EDLIB_EDOP_MISMATCH) {
            Cop = 8;
            qi += Ccnt;
            si += Ccnt;
        } else if (opuc == EDLIB_EDOP_INSERT) {
            Cop = 1;
            qi += Ccnt;
        } else if (opuc == EDLIB_EDOP_DELETE) {
            Cop = 2;
            si += Ccnt;
        }
        cigar.push_back((Ccnt<<4)|Cop);
        i = j;
    }
    hbn_assert(qi == qryEnd);
    hbn_assert(si == tgtEnd);
}


int
edlib_nw_cigar(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    std::vector<u64>& cigar)
{
    if (query_size == 0 || target_size == 0) return 1;

    EdlibAlignTask task = EDLIB_TASK_PATH;
    int tolerance = abs(query_size - target_size);
    tolerance = max(tolerance, 32);
    EdlibAlignResult align = edlibAlign(data->km, query, query_size, target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task));

    if (align.numLocations == 0) {
        tolerance = hbn_max(query_size, target_size);
        align = edlibAlign(data->km, query, query_size, target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task));
    }
    hbn_assert(align.numLocations);
    if (align.numLocations == 0) return 0;

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    x_edlibAlignmentToCigar(align.alignment, align.alignmentLength, tBgn, tEnd, qBgn, qEnd, cigar);
    edlibFreeAlignResult(data->km, align);

    return 1;
}

int
edlib_nw(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    std::string& qaln,
    std::string& taln)
{
    qaln.clear();
    taln.clear();
    if (query_size == 0 || target_size == 0) return 1;

    EdlibAlignTask task = EDLIB_TASK_PATH;
    int tolerance = abs(query_size - target_size);
    tolerance = max(tolerance, 32);
    EdlibAlignResult align = edlibAlign(data->km, query, query_size, target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task));

    if (align.numLocations == 0) {
        tolerance = hbn_max(query_size, target_size);
        align = edlibAlign(data->km, query, query_size, target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task));
    }
    hbn_assert(align.numLocations);
    if (align.numLocations == 0) return 0;

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    data->vqaln.resize(align.alignmentLength+1);
    data->vtaln.resize(align.alignmentLength+1);
    edlibAlignmentToStrings(align.alignment,
        align.alignmentLength,
        tBgn,
        tEnd,
        qBgn,
        qEnd,
        target,
        query,
        data->vtaln.data(),
        data->vqaln.data());
    edlibFreeAlignResult(data->km, align);
    data->vqaln.pop_back();
    data->vtaln.pop_back();
    qaln.assign(data->vqaln.begin(), data->vqaln.end());
    taln.assign(data->vtaln.begin(), data->vtaln.end());
    hbn_assert(qaln.size() == taln.size());
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        query,
        qBgn,
        qEnd,
        qaln.c_str(),
        0,
        target,
        tBgn,
        tEnd,
        taln.c_str(),
        qaln.size(),
        TRUE);

    return 1;
}

int
edlib_nw_dist(EdlibAlignData* data,
	const u8* query,
    int query_size,
	const u8* target,
    int target_size,
    int tolerance,
    int* score)
{
    if (query_size == 0 || target_size == 0) return 1;

    EdlibAlignTask task = EDLIB_TASK_DISTANCE;
    if (tolerance < 0) tolerance = abs(query_size - target_size);
    tolerance = max(tolerance, 32);
    EdlibAlignResult align = edlibAlign(data->km, query, query_size, target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task));

    if (align.numLocations == 0) {
        tolerance = hbn_max(query_size, target_size);
        align = edlibAlign(data->km, query, query_size, target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task));      
    }
    hbn_assert(align.numLocations);
    edlibFreeAlignResult(data->km, align);
    if (align.numLocations == 0) return -1;
    if (score) {
        *score = ((query_size + target_size) / 2 - align.editDistance) * 2 - align.editDistance * 5;
    }
    return align.editDistance;
}

int
edlib_shw(EdlibAlignData* data,
    const u8* _query,
    int query_size,
    const u8* _target,
    int target_size,
    int* qend,
    int* tend,
    std::string* qaln,
    std::string* taln)
{
    if (qaln) qaln->clear();
    if (taln) taln->clear();
    *qend = *tend = 0;
    if (query_size == 0 || target_size == 0) return 0;

    EdlibAlignTask task = EDLIB_TASK_PATH;
    int tolerance = abs(query_size - target_size);// hbn_max(query_size, target_size) * 0.35;
    tolerance = max(tolerance, 32);
    EdlibAlignResult align = edlibAlign(data->km, _query, query_size, _target, target_size,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_SHW, task));
    if (align.numLocations == 0) return 0;
    if (qaln == nullptr || taln == nullptr) {
        *qend = query_size;
        *tend = align.endLocations[0] + 1;
        edlibFreeAlignResult(data->km, align);
        return 1;
    }

    int qBgn = 0;
    int qEnd = query_size;
    int tBgn = align.startLocations[0];
    int tEnd = align.endLocations[0] + 1;
    hbn_assert(tBgn == 0);
    data->vqaln.resize(align.alignmentLength+1);
    data->vtaln.resize(align.alignmentLength+1);
    edlibAlignmentToStrings(align.alignment,
        align.alignmentLength,
        tBgn,
        tEnd,
        qBgn,
        qEnd,
        _target,
        _query,
        data->vtaln.data(),
        data->vqaln.data());
    edlibFreeAlignResult(data->km, align);
    data->vqaln.pop_back();
    data->vtaln.pop_back();
    qaln->assign(data->vqaln.begin(), data->vqaln.end());
    taln->assign(data->vtaln.begin(), data->vtaln.end());
    hbn_assert(qaln->size() == taln->size());
    *qend = qEnd;
    *tend = tEnd;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0,
        _query,
        qBgn,
        qEnd,
        qaln->c_str(),
        0,
        _target,
        tBgn,
        tEnd,
        taln->c_str(),
        qaln->size(),
        TRUE);

    return 1;    
}