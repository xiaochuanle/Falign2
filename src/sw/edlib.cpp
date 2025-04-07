#include "edlib.h"

#include "kalloc.h"

#include <cassert>
#include <stdint.h>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <cstring>
#include <string>

using namespace std;

typedef uint64_t Word;
static const int WORD_SIZE = sizeof(Word) * 8; // Size of Word in bits
static const Word WORD_1 = static_cast<Word>(1);
static const Word HIGH_BIT_MASK = WORD_1 << (WORD_SIZE - 1);  // 100..00
static const int MAX_UCHAR = 3;
static const int kAlphabetLength = 4;

typedef struct {
    Word* Ps;
    Word* Ms;
    int* scores;
    int* firstBlocks;
    int* lastBlocks;
} AlignmentData;

AlignmentData* AlignmentDataNew(void* km, int maxNumBlocks, int targetLength)
{
    AlignmentData* data = Kmalloc(km, AlignmentData, 1);
    data->Ps = Kmalloc(km, Word, maxNumBlocks * targetLength);
    data->Ms = Kmalloc(km, Word, maxNumBlocks * targetLength);
    data->scores = Kmalloc(km, int, maxNumBlocks * targetLength);
    data->firstBlocks = Kmalloc(km, int, targetLength);
    data->lastBlocks = Kmalloc(km, int, targetLength);    
    return data;
}

AlignmentData* AlignmentDataFree(void* km, AlignmentData* data)
{
    if (!data) return nullptr;
    kfree(km, data->Ps);
    kfree(km, data->Ms);
    kfree(km, data->scores);
    kfree(km, data->firstBlocks);
    kfree(km, data->lastBlocks); 
    kfree(km, data);
    return nullptr;
}

struct Block {
    Word P;  // Pvin
    Word M;  // Mvin
    int score; // score of last cell in block;

    Block() {}
    Block(Word p, Word m, int s) :P(p), M(m), score(s) {}
};

static int myersCalcEditDistanceSemiGlobal(void* km, const Word* Peq, int W, int maxNumBlocks,
                                           int queryLength,
                                           const unsigned char* target, int targetLength,
                                           int k, EdlibAlignMode mode,
                                           int* bestScore_, int** positions_, int* numPositions_);

static int myersCalcEditDistanceNW(void* km, const Word* Peq, int W, int maxNumBlocks,
                                   int queryLength,
                                   const unsigned char* target, int targetLength,
                                   int k, int* bestScore_,
                                   int* position_, bool findAlignment,
                                   AlignmentData** alignData, int targetStopPosition);


static int obtainAlignment(void* km,
        const unsigned char* const query, const unsigned char* const rQuery, const int queryLength,
        const unsigned char* const target, const unsigned char* const rTarget, const int targetLength,
        const int bestScore, unsigned char** const alignment, int* const alignmentLength);

static int obtainAlignmentHirschberg(void* km,
        const unsigned char* const query, const unsigned char* const rQuery, const int queryLength,
        const unsigned char* const target, const unsigned char* const rTarget, const int targetLength,
        const int bestScore, unsigned char** const alignment, int* const alignmentLength);

static int obtainAlignmentTraceback(void* km, const int queryLength, const int targetLength,
                                    const int bestScore, const AlignmentData* const alignData,
                                    unsigned char** const alignment, int* const alignmentLength);

static inline int ceilDiv(int x, int y);

static inline unsigned char* createReverseCopy(void* km, const unsigned char* seq, int length);

static inline Word* buildPeq(void* km, const unsigned char* query, const int queryLength);


/**
 * Main edlib method.
 */
extern "C" EdlibAlignResult edlibAlign(void* km, const uint8_t* const query, const int queryLength,
                                       const uint8_t* const target, const int targetLength,
                                       const EdlibAlignConfig config) {
    EdlibAlignResult result;
    result.status = EDLIB_STATUS_OK;
    result.editDistance = -1;
    result.endLocations = result.startLocations = NULL;
    result.numLocations = 0;
    result.alignment = NULL;
    result.alignmentLength = 0;
    result.alphabetLength = kAlphabetLength;

    // Handle special situation when at least one of the sequences has length 0.
    if (queryLength == 0 || targetLength == 0) {
        if (config.mode == EDLIB_MODE_NW) {
            result.editDistance = std::max(queryLength, targetLength);
            result.endLocations = Kmalloc(km, int, 1);
            result.endLocations[0] = targetLength - 1;
            result.numLocations = 1;
        } else if (config.mode == EDLIB_MODE_SHW || config.mode == EDLIB_MODE_HW) {
            result.editDistance = queryLength;
            result.endLocations = Kmalloc(km, int, 1);
            result.endLocations[0] = -1;
            result.numLocations = 1;
        } else {
            result.status = EDLIB_STATUS_ERROR;
        }
        return result;
    }

    /*--------------------- INITIALIZATION ------------------*/
    int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE); // bmax in Myers
    int W = maxNumBlocks * WORD_SIZE - queryLength; // number of redundant cells in last level blocks
    Word* Peq = buildPeq(km, query, queryLength);
    /*-------------------------------------------------------*/

    /*------------------ MAIN CALCULATION -------------------*/
    // TODO: Store alignment data only after k is determined? That could make things faster.
    int positionNW; // Used only when mode is NW.
    AlignmentData* alignData = NULL;
    bool dynamicK = false;
    int k = config.k;
    if (k < 0) { // If valid k is not given, auto-adjust k until solution is found.
        dynamicK = true;
        k = WORD_SIZE; // Gives better results than smaller k.
    }

    do {
        if (config.mode == EDLIB_MODE_HW || config.mode == EDLIB_MODE_SHW) {
            myersCalcEditDistanceSemiGlobal(km, Peq, W, maxNumBlocks,
                                            queryLength, target, targetLength,
                                            k, config.mode, &(result.editDistance),
                                            &(result.endLocations), &(result.numLocations));
        } else {  // mode == EDLIB_MODE_NW
            myersCalcEditDistanceNW(km, Peq, W, maxNumBlocks,
                                    queryLength, target, targetLength,
                                    k, &(result.editDistance), &positionNW,
                                    false, &alignData, -1);
        }
        k *= 2;
    } while(dynamicK && result.editDistance == -1);

    if (result.editDistance >= 0) {  // If there is solution.
        // If NW mode, set end location explicitly.
        if (config.mode == EDLIB_MODE_NW) {
            result.endLocations = Kmalloc(km, int, 1);
            result.endLocations[0] = targetLength - 1;
            result.numLocations = 1;
        }

        // Find starting locations.
        if (config.task == EDLIB_TASK_LOC || config.task == EDLIB_TASK_PATH) {
            result.startLocations = Kmalloc(km, int, result.numLocations);
            if (config.mode == EDLIB_MODE_HW) {  // If HW, I need to calculate start locations.
                unsigned char* rTarget = createReverseCopy(km, target, targetLength);
                unsigned char* rQuery  = createReverseCopy(km, query, queryLength);
                // Peq for reversed query.
                Word* rPeq = buildPeq(km, rQuery, queryLength);
                for (int i = 0; i < result.numLocations; i++) {
                    int endLocation = result.endLocations[i];
                    if (endLocation == -1) {
                        // NOTE: Sometimes one of optimal solutions is that query starts before target, like this:
                        //                       AAGG <- target
                        //                   CCTT     <- query
                        //   It will never be only optimal solution and it does not happen often, however it is
                        //   possible and in that case end location will be -1. What should we do with that?
                        //   Should we just skip reporting such end location, although it is a solution?
                        //   If we do report it, what is the start location? -4? -1? Nothing?
                        // TODO: Figure this out. This has to do in general with how we think about start
                        //   and end locations.
                        //   Also, we have alignment later relying on this locations to limit the space of it's
                        //   search -> how can it do it right if these locations are negative or incorrect?
                        result.startLocations[i] = 0;  // I put 0 for now, but it does not make much sense.
                    } else {
                        int bestScoreSHW, numPositionsSHW;
                        int* positionsSHW;
                        myersCalcEditDistanceSemiGlobal(km,
                                rPeq, W, maxNumBlocks,
                                queryLength, rTarget + targetLength - endLocation - 1, endLocation + 1,
                                result.editDistance, EDLIB_MODE_SHW,
                                &bestScoreSHW, &positionsSHW, &numPositionsSHW);
                        // Taking last location as start ensures that alignment will not start with insertions
                        // if it can start with mismatches instead.
                        result.startLocations[i] = endLocation - positionsSHW[numPositionsSHW - 1];
                        kfree(km, positionsSHW); positionsSHW = nullptr;
                    }
                }
                kfree(km, rTarget); rTarget = nullptr;
                kfree(km, rQuery); rQuery = nullptr;
                kfree(km, rPeq); rPeq = nullptr;
            } else {  // If mode is SHW or NW
                for (int i = 0; i < result.numLocations; i++) {
                    result.startLocations[i] = 0;
                }
            }
        }

        // Find alignment -> all comes down to finding alignment for NW.
        // Currently we return alignment only for first pair of locations.
        if (config.task == EDLIB_TASK_PATH) {
            int alnStartLocation = result.startLocations[0];
            int alnEndLocation = result.endLocations[0];
            const unsigned char* alnTarget = target + alnStartLocation;
            const int alnTargetLength = alnEndLocation - alnStartLocation + 1;
            unsigned char* rAlnTarget = createReverseCopy(km, alnTarget, alnTargetLength);
            unsigned char* rQuery  = createReverseCopy(km, query, queryLength);
            obtainAlignment(km, query, rQuery, queryLength,
                            alnTarget, rAlnTarget, alnTargetLength,
                            result.editDistance, &(result.alignment), &(result.alignmentLength));
            kfree(km, rAlnTarget); rAlnTarget = nullptr;
            kfree(km, rQuery); rQuery = nullptr;
        }
    }
    /*-------------------------------------------------------*/

    //--- Free memory ---//
    kfree(km, Peq); Peq = nullptr;
    if (alignData) AlignmentDataFree(km, alignData); alignData = nullptr;
    //-------------------//

    return result;
}

extern "C" char* edlibAlignmentToCigar(void* km, const unsigned char* const alignment, const int alignmentLength,
                                       const EdlibCigarFormat cigarFormat) {
    if (cigarFormat != EDLIB_CIGAR_EXTENDED && cigarFormat != EDLIB_CIGAR_STANDARD) {
        return 0;
    }

    // Maps move code from alignment to char in cigar.
    //                        0    1    2    3
    char moveCodeToChar[] = {'=', 'I', 'D', 'X'};
    if (cigarFormat == EDLIB_CIGAR_STANDARD) {
        moveCodeToChar[0] = moveCodeToChar[3] = 'M';
    }

    int cigar_length = 0;
    char lastMove = 0;  // Char of last move. 0 if there was no previous move.
    int numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                ++cigar_length;
                numDigits++;
            }
            // Write code of move to cigar string.
            ++cigar_length;
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    return 0;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }

    char* cigar = Kmalloc(km, char, cigar_length+1);
    int ci = 0;
    lastMove = 0;  // Char of last move. 0 if there was no previous move.
    numOfSameMoves = 0;
    for (int i = 0; i <= alignmentLength; i++) {
        // if new sequence of same moves started
        if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
            // Write number of moves to cigar string.
            int numDigits = 0;
            for (; numOfSameMoves; numOfSameMoves /= 10) {
                cigar[ci++] = '0' + numOfSameMoves % 10;
                numDigits++;
            }
            reverse(cigar + ci - numDigits, cigar + ci);
            // Write code of move to cigar string.
            cigar[ci++] = lastMove;
            // If not at the end, start new sequence of moves.
            if (i < alignmentLength) {
                // Check if alignment has valid values.
                if (alignment[i] > 3) {
                    kfree(km, cigar);
                    return 0;
                }
                numOfSameMoves = 0;
            }
        }
        if (i < alignmentLength) {
            lastMove = moveCodeToChar[alignment[i]];
            numOfSameMoves++;
        }
    }
    cigar[ci++] = 0;
    if (ci != cigar_length + 1) {
        fprintf(stderr, "CIGAR length error: %d (should be %d)\n", ci, cigar_length+1);
        abort();
    }

    return cigar;
}

/**
 * Build Peq table for given query and alphabet.
 * Peq is table of dimensions alphabetLength+1 x maxNumBlocks.
 * Bit i of Peq[s * maxNumBlocks + b] is 1 if i-th symbol from block b of query equals symbol s, otherwise it is 0.
 * NOTICE: free returned array with delete[]!
 */
static inline Word* buildPeq(void* km, const unsigned char* const query, const int queryLength) {
    int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
    // table of dimensions alphabetLength+1 x maxNumBlocks. Last symbol is wildcard.
    size_t num_words = (kAlphabetLength+1) * maxNumBlocks;
    Word* Peq = Kmalloc(km, Word, num_words);

    // Build Peq (1 is match, 0 is mismatch). NOTE: last column is wildcard(symbol that matches anything) with just 1s
    for (int symbol = 0; symbol <= kAlphabetLength; symbol++) {
        for (int b = 0; b < maxNumBlocks; b++) {
            if (symbol < kAlphabetLength) {
                Peq[symbol * maxNumBlocks + b] = 0;
                for (int r = (b+1) * WORD_SIZE - 1; r >= b * WORD_SIZE; r--) {
                    Peq[symbol * maxNumBlocks + b] <<= 1;
                    // NOTE: We pretend like query is padded at the end with W wildcard symbols
                    if (r >= queryLength || query[r] == symbol)
                        Peq[symbol * maxNumBlocks + b] += 1;
                }
            } else { // Last symbol is wildcard, so it is all 1s
                Peq[symbol * maxNumBlocks + b] = static_cast<Word>(-1);
            }
        }
    }

    return Peq;
}


/**
 * Returns new sequence that is reverse of given sequence.
 * Free returned array with delete[].
 */
static inline unsigned char* createReverseCopy(void* km, const unsigned char* const seq, const int length) {
    unsigned char* rSeq = Kmalloc(km, unsigned char, length);
    copy(seq, seq + length, rSeq);
    reverse(rSeq, rSeq + length);
    return rSeq;
}

/**
 * Corresponds to Advance_Block function from Myers.
 * Calculates one word(block), which is part of a column.
 * Highest bit of word (one most to the left) is most bottom cell of block from column.
 * Pv[i] and Mv[i] define vin of cell[i]: vin = cell[i] - cell[i-1].
 * @param [in] Pv  Bitset, Pv[i] == 1 if vin is +1, otherwise Pv[i] == 0.
 * @param [in] Mv  Bitset, Mv[i] == 1 if vin is -1, otherwise Mv[i] == 0.
 * @param [in] Eq  Bitset, Eq[i] == 1 if match, 0 if mismatch.
 * @param [in] hin  Will be +1, 0 or -1.
 * @param [out] PvOut  Bitset, PvOut[i] == 1 if vout is +1, otherwise PvOut[i] == 0.
 * @param [out] MvOut  Bitset, MvOut[i] == 1 if vout is -1, otherwise MvOut[i] == 0.
 * @param [out] hout  Will be +1, 0 or -1.
 */
static inline int calculateBlock(Word Pv, Word Mv, Word Eq, const int hin,
                                 Word &PvOut, Word &MvOut) {
    // hin can be 1, -1 or 0.
    // 1  -> 00...01
    // 0  -> 00...00
    // -1 -> 11...11 (2-complement)

    Word hinIsNeg = static_cast<Word>(hin >> 2) & WORD_1; // 00...001 if hin is -1, 00...000 if 0 or 1

    Word Xv = Eq | Mv;
    // This is instruction below written using 'if': if (hin < 0) Eq |= (Word)1;
    Eq |= hinIsNeg;
    Word Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq;

    Word Ph = Mv | ~(Xh | Pv);
    Word Mh = Pv & Xh;

    int hout = 0;
    // This is instruction below written using 'if': if (Ph & HIGH_BIT_MASK) hout = 1;
    hout = (Ph & HIGH_BIT_MASK) >> (WORD_SIZE - 1);
    // This is instruction below written using 'if': if (Mh & HIGH_BIT_MASK) hout = -1;
    hout -= (Mh & HIGH_BIT_MASK) >> (WORD_SIZE - 1);

    Ph <<= 1;
    Mh <<= 1;

    // This is instruction below written using 'if': if (hin < 0) Mh |= (Word)1;
    Mh |= hinIsNeg;
    // This is instruction below written using 'if': if (hin > 0) Ph |= (Word)1;
    Ph |= static_cast<Word>((hin + 1) >> 1);

    PvOut = Mh | ~(Xv | Ph);
    MvOut = Ph & Xv;

    return hout;
}

/**
 * Does ceiling division x / y.
 * Note: x and y must be non-negative and x + y must not overflow.
 */
static inline int ceilDiv(const int x, const int y) {
    return x % y ? x / y + 1 : x / y;
}

static inline int min(const int x, const int y) {
    return x < y ? x : y;
}

static inline int max(const int x, const int y) {
    return x > y ? x : y;
}


/**
 * @param [in] block
 * @return Values of cells in block, starting with bottom cell in block.
 */
static inline vector<int> getBlockCellValues(const Block block) {
    vector<int> scores(WORD_SIZE);
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; i++) {
        scores[i] = score;
        if (block.P & mask) score--;
        if (block.M & mask) score++;
        mask >>= 1;
    }
    scores[WORD_SIZE - 1] = score;
    return scores;
}

/**
 * Writes values of cells in block into given array, starting with first/top cell.
 * @param [in] block
 * @param [out] dest  Array into which cell values are written. Must have size of at least WORD_SIZE.
 */
static inline void readBlock(const Block block, int* const dest) {
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; i++) {
        dest[WORD_SIZE - 1 - i] = score;
        if (block.P & mask) score--;
        if (block.M & mask) score++;
        mask >>= 1;
    }
    dest[0] = score;
}

/**
 * Writes values of cells in block into given array, starting with last/bottom cell.
 * @param [in] block
 * @param [out] dest  Array into which cell values are written. Must have size of at least WORD_SIZE.
 */
static inline void readBlockReverse(const Block block, int* const dest) {
    int score = block.score;
    Word mask = HIGH_BIT_MASK;
    for (int i = 0; i < WORD_SIZE - 1; i++) {
        dest[i] = score;
        if (block.P & mask) score--;
        if (block.M & mask) score++;
        mask >>= 1;
    }
    dest[WORD_SIZE - 1] = score;
}

/**
 * @param [in] block
 * @param [in] k
 * @return True if all cells in block have value larger than k, otherwise false.
 */
static inline bool allBlockCellsLarger(const Block block, const int k) {
    vector<int> scores = getBlockCellValues(block);
    for (int i = 0; i < WORD_SIZE; i++) {
        if (scores[i] <= k) return false;
    }
    return true;
}


/**
 * Uses Myers' bit-vector algorithm to find edit distance for one of semi-global alignment methods.
 * @param [in] Peq  Query profile.
 * @param [in] W  Size of padding in last block.
 *                TODO: Calculate this directly from query, instead of passing it.
 * @param [in] maxNumBlocks  Number of blocks needed to cover the whole query.
 *                           TODO: Calculate this directly from query, instead of passing it.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] targetLength
 * @param [in] k
 * @param [in] mode  EDLIB_MODE_HW or EDLIB_MODE_SHW
 * @param [out] bestScore_  Edit distance.
 * @param [out] positions_  Array of 0-indexed positions in target at which best score was found.
                            Make sure to free this array with free().
 * @param [out] numPositions_  Number of positions in the positions_ array.
 * @return Status.
 */
static int myersCalcEditDistanceSemiGlobal(void* km,
        const Word* const Peq, const int W, const int maxNumBlocks,
        const int queryLength,
        const unsigned char* const target, const int targetLength,
        int k, const EdlibAlignMode mode,
        int* const bestScore_, int** const positions_, int* const numPositions_) {
    *positions_ = NULL;
    *numPositions_ = 0;

    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of last block in Ukkonen band.
    int firstBlock = 0;
    int lastBlock = min(ceilDiv(k + 1, WORD_SIZE), maxNumBlocks) - 1; // y in Myers
    Block *bl; // Current block

    Block* blocks = Kmalloc(km, Block, maxNumBlocks);

    // For HW, solution will never be larger then queryLength.
    if (mode == EDLIB_MODE_HW) {
        k = min(queryLength, k);
    }

    // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
    // This gives speed up of about 2 times for small k.
    const int STRONG_REDUCE_NUM = 2048;

    // Initialize P, M and score
    bl = blocks;
    for (int b = 0; b <= lastBlock; b++) {
        bl->score = (b + 1) * WORD_SIZE;
        bl->P = static_cast<Word>(-1); // All 1s
        bl->M = static_cast<Word>(0);
        bl++;
    }

    int bestScore = -1;
    vector<int> positions; // TODO: Maybe put this on heap?
    const int startHout = mode == EDLIB_MODE_HW ? 0 : 1; // If 0 then gap before query is not penalized;
    const unsigned char* targetChar = target;
    for (int c = 0; c < targetLength; c++) { // for each column
        const Word* Peq_c = Peq + (*targetChar) * maxNumBlocks;

        //----------------------- Calculate column -------------------------//
        int hout = startHout;
        bl = blocks + firstBlock;
        Peq_c += firstBlock;
        for (int b = firstBlock; b <= lastBlock; b++) {
            hout = calculateBlock(bl->P, bl->M, *Peq_c, hout, bl->P, bl->M);
            bl->score += hout;
            bl++; Peq_c++;
        }
        bl--; Peq_c--;
        //------------------------------------------------------------------//

        //---------- Adjust number of blocks according to Ukkonen ----------//
        if ((lastBlock < maxNumBlocks - 1) && (bl->score - hout <= k) // bl is pointing to last block
            && ((*(Peq_c + 1) & WORD_1) || hout < 0)) { // Peq_c is pointing to last block
            // If score of left block is not too big, calculate one more block
            lastBlock++; bl++; Peq_c++;
            bl->P = static_cast<Word>(-1); // All 1s
            bl->M = static_cast<Word>(0);
            bl->score = (bl - 1)->score - hout + WORD_SIZE + calculateBlock(bl->P, bl->M, *Peq_c, hout, bl->P, bl->M);
        } else {
            while (lastBlock >= firstBlock && bl->score >= k + WORD_SIZE) {
                lastBlock--; bl--; Peq_c--;
            }
        }

        // Every some columns, do some expensive but also more efficient block reducing.
        // This is important!
        //
        // Reduce the band by decreasing last block if possible.
        if (c % STRONG_REDUCE_NUM == 0) {
            while (lastBlock >= 0 && lastBlock >= firstBlock && allBlockCellsLarger(*bl, k)) {
                lastBlock--; bl--; Peq_c--;
            }
        }
        // For HW, even if all cells are > k, there still may be solution in next
        // column because starting conditions at upper boundary are 0.
        // That means that first block is always candidate for solution,
        // and we can never end calculation before last column.
        if (mode == EDLIB_MODE_HW && lastBlock == -1) {
            lastBlock++; bl++; Peq_c++;
        }

        // Reduce band by increasing first block if possible. Not applicable to HW.
        if (mode != EDLIB_MODE_HW) {
            while (firstBlock <= lastBlock && blocks[firstBlock].score >= k + WORD_SIZE) {
                firstBlock++;
            }
            if (c % STRONG_REDUCE_NUM == 0) { // Do strong reduction every some blocks
                while (firstBlock <= lastBlock && allBlockCellsLarger(blocks[firstBlock], k)) {
                    firstBlock++;
                }
            }
        }

        // If band stops to exist finish
        if (lastBlock < firstBlock) {
            *bestScore_ = bestScore;
            if (bestScore != -1) {
                *numPositions_ = static_cast<int>(positions.size());
                *positions_ = Kmalloc(km, int, *numPositions_);
                copy(positions.begin(), positions.end(), *positions_);
            }
            kfree(km, blocks);
            return EDLIB_STATUS_OK;
        }
        //------------------------------------------------------------------//

        //------------------------- Update best score ----------------------//
        if (lastBlock == maxNumBlocks - 1) {
            int colScore = bl->score;
            if (colScore <= k) { // Scores > k dont have correct values (so we cannot use them), but are certainly > k.
                // NOTE: Score that I find in column c is actually score from column c-W
                if (bestScore == -1 || colScore <= bestScore) {
                    if (colScore != bestScore) {
                        positions.clear();
                        bestScore = colScore;
                        // Change k so we will look only for equal or better
                        // scores then the best found so far.
                        k = bestScore;
                    }
                    positions.push_back(c - W);
                }
            }
        }
        //------------------------------------------------------------------//

        targetChar++;
    }


    // Obtain results for last W columns from last column.
    if (lastBlock == maxNumBlocks - 1) {
        vector<int> blockScores = getBlockCellValues(*bl);
        for (int i = 0; i < W; i++) {
            int colScore = blockScores[i + 1];
            if (colScore <= k && (bestScore == -1 || colScore <= bestScore)) {
                if (colScore != bestScore) {
                    positions.clear();
                    k = bestScore = colScore;
                }
                positions.push_back(targetLength - W + i);
            }
        }
    }

    *bestScore_ = bestScore;
    if (bestScore != -1) {
        *numPositions_ = static_cast<int>(positions.size());
        *positions_ = Kmalloc(km, int, *numPositions_);
        copy(positions.begin(), positions.end(), *positions_);
    }

    kfree(km, blocks);
    return EDLIB_STATUS_OK;
}


/**
 * Uses Myers' bit-vector algorithm to find edit distance for global(NW) alignment method.
 * @param [in] Peq  Query profile.
 * @param [in] W  Size of padding in last block.
 *                TODO: Calculate this directly from query, instead of passing it.
 * @param [in] maxNumBlocks  Number of blocks needed to cover the whole query.
 *                           TODO: Calculate this directly from query, instead of passing it.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] targetLength
 * @param [in] k
 * @param [out] bestScore_  Edit distance.
 * @param [out] position_  0-indexed position in target at which best score was found.
 * @param [in] findAlignment  If true, whole matrix is remembered and alignment data is returned.
 *                            Quadratic amount of memory is consumed.
 * @param [out] alignData  Data needed for alignment traceback (for reconstruction of alignment).
 *                         Set only if findAlignment is set to true, otherwise it is NULL.
 *                         Make sure to free this array using delete[].
 * @param [out] targetStopPosition  If set to -1, whole calculation is performed normally, as expected.
 *         If set to p, calculation is performed up to position p in target (inclusive)
 *         and column p is returned as the only column in alignData.
 * @return Status.
 */
static int myersCalcEditDistanceNW(void* km, const Word* const Peq, const int W, const int maxNumBlocks,
                                   const int queryLength,
                                   const unsigned char* const target, const int targetLength,
                                   int k, int* const bestScore_,
                                   int* const position_, const bool findAlignment,
                                   AlignmentData** const alignData, const int targetStopPosition) {
    if (targetStopPosition > -1 && findAlignment) {
        // They can not be both set at the same time!
        return EDLIB_STATUS_ERROR;
    }

    // Each STRONG_REDUCE_NUM column is reduced in more expensive way.
    const int STRONG_REDUCE_NUM = 2048; // TODO: Choose this number dinamically (based on query and target lengths?), so it does not affect speed of computation

    if (k < abs(targetLength - queryLength)) {
        *bestScore_ = *position_ = -1;
        return EDLIB_STATUS_OK;
    }

    k = min(k, max(queryLength, targetLength));  // Upper bound for k

    // firstBlock is 0-based index of first block in Ukkonen band.
    // lastBlock is 0-based index of last block in Ukkonen band.
    int firstBlock = 0;
    // This is optimal now, by my formula.
    int lastBlock = min(maxNumBlocks, ceilDiv(min(k, (k + queryLength - targetLength) / 2) + 1, WORD_SIZE)) - 1;
    Block* bl; // Current block

    Block* blocks = Kmalloc(km, Block, maxNumBlocks); 

    // Initialize P, M and score
    bl = blocks;
    for (int b = 0; b <= lastBlock; b++) {
        bl->score = (b + 1) * WORD_SIZE;
        bl->P = static_cast<Word>(-1); // All 1s
        bl->M = static_cast<Word>(0);
        bl++;
    }

    // If we want to find alignment, we have to store needed data.
    if (findAlignment)
        *alignData = AlignmentDataNew(km, maxNumBlocks, targetLength);
    else if (targetStopPosition > -1)
        *alignData = AlignmentDataNew(km, maxNumBlocks, 1);
    else
        *alignData = NULL;

    const unsigned char* targetChar = target;
    for (int c = 0; c < targetLength; c++) { // for each column
        const Word* Peq_c = Peq + *targetChar * maxNumBlocks;

        //----------------------- Calculate column -------------------------//
        int hout = 1;
        bl = blocks + firstBlock;
        for (int b = firstBlock; b <= lastBlock; b++) {
            hout = calculateBlock(bl->P, bl->M, Peq_c[b], hout, bl->P, bl->M);
            bl->score += hout;
            bl++;
        }
        bl--;
        //------------------------------------------------------------------//
        // bl now points to last block

        // Update k. I do it only on end of column because it would slow calculation too much otherwise.
        // NOTICE: I add W when in last block because it is actually result from W cells to the left and W cells up.
        k = min(k, bl->score
                + max(targetLength - c - 1, queryLength - ((1 + lastBlock) * WORD_SIZE - 1) - 1)
                + (lastBlock == maxNumBlocks - 1 ? W : 0));

        //---------- Adjust number of blocks according to Ukkonen ----------//
        //--- Adjust last block ---//
        // If block is not beneath band, calculate next block. Only next because others are certainly beneath band.
        if (lastBlock + 1 < maxNumBlocks
            && !(//score[lastBlock] >= k + WORD_SIZE ||  // NOTICE: this condition could be satisfied if above block also!
                 ((lastBlock + 1) * WORD_SIZE - 1
                  > k - bl->score + 2 * WORD_SIZE - 2 - targetLength + c + queryLength))) {
            lastBlock++; bl++;
            bl->P = static_cast<Word>(-1); // All 1s
            bl->M = static_cast<Word>(0);
            int newHout = calculateBlock(bl->P, bl->M, Peq_c[lastBlock], hout, bl->P, bl->M);
            bl->score = (bl - 1)->score - hout + WORD_SIZE + newHout;
            hout = newHout;
        }

        // While block is out of band, move one block up.
        // NOTE: Condition used here is more loose than the one from the article, since I simplified the max() part of it.
        // I could consider adding that max part, for optimal performance.
        while (lastBlock >= firstBlock
               && (bl->score >= k + WORD_SIZE
                   || ((lastBlock + 1) * WORD_SIZE - 1 >
                       // TODO: Does not work if do not put +1! Why???
                       k - bl->score + 2 * WORD_SIZE - 2 - targetLength + c + queryLength + 1))) {
            lastBlock--; bl--;
        }
        //-------------------------//

        //--- Adjust first block ---//
        // While outside of band, advance block
        while (firstBlock <= lastBlock
               && (blocks[firstBlock].score >= k + WORD_SIZE
                   || ((firstBlock + 1) * WORD_SIZE - 1 <
                       blocks[firstBlock].score - k - targetLength + queryLength + c))) {
            firstBlock++;
        }
        //--------------------------/


        // TODO: consider if this part is useful, it does not seem to help much
        if (c % STRONG_REDUCE_NUM == 0) { // Every some columns do more expensive but more efficient reduction
            while (lastBlock >= firstBlock) {
                // If all cells outside of band, remove block
                vector<int> scores = getBlockCellValues(*bl);
                int numCells = lastBlock == maxNumBlocks - 1 ? WORD_SIZE - W : WORD_SIZE;
                int r = lastBlock * WORD_SIZE + numCells - 1;
                bool reduce = true;
                for (int i = WORD_SIZE - numCells; i < WORD_SIZE; i++) {
                    // TODO: Does not work if do not put +1! Why???
                    if (scores[i] <= k && r <= k - scores[i] - targetLength + c + queryLength + 1) {
                        reduce = false;
                        break;
                    }
                    r--;
                }
                if (!reduce) break;
                lastBlock--; bl--;
            }

            while (firstBlock <= lastBlock) {
                // If all cells outside of band, remove block
                vector<int> scores = getBlockCellValues(blocks[firstBlock]);
                int numCells = firstBlock == maxNumBlocks - 1 ? WORD_SIZE - W : WORD_SIZE;
                int r = firstBlock * WORD_SIZE + numCells - 1;
                bool reduce = true;
                for (int i = WORD_SIZE - numCells; i < WORD_SIZE; i++) {
                    if (scores[i] <= k && r >= scores[i] - k - targetLength + c + queryLength) {
                        reduce = false;
                        break;
                    }
                    r--;
                }
                if (!reduce) break;
                firstBlock++;
            }
        }


        // If band stops to exist finish
        if (lastBlock < firstBlock) {
            *bestScore_ = *position_ = -1;
            kfree(km, blocks);
            return EDLIB_STATUS_OK;
        }
        //------------------------------------------------------------------//


        //---- Save column so it can be used for reconstruction ----//
        if (findAlignment && c < targetLength) {
            bl = blocks + firstBlock;
            for (int b = firstBlock; b <= lastBlock; b++) {
                (*alignData)->Ps[maxNumBlocks * c + b] = bl->P;
                (*alignData)->Ms[maxNumBlocks * c + b] = bl->M;
                (*alignData)->scores[maxNumBlocks * c + b] = bl->score;
                bl++;
            }
            (*alignData)->firstBlocks[c] = firstBlock;
            (*alignData)->lastBlocks[c] = lastBlock;
        }
        //----------------------------------------------------------//
        //---- If this is stop column, save it and finish ----//
        if (c == targetStopPosition) {
            for (int b = firstBlock; b <= lastBlock; b++) {
                (*alignData)->Ps[b] = (blocks + b)->P;
                (*alignData)->Ms[b] = (blocks + b)->M;
                (*alignData)->scores[b] = (blocks + b)->score;
            }
            (*alignData)->firstBlocks[0] = firstBlock;
            (*alignData)->lastBlocks[0] = lastBlock;
            *bestScore_ = -1;
            *position_ = targetStopPosition;
            kfree(km, blocks);
            return EDLIB_STATUS_OK;
        }
        //----------------------------------------------------//

        targetChar++;
    }

    if (lastBlock == maxNumBlocks - 1) { // If last block of last column was calculated
        // Obtain best score from block -> it is complicated because query is padded with W cells
        int bestScore = getBlockCellValues(blocks[lastBlock])[W];
        if (bestScore <= k) {
            *bestScore_ = bestScore;
            *position_ = targetLength - 1;
            kfree(km, blocks);
            return EDLIB_STATUS_OK;
        }
    }

    *bestScore_ = *position_ = -1;
    kfree(km, blocks);
    return EDLIB_STATUS_OK;
}


/**
 * Finds one possible alignment that gives optimal score by moving back through the dynamic programming matrix,
 * that is stored in alignData. Consumes large amount of memory: O(queryLength * targetLength).
 * @param [in] queryLength  Normal length, without W.
 * @param [in] targetLength  Normal length, without W.
 * @param [in] bestScore  Best score.
 * @param [in] alignData  Data obtained during finding best score that is useful for finding alignment.
 * @param [out] alignment  Alignment.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentTraceback(void* km, const int queryLength, const int targetLength,
                                    const int bestScore, const AlignmentData* const alignData,
                                    unsigned char** const alignment, int* const alignmentLength) {
    const int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
    const int W = maxNumBlocks * WORD_SIZE - queryLength;

    int maxAlignmentLength = queryLength + targetLength - 1;
    *alignment = Kmalloc(km, unsigned char, maxAlignmentLength);
    *alignmentLength = 0;
    int c = targetLength - 1; // index of column
    int b = maxNumBlocks - 1; // index of block in column
    int currScore = bestScore; // Score of current cell
    int lScore  = -1; // Score of left cell
    int uScore  = -1; // Score of upper cell
    int ulScore = -1; // Score of upper left cell
    Word currP = alignData->Ps[c * maxNumBlocks + b]; // P of current block
    Word currM = alignData->Ms[c * maxNumBlocks + b]; // M of current block
    // True if block to left exists and is in band
    bool thereIsLeftBlock = c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1];
    // We set initial values of lP and lM to 0 only to avoid compiler warnings, they should not affect the
    // calculation as both lP and lM should be initialized at some moment later (but compiler can not
    // detect it since this initialization is guaranteed by "business" logic).
    Word lP = 0, lM = 0;
    if (thereIsLeftBlock) {
        lP = alignData->Ps[(c - 1) * maxNumBlocks + b]; // P of block to the left
        lM = alignData->Ms[(c - 1) * maxNumBlocks + b]; // M of block to the left
    }
    currP <<= W;
    currM <<= W;
    int blockPos = WORD_SIZE - W - 1; // 0 based index of current cell in blockPos

    // TODO(martin): refactor this whole piece of code. There are too many if-else statements,
    // it is too easy for a bug to hide and to hard to effectively cover all the edge-cases.
    // We need better separation of logic and responsibilities.
    while (true) {
        if (c == 0) {
            thereIsLeftBlock = true;
            lScore = b * WORD_SIZE + blockPos + 1;
            ulScore = lScore - 1;
        }

        // TODO: improvement: calculate only those cells that are needed,
        //       for example if I calculate upper cell and can move up,
        //       there is no need to calculate left and upper left cell
        //---------- Calculate scores ---------//
        if (lScore == -1 && thereIsLeftBlock) {
            lScore = alignData->scores[(c - 1) * maxNumBlocks + b]; // score of block to the left
            for (int i = 0; i < WORD_SIZE - blockPos - 1; i++) {
                if (lP & HIGH_BIT_MASK) lScore--;
                if (lM & HIGH_BIT_MASK) lScore++;
                lP <<= 1;
                lM <<= 1;
            }
        }
        if (ulScore == -1) {
            if (lScore != -1) {
                ulScore = lScore;
                if (lP & HIGH_BIT_MASK) ulScore--;
                if (lM & HIGH_BIT_MASK) ulScore++;
            }
            else if (c > 0 && b-1 >= alignData->firstBlocks[c-1] && b-1 <= alignData->lastBlocks[c-1]) {
                // This is the case when upper left cell is last cell in block,
                // and block to left is not in band so lScore is -1.
                ulScore = alignData->scores[(c - 1) * maxNumBlocks + b - 1];
            }
        }
        if (uScore == -1) {
            uScore = currScore;
            if (currP & HIGH_BIT_MASK) uScore--;
            if (currM & HIGH_BIT_MASK) uScore++;
            currP <<= 1;
            currM <<= 1;
        }
        //-------------------------------------//

        // TODO: should I check if there is upper block?

        //-------------- Move --------------//
        // Move up - insertion to target - deletion from query
        if (uScore != -1 && uScore + 1 == currScore) {
            currScore = uScore;
            lScore = ulScore;
            uScore = ulScore = -1;
            if (blockPos == 0) { // If entering new (upper) block
                if (b == 0) { // If there are no cells above (only boundary cells)
                    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT; // Move up
                    for (int i = 0; i < c + 1; i++) // Move left until end
                        (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
                    break;
                } else {
                    blockPos = WORD_SIZE - 1;
                    b--;
                    currP = alignData->Ps[c * maxNumBlocks + b];
                    currM = alignData->Ms[c * maxNumBlocks + b];
                    if (c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1]) {
                        thereIsLeftBlock = true;
                        lP = alignData->Ps[(c - 1) * maxNumBlocks + b]; // TODO: improve this, too many operations
                        lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
                    } else {
                        thereIsLeftBlock = false;
                        // TODO(martin): There may not be left block, but there can be left boundary - do we
                        // handle this correctly then? Are l and ul score set correctly? I should check that / refactor this.
                    }
                }
            } else {
                blockPos--;
                lP <<= 1;
                lM <<= 1;
            }
            // Mark move
            (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
        }
        // Move left - deletion from target - insertion to query
        else if (lScore != -1 && lScore + 1 == currScore) {
            currScore = lScore;
            uScore = ulScore;
            lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE; // Move left
                int numUp = b * WORD_SIZE + blockPos + 1;
                for (int i = 0; i < numUp; i++) // Move up until end
                    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                break;
            }
            currP = lP;
            currM = lM;
            if (c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1]) {
                thereIsLeftBlock = true;
                lP = alignData->Ps[(c - 1) * maxNumBlocks + b];
                lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = true;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = false;
                }
            }
            // Mark move
            (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
        }
        // Move up left - (mis)match
        else if (ulScore != -1) {
            unsigned char moveCode = ulScore == currScore ? EDLIB_EDOP_MATCH : EDLIB_EDOP_MISMATCH;
            currScore = ulScore;
            uScore = lScore = ulScore = -1;
            c--;
            if (c == -1) { // If there are no cells to the left (only boundary cells)
                (*alignment)[(*alignmentLength)++] = moveCode; // Move left
                int numUp = b * WORD_SIZE + blockPos;
                for (int i = 0; i < numUp; i++) // Move up until end
                    (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_INSERT;
                break;
            }
            if (blockPos == 0) { // If entering upper left block
                if (b == 0) { // If there are no more cells above (only boundary cells)
                    (*alignment)[(*alignmentLength)++] = moveCode; // Move up left
                    for (int i = 0; i < c + 1; i++) // Move left until end
                        (*alignment)[(*alignmentLength)++] = EDLIB_EDOP_DELETE;
                    break;
                }
                blockPos = WORD_SIZE - 1;
                b--;
                currP = alignData->Ps[c * maxNumBlocks + b];
                currM = alignData->Ms[c * maxNumBlocks + b];
            } else { // If entering left block
                blockPos--;
                currP = lP;
                currM = lM;
                currP <<= 1;
                currM <<= 1;
            }
            // Set new left block
            if (c > 0 && b >= alignData->firstBlocks[c-1] && b <= alignData->lastBlocks[c-1]) {
                thereIsLeftBlock = true;
                lP = alignData->Ps[(c - 1) * maxNumBlocks + b];
                lM = alignData->Ms[(c - 1) * maxNumBlocks + b];
            } else {
                if (c == 0) { // If there are no cells to the left (only boundary cells)
                    thereIsLeftBlock = true;
                    lScore = b * WORD_SIZE + blockPos + 1;
                    ulScore = lScore - 1;
                } else {
                    thereIsLeftBlock = false;
                }
            }
            // Mark move
            (*alignment)[(*alignmentLength)++] = moveCode;
        } else {
            // Reached end - finished!
            break;
        }
        //----------------------------------//
    }

    reverse(*alignment, *alignment + (*alignmentLength));
    return EDLIB_STATUS_OK;
}


/**
 * Finds one possible alignment that gives optimal score (bestScore).
 * It will split problem into smaller problems using Hirschberg's algorithm and when they are small enough,
 * it will solve them using traceback algorithm.
 * @param [in] query
 * @param [in] rQuery  Reversed query.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] rTarget  Reversed target.
 * @param [in] targetLength
 * @param [in] equalityDefinition
 * @param [in] alphabetLength
 * @param [in] bestScore  Best(optimal) score.
 * @param [out] alignment  Sequence of edit operations that make target equal to query.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignment(void* km,
        const unsigned char* const query, const unsigned char* const rQuery, const int queryLength,
        const unsigned char* const target, const unsigned char* const rTarget, const int targetLength,
        const int bestScore, unsigned char** const alignment, int* const alignmentLength) {

    // Handle special case when one of sequences has length of 0.
    if (queryLength == 0 || targetLength == 0) {
        *alignmentLength = targetLength + queryLength;
        *alignment = Kmalloc(km, unsigned char, *alignmentLength);
        for (int i = 0; i < *alignmentLength; i++) {
            (*alignment)[i] = queryLength == 0 ? EDLIB_EDOP_DELETE : EDLIB_EDOP_INSERT;
        }
        return EDLIB_STATUS_OK;
    }

    const int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
    const int W = maxNumBlocks * WORD_SIZE - queryLength;
    int statusCode;

    // TODO: think about reducing number of memory allocations in alignment functions, probably
    // by sharing some memory that is allocated only once. That refers to: Peq, columns in Hirschberg,
    // and it could also be done for alignments - we could have one big array for alignment that would be
    // sparsely populated by each of steps in recursion, and at the end we would just consolidate those results.

    // If estimated memory consumption for traceback algorithm is smaller than 1MB use it,
    // otherwise use Hirschberg's algorithm. By running few tests I choose boundary of 1MB as optimal.
    long long alignmentDataSize = (2ll * sizeof(Word) + sizeof(int)) * maxNumBlocks * targetLength
        + 2ll * sizeof(int) * targetLength;
    if (alignmentDataSize < 1024 * 1024) {
        int score_, endLocation_;  // Used only to call function.
        AlignmentData* alignData = NULL;
        Word* Peq = buildPeq(km, query, queryLength);
        myersCalcEditDistanceNW(km, Peq, W, maxNumBlocks,
                                queryLength,
                                target, targetLength,
                                bestScore,
                                &score_, &endLocation_, true, &alignData, -1);
        //assert(score_ == bestScore);
        //assert(endLocation_ == targetLength - 1);

        statusCode = obtainAlignmentTraceback(km, queryLength, targetLength,
                                              bestScore, alignData, alignment, alignmentLength);
        AlignmentDataFree(km, alignData);
        kfree(km, Peq);
    } else {
        statusCode = obtainAlignmentHirschberg(km, query, rQuery, queryLength,
                                               target, rTarget, targetLength,
                                               bestScore, alignment, alignmentLength);
    }
    return statusCode;
}


/**
 * Finds one possible alignment that gives optimal score (bestScore).
 * Uses Hirschberg's algorithm to split problem into two sub-problems, solve them and combine them together.
 * @param [in] query
 * @param [in] rQuery  Reversed query.
 * @param [in] queryLength
 * @param [in] target
 * @param [in] rTarget  Reversed target.
 * @param [in] targetLength
 * @param [in] alphabetLength
 * @param [in] bestScore  Best(optimal) score.
 * @param [out] alignment  Sequence of edit operations that make target equal to query.
 * @param [out] alignmentLength  Length of alignment.
 * @return Status code.
 */
static int obtainAlignmentHirschberg(void* km,
        const unsigned char* const query, const unsigned char* const rQuery, const int queryLength,
        const unsigned char* const target, const unsigned char* const rTarget, const int targetLength,
        const int bestScore, unsigned char** const alignment, int* const alignmentLength) {

    const int maxNumBlocks = ceilDiv(queryLength, WORD_SIZE);
    const int W = maxNumBlocks * WORD_SIZE - queryLength;

    Word* Peq = buildPeq(km, query, queryLength);
    Word* rPeq = buildPeq(km, rQuery, queryLength);

    // Used only to call functions.
    int score_, endLocation_;

    // Divide dynamic matrix into two halfs, left and right.
    const int leftHalfWidth = targetLength / 2;
    const int rightHalfWidth = targetLength - leftHalfWidth;

    // Calculate left half.
    AlignmentData* alignDataLeftHalf = NULL;
    int leftHalfCalcStatus = myersCalcEditDistanceNW(km,
            Peq, W, maxNumBlocks, queryLength, target, targetLength, bestScore,
            &score_, &endLocation_, false, &alignDataLeftHalf, leftHalfWidth - 1);

    // Calculate right half.
    AlignmentData* alignDataRightHalf = NULL;
    int rightHalfCalcStatus = myersCalcEditDistanceNW(km,
            rPeq, W, maxNumBlocks, queryLength, rTarget, targetLength, bestScore,
            &score_, &endLocation_, false, &alignDataRightHalf, rightHalfWidth - 1);

    kfree(km, Peq);
    kfree(km, rPeq);

    if (leftHalfCalcStatus == EDLIB_STATUS_ERROR || rightHalfCalcStatus == EDLIB_STATUS_ERROR) {
        if (alignDataLeftHalf) AlignmentDataFree(km, alignDataLeftHalf);
        if (alignDataRightHalf) AlignmentDataFree(km, alignDataRightHalf);
        return EDLIB_STATUS_ERROR;
    }

    // Unwrap the left half.
    int firstBlockIdxLeft = alignDataLeftHalf->firstBlocks[0];
    int lastBlockIdxLeft = alignDataLeftHalf->lastBlocks[0];
    // TODO: avoid this allocation by using some shared array?
    // scoresLeft contains scores from left column, starting with scoresLeftStartIdx row (query index)
    // and ending with scoresLeftEndIdx row (0-indexed).
    int scoresLeftLength = (lastBlockIdxLeft - firstBlockIdxLeft + 1) * WORD_SIZE;
    int* scoresLeft = Kmalloc(km, int, scoresLeftLength);
    for (int blockIdx = firstBlockIdxLeft; blockIdx <= lastBlockIdxLeft; blockIdx++) {
        Block block(alignDataLeftHalf->Ps[blockIdx], alignDataLeftHalf->Ms[blockIdx],
                    alignDataLeftHalf->scores[blockIdx]);
        readBlock(block, scoresLeft + (blockIdx - firstBlockIdxLeft) * WORD_SIZE);
    }
    int scoresLeftStartIdx = firstBlockIdxLeft * WORD_SIZE;
    // If last block contains padding, shorten the length of scores for the length of padding.
    if (lastBlockIdxLeft == maxNumBlocks - 1) {
        scoresLeftLength -= W;
    }

    // Unwrap the right half (I also reverse it while unwraping).
    int firstBlockIdxRight = alignDataRightHalf->firstBlocks[0];
    int lastBlockIdxRight = alignDataRightHalf->lastBlocks[0];
    int scoresRightLength = (lastBlockIdxRight - firstBlockIdxRight + 1) * WORD_SIZE;
    int* scoresRight = Kmalloc(km, int, scoresRightLength);
    int* scoresRightOriginalStart = scoresRight;
    for (int blockIdx = firstBlockIdxRight; blockIdx <= lastBlockIdxRight; blockIdx++) {
        Block block(alignDataRightHalf->Ps[blockIdx], alignDataRightHalf->Ms[blockIdx],
                    alignDataRightHalf->scores[blockIdx]);
        readBlockReverse(block, scoresRight + (lastBlockIdxRight - blockIdx) * WORD_SIZE);
    }
    int scoresRightStartIdx = queryLength - (lastBlockIdxRight + 1) * WORD_SIZE;
    // If there is padding at the beginning of scoresRight (that can happen because of reversing that we do),
    // move pointer forward to remove the padding (that is why we remember originalStart).
    if (scoresRightStartIdx < 0) {
        //assert(scoresRightStartIdx == -1 * W);
        scoresRight += W;
        scoresRightStartIdx += W;
        scoresRightLength -= W;
    }

    AlignmentDataFree(km, alignDataLeftHalf);
    AlignmentDataFree(km, alignDataRightHalf);

    //--------------------- Find the best move ----------------//
    // Find the query/row index of cell in left column which together with its lower right neighbour
    // from right column gives the best score (when summed). We also have to consider boundary cells
    // (those cells at -1 indexes).
    //  x|
    //  -+-
    //   |x
    int queryIdxLeftStart = max(scoresLeftStartIdx, scoresRightStartIdx - 1);
    int queryIdxLeftEnd = min(scoresLeftStartIdx + scoresLeftLength - 1,
                          scoresRightStartIdx + scoresRightLength - 2);
    int leftScore = -1, rightScore = -1;
    int queryIdxLeftAlignment = -1;  // Query/row index of cell in left column where alignment is passing through.
    bool queryIdxLeftAlignmentFound = false;
    for (int queryIdx = queryIdxLeftStart; queryIdx <= queryIdxLeftEnd; queryIdx++) {
        leftScore = scoresLeft[queryIdx - scoresLeftStartIdx];
        rightScore = scoresRight[queryIdx + 1 - scoresRightStartIdx];
        if (leftScore + rightScore == bestScore) {
            queryIdxLeftAlignment = queryIdx;
            queryIdxLeftAlignmentFound = true;
            break;
        }
    }
    // Check boundary cells.
    if (!queryIdxLeftAlignmentFound && scoresLeftStartIdx == 0 && scoresRightStartIdx == 0) {
        leftScore = leftHalfWidth;
        rightScore = scoresRight[0];
        if (leftScore + rightScore == bestScore) {
            queryIdxLeftAlignment = -1;
            queryIdxLeftAlignmentFound = true;
        }
    }
    if (!queryIdxLeftAlignmentFound && scoresLeftStartIdx + scoresLeftLength == queryLength
        && scoresRightStartIdx + scoresRightLength == queryLength) {
        leftScore = scoresLeft[scoresLeftLength - 1];
        rightScore = rightHalfWidth;
        if (leftScore + rightScore == bestScore) {
            queryIdxLeftAlignment = queryLength - 1;
            queryIdxLeftAlignmentFound = true;
        }
    }

    kfree(km, scoresLeft);
    kfree(km, scoresRightOriginalStart);


    if (queryIdxLeftAlignmentFound == false) {
        // If there was no move that is part of optimal alignment, then there is no such alignment
        // or given bestScore is not correct!
        return EDLIB_STATUS_ERROR;
    }
    //----------------------------------------------------------//

    // Calculate alignments for upper half of left half (upper left - ul)
    // and lower half of right half (lower right - lr).
    const int ulHeight = queryIdxLeftAlignment + 1;
    const int lrHeight = queryLength - ulHeight;
    const int ulWidth = leftHalfWidth;
    const int lrWidth = rightHalfWidth;
    unsigned char* ulAlignment = NULL; int ulAlignmentLength;
    int ulStatusCode = obtainAlignment(km, query, rQuery + lrHeight, ulHeight,
                                       target, rTarget + lrWidth, ulWidth,
                                       leftScore, &ulAlignment, &ulAlignmentLength);
    unsigned char* lrAlignment = NULL; int lrAlignmentLength;
    int lrStatusCode = obtainAlignment(km, query + ulHeight, rQuery, lrHeight,
                                       target + ulWidth, rTarget, lrWidth,
                                       rightScore, &lrAlignment, &lrAlignmentLength);
    if (ulStatusCode == EDLIB_STATUS_ERROR || lrStatusCode == EDLIB_STATUS_ERROR) {
        if (ulAlignment) kfree(km, ulAlignment);
        if (lrAlignment) kfree(km, lrAlignment);
        return EDLIB_STATUS_ERROR;
    }

    // Build alignment by concatenating upper left alignment with lower right alignment.
    *alignmentLength = ulAlignmentLength + lrAlignmentLength;
    *alignment = Kmalloc(km, unsigned char, *alignmentLength);
    memcpy(*alignment, ulAlignment, ulAlignmentLength);
    memcpy(*alignment + ulAlignmentLength, lrAlignment, lrAlignmentLength);

    kfree(km, ulAlignment);
    kfree(km, lrAlignment);
    return EDLIB_STATUS_OK;
}

extern "C" EdlibAlignConfig edlibNewAlignConfig(int k, EdlibAlignMode mode, EdlibAlignTask task) {
    EdlibAlignConfig config;
    config.k = k;
    config.mode = mode;
    config.task = task;
    return config;
}

extern "C" EdlibAlignConfig edlibDefaultAlignConfig(void) {
    return edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE);
}

extern "C" void edlibFreeAlignResult(void* km, EdlibAlignResult result) {
    if (result.endLocations) kfree(km, result.endLocations);
    if (result.startLocations) kfree(km, result.startLocations);
    if (result.alignment) kfree(km, result.alignment);
}

extern "C" void edlibAlignmentToStrings(const unsigned char* alignment, 
							 int alignmentLength, 
							 int tgtStart, int tgtEnd, 
							 int qryStart, int qryEnd, 
							 const uint8_t *tgt, 
							 const uint8_t *qry, char *tgt_aln_str, 
							 char *qry_aln_str) {
	for (int a = 0, qryPos=qryStart, tgtPos=tgtStart; a < alignmentLength; a++) {
		assert(qryPos <= qryEnd && tgtPos <= tgtEnd);
		if (alignment[a] == EDLIB_EDOP_MATCH || alignment[a] == EDLIB_EDOP_MISMATCH) {  // match or mismatch
            int qc = qry[qryPos]; qc = "ACGT"[qc]; 
			qry_aln_str[a] = qc;
            int tc = tgt[tgtPos]; tc = "ACGT"[tc];
			tgt_aln_str[a] = tc;
            assert(qry_aln_str[a] != '\0');
            assert(tgt_aln_str[a] != '\0');
			qryPos++;
			tgtPos++;
		} else if (alignment[a] == EDLIB_EDOP_INSERT) { // insertion in target
			tgt_aln_str[a] = '-';
            int qc = qry[qryPos]; qc = "ACGT"[qc]; 
			qry_aln_str[a] = qc;
            assert(qry_aln_str[a] != '\0');
			qryPos++;
		} else if (alignment[a] == EDLIB_EDOP_DELETE) { // insertion in query
            int tc = tgt[tgtPos]; tc = "ACGT"[tc];
			tgt_aln_str[a] = tc;
            assert(tgt_aln_str[a] != '\0');
			qry_aln_str[a] = '-';
			tgtPos++;
		}
	}
	qry_aln_str[alignmentLength] = tgt_aln_str[alignmentLength] = '\0';
    //fprintf(stderr, "qlen = %zu\n", strlen(qry_aln_str));
    //fprintf(stderr, "tlen = %zu\n", strlen(tgt_aln_str));
	assert(strlen(qry_aln_str) == (size_t)alignmentLength && strlen(tgt_aln_str) == (size_t)alignmentLength);
}
