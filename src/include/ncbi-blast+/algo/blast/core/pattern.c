/* $Id: pattern.c 134303 2008-07-17 17:42:49Z camacho $
 * ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================
 *
 * Author: Ilya Dondoshansky
 *
 */

/** @file pattern.c
 * Functions for finding pattern matches in sequence.
 * The following functions are defined here. @sa phi_lookup.h
 *
 * <pre>
 * SPHIQueryInfoNew, SPHIQueryInfoFree, SPHIQueryInfoCopy - life cycle functions
 * for the SPHIQueryInfo structure for saving pattern occurrences in query. 
 * 
 * Main API function to find and save pattern occurrences in query, and functions 
 * called from it:
 *
 * PHIGetPatternOccurrences
 *     FindPatternHits
 *         if ( pattern fits into a single word)
 *             s_FindHitsShortHead
 *         else if ( pattern fits into several words )
 *             s_FindHitsLong
 *         else if ( pattern contains parts longer than a word )
 *             s_FindHitsVeryLong
 *                 calls s_FindHitsShortHead for every word and extends them
 * 
 * For pattern occurrences in subject (database), 
 * FindPatternHits is called from PHIBlastScanSubject.
 * </pre>
 *
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: pattern.c 134303 2008-07-17 17:42:49Z camacho $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include <algo/blast/core/pattern.h>
#include "pattern_priv.h"

void 
_PHIGetRightOneBits(Int4 s, Int4 mask, Int4* rightOne, Int4* rightMaskOnly)
{
    const int kCheckingMatches = s & mask;  /*look for 1 bits in same position*/
    Int4 right_index; /*loop index looking for 1 in kCheckingMatches*/
    Int4 left_index; /*rightmost bit that is 1 in the mask only*/

    left_index = -1;

    for (right_index = 0; right_index < PHI_BITS_PACKED_PER_WORD; right_index++) {
       if ((kCheckingMatches >> right_index) % 2  == 1) 
           break;
       if ((mask >> right_index) %2  == 1) 
          left_index = right_index;
    }
    /* If there was no break from the loop, s and mask have no 1 bits in same
       position, so rightOne should be 0. */
    if (right_index == PHI_BITS_PACKED_PER_WORD)
        right_index = 0;

    *rightOne = right_index;
    *rightMaskOnly = left_index;
}

/** Looks for 1 bits in the same position of s and mask
 * Let R be the rightmost position where s and mask both have a 1.
 * Let L < R be the rightmost position where mask has a 1, if any, 
 * or -1 otherwise.
 * @param s Some number [in]
 * @param mask Mask [in]
 * @return (R - L). 
 */
static Int4 
s_LenOf(Int4 s, Int4 mask)
{
    Int4 rightOne; /*loop index looking for 1 in kCheckingMatches*/
    Int4 rightMaskOnly; /*rightmost bit that is 1 in the mask only*/

    _PHIGetRightOneBits(s, mask, &rightOne, &rightMaskOnly);

    return (rightOne - rightMaskOnly);
}

Int4 
_PHIBlastFindHitsShort(Int4 *hitArray, const Uint1* seq, Int4 len1, 
                      const SPHIPatternSearchBlk *pattern_blk)
{
    Int4 i; /*loop index on sequence*/
    Int4 prefixMatchedBitPattern = 0; /*indicates where pattern aligns
                 with seq; e.g., if value is 9 = 0101 then 
                 last 3 chars of seq match first 3 positions in pattern
                 and last 1 char of seq matches 1 position of pattern*/
    Int4 numMatches = 0; /*number of matches found*/
    Int4 mask;  /*mask of input pattern positions after which
                  a match can be declared*/
    Int4 maskShiftPlus1; /*mask shifted left 1 plus 1 */
    const SShortPatternItems* pattern_items = pattern_blk->one_word_items;

    mask = pattern_items->match_mask; 
    maskShiftPlus1 = (mask << 1) + 1;
    for (i = 0; i < len1; i++) {
      /*shift the positions matched by 1 and try to match up against
        the next character, also allow next character to match the
        first position*/
      prefixMatchedBitPattern =  
         ((prefixMatchedBitPattern << 1) | maskShiftPlus1) & 
         pattern_items->whichPositionPtr[seq[i]];
      if (prefixMatchedBitPattern & mask) { 
         /*first part of pair is index of place in seq where match
           ends; second part is where match starts*/
         hitArray[numMatches++] = i;
         hitArray[numMatches++] = i - s_LenOf(prefixMatchedBitPattern, mask)+1;
         if (numMatches == PHI_MAX_HIT)
         {
             /** @todo FIXME: Pass back an error message saying that
                 numMatches matches are saved, others discarded. */
             break;
         }
      }
    }
    return numMatches;
}

/** Find hits when sequence is DNA and pattern is short returns twice the number
 * of hits.
 * @param hitArray Array of hits to pass back [out]
 * @param seq The input sequence [in]
 * @param pos Starting position [in]
 * @param len Length of sequence seq [in]
 * @param pattern_blk Carries variables that keep track of search 
 *                      parameters. [in]
 * @return Number of hits found.
 */
static Int4 
s_FindHitsShortDNA(Int4* hitArray, const Uint1* seq, Int4 pos, Int4 len,
                   const SPHIPatternSearchBlk *pattern_blk)
{
    Uint4 prefixMatchedBitPattern; /*indicates where pattern aligns
                                     with sequence*/
    Uint4 tmp; /*intermediate result of masked comparisons*/
    Int4 i; /*index on seq*/
    Int4 end; /*count of number of 4-mer iterations needed*/
    Int4 remain; /*0,1,2,3 DNA letters left over*/
    Int4 j; /*index on suffixRemnant*/
    Int4 twiceNumHits = 0; /*twice the number of hits*/
    const SShortPatternItems* pattern_items = pattern_blk->one_word_items;
    const Int4 kMatchMask = pattern_items->match_mask;
    /* Mask to match agaist */
    const Uint4 kMask2 = kMatchMask*PHI_BITS_PACKED_PER_WORD+15; 
    const Int4 kMaskShiftPlus1 = (kMatchMask << 1)+1; /* kMask2 shifted plus 1*/
    
    if (pos != 0) {
        pos = 4 - pos;
        prefixMatchedBitPattern = 
            ((kMatchMask * ((1 << (pos+1))-1)*2) + (1 << (pos+1))-1) & 
            pattern_items->dna_items->DNAwhichSuffixPosPtr[seq[0]];
        seq++;
        end = (len-pos)/4; 
        remain = (len-pos) % 4;
    } 
    else {
        prefixMatchedBitPattern = kMaskShiftPlus1;
        end = len/4; 
        remain = len % 4;
    }
    for (i = 0; i < end; i++) {
        if ( (tmp = (prefixMatchedBitPattern &
                     pattern_items->dna_items->DNAwhichPrefixPosPtr[seq[i]]))) {
            for (j = 0; j < 4; j++) {
                if (tmp & kMatchMask) {
                    hitArray[twiceNumHits++] = i*4 + j + pos;
                    hitArray[twiceNumHits++] = i*4 + j + pos - 
                        s_LenOf(tmp & kMatchMask, kMatchMask) + 1;
                }
                tmp = (tmp << 1);
            }
        }
        prefixMatchedBitPattern = 
            (((prefixMatchedBitPattern << 4) | kMask2) & 
             pattern_items->dna_items->DNAwhichSuffixPosPtr[seq[i]]);
    }
    /* In the last byte check bits only up to 'remain' */
    if ( (tmp = (prefixMatchedBitPattern &
                 pattern_items->dna_items->DNAwhichPrefixPosPtr[seq[i]]))) {
        for (j = 0; j < remain; j++) {
            if (tmp & kMatchMask) {
                hitArray[twiceNumHits++] = i*4+j + pos;
                hitArray[twiceNumHits++] = i*4+j + pos - 
                    s_LenOf(tmp & kMatchMask, kMatchMask) + 1;
            }
            tmp = (tmp << 1);
        }
    }
    return twiceNumHits;
}

/** Top level routine to find hits when pattern has a short description.
 * @param hitArray Array of hits to pass back [out]
 * @param seq Input sequence [in]
 * @param start Position to start at in seq [in]
 * @param len Length of seq [in]
 * @param is_dna 1 if and only if seq is a DNA sequence [in]
 * @param pattern_blk Carries variables that keep track of search 
 *                      parameters. [in]
 * @return Number of matches found.
 */
static Int4  
s_FindHitsShortHead(Int4* hitArray, const Uint1* seq, Int4 start, Int4 len, 
                Uint1 is_dna, const SPHIPatternSearchBlk *pattern_blk)
{
  if (is_dna) 
    return s_FindHitsShortDNA(hitArray, &seq[start/4], start % 4, len, pattern_blk);
  return _PHIBlastFindHitsShort(hitArray, &seq[start], len, pattern_blk);
}

void 
_PHIPatternWordsLeftShift(Int4 *a, Uint1 b, Int4 numWords)
{
    Int4 x;
    Int4 i; /*index on words*/
    /* Overflow threshold */
    const Int4 kOverflowThreshold = (1 << PHI_BITS_PACKED_PER_WORD);

    for (i = 0; i < numWords; i++) {
        x = (a[i] << 1) + b;
        if (x >= kOverflowThreshold) {
            a[i] = x - kOverflowThreshold; 
            b = 1;
        }
        else { 
            a[i] = x; 
            b = 0;
        }
    }
}  

void 
_PHIPatternWordsBitwiseOr(Int4 *a, Int4 *b, Int4 numWords)
{
    Int4 i; /*index over words*/
    for (i = 0; i < numWords; i++) 
        a[i] = (a[i] | b[i]);
}

Int4
_PHIPatternWordsBitwiseAnd(Int4 *result, Int4 *a, Int4 *b, Int4 numWords) 
{
    Int4 i; /*index over words*/
    Int4 returnValue = 0;
    
    for (i = 0; i < numWords; i++) {
        if ((result[i] = (a[i] & b[i])))
            returnValue = 1;
    }
    return returnValue;
}

/** Returns the difference between the offset F of a first 1-bit in a word
 * sequence and the first offset G < F of a 1-bit in the pattern mask. If 
 * such G does not exist, it is set to -1.
 * @param s Input sequence [in]
 * @param mask Array of word masks [in]
 * @param numWords Number of words in s. [in]
 * @return F - G, see explanation above.
 */
static Int4 
s_LenOfL(Int4 *s, Int4 *mask, Int4 numWords)
{
    Int4 bitIndex; /*loop index over bits in a word*/
    Int4 wordIndex;  /*loop index over words*/
    Int4 firstOneInMask;

    firstOneInMask = -1;
    for (wordIndex = 0; wordIndex < numWords; wordIndex++) {
        for (bitIndex = 0; bitIndex < PHI_BITS_PACKED_PER_WORD; bitIndex++) { 
            if ((s[wordIndex] >> bitIndex) % 2  == 1) 
                return wordIndex*PHI_BITS_PACKED_PER_WORD+bitIndex-firstOneInMask;
            if ((mask[wordIndex] >> bitIndex) %2  == 1) 
                firstOneInMask = wordIndex*PHI_BITS_PACKED_PER_WORD+bitIndex;
        }
    }
    /* This point should never be reached. */
    return -1;
}

/** Finds places where pattern matches seq and returns them as
 * pairs of positions in consecutive entries of hitArray;
 * similar to _PHIBlastFindHitsShort
 * @param hitArray Array of hits to return [out]
 * @param seq Input sequence [in]
 * @param len1 Length of seq [in]
 * @param pattern_blk carries all the pattern variables
 * @return twice the number of hits.
 */
static Int4 
s_FindHitsLong(Int4 *hitArray, const Uint1* seq, Int4 len1, 
               const SPHIPatternSearchBlk *pattern_blk)
{
    Int4 wordIndex; /*index on words in mask*/
    Int4 i; /*loop index on seq */
    Int4  *prefixMatchedBitPattern; /*see similar variable in
                                      _PHIBlastFindHitsShort*/
    Int4 twiceNumHits = 0; /*counter for hitArray*/
    Int4 *mask; /*local copy of match_maskL version of pattern
                  indicates after which positions a match can be declared*/
    Int4 *matchResult; /*Array of words to hold the result of the
                         final test for a match*/
    SLongPatternItems* pattern_items = pattern_blk->multi_word_items;
    Int4 num_words = pattern_items->numWords;

    matchResult = (Int4 *) calloc(num_words, sizeof(Int4));
    mask = (Int4 *) calloc(num_words, sizeof(Int4));
    prefixMatchedBitPattern = (Int4 *) calloc(num_words, sizeof(Int4));
    for (wordIndex = 0; wordIndex < num_words; wordIndex++) {
      mask[wordIndex] = pattern_items->match_maskL[wordIndex];
      prefixMatchedBitPattern[wordIndex] = 0;
    }
    /* This is a multiword version of the algorithm in _PHIBlastFindHitsShort */
    _PHIPatternWordsLeftShift(mask, 1, num_words);
    for (i = 0; i < len1; i++) {
      _PHIPatternWordsLeftShift(prefixMatchedBitPattern, 0, num_words);
      _PHIPatternWordsBitwiseOr(prefixMatchedBitPattern, mask, num_words); 
      _PHIPatternWordsBitwiseAnd(prefixMatchedBitPattern, prefixMatchedBitPattern, 
                                pattern_items->bitPatternByLetter[seq[i]], 
                                num_words);
      if (_PHIPatternWordsBitwiseAnd(matchResult, prefixMatchedBitPattern, 
                                    pattern_items->match_maskL, num_words)) { 
          hitArray[twiceNumHits++] = i; 
          hitArray[twiceNumHits++] = i - 
              s_LenOfL(matchResult, pattern_items->match_maskL, num_words) + 1;
      }
    }
    sfree(prefixMatchedBitPattern); 
    sfree(matchResult); 
    sfree(mask);
    return twiceNumHits;
}

/** Find matches when pattern is very long,
 * @param hitArray Array to pass back pairs of start and end positions for 
 *                 hits [out]
 * @param seq Input sequence [in]
 * @param len Length of seq [in]
 * @param is_dna Is sequence DNA or protein? [in]
 * @param pattern_blk carries all the pattern variables [in]
 * @return Twice the number of hits found.
 */
static Int4 
s_FindHitsVeryLong(Int4 *hitArray, const Uint1* seq, Int4 len, Boolean is_dna,
                   const SPHIPatternSearchBlk *pattern_blk)
{
    Int4 twiceNumHits; /*twice the number of matches*/
    Int4 twiceHitsOneCall; /*twice the number of hits in one call to 
                                 _PHIBlastFindHitsShort */
    Int4 wordIndex;  /*index over words in pattern*/
    Int4 start; /*start position in sequence for calls to _PHIBlastFindHitsShort */
    Int4 hitArray1[PHI_MAX_HIT]; /*used to get hits against different words*/
    Int4 nextPosInHitArray; /*next available position in hitArray1 */
    Int4 hitIndex1, hitIndex2;  /*indices over hitArray1*/
    SLongPatternItems* multiword_items = pattern_blk->multi_word_items;
    SShortPatternItems* word_items = pattern_blk->one_word_items;
    SExtraLongPatternItems* extra_items = multiword_items->extra_long_items;
    Int4 most_specific_word = extra_items->whichMostSpecific;
    
    word_items->whichPositionPtr = multiword_items->SLL[most_specific_word]; 
    word_items->match_mask = multiword_items->match_maskL[most_specific_word];
    if (is_dna) {
      word_items->dna_items->DNAwhichPrefixPosPtr = 
          multiword_items->dna_items->DNAprefixSLL[most_specific_word];
      word_items->dna_items->DNAwhichSuffixPosPtr = 
          multiword_items->dna_items->DNAsuffixSLL[most_specific_word];
    }
    /*find matches to most specific word of pattern*/
    twiceNumHits = 
        s_FindHitsShortHead(hitArray, seq, 0, len, is_dna, pattern_blk);
    if (twiceNumHits < 2) 
      return 0;
    /*extend matches word by word*/
    for (wordIndex = most_specific_word+1; 
         wordIndex < multiword_items->numWords; wordIndex++) {
        word_items->whichPositionPtr = multiword_items->SLL[wordIndex]; 
        word_items->match_mask = multiword_items->match_maskL[wordIndex];
        if (is_dna) {
            word_items->dna_items->DNAwhichPrefixPosPtr = 
                multiword_items->dna_items->DNAprefixSLL[wordIndex]; 
            word_items->dna_items->DNAwhichSuffixPosPtr = 
                multiword_items->dna_items->DNAsuffixSLL[wordIndex];
        }
        nextPosInHitArray = 0;
        for (hitIndex2 = 0; hitIndex2 < twiceNumHits; hitIndex2 += 2) {
            twiceHitsOneCall = 
                s_FindHitsShortHead(&hitArray1[nextPosInHitArray], seq, 
                                    hitArray[hitIndex2]+1, 
                                    MIN(len-hitArray[hitIndex2]-1, 
                                        extra_items->spacing[wordIndex-1] +
                                        extra_items->numPlacesInWord[wordIndex]),
                                    is_dna, pattern_blk);
	  for (hitIndex1 = 0; hitIndex1 < twiceHitsOneCall; hitIndex1+= 2) {
	    hitArray1[nextPosInHitArray+hitIndex1] = 
	      hitArray[hitIndex2]+hitArray1[nextPosInHitArray+hitIndex1]+1;
	    hitArray1[nextPosInHitArray+hitIndex1+1] = hitArray[hitIndex2+1];
	  }
	  nextPosInHitArray += twiceHitsOneCall;
	}
	twiceNumHits = nextPosInHitArray;
	if (twiceNumHits < 2) 
	  return 0;
        /*copy back matches that extend */
	for (hitIndex2 = 0; hitIndex2 < nextPosInHitArray; hitIndex2++) 
	  hitArray[hitIndex2] = hitArray1[hitIndex2];
    }
    /*extend each match back one word at a time*/
    for (wordIndex = most_specific_word-1; wordIndex >=0; 
	 wordIndex--) {
        word_items->whichPositionPtr = multiword_items->SLL[wordIndex]; 
        word_items->match_mask = multiword_items->match_maskL[wordIndex];
        if (is_dna) {
            word_items->dna_items->DNAwhichPrefixPosPtr = 
                multiword_items->dna_items->DNAprefixSLL[wordIndex]; 
            word_items->dna_items->DNAwhichSuffixPosPtr = 
                multiword_items->dna_items->DNAsuffixSLL[wordIndex];
      }
      nextPosInHitArray = 0;
      for (hitIndex2 = 0; hitIndex2 < twiceNumHits; hitIndex2 += 2) {
          start = hitArray[hitIndex2+1] - extra_items->spacing[wordIndex] - 
              extra_items->numPlacesInWord[wordIndex];
          if (start < 0) 
              start = 0;
          twiceHitsOneCall = 
              s_FindHitsShortHead(&hitArray1[nextPosInHitArray], seq, start, 
                                  hitArray[hitIndex2+1]-start, is_dna, pattern_blk);
          for (hitIndex1 = 0; hitIndex1 < twiceHitsOneCall; hitIndex1+= 2) {
              hitArray1[nextPosInHitArray+hitIndex1] = hitArray[hitIndex2];
              hitArray1[nextPosInHitArray+hitIndex1+1] = start + 
                  hitArray1[nextPosInHitArray+hitIndex1+1];
          }
          nextPosInHitArray += twiceHitsOneCall;
      }
      twiceNumHits = nextPosInHitArray;
      if (twiceNumHits < 2) 
          return 0;
      /*copy back matches that extend*/
      for (hitIndex2 = 0; hitIndex2 < nextPosInHitArray; hitIndex2++) 
          hitArray[hitIndex2] = hitArray1[hitIndex2];
    }
    return twiceNumHits;
}

Int4 FindPatternHits(Int4 * hitArray, const Uint1* seq, Int4 len, 
               Boolean is_dna, const SPHIPatternSearchBlk * pattern_blk)
{
    if (pattern_blk->flagPatternLength == eOneWord) 
      return s_FindHitsShortHead(hitArray, seq, 0, len, is_dna, pattern_blk);
    if (pattern_blk->flagPatternLength == eMultiWord) 
      return s_FindHitsLong(hitArray, seq, len, pattern_blk);
    return s_FindHitsVeryLong(hitArray, seq, len, is_dna, pattern_blk);
}

SPHIQueryInfo* SPHIQueryInfoNew()
{
    SPHIQueryInfo* pattern_info;
    const Int4 kMinPhiLookupSize = 8; /* Minimal allocation size for the array
                                          of pattern occurrences. */

    if ((pattern_info = 
         (SPHIQueryInfo*) calloc(1, sizeof(SPHIQueryInfo))) == NULL)
        return NULL;
    pattern_info->allocated_size = kMinPhiLookupSize;
    if ((pattern_info->occurrences = 
         (SPHIPatternInfo*) calloc(kMinPhiLookupSize, sizeof(SPHIPatternInfo)))
        == NULL)
        return NULL;
    return pattern_info;
}

SPHIQueryInfo*
SPHIQueryInfoFree(SPHIQueryInfo* pat_info)
{
    if (pat_info) {
        sfree(pat_info->occurrences);
        sfree(pat_info->pattern);
        sfree(pat_info);
    }
    return NULL;
}

SPHIQueryInfo* 
SPHIQueryInfoCopy(const SPHIQueryInfo* pat_info)
{
    SPHIQueryInfo* retval = NULL;
    
    if (!pat_info)
        return retval;

    retval = 
        (SPHIQueryInfo*) BlastMemDup(pat_info, sizeof(SPHIQueryInfo));
    retval->pattern = 
        (char *) BlastMemDup(pat_info->pattern, 1+strlen(pat_info->pattern));
    retval->occurrences = (SPHIPatternInfo*)
        BlastMemDup(pat_info->occurrences, 
                    pat_info->num_patterns*sizeof(SPHIPatternInfo));
    return retval;
}

/** Adds a new pattern hit to the PHI BLAST pseudo lookup table.
 * @param pattern_info The query pattern information structure. [in] [out]
 * @param offset Offset in query at which pattern was found. [in]
 * @param length Length of the pattern at this offset. [in] 
 */
static Int2 
s_PHIBlastAddPatternHit(SPHIQueryInfo* pattern_info, 
                        Int4 offset, Int4 length)
{
    SPHIPatternInfo* occurrence_array;
    Int4 pat_index = pattern_info->num_patterns;
    
   if (pat_index >= pattern_info->allocated_size) {
       if ((occurrence_array = (SPHIPatternInfo*) 
            realloc(pattern_info->occurrences, 
                    2*pattern_info->allocated_size*sizeof(SPHIPatternInfo)))
           == NULL)
           return -1;
       pattern_info->occurrences = occurrence_array;
       pattern_info->allocated_size *= 2;
   }      
   
   pattern_info->occurrences[pat_index].offset = offset;
   pattern_info->occurrences[pat_index].length = length;
   ++pattern_info->num_patterns;

   return 0;
}

Int4 PHIGetPatternOccurrences(const SPHIPatternSearchBlk * pattern_blk, 
                              const BLAST_SequenceBlk    * query,
                              const BlastSeqLoc          * location, 
                              Boolean                      is_dna,
                              BlastQueryInfo             * query_info)
{
   const BlastSeqLoc* loc;
   Int4* hitArray;
   EBlastProgramType program = (is_dna ? eBlastTypePhiBlastn : eBlastTypePhiBlastp);
   SPHIQueryInfo* pattern_info = query_info->pattern_info;

   ASSERT(pattern_info);
   
   hitArray = (Int4 *) calloc(2*query->length, sizeof(Int4));

   for(loc=location; loc; loc=loc->next) {
      Int4 i, twiceNumHits;
      Int4 from, to;
      Int4 loc_length;
      Uint1* sequence;

      from = loc->ssr->left;
      to = loc->ssr->right;
      loc_length = to - from + 1;
      sequence = query->sequence + from;
      
      twiceNumHits = FindPatternHits(hitArray, sequence, loc_length, is_dna,
                                     pattern_blk);
      
      for (i = 0; i < twiceNumHits; i += 2) {
         if (hitArray[i+1]+from == 0 &&
              hitArray[i]-hitArray[i+1]+1 == BlastQueryInfoGetQueryLength(query_info, program, 0))
         {
            return INT4_MAX;
         }
         s_PHIBlastAddPatternHit(pattern_info, hitArray[i+1]+from, 
                                 hitArray[i]-hitArray[i+1]+1);
      }
   }

   sfree(hitArray);

   return pattern_info->num_patterns;
}

