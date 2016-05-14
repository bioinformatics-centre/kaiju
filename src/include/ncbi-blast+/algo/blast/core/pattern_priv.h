/* $Id: pattern_priv.h 103491 2007-05-04 17:18:18Z kazimird $
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

/** @file pattern_priv.h
 * Auxiliary functions for finding pattern matches in sequence (PHI-BLAST), that
 * are used in multiple source files.
 */

#ifndef ALGO_BLAST_CORE__PATTERN_PRIV_H
#define ALGO_BLAST_CORE__PATTERN_PRIV_H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/pattern.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Routine to find hits of pattern to sequence when sequence is proteins
 * @param hitArray An array of matches to pass back [out]
 * @param seq The input sequence [in]
 * @param len1 Length of the input sequence. [in]
 * @param pattern_blk Carries variables that keep track of search 
 *                      parameters. [in]
 * @return the number of matches found.
 */
Int4 
_PHIBlastFindHitsShort(Int4 *hitArray, const Uint1* seq, Int4 len1, 
                       const SPHIPatternSearchBlk *pattern_blk);

/** Shift each word in the array left by 1 bit and add bit b.
 * If the new values is bigger than an overflow threshold, then subtract the
 * overflow threshold.
 * @param a Array of integers, representing words in a pattern [in] [out]
 * @param b bit to add [in]
 * @param numWords Number of words to process [in]
 */
void 
_PHIPatternWordsLeftShift(Int4 *a, Uint1 b, Int4 numWords);

/** Do a word-by-word bit-wise or of two integer arrays and put the result back
 * in the first array.
 * @param a First array [in] [out]
 * @param b Second array [in]
 * @param numWords Number of words in a and b [in]
 */
void 
_PHIPatternWordsBitwiseOr(Int4 *a, Int4 *b, Int4 numWords);

/** Do a word-by-word bit-wise and of two integer arrays and put the result in
 * a new array.
 * @param result Result of the operation [out]
 * @param a First array [in]
 * @param b Second array [in]
 * @param numWords Size of the two input arrays [in]
 * @return 1 if there are any non-zero words, otherwize 0. 
 */
Int4
_PHIPatternWordsBitwiseAnd(Int4 *result, Int4 *a, Int4 *b, Int4 numWords);

/** Masks all bits corresponding to the aminoacid alphabet, i.e. the first 26
 * bits of an integer number.
 */
extern const int kMaskAaAlphabetBits;

/** Looks for 1 bits in the same position of s and mask
 * Let R be the rightmost position where s and mask both have a 1.
 * Let L < R be the rightmost position where mask has a 1, if any, 
 * or -1 otherwise.
 * @param s Number to check bits in [in]
 * @param mask Mask to apply [in]
 * @param rightOne The rightmost position where s and mask both have a 1 [out]
 * @param rightMaskOnly The rightmost position < rightOne, where mask has a 1,
 *                       if any, or -1 otherwise [out]
 */
void
_PHIGetRightOneBits(Int4 s, Int4 mask, Int4* rightOne, Int4* rightMaskOnly);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__PATTERN_PRIV_H */
