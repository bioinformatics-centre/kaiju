/* $Id: pattern.h 124526 2008-04-15 15:27:44Z camacho $
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

/** @file pattern.h
 * Functions for finding pattern matches in sequence (PHI-BLAST).
 * @todo FIXME The structures used for finding pattern contain a number of
 *             arrays with static fixed sizes. It remains to be determined
 *             whether some of these arrays may be dynamically allocated 
 *             instead.
 */

#ifndef ALGO_BLAST_CORE__PATTERN_H
#define ALGO_BLAST_CORE__PATTERN_H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_export.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_query_info.h>

#ifdef __cplusplus
extern "C" {
#endif

#define PHI_BUF_SIZE 100       /**< Default size for buffers */
#define PHI_ASCII_SIZE 256     /**< Size of ASCII alphabet */ 
#define PHI_BITS_PACKED_PER_WORD 30 /**< Number of bits packed in a word */
#define PHI_MAX_WORD_SIZE   11               /**< Maximal word size */ 
/** Threshold pattern length. */
#define PHI_MAX_PATTERN_LENGTH   (PHI_BITS_PACKED_PER_WORD * PHI_MAX_WORD_SIZE) 
#define PHI_MAX_WORDS_IN_PATTERN 100  /**< Maximal number of words in pattern */
#define PHI_MAX_HIT 20000 /**< Maximal size of an array of pattern hits. 
                             @todo FIXME This limit can be avoided by using
                             dynamically allocated arrays! */

/** Options for running the pattern search. */
typedef enum EPatternProgram {
    eSeed = 1,/**< Use only those query occurrences that are specified in the input
                 pattern file. @todo Not implemented. */
    ePattern, /**< Only find pattern occurrences in database, but do not perform
                 alignments. @todo Not implemented. */
    ePatSeed, /**< Search a BLAST database using pattern occurrences as seeds */
    ePatMatch /**< Only find pattern occurrences in query, but do not search the
                 database. @todo Not implemented. */
} EPatternProgram;

/** Type of pattern: fits in single word, several words, or is very long. */
typedef enum EPatternType {
    eOneWord = 0, /**< Does pattern consist of a single word? */
    eMultiWord,   /**< Does pattern consist of a multiple words? */
    eVeryLong     /**< Is pattern too long for a simple multi-word processing? */
} EPatternType;

/** Structure containing auxiliary items needed for a DNA search with a pattern
 * that fits in a single word. 
 */
typedef struct SDNAShortPatternItems {
    Uint4 *DNAwhichPrefixPosPtr; /**< Prefix position array for DNA patterns */
    Uint4 *DNAwhichSuffixPosPtr; /**< Suffix position array for DNA patterns*/
    /** Where prefix of DNA 4-mer matches pattern */
    Uint4 DNAwhichPrefixPositions[PHI_ASCII_SIZE]; 
    /** Where suffix of DNA 4-mer matches pattern */
    Uint4 DNAwhichSuffixPositions[PHI_ASCII_SIZE]; 
} SDNAShortPatternItems;

/** Auxiliary items needed for a PHI BLAST search with a pattern that fits in a
 * single word.
 */
typedef struct SShortPatternItems {
    Int4 match_mask;/**< Bit mask representation of input pattern
                       for patterns that fit in a word*/
    Int4 *whichPositionPtr; /**< Array of positions where pattern lettern should
                               match, for a single word of the pattern. */
    SDNAShortPatternItems* dna_items; /**< Additional items for a DNA search. */
} SShortPatternItems;

/** Auxiliary items needed for a PHI BLAST search with pattern that contains
 * pieces longer than a word. 
 */
typedef struct SExtraLongPatternItems {
    /** When pattern has more than 7 words, keep track of how many places of 
       pattern in each word of the representation; */
    Int4   numPlacesInWord[PHI_MAX_WORDS_IN_PATTERN]; 
    /** Spaces until next word due to wildcard*/
    Int4   spacing[PHI_MAX_WORDS_IN_PATTERN]; 
    Int4   highestPlace; /**< Number of places in pattern representation
                            as computed in input_pattern */
    Int4   whichMostSpecific; /**< Which word in an extra long pattern
                                 has the lowest probability of a match*/
} SExtraLongPatternItems;

/** Auxiliary items needed for a DNA pattern search with pattern containing 
 * multiple words. 
 */
typedef struct SDNALongPatternItems {
    /** Where prefix of DNA 4-mer matches pattern, for multiple-word patterns */
    Uint4   DNAprefixSLL[PHI_MAX_WORDS_IN_PATTERN][PHI_ASCII_SIZE];
    /** Where suffix of DNA 4-mer matches pattern, for multiple-word patterns */
    Uint4   DNAsuffixSLL[PHI_MAX_WORDS_IN_PATTERN][PHI_ASCII_SIZE];
} SDNALongPatternItems;


/** Auxiliary items needed for a PHI BLAST search with pattern containing 
 * multiple words. 
 */
typedef struct SLongPatternItems {
    Int4 numWords;  /**< Number of words need to hold bit representation
                       of pattern*/
    Int4 match_maskL[PHI_BUF_SIZE]; /**< Bit mask representation of input pattern
                                   for long patterns*/
    /** Which positions can a character occur in for long patterns*/
    Int4 bitPatternByLetter[PHI_ASCII_SIZE][PHI_MAX_WORD_SIZE]; 
    /** For each letter in the alphabet and each word in the masked
        pattern representation, holds a bit pattern saying for which
        positions the letter will match. Similar to whichPositionsByCharacter 
        for many-word patterns. */
    Int4   SLL[PHI_MAX_WORDS_IN_PATTERN][PHI_ASCII_SIZE]; 
    Int4   inputPatternMasked[PHI_MAX_PATTERN_LENGTH]; /**< Masked input pattern */

    SDNALongPatternItems* dna_items; /**< Additional items necessary for a DNA
                                        pattern. */
    SExtraLongPatternItems* extra_long_items; /**< Additional items necessary if
                                                 pattern contains pieces longer
                                                 than a word. */
} SLongPatternItems;

/** Structure containing all auxiliary information needed in a pattern 
 * search.
 */
typedef struct SPHIPatternSearchBlk {
    /** Indicates if the whole pattern fits in 1 word, each of several parts of
        the pattern fit in a word, or some parts of the pattern are too long to
        fit in a word. */
    EPatternType flagPatternLength; 
    double  patternProbability; /**< Probability of this letter combination */
    Int4   minPatternMatchLength; /**< Minimum length of string to match this 
                                     pattern*/
    SShortPatternItems* one_word_items; /**< Items necessary when pattern fits
                                           in one word. */
                                           
    SLongPatternItems* multi_word_items; /**< Additional items, when pattern
                                            requires multiple words. */
    Int4 num_patterns_db; /**< Number of patterns actually found during the 
                                            database search. */
    char* pattern; /**< Pattern used, saved here for error reporting. */
} SPHIPatternSearchBlk;

    
/** Find the places where the pattern matches seq;
 * 3 different methods are used depending on the length of the pattern.
 * @param hitArray Stores the results as pairs of positions in consecutive
 *                 entries [out]
 * @param seq Sequence [in]
 * @param len Length of the sequence [in]
 * @param is_dna Indicates whether seq is made of DNA or protein letters [in]
 * @param patternSearch Pattern information [in]
 * @return Twice the number of hits (length of hitArray filled in)
*/
NCBI_XBLAST_EXPORT
Int4 FindPatternHits(Int4 * hitArray, const Uint1 * seq, Int4 len,
                     Boolean is_dna,
                     const SPHIPatternSearchBlk * patternSearch);

/** Allocates the pattern occurrences structure. */
NCBI_XBLAST_EXPORT
SPHIQueryInfo* SPHIQueryInfoNew(void);

/** Frees the pattern information structure. 
 * @param pat_info Structure to free. [in]
 * @return NULL.
 */
NCBI_XBLAST_EXPORT
SPHIQueryInfo*
SPHIQueryInfoFree(SPHIQueryInfo* pat_info);

/** Copies the SPHIQueryInfo structure.
 * @param pat_info Structure to copy [in]
 * @return New structure.
 */
NCBI_XBLAST_EXPORT
SPHIQueryInfo* 
SPHIQueryInfoCopy(const SPHIQueryInfo* pat_info);

/** Finds all pattern hits in a given query and saves them in the 
 * previously allocated SPHIQueryInfo structure.
 * @param pattern_blk Structure containing pattern structure. [in]
 * @param query Query sequence(s) [in]
 * @param location Segments in the query sequence where to look for 
 *                 pattern [in]
 * @param is_dna Is this a nucleotide sequence? [in]
 * @param query_info Used to store pattern occurrences and get length 
 *                     of query (for error checking) [out]
 * @return a negative number is an unknown error, INT4_MAX indicates the
 *       pattern (illegally) covered the entire query, other non-negative numbers
 *       indicate the nubmer of pattern occurrences found.
 */
NCBI_XBLAST_EXPORT
Int4 PHIGetPatternOccurrences(const SPHIPatternSearchBlk * pattern_blk,
                              const BLAST_SequenceBlk    * query,
                              const BlastSeqLoc          * location, 
                              Boolean                      is_dna,
                              BlastQueryInfo*             query_info);


#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__PATTERN_H */
