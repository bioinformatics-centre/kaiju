/* $Id: blast_rps.h 374315 2012-09-10 13:46:47Z fongah2 $
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
 * Author:  Jason Papadopoulos
 *
 */

/** @file blast_rps.h
 * RPS BLAST structure definitions.
 */

#ifndef ALGO_BLAST_CORE__BLAST_RPS__H
#define ALGO_BLAST_CORE__BLAST_RPS__H

#include <algo/blast/core/ncbi_std.h>

#ifdef __cplusplus
extern "C" {
#endif

#define RPS_MAGIC_NUM 0x1e16    /**< RPS blast version number */
#define RPS_MAGIC_NUM_28 0x1e17    /**< Version number for 28-letter alphabet */
#define NUM_EXPANSION_WORDS 3   /**< Intentionally unused words in .loo file */

/** header of RPS blast '.loo' file */

typedef struct BlastRPSLookupFileHeader {
    Int4 magic_number;               /**< value should be RPS_MAGIC_NUM */
    Int4 num_lookup_tables;          /**< hardwired to 1 at present */
    Int4 num_hits;                   /**< number of hits in the lookup table */
    Int4 num_filled_backbone_cells;  /**< backbone cells that contain hits */
    Int4 overflow_hits;              /**< number of hits in overflow array */
    Int4 unused[NUM_EXPANSION_WORDS];/**< empty space in the on-disk format */
    Int4 start_of_backbone;          /**< byte offset of start of backbone */
    Int4 end_of_overflow;            /**< byte offset to end of overflow array */
} BlastRPSLookupFileHeader;

/** header of RPS blast '.rps' file */

typedef struct BlastRPSProfileHeader {
    Int4 magic_number;     /**< value should be RPS_MAGIC_NUM */
    Int4 num_profiles;     /**< number of PSSMs in the file */
    Int4 start_offsets[1]; /**< start of an Int4 array that gives the starting 
                              byte offset of each RPS DB sequence. There
                              are num_profiles+1 entries in the list, and
                              the last entry effectively contains the length
                              of all protein sequences combined. Note that
                              the length of each sequence includes one byte
                              at the end for an end-of-sequence sentinel */

    /* After the list of sequence start offsets comes the list
       of PSSM rows. There is one row for each letter in the RPS
       sequence database, and each row has BLASTAA_SIZE entries.
       Because there is a sentinel byte at the end of each sequence,
       there is also a PSSM row for each sentinel byte */

} BlastRPSProfileHeader;

/** header for RPS blast frequency ratios ('.freq') file */

#define FREQ_RATIO_SCALE 1000000000

typedef struct BlastRPSFreqRatiosHeader {
	Int4 magic_number;
    Int4 num_profiles;     /**< number of PSSMs in the file */
    Int4 start_offsets[1];   /**< start of an Int4 array that gives the starting
                              byte offset of each RPS DB sequence. There
                              are num_profiles+1 entries in the list, and
                              the last entry effectively contains the length
                              of all protein sequences combined. Note that
                              the length of each sequence includes one byte
                              at the end for an end-of-sequence sentinel */

} BlastRPSFreqRatiosHeader;

/** information derived from RPS blast '.aux' file */

typedef struct BlastRPSAuxInfo {
    char* orig_score_matrix; /**< score matrix used to derive PSSMs */
    Int4 gap_open_penalty;   /**< gap open penalty used in deriving PSSMs */
    Int4 gap_extend_penalty; /**< gap extend penalty used in deriving PSSMs */
    double ungapped_k;       /**< ungapped Karlin value for orig_score_matrix
                                  (not used) */
    double ungapped_h;       /**< ungapped Karlin value for orig_score_matrix
                                  (not used) */
    Int4 max_db_seq_length;  /**< maximum DB sequence length (not used) */
    Int4 db_length;          /**< RPS DB search space (not used) */
    double scale_factor;     /**< the PSSMs are scaled by this amount, and so
                                  all scores and all cutoff values must be
                                  similarly scaled during the search */
    double *karlin_k;        /**< one Karlin value for each DB sequence */
} BlastRPSAuxInfo;

/** The RPS engine uses this structure to access all of the
 *  RPS blast related data (assumed to be collected in an 
 *  implementation-specific manner). 
 */
typedef struct BlastRPSInfo {
    BlastRPSLookupFileHeader *lookup_header; /**< for '.loo' file */
    BlastRPSProfileHeader *profile_header;   /**< for '.rps' file */
    BlastRPSAuxInfo aux_info;                /**< for '.aux' file */

    BlastRPSProfileHeader *freq_header;      /**< for '.wcounts' file */
    BlastRPSProfileHeader *obsr_header;      /**< for '.obsr' file */
    BlastRPSFreqRatiosHeader *freq_ratios_header; /**< for '.freq' file */
} BlastRPSInfo;

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLAST_RPS__H */
