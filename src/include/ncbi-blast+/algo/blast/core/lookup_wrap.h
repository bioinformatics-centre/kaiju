/* $Id: lookup_wrap.h 369355 2012-07-18 17:07:15Z morgulis $
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

/** @file lookup_wrap.h
 * Wrapper for all lookup tables used in BLAST
 */

#ifndef ALGO_BLAST_CORE__LOOKUP_WRAP__H
#define ALGO_BLAST_CORE__LOOKUP_WRAP__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_rps.h>
#include <algo/blast/core/blast_stat.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Wrapper structure for different types of BLAST lookup tables */
typedef struct LookupTableWrap {
   ELookupTableType lut_type; /**< What kind of a lookup table it is? */
   void* lut; /**< Pointer to the actual lookup table structure */
   void* read_indexed_db; /**< function used to retrieve hits
                              from an indexed database */
   void* check_index_oid; /**< function used to check if seeds
                               for a given oid are present */
   void * end_search_indication; /**< function used to report that
                                      a thread is done iterating over
                                      the database in preliminary
                                      search */
   void* lookup_callback;    /**< function used to look up an
                                  index->q_off pair */
} LookupTableWrap;

/** Function pointer type to check the presence of index->q_off pair */
typedef Boolean (*T_Lookup_Callback)(const LookupTableWrap *, Int4, Int4);

/** Create the lookup table for all query words.
 * @param query The query sequence [in]
 * @param lookup_options What kind of lookup table to build? [in]
 * @param query_options options for query setup [in]
 * @param lookup_segments Locations on query to be used for lookup table
 *                        construction [in]
 * @param sbp Scoring block containing matrix [in]
 * @param lookup_wrap_ptr The initialized lookup table [out]
 * @param rps_info Structure containing RPS blast setup information [in]
 * @param error_msg message with warning or errors [in|out]
 */
NCBI_XBLAST_EXPORT
Int2 LookupTableWrapInit(BLAST_SequenceBlk* query, 
        const LookupTableOptions* lookup_options,	
        const QuerySetUpOptions* query_options,
        BlastSeqLoc* lookup_segments, BlastScoreBlk* sbp, 
        LookupTableWrap** lookup_wrap_ptr, const BlastRPSInfo *rps_info,
        Blast_Message* *error_msg);

/** Deallocate memory for the lookup table */
NCBI_XBLAST_EXPORT
LookupTableWrap* LookupTableWrapFree(LookupTableWrap* lookup);

/** Default size of offset arrays filled in a single ScanSubject call. */
#define OFFSET_ARRAY_SIZE 4096

/** Determine the size of the offsets arrays to be filled by
 * the ScanSubject function.
 */
NCBI_XBLAST_EXPORT
Int4 GetOffsetArraySize(LookupTableWrap* lookup);

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__LOOKUP_WRAP__H */
