/*  $Id: hspfilter_collector.h 161402 2009-05-27 17:35:47Z camacho $
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
 * Author:  Ning Ma
 *
 */

/** @file hspfilter_collector.h
 * Implementation of a number of BlastHSPWriters to save hits from
 * a BLAST search, and subsequently return them in sorted order.
 */

#ifndef ALGO_BLAST_CORE__HSPFILTER_COLLECTOR__H
#define ALGO_BLAST_CORE__HSPFILTER_COLLECTOR__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hspfilter.h>
#include <algo/blast/core/blast_hits.h>
#include <connect/ncbi_core.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Keeps prelim_hitlist_size and HitSavingOptions together. */
typedef struct BlastHSPCollectorParams {
   EBlastProgramType program;/**< program type */
   Int4 prelim_hitlist_size; /**< number of hits saved during preliminary
                                  part of search. */
   Int4 hsp_num_max;         /**< number of HSPs to save per db sequence.*/
} BlastHSPCollectorParams;

/** Sets up parameter set for use by collector.
 * @param program Blast program type.[in]
 * @param hit_options field hitlist_size and hsp_num_max needed, a pointer to 
 *      this structure will be stored on resulting structure.[in]
 * @param ext_options field compositionBasedStats needed here. [in]
 * @param scoring_options gapped_calculation needed here. [in]
 * @return the pointer to the allocated parameter
 */
NCBI_XBLAST_EXPORT
BlastHSPCollectorParams*
BlastHSPCollectorParamsNew(const BlastHitSavingOptions* hit_options,
                           Int4 compositionBasedStats,
                           Boolean gapped_calculation);

/** Deallocates the BlastHSPCollectorParams structure passed in
 * @param opts structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
BlastHSPCollectorParams*
BlastHSPCollectorParamsFree(BlastHSPCollectorParams* opts);

/** WriterInfo to create a default writer: the collecter
 * @param params The collector parameters.
 * @return pointer to WriterInfo
 */
NCBI_XBLAST_EXPORT
BlastHSPWriterInfo* 
BlastHSPCollectorInfoNew(BlastHSPCollectorParams* params);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__HSPFILTER_COLLECTOR__H */
