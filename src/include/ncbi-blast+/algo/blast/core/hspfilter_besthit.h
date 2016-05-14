/*  $Id: hspfilter_besthit.h 274418 2011-04-14 13:40:26Z maning $
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

/** @file hspfilter_besthit.h
 * Implementation of a number of BlastHSPWriters to save hits from
 * a BLAST search, and subsequently return them in sorted order.
 */

#ifndef ALGO_BLAST_CORE__HSPFILTER_BESTHIT__H
#define ALGO_BLAST_CORE__HSPFILTER_BESTHIT__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_hspfilter.h>
#include <algo/blast/core/blast_hits.h>
#include <connect/ncbi_core.h>

#ifdef __cplusplus
extern "C" {
#endif

/************************************************************************/
/** The "best hit" writer

   Prune the hsp_list for each query and keeps only the best ones.
   1. For a pair of hits A and B, check based on 10% overhangs whether
       A can be dropped because of B due to end points of A being within
       10% extension of B and vice versa. Note that this would allow A
       to be dropped even if it is at most 20% longer than B.

   2. If A can be dropped because of B, check if Evalue(A) >= Evalue(B);
       that is A has the same or worse evalue than B. Do the same check for
       whether B can be dropped because of A.

   3. If A can still be dropped because of B, check if density(A) <= density(B).
       Do the same check for whether B can be dropped because of A.

   4. If only one can be dropped, then drop that one. If both are mutually
       replaceable, use length criteria and drop the shorter one
       only if it is at least 10% shorter (90% coverage).

   So, essentially length coverage is being used a tie-breaker and if the 
   tie-breaker does not break the tie, both alignments are kept. Above is not 
   very different than what you have now, just rearranged in conditions 
   so that we do not have non-deterministic behavior between a pair of 
   alignments. We could still have issues with cascades where A was dropped 
   because of B and then B gets dropped because of C, but A would not have 
   been dropped because of C becuase of condition 4. However, I think this 
   will be extremely rare.
  */

/// Default value for overhang
#define kBestHit_OverhangDflt 0.1
/// Minimum value for overhang
#define kBestHit_OverhangMin 0.0
/// Maximum value for overhang
#define kBestHit_OverhangMax 0.5

/// Default value for score_edge
#define kBestHit_ScoreEdgeDflt 0.1
/// Minimum value for score_edge
#define kBestHit_ScoreEdgeMin  0.0
/// Maximum value for score_edge
#define kBestHit_ScoreEdgeMax  0.5

/** Keeps parameters used in best hit algorithm.*/
typedef struct BlastHSPBestHitParams {
   EBlastProgramType program;/**< program type. */
   Int4 prelim_hitlist_size; /**< number of hits saved during preliminary
                                  part of search. */
   Int4 hsp_num_max;         /**< number of HSPs to save per db sequence. */
   double overhang;          /**< overhang used in condition 1. */
   double score_edge;        /**< fraction of score margin in condition 4*/
} BlastHSPBestHitParams;

/** create a set of parameters 
 * @param program Blast program type.[in]
 * @param hit_options field hitlist_size and hsp_num_max needed, a pointer to 
 *      this structure will be stored on resulting structure.[in]
 * @param best_hit_opts Specifies the ratio of overhang to length, which is used to
        determine if hit A is contained in hit B
 * @param compostionBasedStats the compsotion based stats needed. [in]
 * @param gapped_calculation if gapped_calculation is needed. [in]
 * @return the pointer to the allocated parameter
 */
NCBI_XBLAST_EXPORT
BlastHSPBestHitParams*
BlastHSPBestHitParamsNew(const BlastHitSavingOptions* hit_options,
                         const BlastHSPBestHitOptions* best_hit_opts,
                         Int4 compositionBasedStats,
                         Boolean gapped_calculation);

/** Deallocates the BlastHSPBestHitParams structure passed in
 * @param opts structure to deallocate [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
BlastHSPBestHitParams*
BlastHSPBestHitParamsFree(BlastHSPBestHitParams* opts);

/** WriterInfo and PipeInfo to create a best hit writer/pipe
 * @param params Specifies writer parameters. [in]
 * @return the newly allocated writer/pipe info
 */
NCBI_XBLAST_EXPORT
BlastHSPWriterInfo* 
BlastHSPBestHitInfoNew(BlastHSPBestHitParams* params);

NCBI_XBLAST_EXPORT
BlastHSPPipeInfo*
BlastHSPBestHitPipeInfoNew(BlastHSPBestHitParams* params);
                 
#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__HSPFILTER_BESTHIT__H */
