/*  $Id: blast_hits_priv.h 103491 2007-05-04 17:18:18Z kazimird $
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
 * Author:  Christiam Camacho
 *
 */

/** @file blast_hits_priv.h
 *  Utilities for dealing with BLAST HSPs in the core of BLAST.
 */

#ifndef ALGO_BLAST_CORE___BLAST_HITS_PRIV__H
#define ALGO_BLAST_CORE___BLAST_HITS_PRIV__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_hits.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Check the gapped alignments for an overlap of two different alignments.
 * A sufficient overlap is when two alignments have the same start values
 * of have the same final values. 
 * @param hsp_array Pointer to an array of BlastHSP structures [in]
 * @param hsp_count The size of the hsp_array [in]
 * @return The number of valid alignments remaining. 
*/
Int4
Blast_CheckHSPsForCommonEndpoints(BlastHSP* *hsp_array, Int4 hsp_count);

/** Comparison callback function for sorting HSPs, first by score in descending
 * order, then by location. Among alignments with equal score, an HSP will 
 * precede any other HSPs that are completely contained within its endpoints.
 *
 * H2 is contained in H1 if and only if                                         
 * H1.query.offset <= H2.query.offset <= H2.query.end <= H1.query.end 
 * H1.sbjct.offset <= H2.sbjct.offset <= H2.sbjct.end <= H1.sbjct.end
 */
int
ScoreCompareHSPs(const void* h1, const void* h2);

/** TRUE if c is between a and b; f between d and e.  Determines if the
 * coordinates are already in an HSP that has been evaluated. 
*/
#define CONTAINED_IN_HSP(a,b,c,d,e,f) \
    (((a <= c && b >= c) && (d <= f && e >= f)) ? TRUE : FALSE)

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_HITS_PRIV__H */
