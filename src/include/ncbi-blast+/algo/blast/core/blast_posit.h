/*  $Id: blast_posit.h 103491 2007-05-04 17:18:18Z kazimird $
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
 * Author:  Alejandro Schaffer (ported by Christiam Camacho)
 *
 */

/** @file blast_posit.h
 * Port of posit.h structures and impalaScaling for implementing composition
 * based statistics for PSI-BLAST.
 */

#ifndef ALGO_BLAST_CORE___BLAST_POSIT__H
#define ALGO_BLAST_CORE___BLAST_POSIT__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_stat.h>
#include "matrix_freq_ratios.h"

#ifdef __cplusplus
extern "C" {
#endif

/** number of real aminoacids (i.e.: does not include U, X, B, etc) */
#define PRO_TRUE_ALPHABET_SIZE 20
/** range of scores in a matrix */
#define kScoreMatrixScoreRange 10000

/** positions of true characters in protein alphabet*/
extern const Int4 trueCharPositions[PRO_TRUE_ALPHABET_SIZE];

/** Structure used to pass data into the scaling routines. All fields marked as
 * alias are not owned by this structure. */
typedef struct Kappa_posSearchItems {
    /** PSSM */
    int**               posMatrix;
    /** Scaled PSSM [alias] */
    int**               posPrivateMatrix;       
    /** PSSM's frequency ratios [alias] */
    double**            posFreqs;               
    /** Frequecy ratios for underlying scoring matrix */
    SFreqRatios*        stdFreqRatios;
    /** Length of the query sequence, specifies the number of columns in the
     * matrices in this structure */
    unsigned int        queryLength;
} Kappa_posSearchItems;

/** Structure used to pass data into the scaling routines. All fields marked as
 * aliases refer to fields in the BlastScoreBlk structure and are NOT owned by
 * this structure */
typedef struct Kappa_compactSearchItems {
    /** Query sequence data in ncbistdaa format without sentinel bytes [alias]*/
    Uint1*              query;                 
    /** Length of the sequence above */
    int                 qlength;
    /** Size of the alphabet @sa BLASTAA_SIZE */
    int                 alphabetSize;
    /** Standard substitution scoring matrix [alias] */
    int**               matrix;
    /** Ungapped Karlin-Altschul parameters [alias] */
    Blast_KarlinBlk**   kbp_std;
    /** Ungapped PSI-BLAST Karlin-Altschul parameters [alias] */
    Blast_KarlinBlk**   kbp_psi;
    /** Gapped Karlin-Altschul parameters [alias] */
    Blast_KarlinBlk**   kbp_gap_std;
    /** Gapped PSI-BLAST Karlin-Altschul parameters [alias] */
    Blast_KarlinBlk**   kbp_gap_psi;
    /** Lambda calculated using standard residue compositions for the query and
     * database sequences */
    double              lambda_ideal;
    /** K calculated using standard residue compositions for the query and
     * database sequences */
    double              K_ideal;
    /** Array of standard residue probabilities, as those returned by
     * BLAST_GetStandardAaProbabilities */
    double*             standardProb;

} Kappa_compactSearchItems;

/** Allocates a new Kappa_posSearchItems structure
 * @param queryLength length of the query sequence [in]
 * @param matrix_name name of the underlying matrix name to use [in]
 * @param posPrivateMatrix scaled pssm, allocated with dimensions queryLength
 * by BLASTAA_SIZE. This is owned by the caller [in|out]
 * @param posFreqs PSSM's frequency ratios, allocated with dimensions
 * queryLength by BLASTAA_SIZE. This is owned by the caller [in|out]
 * @return newly allocated structure or NULL if out of memory
 */
Kappa_posSearchItems*
Kappa_posSearchItemsNew(unsigned int queryLength,
                        const char* matrix_name,
                        int** posPrivateMatrix,
                        double** posFreqs);

/** Deallocates the Kappa_posSearchItems structure.
 * @param posSearchItems data structure to deallocate [in] 
 * @return NULL
 */
Kappa_posSearchItems*
Kappa_posSearchItemsFree(Kappa_posSearchItems* posSearchItems);

/** Creates a new Kappa_compactSearchItems structure
 * @param query query sequence data in ncbistdaa format without sentinel 
 * bytes [in]
 * @param queryLength length of the sequence above [in]
 * @param sbp BLAST scoring block structure [in]
 * @return newly allocated structure or NULL if out of memory
 */
Kappa_compactSearchItems*
Kappa_compactSearchItemsNew(const Uint1* query, unsigned int queryLength, 
                            BlastScoreBlk* sbp);

/** Deallocates the Kappa_compactSearchItems structure.
 * @param compactSearchItems data structure to deallocate [in] 
 * @return NULL
 */
Kappa_compactSearchItems*
Kappa_compactSearchItemsFree(Kappa_compactSearchItems* compactSearchItems);

/** Copied from posit2.c 
 * @return 0 on success, 1 on failure
 */
int Kappa_impalaScaling(Kappa_posSearchItems* posSearch,
                        Kappa_compactSearchItems* compactSearch,
                        double scalingFactor,
                        Boolean doBinarySearch,
                        BlastScoreBlk* sbp);

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_POSIT__H */
