#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_posit.c 99676 2007-03-05 20:41:55Z kazimird $";
#endif /* SKIP_DOXYGEN_PROCESSING */
/* ===========================================================================
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

/** @file blast_posit.c
 *  Port of posit.c/posit2.c to implement composition based statistics for
 *  PSI-BLAST
 */
    
#include <algo/blast/core/ncbi_math.h>
#include <algo/blast/core/blast_util.h>
#include "blast_posit.h"
#include "blast_psi_priv.h"

const Int4 trueCharPositions[PRO_TRUE_ALPHABET_SIZE] =
  {1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22};


Kappa_posSearchItems*
Kappa_posSearchItemsNew(unsigned int queryLength,
                  const char* matrix_name,
                  int** posPrivateMatrix,
                  double** posFreqs)
{
    Kappa_posSearchItems* retval = NULL;

    retval = (Kappa_posSearchItems*) calloc(1, sizeof(Kappa_posSearchItems));
    if ( !retval ) {
        return NULL;
    }

    retval->posMatrix = (int**) _PSIAllocateMatrix(queryLength, BLASTAA_SIZE,
                                                   sizeof(int));
    if ( !retval->posMatrix ) {
        return Kappa_posSearchItemsFree(retval);
    }

    retval->stdFreqRatios = _PSIMatrixFrequencyRatiosNew(matrix_name);
    if ( !retval->stdFreqRatios ) {
        return Kappa_posSearchItemsFree(retval);
    }

    retval->queryLength = queryLength;
    retval->posPrivateMatrix = posPrivateMatrix;
    retval->posFreqs = posFreqs;

    return retval;
}

Kappa_posSearchItems*
Kappa_posSearchItemsFree(Kappa_posSearchItems* posSearch)
{
    if ( !posSearch ) {
        return NULL;
    }

    if (posSearch->posMatrix) {
        posSearch->posMatrix = 
            (int**) _PSIDeallocateMatrix((void**) posSearch->posMatrix,
                                         posSearch->queryLength);
    }

    if (posSearch->stdFreqRatios) {
        posSearch->stdFreqRatios =
            _PSIMatrixFrequencyRatiosFree(posSearch->stdFreqRatios);
    }

    posSearch->posPrivateMatrix = NULL;
    posSearch->posFreqs = NULL;
    
    sfree(posSearch);

    return NULL;
}

Kappa_compactSearchItems*
Kappa_compactSearchItemsNew(const Uint1* query, unsigned int queryLength, 
                            BlastScoreBlk* sbp)
{
    Kappa_compactSearchItems* retval = NULL;

    ASSERT(sbp);
    ASSERT(query);

    retval = (Kappa_compactSearchItems*) 
        calloc(1, sizeof(Kappa_compactSearchItems));
    if ( !retval ) {
        return NULL;
    }

    retval->standardProb = BLAST_GetStandardAaProbabilities();
    if ( !retval->standardProb ) {
        return Kappa_compactSearchItemsFree(retval);
    }

    ASSERT(sbp->alphabet_code == BLASTAA_SEQ_CODE);
    ASSERT(sbp->protein_alphabet == TRUE);
    ASSERT(sbp->alphabet_size == BLASTAA_SIZE);
    ASSERT(sbp->matrix);

    retval->query = (Uint1*) query;
    retval->qlength = queryLength;
    retval->alphabetSize = BLASTAA_SIZE;
    retval->matrix = sbp->matrix->data;
    retval->kbp_std = sbp->kbp_std;
    retval->kbp_psi = sbp->kbp_psi;
    retval->kbp_gap_std = sbp->kbp_gap_std;
    retval->kbp_gap_psi = sbp->kbp_gap_psi;
    retval->lambda_ideal = sbp->kbp_ideal->Lambda;
    retval->K_ideal = sbp->kbp_ideal->K;

    return retval;
}

Kappa_compactSearchItems*
Kappa_compactSearchItemsFree(Kappa_compactSearchItems* compactSearch)
{
    if ( !compactSearch ) {
        return NULL;
    } 

    if (compactSearch->standardProb) {
        sfree(compactSearch->standardProb);
    }

    compactSearch->query = NULL;
    compactSearch->qlength = 0;
    compactSearch->alphabetSize = 0;
    compactSearch->matrix = NULL;
    compactSearch->kbp_std = NULL;
    compactSearch->kbp_psi = NULL;
    compactSearch->kbp_gap_std = NULL;
    compactSearch->kbp_gap_psi = NULL;
    compactSearch->lambda_ideal = 0.0;
    compactSearch->K_ideal = 0.0;

    sfree(compactSearch);

    return NULL;
}

/** Replace/refactor with either one of the already existing variants of this
 * function: posfillSfp (blast_kappa.c), _PSIComputeScoreProbabilities
 * (blast_psi_priv.c) or BlastScoreFreqCalc (blast_stat.c) or RPSFillScores
 * (blast_stat.c) */
static Blast_ScoreFreq* fillSfp(int ** matrix, int matrixLength,
                                double * queryProbArray,
                                double * scoreArray,
                                Blast_ScoreFreq* return_sfp)
{
    int minScore, maxScore;    /* observed minimum and maximum scores */
    int i, j, k;               /* indices */
    double onePosFrac;     /* 1/matrix length as a double */

    minScore = BLAST_SCORE_MAX;
    maxScore = BLAST_SCORE_MIN;

    for (i = 0; i < matrixLength; i++) {
        for (j = 0; j < PRO_TRUE_ALPHABET_SIZE; j++) {
            k = trueCharPositions[j];
            if ((matrix[i][k] != BLAST_SCORE_MIN)
                && (matrix[i][k] < minScore))
                minScore = matrix[i][k];
            if (matrix[i][k] > maxScore)
                maxScore = matrix[i][k];
        }
    }
    ASSERT(minScore != BLAST_SCORE_MAX);
    ASSERT(maxScore != BLAST_SCORE_MIN);

    return_sfp->obs_min = minScore;
    return_sfp->obs_max = maxScore;
    if ((maxScore - minScore) >= kScoreMatrixScoreRange) {
        return NULL;
    }
    for (i = 0; i < kScoreMatrixScoreRange; i++)
        scoreArray[i] = 0.0;
    return_sfp->sprob = &(scoreArray[-minScore]);       /* center around 0 */
    onePosFrac = 1.0 / ((double) matrixLength);
    for (i = 0; i < matrixLength; i++) {
        for (j = 0; j < PRO_TRUE_ALPHABET_SIZE; j++) {
            k = trueCharPositions[j];
            if (matrix[i][k] >= minScore) {
                return_sfp->sprob[matrix[i][k]] +=
                    (onePosFrac * queryProbArray[k]);
            }
        }
    }
    return_sfp->score_avg = 0;
    for (i = minScore; i <= maxScore; i++)
        return_sfp->score_avg += i * return_sfp->sprob[i];
    return (return_sfp);
}

/** Copy of posit2.c's impalaScaleMatrix
 * @todo refactor this function along with all scaling code!
 * @param compactSearch [in]
 * @param posMatrix PSSM [in|out]
 * @param posPrivateMatrix scaled PSSM [in|out]
 * @param scalingFactor impala scaling factor [in]
 * @param doBinarySearch perform binary search? [in]
 * @param sbp BLAST scoring block structure, Karlin-Altschul parameters are
 * updated [in|out]
 */
static Boolean
impalaScaleMatrix(Kappa_compactSearchItems* compactSearch, 
                  int** posMatrix,
                  int** posPrivateMatrix,
                  double scalingFactor, 
                  Boolean doBinarySearch,
                  BlastScoreBlk* sbp)
{
    int dim1, dim2;            /* number of rows and number of columns */
    int a, c;                  /* loop indices */
    Boolean too_high = TRUE, done, first_time;  /* control variables for
                                                   binary search */
    /* multiplicative factors in binary search */
    double factor, factor_low = 1.0, factor_high = 1.0;    
    double lambda, new_lambda;     /* Karlin-Altschul parameter */
    unsigned int index;                 /* loop index for binary search */
    int** private_matrix;    /* pointer to locally manipulated version of
                                   the matrix */
    int** matrix;
    double scalefactor;    /* local version of amount of scaling */
    double divFactor;      /* 1/scalefacor */
    Blast_ScoreFreq *this_sfp, *return_sfp;    /* score frequency pointers
                                                   to compute lambda */
    double *scoreArray;    /* array of score probabilities */

    scoreArray = (double *) calloc(kScoreMatrixScoreRange, sizeof(double));
    return_sfp = (Blast_ScoreFreq*) calloc(1, sizeof(Blast_ScoreFreq));


    private_matrix = posPrivateMatrix;
    matrix = posMatrix;

    /* Bracket the values. */
    dim1 = compactSearch->qlength;
    dim2 = compactSearch->alphabetSize;

    lambda = compactSearch->lambda_ideal / scalingFactor;
    divFactor = ((double) kPSIScaleFactor) / scalingFactor;
    factor = 1.0;

    if (doBinarySearch) {
        done = FALSE;
        first_time = TRUE;
        while (done != TRUE) {
            for (c = 0; c < dim1; c++) {
                for (a = 0; a < dim2; a++) {
                    if (private_matrix[c][a] == BLAST_SCORE_MIN) {
                        matrix[c][a] = BLAST_SCORE_MIN;
                    } else {
                        matrix[c][a] = (int)
                             ((factor * private_matrix[c][a]) / divFactor);
                    }
                }
            }

            this_sfp =
                fillSfp(matrix, dim1, compactSearch->standardProb,
                        scoreArray, return_sfp);
            if (!this_sfp) {
                sfree(scoreArray);
                Blast_ScoreFreqFree(return_sfp);
                return FALSE;
            }
            new_lambda =
                Blast_KarlinLambdaNR(this_sfp,
                                     compactSearch->kbp_psi[0]->Lambda /
                                     scalingFactor);

            if (new_lambda > lambda) {
                if (first_time) {
                    factor_high = 1.0 + kPositScalingPercent;
                    factor = factor_high;
                    factor_low = 1.0;
                    too_high = TRUE;
                    first_time = FALSE;
                } else {
                    if (too_high == FALSE)
                        break;
                    factor_high += (factor_high - 1.0);
                    factor = factor_high;
                }
            } else {
                if (first_time) {
                    factor_high = 1.0;
                    factor_low = 1.0 - kPositScalingPercent;
                    factor = factor_low;
                    too_high = FALSE;
                    first_time = FALSE;
                } else {
                    if (too_high == TRUE)
                        break;
                    factor_low += (factor_low - 1.0);
                    factor = factor_low;
                }
            }
        }

        /* binary search for ten times. */
        for (index = 0; index < kPositScalingNumIterations; index++) {
            factor = 0.5 * (factor_high + factor_low);
            for (c = 0; c < dim1; c++) {
                for (a = 0; a < dim2; a++) {
                    if (private_matrix[c][a] == BLAST_SCORE_MIN) {
                        matrix[c][a] = BLAST_SCORE_MIN;
                    } else {
                        matrix[c][a] = (int)
                            ((factor * private_matrix[c][a]) / divFactor);
                    }
                }
            }

            this_sfp =
                fillSfp(matrix, dim1, compactSearch->standardProb,
                        scoreArray, return_sfp);
            if (!this_sfp) {
                sfree(scoreArray);
                Blast_ScoreFreqFree(return_sfp);
                return FALSE;
            }
            new_lambda =
                Blast_KarlinLambdaNR(this_sfp,
                                     compactSearch->kbp_psi[0]->Lambda /
                                     scalingFactor);

            if (new_lambda > lambda) {
                factor_low = factor;
            } else {
                factor_high = factor;
            }
        }
    }

    for (c = 0; c < dim1; c++) {
        for (a = 0; a < dim2; a++) {
            if (BLAST_SCORE_MIN != private_matrix[c][a]) {
                matrix[c][a] = BLAST_Nint((double) private_matrix[c][a] *
                                        factor / kPSIScaleFactor);
            }
        }
    }

    _PSIUpdateLambdaK((const int**) matrix, 
                      compactSearch->query, 
                      compactSearch->qlength, 
                      compactSearch->standardProb, 
                      sbp);

    scalefactor = ((double) scalingFactor) / kPSIScaleFactor;
    for (c = 0; c < dim1; c++) {
        for (a = 0; a < dim2; a++) {
            if (BLAST_SCORE_MIN != private_matrix[c][a]) {
                private_matrix[c][a] =
                    BLAST_Nint((double) private_matrix[c][a] * factor *
                             scalefactor);
            }
        }
    }

    sfree(scoreArray);
    Blast_ScoreFreqFree(return_sfp);
    return TRUE;
}

int
Kappa_impalaScaling(Kappa_posSearchItems* posSearch,
                    Kappa_compactSearchItems* compactSearch,
                    double scalingFactor,
                    Boolean doBinarySearch,
                    BlastScoreBlk* sbp)
{
    /* Sorry about the aliasing, this needs refactoring! */
    ASSERT(sbp->kbp_std == compactSearch->kbp_std);
    ASSERT(sbp->kbp_psi == compactSearch->kbp_psi);
    ASSERT(sbp->kbp_gap_std == compactSearch->kbp_gap_std);
    ASSERT(sbp->kbp_gap_psi == compactSearch->kbp_gap_psi);
    ASSERT(sbp->kbp_ideal->Lambda == compactSearch->lambda_ideal);
    ASSERT(sbp->kbp_ideal->K == compactSearch->K_ideal);

    return ((impalaScaleMatrix(compactSearch, 
                               posSearch->posMatrix,
                               posSearch->posPrivateMatrix,
                               scalingFactor, 
                               doBinarySearch,
                               sbp) == TRUE)? 0 : 1);
}
