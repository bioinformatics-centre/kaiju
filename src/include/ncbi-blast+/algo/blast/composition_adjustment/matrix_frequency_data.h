/* $Id: matrix_frequency_data.h 103491 2007-05-04 17:18:18Z kazimird $
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
 * ===========================================================================*/
/**
 * @file matrix_frequency_data.h
 * Definitions used to get joint probabilities for a scoring matrix
 *
 * @author Alejandro Schaffer, E. Michael Gertz
 */
#ifndef __MATRIX_FREQUENCY_DATA__
#define __MATRIX_FREQUENCY_DATA__

#include <algo/blast/core/blast_export.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Get joint probabilities for the named matrix.
 *
 * @param probs        the joint probabilities [out]
 * @param row_sums     sum of the values in each row of probs [out]
 * @param col_sums     sum of the values in each column of probs [out]
 * @param matrixName   the name of the matrix sought [in]
 * @returns 0 if successful; -1 if the named matrix is not known.
 */
NCBI_XBLAST_EXPORT
int Blast_GetJointProbsForMatrix(double ** probs, double row_sums[],
                                 double col_sums[], const char *matrixName);


/** Return true if frequency data is available for the given matrix name. */
NCBI_XBLAST_EXPORT
const double * Blast_GetMatrixBackgroundFreq(const char *matrix_name);


/** Retrieve the background letter probabilities implicitly used in
 * constructing the score matrix matrix_name. */
NCBI_XBLAST_EXPORT
int Blast_FrequencyDataIsAvailable(const char *matrix_name);

#ifdef __cplusplus
}
#endif

#endif
