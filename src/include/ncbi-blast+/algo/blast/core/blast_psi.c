#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] =
    "$Id: blast_psi.c 459883 2015-02-23 16:45:34Z camacho $";
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
 * Author:  Christiam Camacho
 *
 */

/** @file blast_psi.c
 * Implementation of the high level functions of PSI-BLAST's PSSM engine.
 */
    
#include <algo/blast/core/blast_encoding.h>
#include "blast_psi_priv.h"

/* needed for BLAST_GetStandardAaProbabilities(); */
#include <algo/blast/core/blast_util.h> 

/****************************************************************************/
/* Function prototypes */

/** Convenience function to deallocate data structures allocated in
 * PSICreatePssmWithDiagnostics.
 * @param pssm PSSM and statistical information [in|out]
 * @param packed_msa compact multiple sequence alignment structure[in]
 * @param msa multiple sequence alignment structure[in]
 * @param aligned_block aligned blocks data structure [in] 
 * @param seq_weights sequence weights data structure [in]
 * @param internal_pssm PSSM being computed [in]
 */
static void
s_PSICreatePssmCleanUp(PSIMatrix** pssm,
                       _PSIPackedMsa* packed_msa,
                       _PSIMsa* msa,
                       _PSIAlignedBlock* aligned_block,
                       _PSISequenceWeights* seq_weights,
                       _PSIInternalPssmData* internal_pssm);

/** Copies pssm data from internal_pssm and sbp into pssm. None of its
 * parameters can be NULL. 
 * @param internal_pssm PSSM being computed [in]
 * @param sbp Score block structure containing the calculated lambda and K
 * which will be saved in the pssm parameter [in]
 * @param pssm PSSM and statistical information [in|out]
 */
static void
s_PSISavePssm(const _PSIInternalPssmData* internal_pssm,
              const BlastScoreBlk* sbp,
              PSIMatrix* pssm);

/** Private function which performs the last 2 stages of the PSSM creation:
 * conversion of PSSM frequecy ratios to PSSM and scaling of the PSSM.
 * @param internal_pssm PSSM being computed, must be already allocated [in|out]
 * @param query query sequence in ncbistdaa encoding. [in]
 * @param query_length length of the query sequence above [in]
 * @param std_prob array containing the standard background residue 
 * probabilities [in]
 * @param sbp Score block structure where the calculated lambda and K will be
 * returned [in|out]
 * @param impala_scaling_factor scaling factor used in IMPALA-style scaling if
 * its value is NOT kPSSM_NoImpalaScaling (otherwise it performs standard
 * PSI-BLAST scaling) [in]
 */
static int
_PSICreateAndScalePssmFromFrequencyRatios(_PSIInternalPssmData* internal_pssm,
                                          const Uint1* query,
                                          Uint4 query_length,
                                          double* std_prob,
                                          BlastScoreBlk* sbp,
                                          double impala_scaling_factor);
      
/****************************************************************************/

int
PSICreatePssm(const PSIMsa* msap,
              const PSIBlastOptions* options,
              BlastScoreBlk* sbp,
              PSIMatrix** pssm)
{
    return PSICreatePssmWithDiagnostics(msap, options, sbp, NULL,
                                        pssm, NULL);
}

int
PSICreatePssmWithDiagnostics(const PSIMsa* msap,                    /* [in] */
                             const PSIBlastOptions* options,        /* [in] */
                             BlastScoreBlk* sbp,                    /* [in] */
                             const PSIDiagnosticsRequest* request,  /* [in] */
                             PSIMatrix** pssm,                      /* [out] */
                             PSIDiagnosticsResponse** diagnostics)  /* [out] */
{
    _PSIMsa* msa = NULL;
    _PSIAlignedBlock* aligned_block = NULL;
    _PSISequenceWeights* seq_weights = NULL; 
    _PSIInternalPssmData* internal_pssm = NULL;
    _PSIPackedMsa* packed_msa = NULL;
    int status = 0;

    if ( !msap || !options || !sbp || !pssm ) {
        return PSIERR_BADPARAM;
    }

    packed_msa = _PSIPackedMsaNew(msap);

    /*** Run the engine's stages ***/

    status = _PSIPurgeBiasedSegments(packed_msa);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return status;
    }

    /*** Allocate data structures ***/
    msa = _PSIMsaNew(packed_msa, (Uint4) sbp->alphabet_size);
    aligned_block = _PSIAlignedBlockNew(msa->dimensions->query_length);
    seq_weights = _PSISequenceWeightsNew(msa->dimensions, sbp);
    internal_pssm = _PSIInternalPssmDataNew(msa->dimensions->query_length,
                                            (Uint4) sbp->alphabet_size);
    *pssm = PSIMatrixNew(msa->dimensions->query_length, 
                         (Uint4) sbp->alphabet_size);
    if ( !msa || ! aligned_block || !seq_weights || !internal_pssm || !*pssm ) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return PSIERR_OUTOFMEM;
    }
    packed_msa = _PSIPackedMsaFree(packed_msa);

    /*** Enable structure group customization if needed and validate the
     * multiple sequence alignment data ***/
    if (options->nsg_compatibility_mode) {
        _PSIStructureGroupCustomization(msa);
        status = _PSIValidateMSA_StructureGroup(msa);
    } else {
        status = _PSIValidateMSA(msa, options->ignore_unaligned_positions);
    }
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return status;
    }

    status = _PSIComputeAlignmentBlocks(msa, aligned_block);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return status;
    }

    status = _PSIComputeSequenceWeights(msa, aligned_block, 
                                        options->nsg_compatibility_mode,
                                        seq_weights);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return status;
    }

    status = _PSIComputeFreqRatios(msa, seq_weights, sbp, aligned_block, 
                                   options->pseudo_count, 
                                   options->nsg_compatibility_mode,
                                   internal_pssm);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return status;
    }

    status = _PSICreateAndScalePssmFromFrequencyRatios
        (internal_pssm, msa->query, msa->dimensions->query_length, 
         seq_weights->std_prob, sbp, options->impala_scaling_factor);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                               seq_weights, internal_pssm);
        return status;
    }
    /*** Save the pssm outgoing parameter ***/
    s_PSISavePssm(internal_pssm, sbp, *pssm);


    /*** Save diagnostics if required ***/
    if (request && diagnostics) {
        *diagnostics = PSIDiagnosticsResponseNew(msa->dimensions->query_length,
                                                 (Uint4) sbp->alphabet_size,
                                                 request);
        if ( !*diagnostics ) {
            /* FIXME: This could be changed to return a warning and not
             * deallocate PSSM data */
            s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                                   seq_weights, internal_pssm);
            return PSIERR_OUTOFMEM;
        }
        status = _PSISaveDiagnostics(msa, aligned_block, seq_weights, 
                                     internal_pssm, *diagnostics);
        if (status != PSI_SUCCESS) {
            *diagnostics = PSIDiagnosticsResponseFree(*diagnostics);
            s_PSICreatePssmCleanUp(pssm, packed_msa, msa, aligned_block, 
                                   seq_weights, internal_pssm);
            return status;
        }
    }
    s_PSICreatePssmCleanUp(NULL, packed_msa, msa, aligned_block, seq_weights, 
                           internal_pssm);

    return PSI_SUCCESS;
}

int
PSICreatePssmFromCDD(const PSICdMsa* cd_msa,                /* [in] */
                     const PSIBlastOptions* options,        /* [in] */
                     BlastScoreBlk* sbp,                    /* [in] */
                     const PSIDiagnosticsRequest* request,  /* [in] */
                     PSIMatrix** pssm,                      /* [out] */
                     PSIDiagnosticsResponse** diagnostics)  /* [out] */
{
    _PSISequenceWeights* seq_weights = NULL; 
    _PSIInternalPssmData* internal_pssm = NULL;
    int status = 0;

    if ( !cd_msa || !options || !sbp || !pssm ) {
        return PSIERR_BADPARAM;
    }

    /*** Run the engine's stages ***/


    /*** Allocate data structures ***/
    seq_weights = _PSISequenceWeightsNew(cd_msa->dimensions, sbp);
    internal_pssm = _PSIInternalPssmDataNew(cd_msa->dimensions->query_length,
                                            (Uint4) sbp->alphabet_size);
    *pssm = PSIMatrixNew(cd_msa->dimensions->query_length, 
                         (Uint4) sbp->alphabet_size);
    if ( !seq_weights || !internal_pssm || !*pssm ) {
        s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                               internal_pssm);
        return PSIERR_OUTOFMEM;
    }

    status = _PSIValidateCdMSA(cd_msa, sbp->alphabet_size);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                               internal_pssm);
        return status;
    }

    status = _PSIComputeFrequenciesFromCDs(cd_msa, sbp, options, seq_weights);

    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                               internal_pssm);
        return status;
    }

    status = _PSIComputeFreqRatiosFromCDs(cd_msa, seq_weights, sbp,  
                                          options->pseudo_count, 
                                          internal_pssm);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                               internal_pssm);
        return status;
    }

    status = _PSICreateAndScalePssmFromFrequencyRatios
        (internal_pssm, cd_msa->query, cd_msa->dimensions->query_length, 
         seq_weights->std_prob, sbp, options->impala_scaling_factor);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                               internal_pssm);
        return status;
    }
    /*** Save the pssm outgoing parameter ***/
    s_PSISavePssm(internal_pssm, sbp, *pssm);


    /*** Save diagnostics if required ***/
    if (request && diagnostics) {

        *diagnostics = PSIDiagnosticsResponseNew(
                                        cd_msa->dimensions->query_length,
                                        (Uint4) sbp->alphabet_size,
                                        request);
        if ( !*diagnostics ) {
            /* FIXME: This could be changed to return a warning and not
              * deallocate PSSM data */
            s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                                   internal_pssm);
            return PSIERR_OUTOFMEM;
        }
        status = _PSISaveCDDiagnostics(cd_msa, seq_weights, internal_pssm,
                                       *diagnostics);
        if (status != PSI_SUCCESS) {
            *diagnostics = PSIDiagnosticsResponseFree(*diagnostics);
            s_PSICreatePssmCleanUp(pssm, NULL, NULL, NULL, seq_weights,
                                   internal_pssm);
            return status;
        }
    }
    s_PSICreatePssmCleanUp(NULL, NULL, NULL, NULL, seq_weights, internal_pssm);

    return PSI_SUCCESS;

}

/** Convenience function to deallocate data structures allocated in
 * PSICreatePssmFromFrequencyRatios
 * @param pssm PSSM and statistical information [in|out]
 * @param internal_pssm PSSM being computed [in]
 * @param std_prob array containing the standard background residue 
 * probabilities [in]
 */
static void
s_PSICreatePssmFromFrequencyRatiosCleanUp(PSIMatrix** pssm,
                                          _PSIInternalPssmData* internal_pssm,
                                          double* std_prob)
{
    if (pssm) {
        *pssm = PSIMatrixFree(*pssm);
    }
    _PSIInternalPssmDataFree(internal_pssm);
    sfree(std_prob);
}

int
PSICreatePssmFromFrequencyRatios(const Uint1* query,
                                 Uint4 query_length,
                                 BlastScoreBlk* sbp,
                                 double** freq_ratios,
                                 double impala_scaling_factor,
                                 PSIMatrix** pssm)
{
    int status = PSI_SUCCESS;
    double* std_prob = NULL;
    _PSIInternalPssmData* internal_pssm = NULL;

    std_prob = BLAST_GetStandardAaProbabilities();
    *pssm = PSIMatrixNew(query_length, (Uint4) sbp->alphabet_size);
    internal_pssm = _PSIInternalPssmDataNew(query_length, sbp->alphabet_size);

    if ( !std_prob || !*pssm || !internal_pssm ) {
        s_PSICreatePssmFromFrequencyRatiosCleanUp(pssm, internal_pssm,
                                                  std_prob);
        return PSIERR_OUTOFMEM;
    }

    _PSICopyMatrix_double(internal_pssm->freq_ratios, freq_ratios, 
                          internal_pssm->ncols, internal_pssm->nrows);

    status = _PSICreateAndScalePssmFromFrequencyRatios(internal_pssm, 
                                                       query, query_length, 
                                                       std_prob, sbp, 
                                                       impala_scaling_factor);
    if (status != PSI_SUCCESS) {
        s_PSICreatePssmFromFrequencyRatiosCleanUp(pssm, internal_pssm,
                                                  std_prob);
        return status;
    }
    /*** Save the pssm outgoing parameter ***/
    s_PSISavePssm(internal_pssm, sbp, *pssm);

    s_PSICreatePssmFromFrequencyRatiosCleanUp(NULL, internal_pssm, std_prob);
    return status;
}

static int
_PSICreateAndScalePssmFromFrequencyRatios(_PSIInternalPssmData* internal_pssm,
                                          const Uint1* query,
                                          Uint4 query_length,
                                          double* std_prob,
                                          BlastScoreBlk* sbp,
                                          double impala_scaling_factor)
{
    int status = PSI_SUCCESS;

    ASSERT(internal_pssm);
    ASSERT(query);
    ASSERT(std_prob);
    ASSERT(sbp);

    status = _PSIConvertFreqRatiosToPSSM(internal_pssm, query, sbp, std_prob);
    if (status != PSI_SUCCESS) {
        /* clean up is done in calling code */
        return status;
    }

    if (impala_scaling_factor == kPSSM_NoImpalaScaling) {
        status = _PSIScaleMatrix(query, std_prob, internal_pssm, sbp);
    } else {
        status = _IMPALAScaleMatrix(query, std_prob, internal_pssm, sbp,
                                    impala_scaling_factor);
    }
    if (status != PSI_SUCCESS) {
        /* clean up is done in calling code */
        return status;
    }

    return status;
}

/****************************************************************************/

static void
s_PSICreatePssmCleanUp(PSIMatrix** pssm,
                       _PSIPackedMsa* packed_msa,
                       _PSIMsa* msa,
                       _PSIAlignedBlock* aligned_block,
                       _PSISequenceWeights* seq_weights,
                       _PSIInternalPssmData* internal_pssm)
{
    if (pssm) {
        *pssm = PSIMatrixFree(*pssm);
    }
    _PSIPackedMsaFree(packed_msa);
    _PSIMsaFree(msa);
    _PSIAlignedBlockFree(aligned_block);
    _PSISequenceWeightsFree(seq_weights);
    _PSIInternalPssmDataFree(internal_pssm);
}

static void
s_PSISavePssm(const _PSIInternalPssmData* internal_pssm,
             const BlastScoreBlk* sbp,
             PSIMatrix* pssm)
{
    ASSERT(internal_pssm);
    ASSERT(sbp);
    ASSERT(pssm);

    _PSICopyMatrix_int(pssm->pssm, internal_pssm->pssm,
                       pssm->ncols, pssm->nrows);

    pssm->lambda = sbp->kbp_gap_psi[0]->Lambda;
    pssm->kappa = sbp->kbp_gap_psi[0]->K;
    pssm->h = sbp->kbp_gap_psi[0]->H;
    pssm->ung_lambda = sbp->kbp_psi[0]->Lambda;
    pssm->ung_kappa = sbp->kbp_psi[0]->K;
    pssm->ung_h = sbp->kbp_psi[0]->H;
}

/****************************************************************************/

PSIMsa*
PSIMsaNew(const PSIMsaDimensions* dimensions)
{
    PSIMsa* retval = NULL;

    if ( !dimensions ) {
        return NULL;
    }

    retval = (PSIMsa*) calloc(1, sizeof(PSIMsa));
    if ( !retval ) {
        return PSIMsaFree(retval);
    }

    retval->dimensions = (PSIMsaDimensions*) malloc(sizeof(PSIMsaDimensions));
    if ( !retval->dimensions ) {
        return PSIMsaFree(retval);
    }
    memcpy((void*) retval->dimensions,
           (void*) dimensions, 
           sizeof(PSIMsaDimensions));

    retval->data = (PSIMsaCell**) _PSIAllocateMatrix(dimensions->num_seqs + 1,
                                                     dimensions->query_length,
                                                     sizeof(PSIMsaCell));
    if ( !retval->data ) {
        return PSIMsaFree(retval);
    }
    {
        Uint4 s = 0;    /* index on sequences */
        Uint4 p = 0;    /* index on positions */

        for (s = 0; s < dimensions->num_seqs + 1; s++) {
            for (p = 0; p < dimensions->query_length; p++) {
                retval->data[s][p].letter = 0;
                retval->data[s][p].is_aligned = FALSE;
            }
        }
    }

#ifdef DEBUG_PSSM_ENGINE
    retval->seqinfo = (PSISeqInfo*) calloc(dimensions->num_seqs + 1,
                                           sizeof(PSISeqInfo));
    if ( !retval->seqinfo ) {
        return PSIMsaFree(retval);
    }
#endif /* DEBUG_PSSM_ENGINE */

    return retval;
}

PSIMsa*
PSIMsaFree(PSIMsa* msa)
{
    if ( !msa ) {
        return NULL;
    }

    if ( msa->data && msa->dimensions ) {
        _PSIDeallocateMatrix((void**) msa->data,
                             msa->dimensions->num_seqs + 1);
        msa->data = NULL;
    }

    if ( msa->dimensions ) {
        sfree(msa->dimensions);
    }

#ifdef DEBUG_PSSM_ENGINE
    if ( msa->seqinfo ) {
        sfree(msa->seqinfo);
    }
#endif /* DEBUG_PSSM_ENGINE */

    sfree(msa);

    return NULL;
}

PSIMatrix*
PSIMatrixNew(Uint4 query_length, Uint4 alphabet_size)
{
    PSIMatrix* retval = NULL;

    retval = (PSIMatrix*) malloc(sizeof(PSIMatrix));
    if ( !retval ) {
        return NULL;
    }
    retval->ncols = query_length;
    retval->nrows = alphabet_size;

    retval->pssm = (int**) _PSIAllocateMatrix(query_length, alphabet_size,
                                              sizeof(int));
    if ( !(retval->pssm) ) {
        return PSIMatrixFree(retval);
    }

    retval->lambda = 0.0;
    retval->kappa = 0.0;
    retval->h = 0.0;
    retval->ung_lambda = 0.0;
    retval->ung_kappa = 0.0;
    retval->ung_h = 0.0;

    return retval;
}

PSIMatrix*
PSIMatrixFree(PSIMatrix* matrix)
{
    if ( !matrix ) {
        return NULL;
    }

    if (matrix->pssm) {
        _PSIDeallocateMatrix((void**) matrix->pssm, matrix->ncols);
    }

    sfree(matrix);

    return NULL;
}

PSIDiagnosticsRequest* 
PSIDiagnosticsRequestNew(void)
{
    return calloc(1, sizeof(PSIDiagnosticsRequest));
}

PSIDiagnosticsRequest* 
PSIDiagnosticsRequestNewEx(Boolean save_ascii_pssm)
{
    PSIDiagnosticsRequest* retval = PSIDiagnosticsRequestNew(); 
    if ( !retval ) {
        return NULL;
    }

    retval->frequency_ratios = TRUE;
    if (save_ascii_pssm) {
        retval->information_content = TRUE;
        retval->weighted_residue_frequencies = TRUE;
        retval->gapless_column_weights = TRUE;
        retval->sigma = TRUE;
        retval->interval_sizes = TRUE;
        retval->num_matching_seqs = TRUE;
    }
    return retval;
}

PSIDiagnosticsRequest* 
PSIDiagnosticsRequestFree(PSIDiagnosticsRequest* diags_request)
{
    sfree(diags_request);
    return NULL;
}

PSIDiagnosticsResponse*
PSIDiagnosticsResponseNew(Uint4 query_length, Uint4 alphabet_size,
                          const PSIDiagnosticsRequest* wants)
{
    PSIDiagnosticsResponse* retval = NULL;

    if ( !wants ) {
        return NULL;
    }

    /* MUST use calloc to allocate structure because code that uses this
     * structure assumes that non-NULL members will require to be populated */
    retval = (PSIDiagnosticsResponse*) calloc(1, 
                                              sizeof(PSIDiagnosticsResponse));
    if ( !retval ) {
        return NULL;
    }

    retval->query_length = query_length;
    retval->alphabet_size = alphabet_size;

    if (wants->information_content) {
        retval->information_content = (double*) 
            calloc(query_length, sizeof(double));
        if ( !(retval->information_content) ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->residue_frequencies) {
        retval->residue_freqs = (Uint4**) _PSIAllocateMatrix(query_length, 
                                                             alphabet_size, 
                                                             sizeof(Uint4));
        if ( !(retval->residue_freqs) ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->weighted_residue_frequencies) {
        retval->weighted_residue_freqs = (double**) 
            _PSIAllocateMatrix(query_length, 
                               alphabet_size, 
                               sizeof(double));
        if ( !(retval->weighted_residue_freqs) ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->frequency_ratios) {
        retval->frequency_ratios = (double**)
            _PSIAllocateMatrix(query_length, 
                               alphabet_size, 
                               sizeof(double));
        if ( !retval->frequency_ratios ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->gapless_column_weights) {
        retval->gapless_column_weights = (double*) 
            calloc(query_length, sizeof(double));
        if ( !(retval->gapless_column_weights) ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->sigma) {
        retval->sigma = (double*) calloc(query_length, sizeof(double));
        if ( !retval->sigma ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->interval_sizes) {
        retval->interval_sizes = (Uint4*) calloc(query_length, sizeof(Uint4));
        if ( !retval->interval_sizes ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->num_matching_seqs) {
        retval->num_matching_seqs = 
            (Uint4*) calloc(query_length, sizeof(Uint4));
        if ( !retval->num_matching_seqs ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    if (wants->independent_observations) {
        retval->independent_observations =
            (double*) calloc(query_length, sizeof(double));
        if ( !retval->independent_observations ) {
            return PSIDiagnosticsResponseFree(retval);
        }
    }

    return retval;
}

PSIDiagnosticsResponse*
PSIDiagnosticsResponseFree(PSIDiagnosticsResponse* diags)
{
    if ( !diags )
        return NULL;

    if (diags->information_content) {
        sfree(diags->information_content);
    }

    if (diags->residue_freqs) {
        _PSIDeallocateMatrix((void**) diags->residue_freqs,
                             diags->query_length);
    }

    if (diags->weighted_residue_freqs) {
        _PSIDeallocateMatrix((void**) diags->weighted_residue_freqs,
                             diags->query_length);
    }

    if (diags->frequency_ratios) {
        _PSIDeallocateMatrix((void**) diags->frequency_ratios,
                             diags->query_length);
    }

    if (diags->gapless_column_weights) {
        sfree(diags->gapless_column_weights);
    }

    if (diags->independent_observations) {
        sfree(diags->independent_observations);
    }

    sfree(diags);

    return NULL;
}

