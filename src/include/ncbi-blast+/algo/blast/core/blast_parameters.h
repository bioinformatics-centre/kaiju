/* $Id: blast_parameters.h 415473 2013-10-17 13:43:48Z madden $
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
 * Author:  Tom Madden
 *
 */

/** @file blast_parameters.h
 * Structure and function definitions for BLAST parameter structures, which are
 * internal to the CORE of BLAST.
 *
 * <pre>
 * These parameters are normally set by:
 *    1.) reading the options in blast_options.[ch] to find user preferences
 *    2.) making intelligent choices based upon the program, user preferences, 
 *        and other data such as the sequence's length.
 *
 * NOTE: These parameters should be set by calls in algo/blast/core, 
 *       preferrably to functions in this file. User preferences should be 
 *       controlled by the structures and functions in blast_options.[ch]. 
 *       The parameter structures belong to algo/blast/core, the options
 *       structures belong to the user.
 * </pre>
 */

#ifndef ALGO_BLAST_CORE__BLASTPARAMETERS__H
#define ALGO_BLAST_CORE__BLASTPARAMETERS__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_export.h>
#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_query_info.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_stat.h>
#include <algo/blast/core/lookup_wrap.h>


#ifdef __cplusplus
extern "C" {
#endif

/** Default parameters for linking HSPs */
#define BLAST_GAP_PROB 0.5         /**< Gap probability for ungapped search */
#define BLAST_GAP_PROB_GAPPED 1.0  /**< Gap probability for gapped search */
#define BLAST_GAP_DECAY_RATE 0.5   /**< Gap decay rate for ungapped search */
#define BLAST_GAP_DECAY_RATE_GAPPED 0.1/**< Gap decay rate for gapped search */
#define BLAST_GAP_SIZE 40          /**< Default gap size */ 
#define BLAST_OVERLAP_SIZE 9      /**< Default overlap size */

/** Expect values corresponding to the default cutoff
 *  scores for all ungapped and gapped blastn alignments.
 */
#define CUTOFF_E_BLASTN 0.05  /**< default evalue (blastn) */
#define CUTOFF_E_BLASTP 1e-300 /**< default evalue (ungapped blastp) */
#define CUTOFF_E_BLASTX 1.0 /**< default evalue (ungapped blastx) */
#define CUTOFF_E_TBLASTN 1.0 /**< default evalue (ungapped tblastn) */
#define CUTOFF_E_TBLASTX 1e-300/**< default evalue (tblastx) */

/** specifies the data structures used for bookkeeping
 *  during computation of ungapped extensions 
 */
typedef enum ESeedContainerType {
    eDiagArray,         /**< use diagonal structures with array of last hits
                           and levels. */
    eDiagHash,          /**< use hash table (blastn only) */
    eMaxContainerType   /**< maximum value for this enumeration */
} ESeedContainerType;

/** All the ungapped cutoff values that can change 
 *  from context to context
 */
typedef struct BlastUngappedCutoffs {
    Int4 x_dropoff_init; /**< Raw X-dropoff value specified by the
                              bit score in BlastInitialWordOptions */
    Int4 x_dropoff;   /**< Raw X-dropoff value used in the ungapped extension */
    Int4 cutoff_score; /**< Cutoff score for saving ungapped hits. */
    Int4 reduced_nucl_cutoff_score; /**< for blastn, a reduced cutoff score
                                      for use with approximate ungapped
                                      alignments */
} BlastUngappedCutoffs;

/** Parameter block that contains a pointer to BlastInitialWordOptions
 * and the values derived from it.
 */
typedef struct BlastInitialWordParameters {
   BlastInitialWordOptions* options; /**< The original (unparsed) options. */

   Int4 x_dropoff_max; /**< largest X-drop cutoff across all contexts */
   Int4 cutoff_score_min; /**< smallest cutoff score across all contexts */
   BlastUngappedCutoffs *cutoffs;   /**< cutoff values (one per context) */
   ESeedContainerType container_type; /**< How to store offset pairs for initial
                                        seeds? */
   Int4 nucl_score_table[256]; /**< the combined score of all match/mismatch
                                    combinations for aligning four bases */
   Boolean matrix_only_scoring; /**< Use the scoring matrix ( not table )
                                     to score ungapped and gapped alignments 
                                     -RMH- */
   Boolean ungapped_extension; /**< Should an ungapped extension be
                                  performed? */
} BlastInitialWordParameters;
    
/** Computed values used as parameters for gapped alignments.
 *  Note that currently the value of the X-dropoff parameter
 *  is fixed for all search contexts
 */
typedef struct BlastExtensionParameters {
   BlastExtensionOptions* options; /**< The original (unparsed) options. */
   Int4 gap_x_dropoff; /**< X-dropoff value for gapped extension (raw) */
   Int4 gap_x_dropoff_final;/**< X-dropoff value for the final gapped 
                               extension (raw) */
} BlastExtensionParameters;

/** Parameter block for linking HSPs with sum statistics. */
typedef struct BlastLinkHSPParameters {
   double gap_prob;       /**< Probability of decay for linking HSPs */
   Int4 gap_size;          /**< Small gap size for linking HSPs */
   Int4 overlap_size;     /**< Maximal overlap allowed in successive linked
                             HSPs */
   double gap_decay_rate; /**< Decay rate for linking HSPs and calculating
                             cutoff scores. */
   Int4 cutoff_small_gap; /**< Cutoff sum score for linked HSPs with small 
                             gaps. Small gap calculations are ignored if 
                             this value is set to 0. */
   Int4 cutoff_big_gap; /**< Cutoff sum score for linked HSPs with big gaps. */
   Int4 longest_intron; /**< Length of a longest intron for uneven gap linking
                           of HSPs. */
} BlastLinkHSPParameters;

/** All the gapped cutoff values that can change 
 *  from context to context
 */
typedef struct BlastGappedCutoffs {
   Int4 cutoff_score; /**< Raw cutoff score corresponding to the e-value 
                         provided by the user if no sum stats, the lowest score
                         to attempt linking on if sum stats are used. */
   Int4 cutoff_score_max; /**< Raw cutoff score corresponding to the e-value 
                         provided by user, cutoff_score must be <= this. */
} BlastGappedCutoffs;

/** Parameter block that contains a pointer to BlastHitSavingOptions
 * and the values derived from it.
 */
typedef struct BlastHitSavingParameters {
   BlastHitSavingOptions* options; /**< The original (unparsed) options. */
   Int4 cutoff_score_min; /**< smallest cutoff score across all contexts */
   BlastGappedCutoffs *cutoffs; /**< per-context gapped cutoff information */
   BlastLinkHSPParameters* link_hsp_params; /**< Parameters for linking HSPs
                                               with sum statistics; linking 
                                               is not done if NULL. */
   Boolean restricted_align; /**< TRUE if approximate score-only gapped
                                  alignment is used */
   Boolean do_sum_stats;  /**< TRUE if sum stats will be used.  Can override the
                                  do_sum_stats Boolean in the options if criteria for 
                                  doing sum stats are not met.  */
   Int4 mask_level; /**< Only keep the highest scoring HSP when more than
                          one HSP overlaps the same region of the query by
                          more than or equal to mask_level %. -RMH- */
   Int4* low_score;             /**< lowest ungapped score that can trigger a
				gapped alignment if the histlist is already full.
                                One value for each query. */
   double prelim_evalue; /**< evalue for preliminary search (may be higher for CBS). */
} BlastHitSavingParameters;

/** Because approximate gapped alignment adds extra overhead,
 *  it should be avoided if there is no performance benefit
 *  to using it. To be benficial, most score-only gapped alignments
 *  should be the restricted kind, and we can only assure that
 *  if the e-value cutoff is not too generous
 */
#define RESTRICTED_ALIGNMENT_WORST_EVALUE 10.0

/** Scoring parameters block
 *  Contains scoring-related information that is actually used
 *  for the blast search
 */
typedef struct BlastScoringParameters {
   BlastScoringOptions *options; /**< User-provided values for these params */
   Int2 reward;      /**< Reward for a match */
   Int2 penalty;     /**< Penalty for a mismatch */
   Int4 gap_open;    /**< Extra penalty for starting a gap (scaled version) */
   Int4 gap_extend;  /**< Penalty for each gap residue  (scaled version) */
   Int4 shift_pen;   /**< Penalty for shifting a frame in out-of-frame 
                        gapping (scaled version) */
   double scale_factor; /**< multiplier for all cutoff scores */
} BlastScoringParameters;

/** Parameters for setting up effective lengths and search spaces.  
 * The real database size values to be used for statistical calculations, if
 * there are no overriding values in options.
 */
typedef struct BlastEffectiveLengthsParameters {
   BlastEffectiveLengthsOptions* options; /**< User provided values for these 
                                             parameters */
   Int8 real_db_length; /**< Total database length to use in search space
                           calculations. */
   Int4 real_num_seqs;  /**< Number of subject sequences to use for search
                           space calculations */
} BlastEffectiveLengthsParameters;

/********************************************************************************

    Functions to create options blocks with default values
    and free them after use.

*********************************************************************************/


/** Deallocate memory for BlastInitialWordParameters.
 * @param parameters Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastInitialWordParameters*
BlastInitialWordParametersFree(BlastInitialWordParameters* parameters);

/** Allocate memory for BlastInitialWordParameters and set x_dropoff.
 * Calling BlastInitialWordParametersNew calculates the
 * raw x_dropoff from the bit x_dropoff and puts it into
 * the x_dropoff field of BlastInitialWordParameters*.
 * The container type is also set.  For blastn queries over a certain
 * length eDiagHash is set, otherwise it's eDiagArray.
 * The extension method is also set via a call to s_GetBestExtensionMethod
 *
 * @param program_number Type of BLAST program [in]
 * @param word_options The initial word options [in]
 * @param hit_params The hit saving options (needed to calculate cutoff score 
 *                    for ungapped extensions) [in]
 * @param lookup_wrap contains lookup table so the proper extension 
 *        method can be set [in]
 * @param sbp Statistical (Karlin-Altschul) information [in]
 * @param query_info Query information [in]
 * @param subject_length Average subject sequence length [in]
 * @param parameters Resulting parameters [out]
*/
NCBI_XBLAST_EXPORT
Int2
BlastInitialWordParametersNew(EBlastProgramType program_number, 
   const BlastInitialWordOptions* word_options, 
   const BlastHitSavingParameters* hit_params, 
   const LookupTableWrap* lookup_wrap,
   const BlastScoreBlk* sbp, 
   BlastQueryInfo* query_info, 
   Uint4 subject_length,
   BlastInitialWordParameters* *parameters);

/** Update cutoff scores in BlastInitialWordParameters structure.
 * @param program_number Type of BLAST program [in]
 * @param hit_params The hit saving parameters, needed to calculate cutoff 
 *                   score for ungapped extensions. The HSP linking cutoff
 *                   might have to be adjusted here. [in] [out]
 * @param sbp Statistical (Karlin-Altschul) information [in]
 * @param query_info Query information [in]
 * @param subject_length Average subject sequence length [in]
 * @param parameters Preallocated parameters [in] [out]
*/
NCBI_XBLAST_EXPORT
Int2
BlastInitialWordParametersUpdate(EBlastProgramType program_number, 
   const BlastHitSavingParameters* hit_params, 
   const BlastScoreBlk* sbp, 
   BlastQueryInfo* query_info, Uint4 subject_length,
   BlastInitialWordParameters* parameters);

/** Calculate the raw values for the X-dropoff parameters 
 * @param blast_program Program number [in]
 * @param options Already allocated extension options [in]
 * @param sbp Structure containing statistical information [in]
 * @param query_info Query information, needed only for determining the first 
 *                   context [in]
 * @param parameters Extension parameters [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastExtensionParametersNew(EBlastProgramType blast_program, 
        const BlastExtensionOptions* options, 
        BlastScoreBlk* sbp, BlastQueryInfo* query_info, 
        BlastExtensionParameters* *parameters);

/** Deallocate memory for BlastExtensionParameters. 
 * @param parameters Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastExtensionParameters*
BlastExtensionParametersFree(BlastExtensionParameters* parameters);

/**  Deallocate memory for BlastScoringParameters.
 * @param parameters Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastScoringParameters* BlastScoringParametersFree(
                                     BlastScoringParameters* parameters);

/** Calculate scaled cutoff scores and gap penalties
 * @param options Already allocated scoring options [in]
 * @param sbp Structure containing scale factor [in]
 * @param parameters Scoring parameters [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastScoringParametersNew(const BlastScoringOptions *options,
                               BlastScoreBlk* sbp, 
                               BlastScoringParameters* *parameters);

/** Deallocate memory for BlastEffectiveLengthsParameters*. 
 * @param parameters Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastEffectiveLengthsParameters* 
BlastEffectiveLengthsParametersFree(BlastEffectiveLengthsParameters* parameters);

/** Allocate memory for BlastEffectiveLengthsParameters 
 * @param options The user provided options [in]
 * @param db_length The database length [in]
 * @param num_seqs Number of sequences in database [in]
 * @param parameters The parameters structure returned [out]
 */
NCBI_XBLAST_EXPORT
Int2 
BlastEffectiveLengthsParametersNew(const BlastEffectiveLengthsOptions* options, 
                               Int8 db_length, Int4 num_seqs,
                               BlastEffectiveLengthsParameters* *parameters);

/** Deallocate memory for BlastLinkHSPParameters;
 * @param parameters Structure to free [in] 
 */
NCBI_XBLAST_EXPORT
BlastLinkHSPParameters* 
BlastLinkHSPParametersFree(BlastLinkHSPParameters* parameters);

/** Initialize the linking HSPs parameters with default values.
 * @param program_number Type of BLAST program [in]
 * @param gapped_calculation Is this a gapped search? [in]
 * @param link_hsp_params Initialized parameters structure [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastLinkHSPParametersNew(EBlastProgramType program_number, 
                               Boolean gapped_calculation,
                               BlastLinkHSPParameters** link_hsp_params);

/** Update BlastLinkHSPParameters, using calculated values of other parameters.
 * @param word_params Initial word parameters [in]
 * @param hit_params Hit saving parameters, including the link HSP 
 *                   parameters [in] [out]
 * @param gapped_calculation Is this a gapped search? [in]
 */
NCBI_XBLAST_EXPORT
Int2 
BlastLinkHSPParametersUpdate(const BlastInitialWordParameters* word_params,
                             const BlastHitSavingParameters* hit_params,
                             Boolean gapped_calculation);

/** Deallocate memory for BlastHitSavingOptions*. 
 * @param parameters Structure to free [in]
 */
NCBI_XBLAST_EXPORT
BlastHitSavingParameters*
BlastHitSavingParametersFree(BlastHitSavingParameters* parameters);

/** Allocate memory and initialize the BlastHitSavingParameters structure. 
 * Calculates the (raw) score cutoff given an expect value and puts
 * it in the "cutoff_score" field of the returned BlastHitSavingParameters*
 *
 * @param program_number Number of the BLAST program [in]
 * @param options The given hit saving options [in]
 * @param sbp Scoring block, needed for calculating score cutoff from 
 *            e-value [in]
 * @param query_info Query information, needed for calculating score cutoff 
 *                   from e-value [in]
 * @param avg_subject_length average length of subject sequence [in]
 * @param compositionBasedStats follows -t option in blastp [in]
 * @param parameters Resulting parameters [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastHitSavingParametersNew(EBlastProgramType program_number, 
        const BlastHitSavingOptions* options, 
        const BlastScoreBlk* sbp, 
        const BlastQueryInfo* query_info, 
        Int4 avg_subject_length,
   	Int4 compositionBasedStats,
        BlastHitSavingParameters* *parameters);

/** Updates cutoff scores in hit saving parameters. 
 * @param program_number Number of the BLAST program [in]
 * @param sbp Scoring block, needed for calculating score cutoff from 
 *            e-value [in]
 * @param query_info Query information, needed for calculating score cutoff 
 *                   from e-value [in]
 * @param avg_subject_length average length of subject sequence, used in sum_stats
 *            mode [in]
 * @param compositionBasedStats follows -t option in blastp [in]
 * @param parameters Preallocated parameters [in] [out]
 */
NCBI_XBLAST_EXPORT
Int2 BlastHitSavingParametersUpdate(EBlastProgramType program_number, 
        const BlastScoreBlk* sbp, const BlastQueryInfo* query_info, 
        Int4 avg_subject_length,
   	Int4 compositionBasedStats,
        BlastHitSavingParameters* parameters);

/** Calculates cutoff scores and returns them.
 *  Equations provided by Stephen Altschul.
 * @param program BLAST program type [in]
 * @param query_info Query(ies) information [in]
 * @param sbp Scoring statistical parameters [in]
 * @param link_hsp_params Parameters for linking HSPs [in] [out]
 * @param word_params InitialWordParameters  (cutoff_score used for small gaps) [in]
 * @param db_length Total length of database (non-database search if 0) [in]
 * @param subject_length Length of the subject sequence. [in]
 * 
*/
NCBI_XBLAST_EXPORT
void
CalculateLinkHSPCutoffs(EBlastProgramType program, BlastQueryInfo* query_info, 
   const BlastScoreBlk* sbp, BlastLinkHSPParameters* link_hsp_params, 
   const BlastInitialWordParameters* word_params,
   Int8 db_length, Int4 subject_length);

// Routines for debugging parameter settings.  -RMH-
void printBlastScoringParameters( BlastScoringParameters* params );
void printBlastInitialWordParamters ( BlastInitialWordParameters *word_params,
                                      BlastQueryInfo *query_info );
void printBlastExtensionParameters ( BlastExtensionParameters *ext_params );
void printBlastHitSavingParameters ( BlastHitSavingParameters * hit_params,
                                     BlastQueryInfo *query_info );
void printAllParameters ( BlastHitSavingParameters * hit_params,
                          BlastExtensionParameters *ext_params,
                          BlastInitialWordParameters *word_params,
                          BlastQueryInfo *query_info  );

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLASTPARAMETERS__H */

