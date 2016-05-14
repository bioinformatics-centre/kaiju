/* $Id: blast_hits.h 458579 2015-02-06 15:12:03Z fongah2 $
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
 * Author:  Ilya Dondoshansky
 *
 */

/** @file blast_hits.h
 * Structures and API used for saving BLAST hits
 */

#ifndef ALGO_BLAST_CORE__BLAST_HITS__H
#define ALGO_BLAST_CORE__BLAST_HITS__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_export.h>
#include <algo/blast/core/blast_program.h>
#include <algo/blast/core/blast_query_info.h>
#include <algo/blast/core/blast_options.h>
#include <algo/blast/core/blast_parameters.h>
#include <algo/blast/core/blast_stat.h>
#include <algo/blast/core/gapinfo.h>
#include <algo/blast/core/blast_seqsrc.h>
#include <algo/blast/core/pattern.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Keeps prelim_hitlist_size and HitSavingOptions
    together, mostly for use by hspstream. */
typedef struct SBlastHitsParameters {
   Int4 prelim_hitlist_size; /**< number of hits saved during preliminary 
                           part of search. */
   Int4 hsp_num_max; /**< number of HSPs to save per db sequence. */
} SBlastHitsParameters; 

/** Sets up small structures used by blast_hit.c for saving HSPs.
 * @param hit_options field hitlist_size and hsp_num_max needed, a pointer to 
 *      this structure will be stored on resulting structure.[in]
 * @param ext_options field compositionBasedStats needed here. [in]
 * @param scoring_options gapped_calculation needed here. [in]
 * @param retval the allocated SBlastHitsParameters*
 * @return zero on success, 1 on NULL parameter, 2 if calloc fails.
 */
NCBI_XBLAST_EXPORT
Int2 SBlastHitsParametersNew(const BlastHitSavingOptions* hit_options,
                             const BlastExtensionOptions* ext_options,
                             const BlastScoringOptions* scoring_options,
                             SBlastHitsParameters* *retval);

/** Make a deep copy of the SBlastHitsParameters structure passed in
 * @param hit_params source hit parameters structure [in]
 * @return NULL if out of memory, otherwise deep copy of first argument
 */
NCBI_XBLAST_EXPORT
SBlastHitsParameters* 
SBlastHitsParametersDup(const SBlastHitsParameters* hit_params);

/** Deallocated SBlastHitsParameters.
 * @param param object to be freed.
 * @return NULL pointer.
 */
NCBI_XBLAST_EXPORT
SBlastHitsParameters* SBlastHitsParametersFree(SBlastHitsParameters* param);
                   



/** One sequence segment within an HSP */
typedef struct BlastSeg {
   Int2 frame;  /**< Translation frame */
   Int4 offset; /**< Start of hsp */
   Int4 end;    /**< End of hsp */
   Int4 gapped_start;/**< Where the gapped extension started. */
} BlastSeg;

/** In PHI BLAST: information about pattern match in a given HSP. */
typedef struct SPHIHspInfo {
    Int4 index; /**< Index of query pattern occurrence for this HSP. */ 
    Int4 length; /**< Length of this pattern occurrence in subject. */
} SPHIHspInfo;

/** Structure holding all information about an HSP */
typedef struct BlastHSP {
   Int4 score;           /**< This HSP's raw score */
   Int4 num_ident;       /**< Number of identical base pairs in this HSP */
   double bit_score;     /**< Bit score, calculated from score */
   double evalue;        /**< This HSP's e-value */
   BlastSeg query;       /**< Query sequence info. */
   BlastSeg subject;     /**< Subject sequence info. */
   Int4     context;     /**< Context number of query */
   GapEditScript* gap_info;/**< ALL gapped alignment is here */
   Int4 num;             /**< How many HSP's are linked together for sum 
                              statistics evaluation? If unset (0), this HSP is
                              not part of a linked set, i.e. value 0 is treated
                              the same way as 1. */
   Int2		comp_adjustment_method;  /**< which mode of composition
                                              adjustment was used; relevant
                                              only for blastp and tblastn */
   SPHIHspInfo* pat_info; /**< In PHI BLAST, information about this pattern
                                 match. */
   Int4 num_positives;
} BlastHSP;

/** The structure to hold all HSPs for a given sequence after the gapped 
 *  alignment.
 */
typedef struct BlastHSPList {
   Int4 oid;/**< The ordinal id of the subject sequence this HSP list is for */
   Int4 query_index; /**< Index of the query which this HSPList corresponds to.
                        Set to 0 if not applicable */
   BlastHSP** hsp_array; /**< Array of pointers to individual HSPs */
   Int4 hspcnt; /**< Number of HSPs saved */
   Int4 allocated; /**< The allocated size of the hsp_array */
   Int4 hsp_max; /**< The maximal number of HSPs allowed to be saved */
   Boolean do_not_reallocate; /**< Is reallocation of the hsp_array allowed? */
   double best_evalue; /**< Smallest e-value for HSPs in this list. Filled after 
                          e-values are calculated. Necessary because HSPs are
                          sorted by score, but highest scoring HSP may not have
                          the lowest e-value if sum statistics is used. */
} BlastHSPList;

/** The structure to contain all BLAST results for one query sequence */
typedef struct BlastHitList {
   Int4 hsplist_count; /**< Filled size of the HSP lists array */
   Int4 hsplist_max; /**< Maximal allowed size of the HSP lists array */
   double worst_evalue; /**< Highest of the best e-values among the HSP 
                           lists */
   Int4 low_score; /**< The lowest of the best scores among the HSP lists */
   Boolean heapified; /**< Is this hit list already heapified? */
   BlastHSPList** hsplist_array; /**< Array of HSP lists for individual
                                          database hits */
   Int4 hsplist_current; /**< Number of allocated HSP list arrays. */
} BlastHitList;

/** The structure to contain all BLAST results, for multiple queries */
typedef struct BlastHSPResults {
   Int4 num_queries; /**< Number of query sequences */
   BlastHitList** hitlist_array; /**< Array of results for individual
                                          query sequences */
} BlastHSPResults;


/** By how much should the chunks of a subject sequence overlap if it is 
    too long and has to be split */
#define DBSEQ_CHUNK_OVERLAP 100

/********************************************************************************

The following section has four sets of functions (or "APIs"), manipulating with
the following structures:
1. BlastHSP, which is the basic unit to record one alignment.  
2. BlastHSPList, which is a list of BlastHSP's for one database sequence. 
3. BlastHitList, which contains all HSPList's for a given query. 
4. BlastHSPResults, which is a list of BlastHitList's for multiple queries.

 The naming conventions for the functions are the following:

1.) All routines start with "Blast_"

2.) After "Blast_" comes the structure being manipulated, that should be either 
    HSP (all capitals all the time!), HSPList (exactly this capitalization), 
    HitList (capital H and L, all others lower-case), or HSPResults.

3.) finally the task being done, e.g., "Free", "New", "Init".

********************************************************************************/
/********************************************************************************
          HSP API
********************************************************************************/

/** Deallocate memory for an HSP structure */
NCBI_XBLAST_EXPORT
BlastHSP* Blast_HSPFree(BlastHSP* hsp);

/** Allocate and zeros out memory for an HSP structure */
NCBI_XBLAST_EXPORT
BlastHSP* Blast_HSPNew(void);

/** Allocates BlastHSP and inits with information from input.
 * structure.
 * @param query_start Start of query alignment [in]
 * @param query_end End of query alignment [in]
 * @param subject_start Start of subject alignment [in]
 * @param subject_end End of subject alignment [in]
 * @param query_gapped_start Where gapped alignment started on query [in]
 * @param subject_gapped_start Where gapped alignment started on subject [in]
 * @param query_context The index of the query containing this HSP [in]
 * @param query_frame Query frame: -3..3 for translated sequence, 
 *        1 or -1 for blastn, 0 for blastp [in]
 * @param subject_frame Subject frame: -3..3 for translated sequence, 
 *        1 for blastn, 0 for blastp [in]
 * @param score score of alignment [in]
 * @param gap_edit Will be transferred to HSP and nulled out 
 *    if a traceback was not calculated may be NULL [in] [out]
 * @param ret_hsp allocated and filled in BlastHSP [out]
 */
NCBI_XBLAST_EXPORT
Int2
Blast_HSPInit(Int4 query_start, Int4 query_end, 
              Int4 subject_start, Int4 subject_end, 
              Int4 query_gapped_start, Int4 subject_gapped_start, 
              Int4 query_context, Int2 query_frame, Int2 subject_frame,
              Int4 score, GapEditScript* *gap_edit, BlastHSP** ret_hsp);

/** Reevaluate the HSP's score and percent identity after taking
 * into account the ambiguity information. Used only for blastn after a greedy
 * gapped extension with traceback. This function can remove part of the 
 * alignment at either end, if its score becomes negative after reevaluation.
 * Traceback is also adjusted in that case.
 * @param hsp The HSP structure [in] [out]
 * @param query_start Pointer to the start of the query sequence [in]
 * @param query_length Length of the query sequence [in]
 * @param subject_start Pointer to the start of the subject sequence [in]
 * @param subject_length Length of the subject sequence [in]
 * @param hit_params Hit saving parameters containing score cut-off [in]
 * @param score_params Scoring parameters [in]
 * @param sbp Score block with Karlin-Altschul parameters [in]
 * @return Should this HSP be deleted after the score reevaluation?
 */
NCBI_XBLAST_EXPORT
Boolean 
Blast_HSPReevaluateWithAmbiguitiesGapped(BlastHSP* hsp, 
   const Uint1* query_start, const Int4 query_length, 
   const Uint1* subject_start, const Int4 subject_length,
   const BlastHitSavingParameters* hit_params, 
   const BlastScoringParameters* score_params, const BlastScoreBlk* sbp);

/** Reevaluate the HSP's score and percent identity after taking into
 * account the ambiguity information. Used for ungapped searches with 
 * nucleotide database (blastn, tblastn, tblastx).
 * @param hsp The HSP structure [in] [out]
 * @param query_start Pointer to the start of the query sequence [in]
 * @param subject_start Pointer to the start of the subject sequence [in]
 * @param word_params Initial word parameters with ungapped cutoff score [in]
 * @param sbp Score block with Karlin-Altschul parameters [in]
 * @param translated Are sequences protein (with a translated subject)? [in]
 * @return Should this HSP be deleted after the score reevaluation?
 */
NCBI_XBLAST_EXPORT
Boolean 
Blast_HSPReevaluateWithAmbiguitiesUngapped(BlastHSP* hsp, 
   const Uint1* query_start, const Uint1* subject_start,
   const BlastInitialWordParameters* word_params, 
   BlastScoreBlk* sbp, Boolean translated);

/** Calculate number of identities in an HSP and set the BlastHSP::num_ident
 * field (unconditionally)
 * @param query The query sequence [in]
 * @param subject The uncompressed subject sequence [in]
 * @param hsp All information about the HSP, the output of this function will
 * be stored in its num_ident field [in|out]
 * @param score_options Scoring options [in]
 * @param align_length_ptr The alignment length, including gaps (optional) [out]
 * @return 0 on success, -1 on invalid parameters or error
 */
NCBI_XBLAST_EXPORT
Int2
Blast_HSPGetNumIdentities(const Uint1* query, 
                          const Uint1* subject, 
                          BlastHSP* hsp, 
                          const BlastScoringOptions* score_options,
                          Int4* align_length_ptr);

/** Calculate number of identities and positives in an HSP and set the
 *  BlastHSP::num_ident  and BlastHSP::num_positives fields
 * @param query The query sequence [in]
 * @param subject The uncompressed subject sequence [in]
 * @param hsp All information about the HSP, the output of this function will
 * be stored in its num_ident field [in|out]
 * @param score_options Scoring options [in]
 * @param align_length_ptr The alignment length, including gaps (optional) [out]
 * @param sbp Score blk containing the matrix for counting positives [in]
 * @return 0 on success, -1 on invalid parameters or error
 */
NCBI_XBLAST_EXPORT
Int2
Blast_HSPGetNumIdentitiesAndPositives(const Uint1* query,
                          	  	  	  const Uint1* subject,
                          			  BlastHSP* hsp,
                          			  const BlastScoringOptions* score_options,
                          			  Int4* align_length_ptr,
                          			  const BlastScoreBlk* sbp);

/** Determines whether this HSP should be kept or
 * deleted.
 * @param hsp An HSP structure [in] [out]
 * @param hit_options Hit saving options containing percent identity and
 *                    HSP length thresholds.
 * @param align_length alignment length including gaps
 * @return FALSE if HSP passes the test, TRUE if it should be deleted.
 */
NCBI_XBLAST_EXPORT
Boolean
Blast_HSPTest(BlastHSP* hsp,
 	 	 	  const BlastHitSavingOptions* hit_options,
 	 	 	  Int4 align_length);

/** Calculates number of identities and alignment lengths of an HSP via
 * Blast_HSPGetNumIdentities and determines whether this HSP should be kept or
 * deleted. 
 * @param program_number Type of BLAST program [in]
 * @param hsp An HSP structure [in] [out]
 * @param query Query sequence [in]
 * @param subject Subject sequence [in]
 * @param score_options Scoring options, needed to distinguish the 
 *                      out-of-frame case. [in]
 * @param hit_options Hit saving options containing percent identity and
 *                    HSP length thresholds.
 * @return FALSE if HSP passes the test, TRUE if it should be deleted.
 */ 
NCBI_XBLAST_EXPORT
Boolean
Blast_HSPTestIdentityAndLength(EBlastProgramType program_number, 
                               BlastHSP* hsp, const Uint1* query, const Uint1* subject, 
                               const BlastScoringOptions* score_options,
                               const BlastHitSavingOptions* hit_options);

/** Calculate query coverage percentage of an hsp
 *  @param hsp An HSP structure [in]
 *  @param query_length	Length of query [in]
 *  @return percentage query coverage of the input hsp
 */
NCBI_XBLAST_EXPORT
double
Blast_HSPGetQueryCoverage(const BlastHSP* hsp, Int4 query_length);

/** Calculate query coverage percentage of an hsp
 *  @param hsp An HSP structure [in]
 *  @param min_query_coverage_pct Min query coverage pct for saving the hsp[in]
 *  @param query_length	Length of query [in]
 *  @return true if hsp's query coverage pct  < min_query_coverage_pct (delete hsp)
 */
NCBI_XBLAST_EXPORT
Boolean Blast_HSPQueryCoverageTest(BlastHSP* hsp,
        						   double min_query_coverage_pct,
        						   Int4 query_length);


NCBI_XBLAST_EXPORT
Int2 Blast_TrimHSPListByMaxHsps(BlastHSPList* hsp_list,
                                const BlastHitSavingOptions* hit_options);

/** Calculated the number of HSPs that should be saved.
 * @param gapped_calculation ungapped if false [in]
 * @param options HitSavingoptions object [in]
 * @return number of HSPs to save. 
 */
NCBI_XBLAST_EXPORT
Int4
BlastHspNumMax(Boolean gapped_calculation, const BlastHitSavingOptions* options);

/** Calculate length of an HSP as length in query plus length of gaps in 
 * query. If gap information is unavailable, return maximum between length in
 * query and in subject.
 * @param hsp An HSP structure [in]
 * @param length Length of this HSP [out]
 * @param gaps Total number of gaps in this HSP [out]
 * @param gap_opens Number of gap openings in this HSP [out] 
 */
NCBI_XBLAST_EXPORT
void Blast_HSPCalcLengthAndGaps(const BlastHSP* hsp, Int4* length,
                                Int4* gaps, Int4* gap_opens);

/** Adjust HSP endpoint offsets according to strand/frame; return values in
 * 1-offset coordinates instead of internal 0-offset.
 * @param program Type of BLAST program [in]
 * @param hsp An HSP structure [in]
 * @param query_length Length of query [in]
 * @param subject_length Length of subject [in]
 * @param q_start Start of alignment in query [out]
 * @param q_end End of alignment in query [out]
 * @param s_start Start of alignment in subject [out]
 * @param s_end End of alignment in subject [out]
 */
NCBI_XBLAST_EXPORT
void 
Blast_HSPGetAdjustedOffsets(EBlastProgramType program, BlastHSP* hsp, 
                            Int4 query_length, Int4 subject_length, 
                            Int4* q_start, Int4* q_end,
                            Int4* s_start, Int4* s_end);

/** Performs the translation and coordinates adjustment, if only part of the 
 * subject sequence is translated for gapped alignment. 
 * @param subject_blk Subject sequence structure [in]
 * @param hsp The HSP information [in] [out]
 * @param is_ooframe Return a mixed-frame sequence if TRUE [in]
 * @param gen_code_string Database genetic code [in]
 * @param translation_buffer_ptr Pointer to buffer holding the translation [out]
 * @param subject_ptr Pointer to sequence to be passed to the gapped 
 *                    alignment [out]
 * @param subject_length_ptr Length of the translated sequence [out]
 * @param start_shift_ptr How far is the partial sequence shifted w.r.t. the 
 *                        full sequence. [out]
 */
NCBI_XBLAST_EXPORT
Int2
Blast_HSPGetPartialSubjectTranslation(BLAST_SequenceBlk* subject_blk, 
   BlastHSP* hsp, Boolean is_ooframe, const Uint1* gen_code_string, 
   Uint1** translation_buffer_ptr, Uint1** subject_ptr, 
   Int4* subject_length_ptr, Int4* start_shift_ptr);

/** Adjusts offsets if partial sequence was used for extension.
 * @param hsp The hit to work on [in][out]
 * @param start_shift amount of database sequence not used for extension. [in]
*/
NCBI_XBLAST_EXPORT
void
Blast_HSPAdjustSubjectOffset(BlastHSP* hsp, Int4 start_shift);


/** Returns a buffer with a protein translated from nucleotide.
 * @param target_t SBlastTargetTranslation* with information about translation [in]
 * @param hsp The hit to work on [in]
 * @param translated_length length of the protein sequence [in]
*/
NCBI_XBLAST_EXPORT
const Uint1*
Blast_HSPGetTargetTranslation(SBlastTargetTranslation* target_t, const BlastHSP* hsp, Int4* translated_length);

/********************************************************************************
          HSPList API
********************************************************************************/

/** Deallocate memory for an HSP list structure 
 *  as well as all it's components.
 * @param hsp_list the BlastHSPList to be freed [in]. 
*/
NCBI_XBLAST_EXPORT
BlastHSPList* Blast_HSPListFree(BlastHSPList* hsp_list);

/** Creates HSP list structure with a default size HSP array 
 * @param hsp_max the maximum number of HSP's that can ever be
 *    saved at once [in].
*/
NCBI_XBLAST_EXPORT
BlastHSPList* Blast_HSPListNew(Int4 hsp_max);

/** Returns true if the BlastHSPList contains no HSPs
 * @param hsp_list list of HSPs to examine [in]
 */
NCBI_XBLAST_EXPORT
Boolean
Blast_HSPList_IsEmpty(const BlastHSPList* hsp_list);

/** Returns a duplicate (deep copy) of the given hsp list. */
NCBI_XBLAST_EXPORT
BlastHSPList* BlastHSPListDup(const BlastHSPList* hsp_list);

/** Swaps the two HSP lists via structure assignment */
NCBI_XBLAST_EXPORT
void Blast_HSPListSwap(BlastHSPList* list1, BlastHSPList* list2);

/** Saves HSP information into a BlastHSPList structure
 * @param hsp_list Structure holding all HSPs with full gapped alignment 
 *        information [in] [out]
 * @param hsp The new HSP to be inserted into the HSPList [in]
*/
NCBI_XBLAST_EXPORT
Int2
Blast_HSPListSaveHSP(BlastHSPList* hsp_list, BlastHSP* hsp);

/** Calculate the expected values for all HSPs in a hit list, without using 
 * the sum statistics. In case of multiple queries, the offsets are assumed 
 * to be already adjusted to individual query coordinates, and the contexts 
 * are set for each HSP.
 * @param program_number Type of BLAST program [in]
 * @param query_info Auxiliary query information - needed only for effective
 *                   search space calculation if it is not provided [in]
 * @param subject_length Subject length - needed for Spouge's new FSC [in]
 * @param hsp_list List of HSPs for one subject sequence [in] [out]
 * @param gapped_calculation Is this for a gapped or ungapped search? [in]
 * @param RPS_prelim Is this for a RPS preliminary search? [in]
 * @param sbp Structure containing statistical information [in]
 * @param gap_decay_rate Adjustment parameter to compensate for the effects of
 * performing multiple tests when linking HSPs. No adjustment is made if 0. [in]
 * @param scaling_factor Scaling factor by which Lambda should be divided. Used in
 *                       RPS BLAST only; should be set to 1.0 in other cases. [in]
 *                       
 */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListGetEvalues(EBlastProgramType program_number,
                             const BlastQueryInfo* query_info,
                             Int4 subject_length,
                             BlastHSPList* hsp_list,
                             Boolean gapped_calculation, 
                             Boolean RPS_prelim,
                             const BlastScoreBlk* sbp, double gap_decay_rate,
                             double scaling_factor);

/** Calculate e-values for a PHI BLAST HSP list.
 * @param hsp_list HSP list found by PHI BLAST [in] [out]
 * @param sbp Scoring block with statistical parameters [in]
 * @param query_info Structure containing information about pattern counts [in]
 * @param pattern_blk Structure containing information about pattern hits in db [in]
 */
NCBI_XBLAST_EXPORT
void Blast_HSPListPHIGetEvalues(BlastHSPList* hsp_list, BlastScoreBlk* sbp, 
                                const BlastQueryInfo* query_info,
                                const SPHIPatternSearchBlk* pattern_blk);

/** Calculate bit scores from raw scores in an HSP list.
 * @param hsp_list List of HSPs [in] [out]
 * @param gapped_calculation Is this a gapped search? [in]
 * @param sbp Scoring block with statistical parameters [in]
 */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListGetBitScores(BlastHSPList* hsp_list, 
                               Boolean gapped_calculation, 
                               const BlastScoreBlk* sbp);

/** Calculate bit scores from raw scores in an HSP list for a PHI BLAST search.
 * @param hsp_list List of HSPs [in] [out]
 * @param sbp Scoring block with statistical parameters [in]
 */
NCBI_XBLAST_EXPORT
void Blast_HSPListPHIGetBitScores(BlastHSPList* hsp_list, BlastScoreBlk* sbp);
    
/** Discard the HSPs above the e-value threshold from the HSP list 
 * @param hsp_list List of HSPs for one subject sequence [in] [out]
 * @param hit_options Options block containing the e-value cut-off [in]
*/
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListReapByEvalue(BlastHSPList* hsp_list, 
                               const BlastHitSavingOptions* hit_options);

/** Discard the HSPs above the raw threshold from the HSP list 
 * @param hsp_list List of HSPs for one subject sequence [in] [out]
 * @param hit_options Options block containing the e-value cut-off [in]
 * -RMH-
 */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListReapByRawScore(BlastHSPList* hsp_list,
                               const BlastHitSavingOptions* hit_options);

/** Discard the HSPs below the min query coverage pct from the HSP list
 * @param hsp_list List of HSPs for one subject sequence [in] [out]
 * @param hit_options Options block containing the min query coverage pct [in]
 * @param query_info Structure containing information about the queries  [in]
 * @param program_number Type of BLAST program.
*/
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListReapByQueryCoverage(BlastHSPList* hsp_list,
                                      const BlastHitSavingOptions* hit_options,
                                      const BlastQueryInfo* query_info,
                                      EBlastProgramType program_number);

/** Cleans out the NULLed out HSP's from the HSP array that
 * is part of the BlastHSPList.
 * @param hsp_list Contains array of pointers to HSP structures [in]
 * @return status of function call.
*/
NCBI_XBLAST_EXPORT
Int2
Blast_HSPListPurgeNullHSPs(BlastHSPList* hsp_list);

/** Check for an overlap of two different alignments and remove redundant HSPs.
 * A sufficient overlap is when two alignments have the same start or end values
 * If an overlap is found the HSP with the lowest score is removed, if both scores
 * are the same then the first is removed.
 * @param program Type of BLAST program. For some programs (PHI BLAST), the
 *                purge should not be performed. [in]
 * @param hsp_list Contains array of pointers to HSPs to purge [in]
 * @param purge Should the hsp be purged? [in]
 * @return The number of valid alignments remaining. 
*/
NCBI_XBLAST_EXPORT
Int4
Blast_HSPListPurgeHSPsWithCommonEndpoints(EBlastProgramType program, 
                                          BlastHSPList* hsp_list,
                                          Boolean purge);

/** Reevaluate all ungapped HSPs in an HSP list.  
 * This is only done for an ungapped search, or if traceback is 
 * already available.
 * Subject sequence is uncompressed and saved here (for nucleotide sequences). 
 * The number of identities is calculated for each HSP along the way,
 * hence this function is called for all programs. 
 * @param program Type of BLAST program [in]
 * @param hsp_list The list of HSPs for one subject sequence [in] [out]
 * @param query_blk The query sequence [in]
 * @param subject_blk The subject sequence [in] [out]
 * @param word_params Initial word parameters, containing ungapped cutoff 
 *                    score [in]
 * @param hit_params Hit saving parameters, including cutoff score [in]
 * @param query_info Auxiliary query information [in]
 * @param sbp The statistical information [in]
 * @param score_params Parameters related to scoring [in]
 * @param seq_src The BLAST database structure (for retrieving uncompressed
 *             sequence) [in]
 * @param gen_code_string Genetic code string in case of a translated 
 *                        database search. [in]
 */
NCBI_XBLAST_EXPORT
Int2 
Blast_HSPListReevaluateUngapped(EBlastProgramType program, 
   BlastHSPList* hsp_list, BLAST_SequenceBlk* query_blk, 
   BLAST_SequenceBlk* subject_blk, 
   const BlastInitialWordParameters* word_params,
   const BlastHitSavingParameters* hit_params, const BlastQueryInfo* query_info, 
   BlastScoreBlk* sbp, const BlastScoringParameters* score_params, 
   const BlastSeqSrc* seq_src, const Uint1* gen_code_string);

/** Append one HSP list to the other. Discard lower scoring HSPs if there is
 * not enough space to keep all.
 * @param old_hsp_list_ptr list of HSPs, will be NULLed out on return [in|out]
 * @param combined_hsp_list_ptr Pointer to the combined list of HSPs, possibly
 *                              containing previously saved HSPs [in] [out]
 * @param hsp_num_max Maximal allowed number of HSPs to save (unlimited if INT4_MAX) [in]
 * @return Status: 0 on success, -1 on failure.
 */ 
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListAppend(BlastHSPList** old_hsp_list_ptr,
        BlastHSPList** combined_hsp_list_ptr, Int4 hsp_num_max);

/** Merge an HSP list from a chunk of the subject sequence into a previously
 * computed HSP list.
 * @param hsp_list Contains HSPs from the new chunk [in]
 * @param combined_hsp_list_ptr Contains HSPs from previous chunks [in] [out]
 * @param hsp_num_max Maximal allowed number of HSPs to save (unlimited if INT4_MAX) [in]
 * @param split_points Offset The sequence offset (query or subject) that is 
 *             the boundary between HSPs in combined_hsp_list and hsp_list. [in]
 * @param contexts_per_query If positive, the number of query contexts
 *                    that hits can contain. If negative, the (one) split
 *                    point occurs on the subject sequence [in]
 * @param chunk_overlap_size The length of the overlap region between the
 *                    sequence region containing hsp_list and that
 *                    containing combined_hsp_list [in]
 * @param allow_gap Allow merging HSPs at different diagonals [in]
 * @return 0 if HSP lists have been merged successfully, -1 otherwise.
 */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPListsMerge(BlastHSPList** hsp_list, 
                   BlastHSPList** combined_hsp_list_ptr, 
                   Int4 hsp_num_max, Int4* split_points, 
                   Int4 contexts_per_query,
                   Int4 chunk_overlap_size,
                   Boolean allow_gap);
                   
/** Adjust subject offsets in an HSP list if only part of the subject sequence
 * was searched. Used when long subject sequence is split into more manageable
 * chunks.
 * @param hsp_list List of HSPs from a chunk of a subject sequence [in]
 * @param offset Offset where the chunk starts [in]
 */
NCBI_XBLAST_EXPORT
void Blast_HSPListAdjustOffsets(BlastHSPList* hsp_list, Int4 offset);

/** For nucleotide BLAST, if the match reward score is equal to 2, 
 * random alignments are dominated by runs of exact matches, which all have even
 * scores. This makes it impossible to estimate statistical parameters correctly
 * for odd scores. Hence the raw score formula is adjusted - all scores are
 * rounded down to the nearest even value in order to provide a conservative estimate.
 * @param hsp_list HSP list structure to adjust scores for. [in] [out]
 * @param gapped_calculation not an ungapped alignment [in]
 * @param sbp used for round_down Boolean
 */
NCBI_XBLAST_EXPORT
void Blast_HSPListAdjustOddBlastnScores(BlastHSPList* hsp_list, 
                                        Boolean gapped_calculation, 
                                        const BlastScoreBlk* sbp);

/** Check if HSP list is sorted by score.
 * @param hsp_list The list to check [in]
 * @return TRUE if sorted, FALSE if not.
 */
NCBI_XBLAST_EXPORT
Boolean Blast_HSPListIsSortedByScore(const BlastHSPList* hsp_list);

/** Sort the HSPs in an HSP list by score. This type of sorting is done before
 * the e-values are calcaulted, and also at the beginning of the traceback stage, 
 * where it is needed to eliminate the effects of wrong score order because of 
 * application of sum statistics. 
 * Checks if the HSP array is already sorted before proceeding with quicksort.
 * @param hsp_list Structure containing array of HSPs to be sorted. [in] [out]
 */
NCBI_XBLAST_EXPORT
void Blast_HSPListSortByScore(BlastHSPList* hsp_list);

/** Sort the HSPs in an HSP list by e-value, with scores and other criteria
 * used to resolve ties. Checks if the HSP array is already sorted before 
 * proceeding with quicksort.
 * @param hsp_list Structure containing array of HSPs to be sorted. [in] [out]
 */
NCBI_XBLAST_EXPORT
void Blast_HSPListSortByEvalue(BlastHSPList* hsp_list);

/********************************************************************************
          HitList API.
********************************************************************************/

/** Allocate memory for a hit list of a given size.
 * @param hitlist_size Size of the hit list (number of HSP lists) [in]
 */
NCBI_XBLAST_EXPORT
BlastHitList* Blast_HitListNew(Int4 hitlist_size);

/** Deallocate memory for the hit list */
NCBI_XBLAST_EXPORT
BlastHitList* Blast_HitListFree(BlastHitList* hitlist);

/** Deallocate memory for every HSP list on BlastHitList,
 *  as well as all their components.
 * @param hitlist contains the BlastHSPList array to be freed [in/out]. 
*/
NCBI_XBLAST_EXPORT
Int2 Blast_HitListHSPListsFree(BlastHitList* hitlist);

/** Insert a new HSP list into the hit list.
 * Before capacity of the hit list is reached, just add to the end;
 * After that, store in a heap, to ensure efficient insertion and deletion.
 * The heap order is reverse, with worst e-value on top, for convenience
 * of deletion.
 * @param hit_list Contains all HSP lists saved so far [in] [out]
 * @param hsp_list A new HSP list to be inserted into the hit list [in]
*/
NCBI_XBLAST_EXPORT
Int2 Blast_HitListUpdate(BlastHitList* hit_list, BlastHSPList* hsp_list);

/** Combine two hitlists; both HitLists must contain HSPs that
 * represent alignments to the same query sequence
 * @param old_hit_list_ptr Pointer to original HitList, will be NULLed 
 *                          out on return [in|out]
 * @param combined_hit_list_ptr Pointer to the combined HitList [in|out]
 t* @param contexts_per_query The number of different contexts that can
 *             occur in hits from old_hit_list and combined_hit_list [in]
 * @param split_offsets the query offset that marks the boundary between
 *             combined_hit_list and old_hit_list. HSPs in old_hit_list
 *             that hit to context i are assumed to lie to the right 
 *             of split_offsets[i] [in]
 * @param chunk_overlap_size The length of the overlap region between the
 *                    sequence region containing hit_list and that
 *                    containing combined_hit_list [in]
 * @param allow_gap Allow merging HSPs at different diagonals [in]
*/
NCBI_XBLAST_EXPORT
Int2 Blast_HitListMerge(BlastHitList** old_hit_list_ptr,
                        BlastHitList** combined_hit_list_ptr,
                        Int4 contexts_per_query, Int4 *split_offsets,
                        Int4 chunk_overlap_size, Boolean allow_gap);

/** Purges a BlastHitList of NULL HSP lists.
 * @param hit_list BLAST hit list to purge. [in] [out]
 */
NCBI_XBLAST_EXPORT
Int2 
Blast_HitListPurgeNullHSPLists(BlastHitList* hit_list);

/** Sort BlastHitLIst bon evalue
 *  @param hit_lsit BLAST hit list to be sorted [in] [out]
 */
NCBI_XBLAST_EXPORT
Int2
Blast_HitListSortByEvalue(BlastHitList* hit_list);

/********************************************************************************
          HSPResults API.
********************************************************************************/

/** Initialize the results structure.
 * @param num_queries Number of query sequences to allocate results structure
 *                    for [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPResults* Blast_HSPResultsNew(Int4 num_queries);

/** Deallocate memory for BLAST results */
NCBI_XBLAST_EXPORT
BlastHSPResults* Blast_HSPResultsFree(BlastHSPResults* results);

/** Sort each hit list in the BLAST results by best e-value */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPResultsSortByEvalue(BlastHSPResults* results);
/** Sort each hit list in the BLAST results by best e-value, in reverse
    order. */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPResultsReverseSort(BlastHSPResults* results);

/** Reverse order of HSP lists in each hit list in the BLAST results. 
 * This allows to return HSP lists from the end of the arrays when reading
 * from a collector HSP stream.
 */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPResultsReverseOrder(BlastHSPResults* results);

/** Blast_HSPResultsInsertHSPList
 * Insert an HSP list to the appropriate place in the results structure.
 * All HSPs in this list must be from the same query and same subject; the oid
 * and query_index fields must be set in the BlastHSPList input structure.
 * @param results The structure holding results for all queries [in] [out]
 * @param hsp_list The results for one query-subject sequence pair. [in]
 * @param hitlist_size Maximal allowed hit list size. [in]
 */
NCBI_XBLAST_EXPORT
Int2 Blast_HSPResultsInsertHSPList(BlastHSPResults* results, 
        BlastHSPList* hsp_list, Int4 hitlist_size);

/* Forward declaration */
struct BlastHSPStream;

/** Move all of the hits within an HSPStream into a BlastHSPResults
 * structure.
 * @param hsp_stream The HSPStream [in][out]
 * @param num_queries Number of queries in the search [in]
 * @param hit_param Hit parameters [in]
 * @return The generated collection of HSP results
 */
NCBI_XBLAST_EXPORT
BlastHSPResults*
Blast_HSPResultsFromHSPStream(struct BlastHSPStream* hsp_stream, 
                              size_t num_queries, 
                              SBlastHitsParameters* hit_param);

/** As Blast_HSPResultsFromHSPStream, except the total number of
 * HSPs kept for each query does not exceed an explicit limit.
 * The database sequences with the smallest number of hits are
 * saved first, and hits are removed from query i if the average
 * number of hits saved threatens to exceed (max_num_hsps / (number
 * of DB sequences with hits to query i))
 * @param hsp_stream The HSPStream [in][out]
 * @param num_queries Number of queries in the search [in]
 * @param hit_param Hit parameters [in]
 * @param max_num_hsps The limit on the number of HSPs to be
 *                     kept for each query sequence [in]
 * @param removed_hsps Set to TRUE if any hits were removed [out]
 * @return The generated collection of HSP results
 */
BlastHSPResults*
Blast_HSPResultsFromHSPStreamWithLimit(struct BlastHSPStream* hsp_stream, 
                                   Uint4 num_queries, 
                                   SBlastHitsParameters* hit_param,
                                   Uint4 max_num_hsps,
                                   Boolean* removed_hsps);

BlastHSPResults*
/** As Blast_HSPResultsFromHSPStreamWithLimit, except accept and return 
 * array of Boolen flags specifying which query exceeded HSP limits.
 */
Blast_HSPResultsFromHSPStreamWithLimitEx(struct BlastHSPStream* hsp_stream, 
                                   Uint4 num_queries, 
                                   SBlastHitsParameters* hit_param,
                                   Uint4 max_num_hsps,
                                   Boolean* removed_hsps);
/** Splits the BlastHSPResults structure for a PHI BLAST search into an array of
 * BlastHSPResults structures, corresponding to different pattern occurrences in
 * query. All HSPs are copied, so it is safe to free the returned 
 * BlastHSPResults structures independently of the input results structure.
 * @param results All results from a PHI BLAST search, with HSPs for 
 *                different query pattern occurrences mixed together. [in]
 * @param pattern_info Information about pattern occurrences in query. [in]
 * @return Array of pointers to BlastHSPResults structures, corresponding to 
 *         different pattern occurrences.
 */
NCBI_XBLAST_EXPORT
BlastHSPResults** 
PHIBlast_HSPResultsSplit(const BlastHSPResults* results, 
                         const SPHIQueryInfo* pattern_info);


/** Count the number of occurrences of pattern in sequence, which
 * do not overlap by more than half the pattern match length. 
 * @param query_info Query information structure, containing pattern info. [in]
 */
NCBI_XBLAST_EXPORT
Int4
PhiBlastGetEffectiveNumberOfPatterns(const BlastQueryInfo *query_info);

/** Apply Cross_match like masklevel to HSP list.  -RMH-
 */
Int2 Blast_HSPResultsApplyMasklevel(BlastHSPResults *results,
                                    const BlastQueryInfo *query_info,
                                    Int4 masklevel, Int4 query_length);

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLAST_HITS__H */
