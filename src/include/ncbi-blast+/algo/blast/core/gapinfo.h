/* $Id: gapinfo.h 149108 2009-01-07 16:27:23Z ivanov $
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

/** @file gapinfo.h
 * Definitions of structures used for saving traceback information.
 */

#ifndef ALGO_BLAST_CORE__GAPINFO__H
#define ALGO_BLAST_CORE__GAPINFO__H

#include <algo/blast/core/ncbi_std.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Operation types within the edit script*/
typedef enum EGapAlignOpType { 
   eGapAlignDel = 0, /**< Deletion: a gap in query */
   eGapAlignDel2 = 1,/**< Frame shift deletion of two nucleotides */
   eGapAlignDel1 = 2,/**< Frame shift deletion of one nucleotide */
   eGapAlignSub = 3, /**< Substitution */
   eGapAlignIns1 = 4,/**< Frame shift insertion of one nucleotide */
   eGapAlignIns2 = 5,/**< Frame shift insertion of two nucleotides */
   eGapAlignIns = 6, /**< Insertion: a gap in subject */
   eGapAlignDecline = 7, /**< Non-aligned region */
   eGapAlignInvalid = 8 /**< Invalid operation */
} EGapAlignOpType;

/** Edit script: linked list of correspondencies between two sequences */
typedef struct GapEditScript {
   EGapAlignOpType* op_type;    /**< Array of type of operation */
   Int4* num;                   /**< Array of number of operations */
   Int4 size;                   /**< Size of above arrays. */
} GapEditScript;

/** A version of GapEditScript used to store initial results
    from the gapped alignment routines */
typedef struct GapPrelimEditScript {
   EGapAlignOpType op_type;    /**< Type of operation */
   Int4 num;                   /**< Number of operations */
} GapPrelimEditScript;

/** Preliminary version of GapEditBlock, used directly by the low-
 * level dynamic programming routines 
 */
typedef struct GapPrelimEditBlock {
    GapPrelimEditScript *edit_ops;  /**< array of edit operations */
    Int4 num_ops_allocated;        /**< size of allocated array */
    Int4 num_ops;                  /**< number of edit ops presently in use */
    EGapAlignOpType last_op;        /**< most recent operation added */
} GapPrelimEditBlock;

/** Structure to keep memory for state structure. */ 
typedef struct GapStateArrayStruct {
	Int4 	length,		/**< length of the state_array. */
		used;		/**< how much of length is used. */
	Uint1* state_array;	/**< array to be used. */
	struct GapStateArrayStruct* next; /**< Next link in the list. */
} GapStateArrayStruct;

/** Initialize the edit script structure. 
 *  @param size number of elements to allocate.
 *  @return Pointer to the new edit script
 */
NCBI_XBLAST_EXPORT
GapEditScript* 
GapEditScriptNew (Int4 size); 

/** Free edit script structure. 
 *  @param esp Pointer to the edit script [in]
 *  @return NULL
 */
NCBI_XBLAST_EXPORT
GapEditScript* 
GapEditScriptDelete (GapEditScript* esp);

/** Duplicates the edit script structure. 
 *  @param old object to be duplicated [in]
 *  @return Pointer to the new edit script
 */
NCBI_XBLAST_EXPORT
GapEditScript* 
GapEditScriptDup (const GapEditScript* old);

/** Copies the portion of the GapEditScript specified by start and stop to a new one
 * the new one should already exist.
 *  @param new_esp edit script to copy to [in|out]
 *  @param offset starting element in new one to copy to [in]
 *  @param old_esp edit script to copy from [in]
 *  @param start first element to copy from (zero-offset) [in]
 *  @param stop last element to copy [in]
 *  @return 0 on success
 */
NCBI_XBLAST_EXPORT 
Int2
GapEditScriptPartialCopy(GapEditScript* new_esp, int offset, const GapEditScript* old_esp, int start, int stop);

/** Frees a preliminary edit block structure 
 *  @param edit_block The edit block to free [in]
 *  @return Always NULL
 */
NCBI_XBLAST_EXPORT
GapPrelimEditBlock *
GapPrelimEditBlockFree(GapPrelimEditBlock *edit_block);

/** Allocates a preliminary edit block structure 
 *  @return Pointer to the allocated preliminary edit block
 */
NCBI_XBLAST_EXPORT
GapPrelimEditBlock *
GapPrelimEditBlockNew(void);

/** Add a new operation to a preliminary edit block, possibly combining
 *  it with the last operation if the two operations are identical
 *
 *  @param edit_block The script to update [in/modified]
 *  @param op_type The operation type to add [in]
 *  @param num_ops The number of the specified type of operation to add [in]
 */
NCBI_XBLAST_EXPORT
void
GapPrelimEditBlockAdd(GapPrelimEditBlock *edit_block, 
                 EGapAlignOpType op_type, Int4 num_ops);

/** Append one GapPrelimEditBlock to the end of the other.
 * @param edit_block1 First traceback block [in]
 * @param edit_block2 Second traceback block, to be appended at the end of
 *                    the first.
 */
NCBI_XBLAST_EXPORT
void
GapPrelimEditBlockAppend(GapPrelimEditBlock *edit_block1,
                         GapPrelimEditBlock *edit_block2);


/** Reset a preliminary edit block without freeing it 
 * @param edit_block The preliminary edit block to reset
 */
NCBI_XBLAST_EXPORT
void
GapPrelimEditBlockReset(GapPrelimEditBlock *edit_block);

/** Free the gap state structure. 
 * @param state_struct The state structure to free
 * @return Always NULL
 */
NCBI_XBLAST_EXPORT
GapStateArrayStruct* 
GapStateFree(GapStateArrayStruct* state_struct);

#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__GAPINFO__H */
