/*  $Id: blast_hspfilter.h 161402 2009-05-27 17:35:47Z camacho $
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

/** @file blast_hspfilter.h
 * Declaration of ADT to filter and/or process HSPs in the BLAST engine.
 */

#ifndef ALGO_BLAST_CORE__BLAST_HSPFILTER_H
#define ALGO_BLAST_CORE__BLAST_HSPFILTER_H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_export.h>
#include <algo/blast/core/blast_hits.h>
#include <connect/ncbi_core.h>

#ifdef __cplusplus
extern "C" {
#endif

/**--------------------------------writer-----------------------------*/
/** forwarded declarations */
typedef struct BlastHSPWriter BlastHSPWriter;

/** Function pointer typedef to implement hsp_writer */
typedef BlastHSPWriter* (*BlastHSPWriterNewFn)  (void*, BlastQueryInfo*);
typedef BlastHSPWriter* (*BlastHSPWriterFreeFn) (BlastHSPWriter*);
typedef int (*BlastHSPWriterInitFn) (void*, BlastHSPResults*);
typedef int (*BlastHSPWriterRunFn)  (void*, BlastHSPList*);
typedef int (*BlastHSPWriterFinalFn)(void*, BlastHSPResults*);


/** ADT definition of BlastHSPWriter */
struct BlastHSPWriter {
   void * data;                     /**< data structure */
   BlastHSPWriterInitFn         InitFnPtr;
   BlastHSPWriterRunFn          RunFnPtr;
   BlastHSPWriterFinalFn        FinalFnPtr;
   BlastHSPWriterFreeFn         FreeFnPtr;
} ;

/** A wrap of data structure used to create a writer */
typedef struct BlastHSPWriterInfo {
   void*                        params;
   BlastHSPWriterNewFn          NewFnPtr;
} BlastHSPWriterInfo;

/** A generic function to create writer 
 * @param writer_info structure containing information to create the
 * BlastHSPWriter. This will be free'd in this function and thus will be NULL
 * on function exit [in|out]
 */
NCBI_XBLAST_EXPORT
BlastHSPWriter*
BlastHSPWriterNew(BlastHSPWriterInfo** writer_info, BlastQueryInfo *query_info);

/**--------------------------------pipe------------------------------*/
/** forwarded declarations */
typedef struct BlastHSPPipe BlastHSPPipe;

/** Function pointer typedef to implement hsp_pipe */
typedef BlastHSPPipe* (*BlastHSPPipeNewFn)  (void*, BlastQueryInfo*);
typedef BlastHSPPipe* (*BlastHSPPipeFreeFn) (BlastHSPPipe*);
typedef int (*BlastHSPPipeRunFn)    (void*, BlastHSPResults*);

/** ADT definition of BlastHSPPipe */
struct BlastHSPPipe {
   void * data;                     /**< data structure */
   BlastHSPPipeRunFn            RunFnPtr;
   BlastHSPPipeFreeFn           FreeFnPtr;
   BlastHSPPipe*                next;      /**< the next pipe in chain*/
};

/** A wrap of data structure used to create a pipe */
typedef struct BlastHSPPipeInfo {
   void*                        params;
   BlastHSPPipeNewFn            NewFnPtr;
   struct BlastHSPPipeInfo*     next;      /**< the next pipe inf in chain*/
} BlastHSPPipeInfo;

/** Adds node to the linked list starting a head, which should be NULL when
 * initializing the linked list, on subsequent calls, new nodes will be
 * appended to the list.
 * @param head head of the linked list [in|out]
 * @param node the node to add, ownership is assumed by the linked list [in]
 * @return node most recently added
 * @note there is no explicit free function for a linked list of
 * BlastHSPPipeInfo structures because this is free'd in BlastHSPPipeNew
 */
NCBI_XBLAST_EXPORT
BlastHSPPipeInfo* BlastHSPPipeInfo_Add(BlastHSPPipeInfo** head,
                                       BlastHSPPipeInfo* node);

/** A generic function to create pipe.
 * @param pipe_info linked list of pipe info structures which specifies how to
 * construct the BlastHSPPipe objects. This object's ownership is assumed by
 * the BlastHSPPipe, and thus this field is NULL on function exit [in]
 * @param query_info argument used in construction for all BlastHSPPipe objects
 * [in]
 */
NCBI_XBLAST_EXPORT
BlastHSPPipe*
BlastHSPPipeNew(BlastHSPPipeInfo **pipe_info, BlastQueryInfo *query_info);



/**--------------------------------docs------------------------------*/
/**
 * @page _impl_blast_hspfilter_howto Implementing the BlastHSPFilter interface
 *
 * BlastHSPFilter interface includes BlastHSPWriter and BlastHSPPipe.   
 * The former is used only in preliminary stage to process BlastHSPList coming 
 * directly from scans, whereas the latter can be used in both preliminary and 
 * traceback stages to process the already-collected BlastHSPResults.  
 *
 * Implementations of both types are similar.  
 * As an example, to implement MyWriter, the following must be declared in
 * hspfilter_mywriter.h and implemented in hspfilter_mywriter.c:
 * @code
 * extern "C" {
 * // Introduce data structure to describe MyWriter filtering parameters
 * typedef struct BlastHSPMyWriterParams {...} BlastHSPMyWriterParams;
 * // Function to create MyWriter filtering parameters.  
 * BlastHSPMyWriterParams* 
 * BlastHSPMyWriterParamsNew(...);
 * // Function to create the BlastWriterInfo structure associated with MyWriter.
 * BlastHSPWriterInfo*
 * BlastHSPMyWriterInfoNew(BlastHSPMyWriterParams*);
 * }
 * @endcode
 *
 * In addition, the following should be implemented in hspfilter_mywriter.c:
 * @code
 * extern "C" {
 * // Any auxillary data structures MyWriter may use to store between-call data.
 * typedef struct MyWriterData {...} MyWriterData;
 * // The following are functions to implement BlastHSPWriter interface.
 * // Function to initiate MyWriterData from BlastHSPResults
 * static int s_MyWriterInit(void*, BlastHSPResults*);
 * // Function to finalize MyWriterData to BlastHSPResults
 * static int s_MyWriterFinal(void*, BlastHSPResults*);
 * // Function to process BlastHSPList and save results to MyWriterData
 * // Must call Blast_HSPListFree() before returning
 * static int s_MyWriterRun(void*, BlastHSPList*);
 * // Function to free MyWriter.
 * // Must free mywriter parameters before returning
 * static BlastHSPWriter* s_MyWriterFree(BlastHSPWriter*);
 * // Function to create MyWriter.
 * static BlastHSPWriter* s_MyWriterNew(void*, BlastQueryInfo*);
 * }
 * @endcode
 *
 */

#ifdef __cplusplus
}
#endif

#endif /* ALGO_BLAST_CORE__BLAST_HSPFILTER_H */
