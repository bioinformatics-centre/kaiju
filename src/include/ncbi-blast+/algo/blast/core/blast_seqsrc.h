/*  $Id: blast_seqsrc.h 426803 2014-02-13 14:19:14Z fongah2 $
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

/** @file blast_seqsrc.h
 * Declaration of ADT to retrieve sequences for the BLAST engine.
 *
 * @todo The following are pending items to be addressed in the BlastSeqSrc
 * interface and implementations
 * - Make "multi-sequence" BlastSeqSrc implementations MT-safe                  
 * - Provide possibility of retrieving an error from BlastSeqSrc (not just 
 *   initialization errors).                                               
 * - Constrain all use of the BlastSeqSrc to CORE BLAST                         
 */

#ifndef ALGO_BLAST_CORE__BLAST_SEQSRC__H
#define ALGO_BLAST_CORE__BLAST_SEQSRC__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_def.h>
#include <algo/blast/core/blast_export.h>
#include <algo/blast/core/blast_message.h>
#include <algo/blast/core/blast_encoding.h>

#ifdef __cplusplus
extern "C" {
#endif

/** The BlastSeqSrc ADT is an opaque data type that defines an interface which
 *  is used by the core BLAST code to retrieve sequences.
 *  The interface currently provides the following services:
 *  - Retrieving number of sequences in a set
 *  - Retrieving the total length (in bases/residues) of sequences in a set.
 *  - Retrieving the length of the longest sequence in a set
 *  - Retrieving an individual sequence in a user-specified encoding by ordinal
 *    id (index into a set)
 *  - Retrieving the length of a given sequence in a set by ordinal id
 *  - Allow MT-safe iteration over sequences in a set through the
 *    BlastSeqSrcIterator abstraction, as well as resetting of any applicable
 *    implementation internal 'bookmarks' which keep track of the iteration
 *    progress, as to allow multiple passes over the set of sequences (@sa
 *    BlastSeqSrcResetChunkIterator).
 *  .
 *
 *  Currently available client implementations of the BlastSeqSrc API include:
 *  - ReaddbBlastSeqSrcInit (C toolkit)
 *  - SeqDbBlastSeqSrcInit (C++ toolkit)
 *  - MultiSeqBlastSeqSrcInit (C/C++ toolkit)
 *  .
 *  For more details, see also @ref _impl_blast_seqsrc_howto
 */
typedef struct BlastSeqSrc BlastSeqSrc;

/** Blast Sequence Source Iterator, designed to be used in conjunction with the
 * BlastSeqSrc to provide MT-safe iteration over the sequences in the
 * BlastSeqSrc 
 * @todo This is still coupled to the BLAST database implementations of the
 * BlastSeqSrc, could be made more generic if needed.
 */
typedef struct BlastSeqSrcIterator BlastSeqSrcIterator;

/** Structure that contains the information needed for BlastSeqSrcNew to fully
 * populate the BlastSeqSrc structure it returns */
typedef struct BlastSeqSrcNewInfo BlastSeqSrcNewInfo;

/** Allocates memory for a BlastSeqSrc structure and then invokes the
 * constructor function defined in its first argument, passing the 
 * ctor_argument member of that same structure. If the constructor function
 * pointer is not set or there is a memory allocation failure, NULL is
 * returned.
 * @note This function need not be called directly by client code as all the
 * implementations of the BlastSeqSrc interface provide a function(s) which
 * call this function (@ref MultiSeqBlastSeqSrcInit, @ref SeqDbBlastSeqSrcInit)
 * @param bssn_info Structure defining constructor and its argument to be
 *        invoked from this function [in]
 * @return a properly initialized BlastSeqSrc structure or NULL.
 */
NCBI_XBLAST_EXPORT
BlastSeqSrc* BlastSeqSrcNew(const BlastSeqSrcNewInfo* bssn_info);

/** Copy function: needed to guarantee thread safety. Copies the contents of an
 * input BlastSeqSrc, then calls a copier function, provided by the 
 * implementation, to achieve multi-thread safety.
 * @todo Is this function really needed? 
 * @param seq_src BlastSeqSrc to copy [in]
 * @return An MT-safe copy of the structure passed in, NULL in case of memory
 *         allocation failure, or, if no copy function was provided by the
 *         implementation, a bitwise copy of the input.
 */
NCBI_XBLAST_EXPORT
BlastSeqSrc* BlastSeqSrcCopy(const BlastSeqSrc* seq_src);

/** Frees the BlastSeqSrc structure by invoking the destructor function set by
 * the user-defined constructor function when the structure is initialized
 * (indirectly, by BlastSeqSrcNew). If the destructor function pointer is not
 * set, a memory leak could occur. Note that it is the implementation's
 * destructor responsibility to free the BlastSeqSrc structure by calling
 * sfree.
 * @param seq_src BlastSeqSrc to free [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
BlastSeqSrc* BlastSeqSrcFree(BlastSeqSrc* seq_src);

/** Get the number of sequences contained in the sequence source.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int4
BlastSeqSrcGetNumSeqs(const BlastSeqSrc* seq_src);

/** Get the number of sequences used for calculation of expect values etc.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int4
BlastSeqSrcGetNumSeqsStats(const BlastSeqSrc* seq_src);

/** Get the length of the longest sequence in the sequence source.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int4
BlastSeqSrcGetMaxSeqLen(const BlastSeqSrc* seq_src);

/** Get the length of the longest sequence in the sequence source.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int4
BlastSeqSrcGetMinSeqLen(const BlastSeqSrc* seq_src);

/** Get the average length of all sequences in the sequence source.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int4
BlastSeqSrcGetAvgSeqLen(const BlastSeqSrc* seq_src);

/** Get the total length of all sequences in the sequence source.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int8
BlastSeqSrcGetTotLen(const BlastSeqSrc* seq_src);

/** Get the total length of all sequences for calculation of expect value etc.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Int8
BlastSeqSrcGetTotLenStats(const BlastSeqSrc* seq_src);

/** Get the Blast Sequence source name (e.g.: BLAST database name).
 * Here the full name (path and file name) should be returned.  
 * If an alias file is present the return value should be the full name of 
 * the alias file.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
const char*
BlastSeqSrcGetName(const BlastSeqSrc* seq_src);

/** Find if the Blast Sequence Source contains protein or nucleotide sequences.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Boolean
BlastSeqSrcGetIsProt(const BlastSeqSrc* seq_src);

/** Find if the Blast Sequence Source supports partial fetching
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
Boolean
BlastSeqSrcGetSupportsPartialFetching(const BlastSeqSrc* seq_src);

#define BLAST_SEQSRC_MINGAP     1024  /**< Minimal gap allowed in range list */
#define BLAST_SEQSRC_OVERHANG   1024  /**< Extension for each new range added */
#define BLAST_SEQSRC_MINLENGTH  10    /**< Default minimal sequence length */

/** Structure used as the argument to function SetRanges */
typedef struct BlastSeqSrcSetRangesArg {
    /** Oid in BLAST database, index in an array of sequences, etc [in] */
    Int4 oid;
    
    /** initial allocation*/
    Int4 capacity;

    /** Number of actual ranges contained */
    Int4 num_ranges;

    /** Ranges in sorted order [in] */
    Int4 *ranges;
} BlastSeqSrcSetRangesArg;

/** new setrangearg */
BlastSeqSrcSetRangesArg *
BlastSeqSrcSetRangesArgNew(Int4 num_ranges);

/** add new range 
 * @return 0 in case of success, otherwise 1
 */
Int2
BlastSeqSrcSetRangesArgAddRange(BlastSeqSrcSetRangesArg *arg, 
                                Int4 begin, Int4 end);

/** free setrangearg */
void 
BlastSeqSrcSetRangesArgFree(BlastSeqSrcSetRangesArg * arg);


/** build BlastSeqSrcSetRangesArg from range list*/
void
BlastSeqSrcSetRangesArgBuild(BlastSeqSrcSetRangesArg *arg);

/** Setting the ranges for partial fetching */
NCBI_XBLAST_EXPORT
void
BlastSeqSrcSetSeqRanges(const BlastSeqSrc* seq_src, 
                        BlastSeqSrcSetRangesArg* setranges_arg);


/** Structure used as the second argument to functions satisfying the 
  GetSeqBlkFnPtr signature, associated with index-based 
  implementations of the BlastSeqSrc interface. Index-based implementations
  include BLAST databases or an array/vector of sequences. */
typedef struct BlastSeqSrcGetSeqArg {
    /** Oid in BLAST database, index in an array of sequences, etc [in] */
    Int4 oid;

    /** Encoding of sequence, i.e.: eBlastEncodingProtein,
     * eBlastEncodingNucleotide, etc [in] */
    EBlastEncoding encoding;

    /** This option allows the BLAST engine to communicate with the BlastSeqSrc
     * that the offset ranges for a given OID should be reset and the entire
     * sequence data should be fetched. The motivation for this option is to
     * exploit CSeqDB's performance feature that allows to retrieve only
     * pre-selected portions of the sequence data for the traceback stage.
     * BlastSeqSrc implementations that do not have this feature can safely
     * ignore this field.
     * By default, this should be set to FALSE. [in] */
    Boolean reset_ranges;

    /** Check whether an OID is excluded due to overlapping filtering.
     * The disease is rare, and the test for it is somewhat expensive,
     * so it is deferred to the traceback stage.
     * TRUE to enable this test. [in] */
    Boolean check_oid_exclusion;

    /** Sequence to return, if NULL, it should allocated by GetSeqBlkFnPtr
     * (using BlastSeqBlkNew or BlastSetUp_SeqBlkNew), else its contents are 
     * freed (using BlastSequenceBlkClean) and the structure is reused [out]*/
    BLAST_SequenceBlk* seq;
} BlastSeqSrcGetSeqArg;

/* Return values from BlastSeqSrcGetSequence */
#define BLAST_SEQSRC_EXCLUDED	-3  /**< Sequence excluded due to filtering */
#define BLAST_SEQSRC_ERROR      -2  /**< Error while retrieving sequence */
#define BLAST_SEQSRC_EOF        -1  /**< No more sequences available */
#define BLAST_SEQSRC_SUCCESS     0  /**< Successful sequence retrieval */

/** Retrieve an individual sequence.
 * @param seq_src the BLAST sequence source [in]
 * @param getseq_arg arguments to aid retrieval of sequence data from the
 * BlastSeqSrc [in|out]
 * @return one of the BLAST_SEQSRC_* defined in blast_seqsrc.h
 */
NCBI_XBLAST_EXPORT
Int2
BlastSeqSrcGetSequence(const BlastSeqSrc* seq_src, 
                       BlastSeqSrcGetSeqArg* getseq_arg);

/** Retrieve sequence length (number of residues/bases)
 * @param seq_src the BLAST sequence source [in]
 * @param oid ordinal id of the sequence desired (should be Uint4) [in]
 */
NCBI_XBLAST_EXPORT
Int4
BlastSeqSrcGetSeqLen(const BlastSeqSrc* seq_src, void* oid);

/** Deallocate individual sequence.
 * @param seq_src the BLAST sequence source [in]
 * @param getseq_arg contains sequence to deallocate [in|out]
 */
NCBI_XBLAST_EXPORT
void
BlastSeqSrcReleaseSequence(const BlastSeqSrc* seq_src,
                           BlastSeqSrcGetSeqArg* getseq_arg);

/** Function to retrieve NULL terminated string containing the description 
 * of an initialization error or NULL. This function MUST ALWAYS be called 
 * after calling one of the client implementation's Init functions. If the
 * return value is not NULL, invoking any other functionality from the
 * BlastSeqSrc will result in undefined behavior. Caller is responsible for 
 * deallocating the return value. 
 * @param seq_src BlastSeqSrc from which to get an error [in]
 * @return error message or NULL
 */
NCBI_XBLAST_EXPORT
char* BlastSeqSrcGetInitError(const BlastSeqSrc* seq_src);

/* This is a *debugging* feature for the code that implements the various
 * composition-based statistics types of processing, added at the request of
 * the authors of that code */
#ifdef KAPPA_PRINT_DIAGNOSTICS

/** Data structure to contain a list of gis */
typedef struct Blast_GiList {
    Int4* data;             /**< Gis are stored here */
    size_t num_allocated;   /**< Size of the data array */
    size_t num_used;        /**< Number of elements populated in data */
} Blast_GiList;

/** Allocate a gi list with default size */
NCBI_XBLAST_EXPORT
Blast_GiList* Blast_GiListNew(void);

/** Allocate a gi list with the requested size 
 * @param list_size initial list size [in]
 */
NCBI_XBLAST_EXPORT
Blast_GiList* Blast_GiListNewEx(size_t list_size);

/** Deallocate memory associated with the gi list
 * @return NULL
 */
NCBI_XBLAST_EXPORT
Blast_GiList* Blast_GiListFree(Blast_GiList* gilist);

/* Return values */

/** Invalid parameter used in a function call */
NCBI_XBLAST_EXPORT
extern const Int2 kBadParameter;
/** Failure due to out-of-memory condition */
NCBI_XBLAST_EXPORT
extern const Int2 kOutOfMemory;

/** Appends to an existing gi list, allocating memory if necessary
 * @return 0 on success or one of the constants defined above
 */
Int2 Blast_GiList_Append(Blast_GiList* gilist, Int4 gi);

/** Retrieve a sequence's gi (if available). Caller is responsible for
 * deallocating return value with Blast_GiListFree
 * @param seq_src the BLAST sequence source [in]
 * @param oid ordinal id of the sequence desired (should be Uint4) [in]
 * @return the gis corresponding to the ordinal id requested, an empty list if
 * there are no gis in the BlastSeqSrc implementation, or NULL in case of
 * error
 */
NCBI_XBLAST_EXPORT
Blast_GiList*
BlastSeqSrcGetGis(const BlastSeqSrc* seq_src, void* oid);

#endif /* KAPPA_PRINT_DIAGNOSTICS */

/******************** BlastSeqSrcIterator API *******************************/

/** Allocate and initialize an iterator over a BlastSeqSrc with a default chunk
 * size for MT-safe iteration.
 * @return pointer to initialized iterator for BlastSeqSrc
 */
NCBI_XBLAST_EXPORT
BlastSeqSrcIterator* BlastSeqSrcIteratorNew(void);

/** How many database sequences to process in one database chunk. 
 * this value is overriden in seqdb implementation, where the number of sequences
 * is determined by the mmaped slice size */
extern const unsigned int kBlastSeqSrcDefaultChunkSize;

/** Allocate and initialize an iterator over a BlastSeqSrc. 
 * @param chunk_sz sets the chunk size of the iterator, if zero 
 *    use kBlastSeqSrcDefaultChunkSize (above). [in]
 * @return pointer to initialized iterator for BlastSeqSrc
 */
NCBI_XBLAST_EXPORT
BlastSeqSrcIterator* BlastSeqSrcIteratorNewEx(unsigned int chunk_sz);

/** Frees the BlastSeqSrcIterator structure. 
 * @param itr BlastSeqSrcIterator to free [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT
BlastSeqSrcIterator* BlastSeqSrcIteratorFree(BlastSeqSrcIterator* itr);

/** Increments the BlastSeqSrcIterator.
 * @param itr the BlastSeqSrcIterator to increment.
 * @param seq_src the underlying BlastSeqSrc
 * @return one of the BLAST_SEQSRC_* defined in blast_seqsrc.h
 */
NCBI_XBLAST_EXPORT
Int4 BlastSeqSrcIteratorNext(const BlastSeqSrc* seq_src, 
                             BlastSeqSrcIterator* itr);

/** Reset the internal "bookmark" of the last chunk for iteration provided by 
 * this object.
 * @param seq_src the BLAST sequence source [in]
 */
NCBI_XBLAST_EXPORT
void BlastSeqSrcResetChunkIterator(BlastSeqSrc* seq_src);

/** Set the number of threads for MT mode 
 * @param nthreads the number of threads [in]
 */
NCBI_XBLAST_EXPORT
void BlastSeqSrcSetNumberOfThreads(BlastSeqSrc* seq_src, int nthreads);

/*****************************************************************************/

#ifdef __cplusplus
}
#endif

#endif /* !ALGO_BLAST_CORE__BLAST_SEQSRC__H */
