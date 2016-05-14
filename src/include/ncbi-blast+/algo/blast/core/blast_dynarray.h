/* $Id: blast_dynarray.h 103491 2007-05-04 17:18:18Z kazimird $
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
 *  Author: Christiam Camacho
 *
 */

/** @file blast_dynarray.h
 * Declarations for various dynamic array types
 */

#ifndef ALGO_BLAST_CORE__BLAST_DYNARRAY__H
#define ALGO_BLAST_CORE__BLAST_DYNARRAY__H

#include <algo/blast/core/ncbi_std.h>
#include <algo/blast/core/blast_export.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Data structure to maintain a dynamically allocated array of Uint4 */
typedef struct SDynamicUint4Array {
    Uint4 num_used;      /**< number of elements used in this array */
    Uint4 num_allocated; /**< size of array below */
    Uint4* data;         /**< array of Uint4 */
} SDynamicUint4Array;

/** Initial number of elements allocated */
#define INIT_NUM_ELEMENTS 8

/** Allocate a dynamic array of Uint4 with the initial size of
 * INIT_NUM_ELEMENTS
 * @return NULL if out of memory otherwise, newly allocated structure
 */
NCBI_XBLAST_EXPORT SDynamicUint4Array* 
DynamicUint4ArrayNew();

/** Allocate a dynamic array of Uint4 with an initial size of specified as an
 * argument to the function
 * @param init_num_elements Number of elements initially allocated [in]
 * @return NULL if out of memory otherwise, newly allocated structure
 */
NCBI_XBLAST_EXPORT SDynamicUint4Array* 
DynamicUint4ArrayNewEx(Uint4 init_num_elements);

/** Deallocates a dynamic array structure 
 * @param arr data structure to free [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT SDynamicUint4Array* 
DynamicUint4ArrayFree(SDynamicUint4Array* arr);

/** Append a Uint4 to the dynamic array structure
 * @param arr data structure to manipulate [in]
 * @param element element to add [in]
 * @return 0 on success or BLASTERR_MEMORY if memory reallocation was needed
 * and this failed
 */
NCBI_XBLAST_EXPORT Int2
DynamicUint4Array_Append(SDynamicUint4Array* arr, Uint4 element);

/** Make a deep copy of the src dynamic array
 * @param src data structure to duplicate [in]
 * @return newly allocated data structure or NULL if out of memory
 */
NCBI_XBLAST_EXPORT SDynamicUint4Array*
DynamicUint4Array_Dup(const SDynamicUint4Array* src);

/** Make a shallow copy of the src dynamic array into the dest dynamic array
 * @param dest data structure to copy to. Its contents will be overwritten and
 * memory might be allocated if necessary [in|out]
 * @param src data structure to copy from [in]
 * @return newly allocated data structure or NULL if out of memory
 */
NCBI_XBLAST_EXPORT Int4
DynamicUint4Array_Copy(SDynamicUint4Array* dest,
                       const SDynamicUint4Array* src);

/** Compares dynamic arrays a and b for equality of its contents
 * @param a dynamic array to compare [in]
 * @param b dynamic array to compare [in]
 * @return TRUE if equal, FALSE otherwise
 */
NCBI_XBLAST_EXPORT Boolean
DynamicUint4Array_AreEqual(const SDynamicUint4Array* a,
                           const SDynamicUint4Array* b);

/* ========================================================================== */
/** Data structure to maintain a dynamically allocated array of Int4 */
typedef struct SDynamicInt4Array {
    Uint4 num_used;      /**< number of elements used in this array */
    Uint4 num_allocated; /**< size of array below */
    Int4* data;          /**< array of Int4 */
} SDynamicInt4Array;

/** Allocate a dynamic array of Int4 with the initial size of
 * INIT_NUM_ELEMENTS
 * @return NULL if out of memory otherwise, newly allocated structure
 */
NCBI_XBLAST_EXPORT SDynamicInt4Array* 
DynamicInt4ArrayNew();

/** Deallocates a dynamic array structure 
 * @param arr data structure to free [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT SDynamicInt4Array* 
DynamicInt4ArrayFree(SDynamicInt4Array* arr);

/** Append a Int4 to the dynamic array structure
 * @param arr data structure to manipulate [in]
 * @param element element to add [in]
 * @return 0 on success or BLASTERR_MEMORY if memory reallocation was needed
 * and this failed
 */
NCBI_XBLAST_EXPORT Int2
DynamicInt4Array_Append(SDynamicInt4Array* arr, Int4 element);

/* ========================================================================== */
/** Elements of the SDynamicSGenCodeNodeArray dynamic array */
typedef struct SGenCodeNode {
    Uint4   gc_id;      /**< Genetic code id */
    Uint1*  gc_str;     /**< Genetic code string */
} SGenCodeNode;

/** Initial number of elements allocated, based on the fact that there are only
 * a handful of genetic codes available in the NCBI toolkits */
#define INIT_NUM_GEN_CODES 30

/** Data structure to maintain a dynamically allocated array of SGenCodeNode */
typedef struct SDynamicSGenCodeNodeArray {
    Uint4 num_used;      /**< number of elements used in this array */
    Uint4 num_allocated; /**< size of array below */
    SGenCodeNode* data;          /**< array of SGenCodeNode */
} SDynamicSGenCodeNodeArray;

/** Allocate a dynamic array of SGenCodeNode with the initial size of
 * INIT_NUM_ELEMENTS
 * @return NULL if out of memory otherwise, newly allocated structure
 */
NCBI_XBLAST_EXPORT SDynamicSGenCodeNodeArray* 
DynamicSGenCodeNodeArrayNew();

/** Deallocates a dynamic array structure 
 * @param arr data structure to free [in]
 * @return NULL
 */
NCBI_XBLAST_EXPORT SDynamicSGenCodeNodeArray* 
DynamicSGenCodeNodeArrayFree(SDynamicSGenCodeNodeArray* arr);

/** Append a SGenCodeNode to the dynamic array structure
 * @note No duplicate elements will be inserted
 * @param arr data structure to manipulate [in]
 * @param element element to add [in]
 * @return 0 on success or BLASTERR_MEMORY if memory reallocation was needed
 * and this failed, or BLASTERR_INVALIDPARAM if the genetic string field of the
 * element is NULL.
 */
NCBI_XBLAST_EXPORT Int2
DynamicSGenCodeNodeArray_Append(SDynamicSGenCodeNodeArray* arr, 
                                SGenCodeNode element);

/** Searches the dynamic array for any element that matches the requested
 * genetic code id
 * @param arr dynamic array to search [in]
 * @param gen_code_id genetic code id to search [in]
 * @return genetic code string (owned by this structure) or NULL if not found
 */
Uint1*
DynamicSGenCodeNodeArray_Find(const SDynamicSGenCodeNodeArray* arr,
                              Uint4 gen_code_id);


#ifdef __cplusplus
}
#endif
#endif /* !ALGO_BLAST_CORE__BLAST_DYNARRAY__H */
