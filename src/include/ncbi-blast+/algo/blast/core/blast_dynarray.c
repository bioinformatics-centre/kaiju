/* $Id: blast_dynarray.c 113776 2007-11-08 22:38:18Z camacho $
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
 * Author: Christiam Camacho
 *
 */

/** @file blast_dynarray.c
 * Definitions for various dynamic array types
 */


#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: blast_dynarray.c 113776 2007-11-08 22:38:18Z camacho $";
#endif /* SKIP_DOXYGEN_PROCESSING */

#include "blast_dynarray.h"
#include <algo/blast/core/blast_util.h>
#include <algo/blast/core/blast_message.h>  /* for error codes */

/** factor by which these arrays are resized */
const size_t kResizeFactor = 2;

SDynamicUint4Array* DynamicUint4ArrayFree(SDynamicUint4Array* arr)
{
    if ( !arr ) {
        return NULL;
    }
    if (arr->data) {
        sfree(arr->data);
    }
    sfree(arr);
    return NULL;
}

SDynamicUint4Array* DynamicUint4ArrayNew()
{
    return DynamicUint4ArrayNewEx(INIT_NUM_ELEMENTS);
}

SDynamicUint4Array* DynamicUint4ArrayNewEx(Uint4 init_num_elements)
{
    SDynamicUint4Array* retval = 
        (SDynamicUint4Array*) calloc(1, sizeof(SDynamicUint4Array));
    if ( !retval ) {
        return NULL;
    }
    retval->data = (Uint4*) calloc(init_num_elements, sizeof(Uint4));
    if ( !retval->data ) {
        return DynamicUint4ArrayFree(retval);
    }
    retval->num_allocated = init_num_elements;

    return retval;
}

/** Grow a dynamic array of Uint4 elements
 * @param arr Structure of array data [in][out]
 * @return zero on success
 */
static Int2
s_DynamicUint4Array_ReallocIfNecessary(SDynamicUint4Array* arr)
{
    ASSERT(arr);

    if (arr->num_used+1 > arr->num_allocated) {
        /* we need more room for elements */
        Uint4* reallocation = 
            (Uint4*) realloc(arr->data, 
                     kResizeFactor * arr->num_allocated * sizeof(*arr->data));
        if ( !reallocation ) {
            return BLASTERR_MEMORY;
        }
        arr->data = reallocation;
        arr->num_allocated *= kResizeFactor;
    }
    return 0;
}

Int2
DynamicUint4Array_Append(SDynamicUint4Array* arr, Uint4 element)
{
    Int2 retval = 0;
    ASSERT(arr);

    if ( (retval = s_DynamicUint4Array_ReallocIfNecessary(arr)) != 0) {
        return retval;
    }
    arr->data[arr->num_used] = element;
    arr->num_used++;
    return retval;
}

SDynamicUint4Array*
DynamicUint4Array_Dup(const SDynamicUint4Array* src)
{
    SDynamicUint4Array* retval = NULL;

    if ( !src ) {
        return retval;
    }

    retval = DynamicUint4ArrayNewEx(src->num_allocated);
    memcpy((void*) retval->data, (void*) src->data, 
           sizeof(*src->data) * src->num_used);
    return retval;
}

Int4
DynamicUint4Array_Copy(SDynamicUint4Array* dest,
                       const SDynamicUint4Array* src)
{
    Uint4 i = 0;        /* index into arrays */

    if (dest->num_allocated < src->num_allocated) {
        /* we need more room for elements */
        Uint4* reallocation = (Uint4*) realloc(dest->data,
                                     src->num_allocated * sizeof(*dest->data));
        if ( !reallocation ) {
            return BLASTERR_MEMORY;
        }
        dest->data = reallocation;
        dest->num_allocated = src->num_allocated;
    }

    for (i = 0; i < src->num_used; i++) {
        dest->data[i] = src->data[i];
    }
    dest->num_used = src->num_used;
    return 0;
}

Boolean
DynamicUint4Array_AreEqual(const SDynamicUint4Array* a,
                           const SDynamicUint4Array* b)
{
    Uint4 i = 0;        /* index into arrays */

    if (a->num_used != b->num_used) {
        return FALSE;
    }

    for (i = 0; i < a->num_used; i++) {
        if (a->data[i] != b->data[i]) {
            return FALSE;
        }
    }

    return TRUE;
}

SDynamicInt4Array* DynamicInt4ArrayFree(SDynamicInt4Array* arr)
{
    if ( !arr ) {
        return NULL;
    }
    if (arr->data) {
        sfree(arr->data);
    }
    sfree(arr);
    return NULL;
}

SDynamicInt4Array* DynamicInt4ArrayNew()
{
    SDynamicInt4Array* retval = 
        (SDynamicInt4Array*) calloc(1, sizeof(SDynamicInt4Array));
    if ( !retval ) {
        return NULL;
    }
    retval->data = (Int4*) calloc(INIT_NUM_ELEMENTS, sizeof(Int4));
    if ( !retval->data ) {
        return DynamicInt4ArrayFree(retval);
    }
    retval->num_allocated = INIT_NUM_ELEMENTS;

    return retval;
}

/** Grow a dynamic array of Int4 elements
 * @param arr Structure of array data [in][out]
 * @return zero on success
 */
static Int2
s_DynamicInt4Array_ReallocIfNecessary(SDynamicInt4Array* arr)
{
    ASSERT(arr);

    if (arr->num_used+1 > arr->num_allocated) {
        /* we need more room for elements */
        Int4* reallocation = (Int4*) realloc(arr->data,
                     kResizeFactor * arr->num_allocated * sizeof(*arr->data));
        if ( !reallocation ) {
            return BLASTERR_MEMORY;
        }
        arr->data = reallocation;
        arr->num_allocated *= kResizeFactor;
    }
    return 0;
}

Int2
DynamicInt4Array_Append(SDynamicInt4Array* arr, Int4 element)
{
    Int2 retval = 0;
    ASSERT(arr);

    if ( (retval = s_DynamicInt4Array_ReallocIfNecessary(arr)) != 0) {
        return retval;
    }
    arr->data[arr->num_used] = element;
    arr->num_used++;
    return retval;
}

SDynamicSGenCodeNodeArray* 
DynamicSGenCodeNodeArrayFree(SDynamicSGenCodeNodeArray* arr)
{
    if ( !arr ) {
        return NULL;
    }
    if (arr->data) {
        Uint4 i = 0;
        for (i = 0; i < arr->num_used; i++) {
            sfree(arr->data[i].gc_str);
        }
        sfree(arr->data);
    }
    sfree(arr);
    return NULL;
}

SDynamicSGenCodeNodeArray* DynamicSGenCodeNodeArrayNew()
{
    SDynamicSGenCodeNodeArray* retval = 
        (SDynamicSGenCodeNodeArray*) calloc(1, 
                                            sizeof(SDynamicSGenCodeNodeArray));
    if ( !retval ) {
        return NULL;
    }
    retval->data = (SGenCodeNode*) calloc(INIT_NUM_GEN_CODES, 
                                          sizeof(SGenCodeNode));
    if ( !retval->data ) {
        return DynamicSGenCodeNodeArrayFree(retval);
    }
    retval->num_allocated = INIT_NUM_GEN_CODES;

    return retval;
}

/** Grow a dynamic array of SGenCodeNode elements
 * @param arr Structure of array data [in][out]
 * @return zero on success
 */
static Int2
s_DynamicSGenCodeNodeArray_ReallocIfNecessary(SDynamicSGenCodeNodeArray* arr)
{
    ASSERT(arr);

    if (arr->num_used+1 > arr->num_allocated) {
        /* we need more room for elements */
        SGenCodeNode* reallocation = (SGenCodeNode*) realloc(arr->data,
                     kResizeFactor * arr->num_allocated * sizeof(*arr->data));
        if ( !reallocation ) {
            return BLASTERR_MEMORY;
        }
        arr->data = reallocation;
        arr->num_allocated *= kResizeFactor;
    }
    return 0;
}

/** Determine if the array is sorted
 * @param arr array to examine [in]
 */
static Boolean
s_DynamicSGenCodeNodeArray_IsSorted(const SDynamicSGenCodeNodeArray* arr)
{
    Int4 index;

    if ( !arr || arr->num_used <= 1) {
        return TRUE;
    }

    for (index = arr->num_used - 1; index > 0; index--) {
        if (arr->data[index].gc_id < arr->data[index-1].gc_id) {
            return FALSE;
        }
    }
    return TRUE;
}

/** Compare function for sorting SGenCodeNode elements
 * @param a first element to compare [in]
 * @param b second element to compare [in]
 */
static int
s_SGenCodeNodeCompare(const void* a, const void* b)
{
    SGenCodeNode* ptr1 = (SGenCodeNode*) a;
    SGenCodeNode* ptr2 = (SGenCodeNode*) b;
    return BLAST_CMP(ptr1->gc_id, ptr2->gc_id);
}

/** Sort the dynamic array structure
 * @note this function should be called prior to invoking
 * DynamicSGenCodeNodeArray_Find
 * @param arr data structure to manipulate [in]
 */
static void
s_DynamicSGenCodeNodeArray_Sort(SDynamicSGenCodeNodeArray* arr)
{
    if ( !arr || arr->num_used <= 1) {
        return;
    }

    if ( !s_DynamicSGenCodeNodeArray_IsSorted(arr) ) {
        qsort(arr->data, arr->num_used, sizeof(SGenCodeNode),
              s_SGenCodeNodeCompare);
    }
}

Int2
DynamicSGenCodeNodeArray_Append(SDynamicSGenCodeNodeArray* arr, 
                                SGenCodeNode element)
{
    Int2 retval = 0;
    ASSERT(arr);

    if (element.gc_str == NULL) {
        return BLASTERR_INVALIDPARAM;
    }

    if (DynamicSGenCodeNodeArray_Find(arr, element.gc_id)) {
        return retval;
    }

    if ( (retval = s_DynamicSGenCodeNodeArray_ReallocIfNecessary(arr)) != 0) {
        return retval;
    }
    arr->data[arr->num_used].gc_str = (Uint1*)
        BlastMemDup(element.gc_str, GENCODE_STRLEN);
    if ( arr->data[arr->num_used].gc_str == NULL) {
        return BLASTERR_MEMORY;    
    }
    arr->data[arr->num_used].gc_id = element.gc_id;
    arr->num_used++;
    s_DynamicSGenCodeNodeArray_Sort(arr);
    return retval;
}

/** Perform a binary search for the element containing gen_code_id 
 * @param arr dynamic array to search [in]
 * @param gen_code_id genetic code to find [in]
 * @return index where the genetic code was found
 */
static Int4
s_DynamicSGenCodeNodeArray_BinSearch(const SDynamicSGenCodeNodeArray* arr,
                                     Uint4 gen_code_id)
{
    Int4 m = 0, b = 0, e = 0, size = 0;

    size = (Int4)arr->num_used;
    b = 0;
    e = size;
    while (b < e - 1) {
        m = (b + e) / 2;
        if (arr->data[m].gc_id > gen_code_id) {
            e = m;
        } else {
            b = m;
        }
    }
    return b;
}

Uint1*
DynamicSGenCodeNodeArray_Find(const SDynamicSGenCodeNodeArray* arr,
                              Uint4 gen_code_id)
{
    Uint4 index = (Uint4)s_DynamicSGenCodeNodeArray_BinSearch(arr, gen_code_id);
    if (index < arr->num_used && arr->data[index].gc_id == gen_code_id) {
        return arr->data[index].gc_str;
    }
    return NULL;
}

