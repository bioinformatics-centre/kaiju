#ifndef CONNECT___NCBI_TYPES__H
#define CONNECT___NCBI_TYPES__H

/* $Id: ncbi_types.h 445734 2014-09-08 13:43:28Z lavr $
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
 * Author:  Anton Lavrentiev
 *
 * @file
 * File Description:
 *   Special types for core library.
 *
 * Timeout:
 *    struct STimeout
 *
 * Switch:
 *    ESwitch    (on/off/default)
 *
 * Fixed-size size_t and time_t equivalents:
 *    TNCBI_Size
 *    TNCBI_Time
 *       these two we need to use when mixing 32/64 bit programs
 *       which make simultaneous access to inter-process communication
 *       data areas, like shared memory segments
 *
 */

#include <connect/connect_export.h>
#ifndef _WIN32
#  ifndef   __STDC_FORMAT_MACROS
#    define __STDC_FORMAT_MACROS
#  endif /*!__STDC_FORMAT_MACROS*/
#  include <inttypes.h>
#  include <stdint.h>
#endif /*_WIN32*/
#include <stddef.h>


/** @addtogroup UtilityFunc
 *
 * @{
 */


#ifdef __cplusplus
extern "C" {
#endif


/** Timeout structure
 *
 * @sa CTimeout, g_CTimeoutToSTimeout, g_STimeoutToCTimeout
 */
typedef struct STimeoutTag {
    unsigned int sec;  /**< seconds                         */
    unsigned int usec; /**< microseconds (modulo 1,000,000) */
} STimeout;

#define kDefaultTimeout   ((const STimeout*)(-1))
#define kInfiniteTimeout  ((const STimeout*)( 0))


extern NCBI_XCONNECT_EXPORT unsigned long NcbiTimeoutToMs
(const STimeout* timeout
);


extern NCBI_XCONNECT_EXPORT STimeout*     NcbiMsToTimeout
(STimeout*       timeout,
 unsigned long   ms
 );


#ifndef NCBI_ESWITCH_DEFINED
#define NCBI_ESWITCH_DEFINED

/*
 * ATTENTION!   Do not change this enumeration!
 *
 * It must always be kept in sync with its C++ counterpart defined in
 * "corelib/ncbimisc.hpp". If you absolutely(sic!) need to alter this
 * type, please apply equivalent changes to both definitions.
 */

/** Aux. enum to set/unset/default various features.
 */
typedef enum ENcbiSwitch {
    eOff = 0,
    eOn,
    eDefault
} ESwitch;

#endif /*!NCBI_ESWITCH_DEFINED*/


#ifndef NCBI_EOWNERSHIP_DEFINED
#define NCBI_EOWNERSHIP_DEFINED

/*
 * ATTENTION!   Do not change this enumeration!
 *
 * It must always be kept in sync with its C++ counterpart defined in
 * "corelib/ncbimisc.hpp". If you absolutely(sic!) need to alter this
 * type, please apply equivalent changes to both definitions.
 */

/** Ownership relations between objects.
 *
 * Can be used to define or transfer ownership of objects.
 * For example, specify if a CSocket object owns its underlying SOCK object.
 */
typedef enum ENcbiOwnership {
    eNoOwnership,       /** No ownership is assumed                 */
    eTakeOwnership      /** An object can take ownership of another */
} EOwnership;

#endif /*!NCBI_EOWNERSHIP_DEFINED*/


/** Fixed-size analogs of size_t and time_t (mainly for IPC)
 */
typedef unsigned int TNCBI_Size;
typedef unsigned int TNCBI_Time;

#define NCBI_TIME_INFINITE ((TNCBI_Time)(-1))


/** Big unsigned integer for file size and position
 */

#if defined(__MINGW32__)  ||  defined(__MINGW64__)
typedef unsigned long long  TNCBI_BigCount;
#  define NCBI_BIGCOUNT_FORMAT_SPEC      "I64u"
#  define NCBI_BIGCOUNT_FORMAT_SPEC_HEX  "I64x"
#elif defined(_WIN32)
typedef unsigned __int64    TNCBI_BigCount;
#  define NCBI_BIGCOUNT_FORMAT_SPEC      "I64u"
#  define NCBI_BIGCOUNT_FORMAT_SPEC_HEX  "I64x"
#else
typedef uint64_t            TNCBI_BigCount;
#  define NCBI_BIGCOUNT_FORMAT_SPEC      PRIu64
#  define NCBI_BIGCOUNT_FORMAT_SPEC_HEX  PRIx64
#endif


#ifdef __cplusplus
}  /* extern "C" */
#endif


/* @} */

#endif /* CONNECT___NCBI_TYPES__H */
