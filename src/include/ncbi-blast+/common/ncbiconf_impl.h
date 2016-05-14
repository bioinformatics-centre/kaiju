#ifndef COMMON___NCBICONF_IMPL__H
#define COMMON___NCBICONF_IMPL__H

/* $Id: ncbiconf_impl.h 457074 2015-01-20 16:19:10Z ucko $
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
 *  Author:  Anton Lavrentiev
 *
 *
 */

/**
 * @file ncbiconf_impl.h
 *
 * Configuration macros.
 */

#ifndef FORWARDING_NCBICONF_H
#  error "The header can be used from <ncbiconf.h> only."
#endif /*!FORWARDING_NCBICONF_H*/


/** @addtogroup Portability
 *
 * @{
 */


/* Threads configuration
 */

#undef NCBI_NO_THREADS
#undef NCBI_THREADS
#undef NCBI_POSIX_THREADS
#undef NCBI_WIN32_THREADS

#if defined(_MT)  &&  !defined(NCBI_WITHOUT_MT)
#  if defined(NCBI_OS_MSWIN)
#    define NCBI_WIN32_THREADS
#  elif defined(NCBI_OS_UNIX)
#    define NCBI_POSIX_THREADS
#  else
#    define NCBI_NO_THREADS
#  endif
#else
#  define NCBI_NO_THREADS
#endif

#if !defined(NCBI_NO_THREADS)
#  define NCBI_THREADS
#endif

/* Sync Windows/Cygwin preprocessor conditionals governing wide
 * character usage. */

#if defined(UNICODE)  &&  !defined(_UNICODE)
#  define _UNICODE 1
#elif defined(_UNICODE)  &&  !defined(UNICODE)
#  define UNICODE 1
#endif

/* New/nonstandard keywords
 */

#if defined(__cplusplus)  &&  defined(NCBI_RESTRICT_CXX)
#  define NCBI_RESTRICT NCBI_RESTRICT_CXX
#elif !defined(__cplusplus)  &&  defined(NCBI_RESTRICT_C)
#  define NCBI_RESTRICT NCBI_RESTRICT_C
#elif __STDC_VERSION__ >= 199901 /* C99 specifies restrict */
#  define NCBI_RESTRICT restrict
#else
#  define NCBI_RESTRICT
#endif

#ifndef NCBI_FORCEINLINE
#  ifdef __cplusplus
#    define NCBI_FORCEINLINE inline
#  else
#    define NCBI_FORCEINLINE
#  endif
#endif

#ifndef NCBI_NORETURN
#  ifdef __GNUC__
#    define NCBI_NORETURN __attribute__((__noreturn__))
#  else
#    define NCBI_NORETURN
#  endif
#endif

/* Definition of packed enum type, to save some memory */
/* enum EMyEnum NCBI_PACKED_ENUM_TYPE(Type) { ... } NCBI_PACKED_ENUM_END(); */
#ifndef NCBI_PACKED_ENUM_TYPE
#  define NCBI_PACKED_ENUM_TYPE(type)
#endif
#ifndef NCBI_PACKED_ENUM_END
#  ifdef NCBI_PACKED
#    define NCBI_PACKED_ENUM_END() NCBI_PACKED
#  else
#    define NCBI_PACKED_ENUM_END()
#  endif
#endif

#ifndef NCBI_WARN_UNUSED_RESULT
#  define NCBI_WARN_UNUSED_RESULT
#endif

#ifdef __cplusplus
#  if __cplusplus >= 201103L || defined(__GXX_EXPERIMENTAL_CXX0X__) \
      || defined(__GXX_EXPERIMENTAL_CPP0X__)
#    define NCBI_HAVE_CXX11 1
#  endif
#  if defined(NCBI_HAVE_CXX11) \
      ||  (defined(NCBI_COMPILER_MSVC)  &&  _MSC_VER >= 1600)
#    define HAVE_IS_SORTED 1
#    define HAVE_NULLPTR 1
#  endif
#  if defined(NCBI_HAVE_CXX11) /* or recent MSVC too? */
#    if !defined(NCBI_COMPILER_ICC)  ||  NCBI_COMPILER_VERSION >= 1400
       /* Exclude ICC 13.x and below, which don't support using "enum class"
        * in conjunction with switch. */
#      define HAVE_ENUM_CLASS 1
#    endif
#  endif
#endif

#include <common/ncbi_skew_guard.h>


/* @} */

#endif  /* COMMON___NCBICONF_IMPL__H */
