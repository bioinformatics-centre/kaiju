/*  $Id: ncbi_skew_guard.h 346326 2011-12-06 15:28:48Z ucko $
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
 * Author:  Aaron Ucko, NCBI
 *
 */

/** @file ncbi_skew_guard.h
  * Implementation header to catch build setups that mix incompatible
  * C and C++ Toolkit installations.
  * 
  * Available as <common/ncbi_skew_guard.h> on the C++ side and
  * <ncbi_skew_guard.h> on the C side, in each case customized upon
  * installation to identify itself appropriately.
  */

/* In-house C++ Toolkit installations define NCBI_INSTALLED_CXX_VER to
 * the corresponding date stamp (also available as NCBI_DEVELOPMENT_VER
 * or NCBI_PRODUCTION_VER). */
/* #undef NCBI_INSTALLED_CXX_VER */

/* Accompanying copies of the C Toolkit define NCBI_EXPECTED_CXX_VER
 * accordingly. */
/* #undef NCBI_EXPECTED_CXX_VER */

#if defined(_NCBILCL_)  &&  defined(FORWARDING_NCBICONF_H) \
    &&  !defined(NCBI_ALLOW_MISMATCHED_VERSIONS)

/* The last change to shared headers before this guard came along occurred
 * on Nov. 30, 2011. */
#define NCBI_MIN_CXX_VER 20111130

#  if defined(NCBI_INSTALLED_CXX_VER)

#    if !defined(NCBI_EXPECTED_CXX_VER) \
        ||  NCBI_INSTALLED_CXX_VER != NCBI_EXPECTED_CXX_VER
#      error Please use the C Toolkit installation accompanying your C++ Toolkit tree.
#    endif

#  else

#    include <common/ncbi_source_ver.h>
#    if NCBI_DEVELOPMENT_VER < NCBI_MIN_CXX_VER
#      error Please use a fresher C++ Toolkit version for C Toolkit compatibility.
#    elif defined(NCBI_EXPECTED_CXX_VER)
#      if (defined(NCBI_PRODUCTION_VER) ? NCBI_PRODUCTION_VER \
           : NCBI_DEVELOPMENT_VER) \
          != NCBI_EXPECTED_CXX_VER
#        error Please use matching C and C++ Toolkit versions.
#      endif
#    endif

#  endif

#endif
