#ifndef FORWARDING_NCBICONF_H
#define FORWARDING_NCBICONF_H

/*  $Id: ncbiconf.h 485953 2015-11-30 17:53:34Z blastadm $
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
 * Authors: Denis Vakatov, Aaron Ucko
 *
 */

/** @file ncbiconf.h
 ** Front end for a platform-specific configuration summary.
 **/

#ifdef _MSC_VER
#  include <common/config/ncbiconf_msvc.h>
#elif defined(NCBI_XCODE_BUILD)
#  include <common/config/ncbiconf_xcode.h>
#else
#  include <ncbiconf_unix.h>
#endif

#ifdef NCBI_UNIVERSAL_BUILD
/* sort out the remaining details */
#  include <common/config/ncbiconf_universal.h>
#endif

#include <common/ncbiconf_impl.h>

#endif  /* FORWARDING_NCBICONF_H */
