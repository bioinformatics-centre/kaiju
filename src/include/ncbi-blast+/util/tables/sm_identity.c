/*  $Id: sm_identity.c 458581 2015-02-06 15:18:12Z boratyng $
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
* Author:  Greg Boratyn
*
* File Description:
*   Protein alignment score matrices; shared between the two toolkits.
*
* ===========================================================================
*/

#include <util/tables/raw_scoremat.h>

/** Entries for the IDENTITY matrix. */

static const TNCBIScore s_IdentityPSM[25 * 25] = {
    /*       A,  R,  N,  D,  C,  Q,  E,  G,  H,  I,  L,  K,  M,
             F,  P,  S,  T,  W,  Y,  V,  B,  J,  Z,  X,  *        */ 
    /*A*/    9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*R*/   -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*N*/   -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*D*/   -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*C*/   -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*Q*/   -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*E*/   -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*G*/   -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*H*/   -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*I*/   -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*L*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*K*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*M*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*F*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
             9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*P*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*S*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /*T*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5, -5,
    /*W*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5, -5,
    /*Y*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5, -5,
    /*V*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5, -5,
    /*B*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5, -5,
    /*J*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5, -5,
    /*Z*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5,  9, -5, -5,
    /*X*/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
    /***/   -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  9
};
const SNCBIPackedScoreMatrix NCBISM_Identity = {
    "ARNDCQEGHILKMFPSTWYVBJZX*",
    s_IdentityPSM,
    -5
};

