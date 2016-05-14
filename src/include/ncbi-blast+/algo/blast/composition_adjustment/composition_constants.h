/* $Id: composition_constants.h 187049 2010-03-26 14:52:29Z satskyse $
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
 * ===========================================================================*/
/**
 * @file composition_constants.h
 * Constants used in compositional score matrix adjustment
 *
 * @author E. Michael Gertz, Alejandro Schaffer, Yi-Kuo Yu
 */


#ifndef __COMPOSITION_CONSTANTS__
#define __COMPOSITION_CONSTANTS__

#include <algo/blast/core/ncbi_std.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Minimum score in a matrix */
#define COMPO_SCORE_MIN INT2_MIN

/** Number of standard amino acids */
#define COMPO_NUM_TRUE_AA 20

/** The largest alphabet supported by this code (the code supports 26
    or 28 character amino acid alphabets). Used to specify the size of
    structures that are statically allocated. */
#define COMPO_LARGEST_ALPHABET 28

/* NOTE: Please keep these comments in sync with argument descriptions in
 * CCompositionBasedStatsArgs::SetArgumentDescriptions()
 */

/** An collection of constants that specify all permissible
 * modes of composition adjustment */
typedef enum ECompoAdjustModes {
    /** Don't use composition based statistics */
    eNoCompositionBasedStats       = 0, 
    /** Composition-based statistics as in NAR 29:2994-3005, 2001 */
    eCompositionBasedStats         = 1, 
    /** Composition-based score adjustment as in Bioinformatics 21:902-911,
     * 2005, conditioned on sequence properties. Cannot be applied to PSSMs. */
    eCompositionMatrixAdjust       = 2, 
    /** Composition-based score adjustment as in Bioinformatics 21:902-911,
     * 2005, unconditionally. Cannot be applied to PSSMs. */
    eCompoForceFullMatrixAdjust    = 3,
    eNumCompoAdjustModes
} ECompoAdjustModes;


/** An collection of constants that specify all rules that may
 *  be used to generate a compositionally adjusted matrix.  */
typedef enum EMatrixAdjustRule {
    eDontAdjustMatrix              = (-1),
    eCompoScaleOldMatrix           = 0,
    eUnconstrainedRelEntropy       = 1,
    eRelEntropyOldMatrixNewContext = 2,
    eRelEntropyOldMatrixOldContext = 3,
    eUserSpecifiedRelEntropy       = 4
} EMatrixAdjustRule;


#ifdef __cplusplus
}
#endif

#endif
