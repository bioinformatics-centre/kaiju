/*  $Id: hspfilter_collector.c 443210 2014-08-12 13:55:24Z fongah2 $
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
 * Author:  Ilya Dondoshansky
 *
 */

/** @file hspfilter_collector.c
 * Default implementation of the BlastHSPWriter interface to save hits from
 * a BLAST search, and subsequently return them in sorted order.
 */

#ifndef SKIP_DOXYGEN_PROCESSING
static char const rcsid[] = 
    "$Id: hspfilter_collector.c 443210 2014-08-12 13:55:24Z fongah2 $";
#endif /* SKIP_DOXYGEN_PROCESSING */


#include <algo/blast/core/hspfilter_collector.h>
#include <algo/blast/core/blast_util.h>
#include "blast_hits_priv.h"

/** Data structure used by the writer */
typedef struct BlastHSPCollectorData {
   BlastHSPCollectorParams* params;         /**< how many hits to save */
   BlastHSPResults* results;                /**< place to store hits */
} BlastHSPCollectorData;

/*************************************************************/
/** The following are implementations for BlastHSPWriter ADT */

/** Perform pre-run stage-specific initialization 
 * @param data The internal data structure [in][out]
 * @param results The HSP results to operate on  [in]
 */ 
static int 
s_BlastHSPCollectorInit(void* data, BlastHSPResults* results)
{
   BlastHSPCollectorData * col_data = data;
   /* grab the results as destination to store collected hsps */
   col_data->results = results;
   return 0;
}

/** Perform post-run clean-ups
 * @param data The buffered data structure [in]
 * @param results The HSP results to propagate [in][out]
 */ 
static int 
s_BlastHSPCollectorFinal(void* data, BlastHSPResults* results)
{
   BlastHSPCollectorData * col_data = data;
   /* results already stored during run, no action needed */
   col_data->results = NULL;
   return 0;
}

/** Perform writing task
 * ownership of the HSP list and sets the dereferenced pointer to NULL.
 * @param data To store results to [in][out]
 * @param hsp_list Pointer to the HSP list to save in the collector. [in]
 */
static int 
s_BlastHSPCollectorRun(void* data, BlastHSPList* hsp_list)
{
   BlastHSPCollectorData * col_data = data;
   BlastHSPResults* results = col_data->results;
   BlastHSPCollectorParams* params = col_data->params;
   EBlastProgramType program = params->program;

   if (!hsp_list)
      return 0;

   if (!results || !params)
      return -1;

   /* The HSP list should already be sorted by score coming into this function.
    * Still check that this assumption is true. Note that HSP list does not need to be
    * sorted after preliminary stage for vdb search becuase of offset adjustment.
    */
#ifdef ERR_POST_EX_DEFINED
   if(!Blast_HSPListIsSortedByScore(hsp_list)) {
             ErrPostEx(SEV_WARNING, 0, 0, "HSP List is not sorted by score");
   }
#endif
   
   /* Rearrange HSPs into multiple hit lists if more than one query */
   if (results->num_queries > 1) {
      BlastHSP* hsp;
      BlastHSPList** hsp_list_array;
      BlastHSPList* tmp_hsp_list;
      Int4 index;

      hsp_list_array = calloc(results->num_queries, sizeof(BlastHSPList*));
      if (hsp_list_array == NULL)
         return -1;

      for (index = 0; index < hsp_list->hspcnt; index++) {
         Int4 query_index;
         hsp = hsp_list->hsp_array[index];
         query_index = Blast_GetQueryIndexFromContext(hsp->context, program);

         if (!(tmp_hsp_list = hsp_list_array[query_index])) {
            hsp_list_array[query_index] = tmp_hsp_list = 
               Blast_HSPListNew(params->hsp_num_max);
            if (tmp_hsp_list == NULL)
            {
                 sfree(hsp_list_array);
                 return -1;
            }
            tmp_hsp_list->oid = hsp_list->oid;
         }

         Blast_HSPListSaveHSP(tmp_hsp_list, hsp);
         hsp_list->hsp_array[index] = NULL;
      }

      /* All HSPs from the hsp_list structure are now moved to the results 
         structure, so set the HSP count back to 0 */
      hsp_list->hspcnt = 0;
      Blast_HSPListFree(hsp_list);

      /* Insert the hit list(s) into the appropriate places in the results 
         structure */
      for (index = 0; index < results->num_queries; index++) {
         if (hsp_list_array[index]) {
            if (!results->hitlist_array[index]) {
               results->hitlist_array[index] = 
                  Blast_HitListNew(params->prelim_hitlist_size);
            }
            Blast_HitListUpdate(results->hitlist_array[index], 
                                hsp_list_array[index]);
         }
      }
      sfree(hsp_list_array);
   } else if (hsp_list->hspcnt > 0) {
      /* Single query; save the HSP list directly into the results 
         structure */
      if (!results->hitlist_array[0]) {
         results->hitlist_array[0] = 
            Blast_HitListNew(params->prelim_hitlist_size);
      }
      Blast_HitListUpdate(results->hitlist_array[0], hsp_list);
   } else {
       /* Empty HSPList - free it. */
       Blast_HSPListFree(hsp_list);
   }
       
   return 0; 
}

/** Callback used for sorting HSPs by score, with HSPs
 *  from different contexts segregated from each other
 */
static int
s_ScoreCompareHSPWithContext(const void* h1, const void* h2)
{
   BlastHSP* hsp1,* hsp2;   /* the HSPs to be compared */
   int result = 0;      /* the result of the comparison */
   
   hsp1 = *((BlastHSP**) h1);
   hsp2 = *((BlastHSP**) h2);

   /* Null HSPs are "greater" than any non-null ones, so they go to the end
      of a sorted list. */
   if (!hsp1 && !hsp2)
       return 0;
   else if (!hsp1)
       return 1;
   else if (!hsp2)
       return -1;

   if ((result = BLAST_CMP(hsp1->context, hsp2->context)) != 0)
       return result;
   return ScoreCompareHSPs(h1, h2);
}

/** Perform writing task for RPS case
 * For RPS BLAST saving procedure is different, because HSPs from different
 * subjects are bundled in one HSP list 
 * ownership of the HSP list and sets the dereferenced pointer to NULL.
 * @param data To store results to [in][out]
 * @param hsp_list Pointer to the HSP list to save in the collector. [in]
 */
static int 
s_BlastHSPCollectorRun_RPS(void* data, BlastHSPList* hsplist_in)
{
   Int4 index, next_index;
   BlastHitList* hit_list;
   BlastHSPCollectorData * col_data = data;
   BlastHSPResults* results = col_data->results;
   BlastHSPCollectorParams* params = col_data->params;

   if (!hsplist_in || hsplist_in->hspcnt == 0)
      return 0;

   /* Check that the query index is in the correct range. */
   ASSERT(hsplist_in->query_index < results->num_queries);

   /* Check that program is indeed RPS Blast */
   ASSERT(Blast_ProgramIsRpsBlast(params->program));

   /* If hit list for this query has not yet been allocated, do it here. */
   hit_list = results->hitlist_array[hsplist_in->query_index];
   if (!hit_list) {
       results->hitlist_array[hsplist_in->query_index] =
           hit_list = Blast_HitListNew(params->prelim_hitlist_size);
   }

   /* Sort the input HSPList with context (i.e. oid) as the first priority, 
      and then score, etc. */
   qsort(hsplist_in->hsp_array, hsplist_in->hspcnt, sizeof(BlastHSP*), 
         s_ScoreCompareHSPWithContext);

   /* Sequentially extract HSPs corresponding to one subject into a new 
      HSPList, and save these new HSPLists in a normal way, as in all other
      BLAST programs. */
   next_index = 0;
   
   for (index = 0; index < hsplist_in->hspcnt; index = next_index) {
       BlastHSPList* hsp_list;
       Int4 oid = hsplist_in->hsp_array[index]->context;
       Int4 hspcnt;
       /* Find the first HSP that corresponds to a different subject. 
          At the same time, set all HSP contexts to 0, since this is what
          traceback code expects. */
       for (next_index = index; next_index < hsplist_in->hspcnt; 
            ++next_index) {
           if (hsplist_in->hsp_array[next_index]->context != oid)
               break;
           hsplist_in->hsp_array[next_index]->context = 0;
       }
       hspcnt = next_index - index;
       hsp_list = Blast_HSPListNew(hspcnt);
       /* Set the oid field for this HSPList. */
       hsp_list->oid = oid;
       hsp_list->query_index = hsplist_in->query_index;
       /* Save all HSPs corresponding to this subject. */
       for ( ; index < next_index; ++index)
           Blast_HSPListSaveHSP(hsp_list, hsplist_in->hsp_array[index]);
       /* Check that HSPs are correctly sorted by score, as they should be. */
       ASSERT(Blast_HSPListIsSortedByScore(hsp_list));
       /* Insert this HSPList into this query's hit list. */
       Blast_HitListUpdate(hit_list, hsp_list);
   }

   /* All HSPs have been moved from the input HSPList to new HSPLists, so
      set the input HSPList's count to 0. */
   hsplist_in->hspcnt = 0;
   Blast_HSPListFree(hsplist_in);
   
   return 0;
}

/** Free the writer 
 * @param writer The writer to free [in]
 * @return NULL.
 */
static
BlastHSPWriter*
s_BlastHSPCollectorFree(BlastHSPWriter* writer) 
{
   BlastHSPCollectorData *data = writer->data;
   sfree(data->params); 
   sfree(writer->data);
   sfree(writer);
   return NULL;
}

/** create the writer
 * @param params Pointer to the hit paramters [in]
 * @param query_info BlastQueryInfo (not used) [in]
 * @return writer
 */
static
BlastHSPWriter* 
s_BlastHSPCollectorNew(void* params, BlastQueryInfo* query_info)
{
   BlastHSPWriter * writer = NULL;
   BlastHSPCollectorData * data = NULL;
   BlastHSPCollectorParams * col_param = params;

   /* allocate space for writer */
   writer = malloc(sizeof(BlastHSPWriter));

   /* fill up the function pointers */
   writer->InitFnPtr   = &s_BlastHSPCollectorInit;
   writer->FinalFnPtr  = &s_BlastHSPCollectorFinal;
   writer->FreeFnPtr   = &s_BlastHSPCollectorFree;
   writer->RunFnPtr    = (Blast_ProgramIsRpsBlast(col_param->program)) 
                       ? &s_BlastHSPCollectorRun_RPS
                       : &s_BlastHSPCollectorRun;

   /* allocate for data structure */
   writer->data = malloc(sizeof(BlastHSPCollectorData));
   data = writer->data;
   data->params = params;
    
   return writer;
}

/*************************************************************/
/** The following are exported functions to be used by APP   */

BlastHSPCollectorParams*
BlastHSPCollectorParamsNew(const BlastHitSavingOptions* hit_options,
                           Int4 compositionBasedStats,
                           Boolean gapped_calculation)
{
       BlastHSPCollectorParams* retval=NULL;
       Int4 prelim_hitlist_size;

       if (hit_options == NULL)
           return NULL;

       retval = (BlastHSPCollectorParams*) malloc(sizeof(BlastHSPCollectorParams));

       prelim_hitlist_size = hit_options->hitlist_size;
       if (compositionBasedStats)
            prelim_hitlist_size = prelim_hitlist_size * 2 + 50;  
       else if (gapped_calculation)
            prelim_hitlist_size = MIN(2 * prelim_hitlist_size, 
                                      prelim_hitlist_size + 50);

       retval->prelim_hitlist_size = MAX(prelim_hitlist_size, 10);
       retval->hsp_num_max = BlastHspNumMax(gapped_calculation, hit_options);
       retval->program = hit_options->program_number;
       return retval;
}

BlastHSPCollectorParams*
BlastHSPCollectorParamsFree(BlastHSPCollectorParams* opts)
{
    if ( !opts )
        return NULL;
    sfree(opts);
    return NULL;
}

BlastHSPWriterInfo*
BlastHSPCollectorInfoNew(BlastHSPCollectorParams* params) {
   BlastHSPWriterInfo * writer_info =
                        malloc(sizeof(BlastHSPWriterInfo)); 
   writer_info->NewFnPtr = &s_BlastHSPCollectorNew;
   writer_info->params = params;
   return writer_info;
}
