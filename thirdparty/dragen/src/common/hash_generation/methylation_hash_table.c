#include "methylation_hash_table.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern const char seqnameSuffixCtoT[];
extern const char seqnameSuffixGtoA[];

//-------------------------------------------------------------------------------
// Comparator for sorting refRec_t by seqNum.
int refRecCompareNum(const void* a, const void* b)
{
  return (((refRec_t*)a)->seqNum) - (((refRec_t*)b)->seqNum);
}

//-------------------------------------------------------------------------------
// Doubles the amount of refRecs pointed to by #refRecs# and updates their name
// based on the parity of their position within the array. Even contigs are C->T,
// odd G->A. Assumes that #numRefSeqs# is the number of refRec_t's pointed to,
// and that there is room allocated for 2*#numRefSeqs# refRecs_t's.
// Updates #totRefLen# to take the new contigs and padding into account..
void UpdateContigsForMethylation(refRec_t* refRecs, uint32_t numRefSeqs, uint64_t* totRefLen)
{
  // Zip the refrecs array with itself.
  uint32_t i = 0;
  for (i = 0; i < numRefSeqs; i++) {
    refRecs[i + numRefSeqs] = refRecs[i];
  }
  qsort(refRecs, 2 * numRefSeqs, sizeof(refRec_t), refRecCompareNum);

  // Now we have refRecs "zipped" with itself. ie each entry appears twice in the original order.
  // The first of each pair will be C->T, the second G->A. But first we need to update some things
  for (i = 0; i < numRefSeqs; i++) {
    const int CT = (2 * i);
    const int GA = (2 * i) + 1;

    // Update the contig numbers
    refRecs[CT].seqNum = CT;
    refRecs[GA].seqNum = GA;

    // Update and append the conversion suffixes to the contig names.
    const size_t newNameAlloc = strlen(refRecs[CT].name) + strlen(seqnameSuffixCtoT) + 1;
    refRecs[CT].name          = realloc(refRecs[CT].name, newNameAlloc);
    refRecs[GA].name          = (char*)malloc(newNameAlloc);

    strcpy(refRecs[GA].name, refRecs[CT].name);

    strcat(refRecs[CT].name, seqnameSuffixCtoT);
    strcat(refRecs[GA].name, seqnameSuffixGtoA);
  }

  // The padding on the last pair now need to be changed. They both have the padding of the original last
  // contig.
  refRec_t* r                = &refRecs[2 * (numRefSeqs - 1)];
  int       padLen           = (r->isPopAlt) ? REF_SEQ_MIN_PAD_BASES_POP_ALT : REF_SEQ_MIN_PAD_BASES;
  int       refSeqAlignBases = (r->isPopAlt) ? REF_SEQ_ALIGN_BASES_POP_ALT : REF_SEQ_ALIGN_BASES;
  padLen += (refSeqAlignBases - ((r->trimLen + padLen) % refSeqAlignBases)) % refSeqAlignBases;

  const uint32_t prevBlockLen = r->blockLen;
  r->endPad                   = padLen;
  r->blockLen                 = r->trimLen + r->endPad;

  (*totRefLen) = 2 * ((*totRefLen) - prevBlockLen) - REF_SEQ_END_PAD_BASES + r->blockLen;

  r      = &refRecs[2 * numRefSeqs - 1];
  padLen = REF_SEQ_MIN_PAD_BASES;
  padLen += (REF_SEQ_ALIGN_BASES - (((*totRefLen) + r->trimLen + padLen) % REF_SEQ_ALIGN_BASES)) %
            REF_SEQ_ALIGN_BASES;
  r->endPad   = padLen;
  r->blockLen = r->trimLen + r->endPad;
  (*totRefLen) += r->blockLen;
}

//-------------------------------------------------------------------------------
// Update the qid and rid of the liftover records for the methylated converted contigs.
// Appends a copy of #liftRecs# to the end of itself. The first half are the liftovers
// for the CT converted contigs (even), and the second half for the GA converted contigs
// (odd). (Unlike the refRecs, the order of #liftRecs# does not mattter).
// Every liftRec_t with qid=n and rid=m, wll be replaced with two liftRec_t's with qid=2n rid=2m and
// qid=2n+1 rid=2m+1.

