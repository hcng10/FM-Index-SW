//
// Created by hcng on 6/29/18.
//

#ifndef FM_SA_FORMAT_H
#define FM_SA_FORMAT_H


#include <cstdlib>
#include <omp.h>

#include "file_op.h"
#include "../fm-sa-def.h"

#define LINE_SIZE 512

void faToFmt(FILE *fasta_fp, FILE * N_fp, char *fmt, uint8_t bit_for_bp, uint64_t & fmt_cnt,
             uint64_t & ref_cnt, uint64_t & N_cluster); //uint64_t fa_len
void packSymbols(char *sym, uint8_t *pck, uint8_t len);

#endif //FM_SA_H
