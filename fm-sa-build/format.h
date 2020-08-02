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

/**
    Format the input fasta file, compress the reference using 
    2 or 3 bits per nucleotide

    @param  *fasta_fp    input fd for the fasta file
    @param  *N_fp        output fd for the file that contains N nucleotide info
    @param  *fmt         store the input reference, without compression
    @param  bit_for_bp   Bits needed to store a nucleotide
    @param  &fmt_cnt     use for counting the length of the reference 
                         with N nucleotide (reconstructed ref)
    @param  &ref_cnt     to count the length of the reference
    @param  &N_cluster   to count the number of N nucleotide cluster
*/
void faToFmt(FILE *fasta_fp, FILE * N_fp, char *fmt, uint8_t bit_for_bp, uint64_t & fmt_cnt,
             uint64_t & ref_cnt, uint64_t & N_cluster); //uint64_t fa_len

/**
    Use 1 byte (8 bit) to pack 4 nucleotides. With Open_MP, we iterate every
    nucleotide in the reconstructed ref. Then we shift and append the value 
    within 1 byte until we put 4 nucleotide in there

    @param  *sym    reconstructed ref
    @param  *pck    array to store compressed reconstructed ref
    @param  *len         len of reconstructed ref
*/
void packSymbols(char *sym, uint8_t *pck, uint8_t len);

#endif //FM_SA_H
