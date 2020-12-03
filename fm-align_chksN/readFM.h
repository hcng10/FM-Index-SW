//
// Created by hcng on 9/19/18.
//

#ifndef FM_SA_ALIGN_READFM_H
#define FM_SA_ALIGN_READFM_H

#define BWT32_LEN_BYTE ((BUCKET_SIZE * 8 - 32 * FM_BP_RANGE) / 8)
#define PAD32_LEN_BYTE ((FM_BP_RANGE * 32 + BWT32_LEN_BYTE * 8) >= (BUCKET_SIZE * 8) ? (0) : (BUCKET_SIZE * 8 - FM_BP_RANGE * 32 - BWT32_LEN_BYTE * 8 )/8)

#define BWT64_LEN_BYTE (BUCKET_SIZE - 64 * FM_BP_RANGE / 8)
#define PAD64_LEN_BYTE ((FM_BP_RANGE * 64 / 8 + BWT64_LEN_BYTE) >= (BUCKET_SIZE) ? (0) : (BUCKET_SIZE - FM_BP_RANGE * 64 / 8 - BWT64_LEN_BYTE))

#include "file_op.h"


// Everything has to be well defined because it has to be consecutive in memory
struct index32_t {
    uint32_t count[FM_BP_RANGE];
    uint8_t bwt[BWT32_LEN_BYTE];
};


/**
    Read the meta data of the index

    @param  *FM_meta_fp         2nd file in the index, contains the ref for backward search
    @param  *c32                (pointer) bool to indicate if 32bit is enough to hold the BWT of reference
    @param  *fmt_len            (pointer) the length of reconstructed fmt
    @param  *bucket_bwt_len     (pointer) the length of BWT in a bucket
    @param  *wpad               (pointer) bool to make the length of BWT within a bucket a multiple of 2. 
                                This requires padding but forgos the use of divider.
    @param  *bucket_pad_size    (pointer) The padding size, if wpad == true
    @param  *end_char_pos       (pointer) The pos of $ in the BWT.
    @param  *end_char           (pointer) Array that contains $X info, only useful for multi-step FM-index
    @param  *N_cluster          (pointer) The number of N char cluster
    @param  *chrs_num           (pointer) The number of Chromosome
*/
void read_meta(FILE * FM_meta_fp,
                bool * c32,
                uint64_t * fmt_len, 
                uint32_t * bucket_bwt_len,
                bool * wpad,
                uint32_t * bucket_pad_size, 
                uint64_t * end_char_pos, 
                uint8_t * end_char,
                uint64_t * N_cluster,
                uint16_t * chrs_num);
                
uint32_t getOcc(uint8_t sym, uint8_t *bwt, uint32_t s_idx, uint32_t e_idx);
void revComp(char *dest, char *src, uint8_t len);

#endif //FM_SA_ALIGN_READFM_H



/*struct index64_t {
    uint64_t count[FM_BP_RANGE];
    uint8_t bwt[BWT64_LEN_BYTE];

#if PAD64_LEN_BYTE != 0
    uint8_t pad[PAD64_LEN_BYTE];
#endif
};*/