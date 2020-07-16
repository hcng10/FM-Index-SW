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

//#if PAD32_LEN_BYTE != 0
//    uint8_t pad[PAD32_LEN_BYTE];
//#endif
};



/*struct index64_t {
    uint64_t count[FM_BP_RANGE];
    uint8_t bwt[BWT64_LEN_BYTE];

#if PAD64_LEN_BYTE != 0
    uint8_t pad[PAD64_LEN_BYTE];
#endif
};*/


void read_meta(FILE * FM_meta_fp,
                bool * c32,
                uint64_t * fmt_len, 
                uint32_t * bucket_bwt_len,
                bool * wpad,
                uint32_t * bucket_pad_size, 
                uint64_t * endCharPos, uint8_t * endChar);
                
uint32_t getOcc(uint8_t sym, uint8_t *bwt, uint32_t s_idx, uint32_t e_idx);
void revComp(char *dest, char *src, uint8_t len);

#endif //FM_SA_ALIGN_READFM_H
