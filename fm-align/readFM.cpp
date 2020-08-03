//
// Created by hcng on 9/19/18.
//

#include "readFM.h"


void read_meta(FILE * FM_meta_fp,
                bool * c32,
                uint64_t * fmt_len, 
                uint32_t * bucket_bwt_len,
                bool * wpad,
                uint32_t * bucket_pad_size, 
                uint64_t * end_char_pos, uint8_t * end_char) {

    //uint64_t fmt_len = 0;

    uint32_t BP_bit;
    uint32_t BP_range;
    uint32_t bucket_size;
    uint8_t fm_step;

    if (FM_STEP > 1){
        readFile(FM_meta_fp, &fm_step, sizeof(uint8_t));

        if (fm_step != FM_STEP) {
            fprintf(stderr, "error: FM_STEP is different from the index!\n");
            exit(1);
        }

        for (uint8_t s = 0; s < FM_STEP; s++){
            readFile(FM_meta_fp, end_char_pos + s, sizeof(uint64_t));
            cout<<"The val64 "<< (int)end_char_pos[s]<<"\n";
        }

        for (uint8_t s = 0; s < FM_STEP * FM_STEP; s++){
            readFile(FM_meta_fp, end_char + s, sizeof(uint8_t));
            cout<<"The val8 "<< (int)end_char[s]<<"\n";
        }
    }else{
        readFile(FM_meta_fp, end_char_pos, sizeof(uint64_t));
    }

    readFile(FM_meta_fp, fmt_len, sizeof(uint64_t));

    readFile(FM_meta_fp, &BP_bit, sizeof(uint32_t));
    //check if the meta-data is the same as the index
    if (FM_BP_BIT != BP_bit) {
        fprintf(stderr, "error: FM_BP_BIT is different from the index!\n");
        exit(1);
    }

    readFile(FM_meta_fp, &BP_range, sizeof(uint32_t));
    //check if the meta-data is the same as the index
    if (FM_BP_RANGE != BP_range) {
        fprintf(stderr, "error: FM_BP_RANGE is different from the index!\n");
        exit(1);
    }

    readFile(FM_meta_fp, &bucket_size, sizeof(uint32_t));
    //check if the meta-data is the same as the index
    if (BUCKET_SIZE != bucket_size) {
        fprintf(stderr, "error: BUCKET_SIZE is different from the index!\n");
        exit(1);
    }

    readFile(FM_meta_fp, c32, sizeof(bool));
    readFile(FM_meta_fp, bucket_bwt_len, sizeof(uint32_t));
    readFile(FM_meta_fp, bucket_pad_size, sizeof(uint32_t));

    readFile(FM_meta_fp, wpad, sizeof(bool));

    printf("The length of the reference %ld:\n", *fmt_len);
    printf("Bit required to store i-array: %d \n", *c32? 32: 64);

    printf("\t Size of bucket %d: \n", bucket_size);
    printf("\t The length of BWT within a bucket: %d \n", *bucket_bwt_len );
    printf("\t The padding size within a bucket: %d \n", * bucket_pad_size );
    
}

// get value in packed bwt
inline uint8_t getVal(uint8_t *bwt, uint32_t idx)
{
    uint8_t mod = idx % (8/FM_BP_BIT);
    uint8_t tmp = bwt[idx/(8/FM_BP_BIT)];
    tmp = (tmp >> (mod*FM_BP_BIT)) & 0x3;

    return tmp;
}


uint32_t getOcc(uint8_t sym, uint8_t *bwt, uint32_t s_idx, uint32_t e_idx)
{
    uint32_t cnt = 0;

    for (uint32_t i = s_idx; i < e_idx; i++) {
        uint8_t bwt_sym = getVal(bwt, i);
        if (sym == bwt_sym) {
            cnt += 1;

        }
    }

    return cnt;
}

// reverse complement read
void revComp(char *dest, char *src, uint8_t len)
{
    for (uint8_t i = 0; i < len; i++) {
        switch (src[len-1-i]) {
            case 'A': dest[i] = 'T'; break;
            case 'C': dest[i] = 'G'; break;
            case 'G': dest[i] = 'C'; break;
            case 'T': dest[i] = 'A'; break;
            default : dest[i] = 'A';
        }
    }
}

