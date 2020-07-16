//
// Created by hcng on 8/9/18.
//

#ifndef FM_SA_H
#define FM_SA_H


#define BUCKET_SIZE 64
#define FM_BP_BIT 2
#define FM_STEP 1
#define FM_BP_RANGE 4

#define SEED_LEN 10
#define TRIM_READ_LEN 16

#if FM_STEP == 1 
#define FM_I_NUM FM_BP_RANGE
#elif FM_STEP == 2 
#define FM_I_NUM FM_BP_RANGE*FM_BP_RANGE
#elif FM_STEP == 3 
#define FM_I_NUM FM_BP_RANGE*FM_BP_RANGE*FM_BP_RANGE
#endif 

#define MAX_READ_LEN 2000
#define BUFF_SIZE 80000000
//800000000

// round up division
#define CEIL(a, b) (((a)+(b)-1)/(b))

#endif //FM_SA_H

//uint32_t bucket_bwt_offset = (i % (bucket_bwt_len)) / (8/bit_for_bp);//this will over count
//uint32_t bucket_bwt_sft = i % (8 / bit_for_bp);

//TODO: The end char $ sign is not considered
//uint8_t val_sft = val << (bucket_bwt_sft * bit_for_bp);//lsb, bwt in smaller position
//bwt[bucket_bwt_offset] |= val_sft;
