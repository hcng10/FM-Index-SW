 //
// Created by hcng on 7/4/18.
//

#ifndef FM_SA_BUILDFM_H
#define FM_SA_BUILDFM_H


//#define CEIL(a, b) (((a) + (b) - 1) / (b))

#include <omp.h>
#include <math.h>
#include <cmath>
#include <math.h>

#include "divsufsort64/include/divsufsort64.h"
#include "../fm-sa-def.h"
#include "file_op.h"


struct BP_t{
    char bp[FM_STEP];
};


// pre-computed intervals
struct ival_t {
    uint64_t low;
    uint64_t high;
};


// index structure
/*struct index_t {
    uint64_t count[4];
   // uint8_t bwt[CEIL(BUCKET_SIZE*5,8)];
    uint64_t pad;
};*/

template <class T_cnt>
class Index_t {

public:
    //T determines if it is 32 or 64 bit
    T_cnt count[FM_I_NUM];
    uint8_t * bwt;

    Index_t(){};

    ~Index_t(){
        if (bwt){
            delete[] bwt;
        }
    }

    bool initalize(uint32_t bwt_length_uint8){
        bwt = new uint8_t[bwt_length_uint8];
        memset(bwt, 0, bwt_length_uint8 * sizeof(uint8_t));
    }

};



//uint64_t sumArray(uint64_t *a, int size);
//uint32_t bitCounting(uint64_t fmt_len);
void fmtToIdx(FILE *FM_fp, FILE * FM_meta_fp, FILE * SA_ref_fp, char *fmt, uint64_t fmt_len, uint32_t bucket_size, bool wpad);



#endif //FM_SA_H
