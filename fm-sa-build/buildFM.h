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

/**
    This function converts the reconsturcted ref into compressed FM-index with steps: 
    (1) Uses divsufsort64 to generate the suffix array (SA). 
    (2) Generate BWT (char format) using SA by going back to the reconsturcted ref. 
    (3) Make BWT char represented by uint8_t int. 
    (4) Record all the $ pos 
    (5) Set i(x) 
    (6) Calculate the number of buckets needed and the bwt_length in each bucket 
    (7) Set the bucket val 

    @param  *FM_fp          fd to write the compressed FM-index
    @param  *FM_meta_fp     fd for the meta file of the index
    @param  *SA_ref_fp      fd for storing the SA
    @param  *fmt            the reconstructed fmt
    @param  *fmt_len        the length of reconstructed fmt
    @param  bucket_size     the size of a bucket
    @param  wpad            bool to make the length of BWT within a bucket a multiple of 2. 
                            This requires padding but forgos the use of divider.

*/
void fmtToIdx(FILE *FM_fp, FILE * FM_meta_fp, FILE * SA_ref_fp, char *fmt, uint64_t fmt_len, uint32_t bucket_size, bool wpad);

#endif //FM_SA_H
