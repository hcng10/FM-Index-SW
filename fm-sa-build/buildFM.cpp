//
// Created by hcng on 7/4/18.
//
#include <iostream>
#include <bitset>


#include "buildFM.h"
//#include "divsufsort64/include/divsufsort64.h"


using namespace std;

template <typename T_width>
inline void bucket_setVal(uint8_t *bwt, uint8_t val, uint32_t bucket_bwt_len, T_width i);

// sum array
uint64_t sumArray(uint64_t *a, int size)
{
    uint64_t cnt = 0;

    for (int i = 0; i < size; i++)
        cnt += a[i];

    return cnt;
}

uint32_t bitCounting(uint64_t myNum){

    if (myNum > (uint64_t) (pow(2, 32)-1)){
        return 64;
    }else{
        return 32;
    }
}


void write_meta(FILE * FM_meta_fp, 
                bool c32, 
                uint64_t fmt_len, 
                BP_t * bwt,
                uint32_t bucket_bwt_len, 
                bool wpad,
                uint32_t pad_size, 
                uint64_t* end_char_pos,
                uint8_t * endChar ){

    if (FM_STEP > 1){
        uint8_t fm_step = FM_STEP;

        writeFile(FM_meta_fp, &fm_step, sizeof(uint8_t));
        
        for (uint8_t s = 0; s < FM_STEP; s++){
            writeFile(FM_meta_fp, end_char_pos + s, sizeof(uint64_t));
        }

        for (uint8_t s = 0; s < FM_STEP * FM_STEP; s++){
            writeFile(FM_meta_fp, endChar + s, sizeof(uint8_t));
        }
    }else{
        writeFile(FM_meta_fp, end_char_pos, sizeof(uint64_t));
    }

    writeFile(FM_meta_fp, &fmt_len, sizeof(uint64_t));

    uint32_t BP_bit = FM_BP_BIT;
    writeFile(FM_meta_fp, &BP_bit, sizeof(uint32_t));

    uint32_t BP_range = FM_BP_RANGE;
    writeFile(FM_meta_fp, &BP_range, sizeof(uint32_t));

    uint32_t bucket_size = BUCKET_SIZE;
    writeFile(FM_meta_fp, &bucket_size, sizeof(uint32_t));

    writeFile(FM_meta_fp, &c32, sizeof(bool));

    writeFile(FM_meta_fp, &bucket_bwt_len, sizeof(uint32_t));


    writeFile(FM_meta_fp, &pad_size, sizeof(uint32_t));
    writeFile(FM_meta_fp, &wpad, sizeof(bool));

}


/**
    First it determines the offset inside the bucket bwt. Then it mod i by 4 bits,
    as the BWT uses uint8_t array to hold the compressed BWT. It needs to know
    how many bits should be shifted


    @param  *bwt    to be assigned compressed BWT within a bucket
    @param  val     current BWT value at i
    @param  bucket_bwt_len    length of a bucket
    @param  i                 current index, 32bit or 64bit
*/
template <typename T_width>
inline void bucket_setVal(uint8_t *bwt, uint8_t val, uint32_t bucket_bwt_len, T_width i) {

    uint32_t bucket_bwt_i = (i % (bucket_bwt_len));
    uint32_t bucket_bwt_offset = bucket_bwt_i * (uint32_t)pow(FM_BP_BIT, FM_STEP) / 8;

    uint32_t bucket_bwt_biti = bucket_bwt_i * (uint32_t)pow(FM_BP_BIT, FM_STEP);

    // 2 cases:  within the 8-bit space or outside the 8-bit space
    uint32_t cur_mod = bucket_bwt_biti % 8;
    uint32_t next_mod = (bucket_bwt_biti + (uint32_t) pow(FM_BP_BIT, FM_STEP)) % 8;

    //uint32_t bucket_bwt_sft;
    uint8_t val_sft;
    uint8_t next_val_sft;

    if (cur_mod < next_mod | 8 % (uint32_t) pow(FM_BP_BIT, FM_STEP) == 0){

        // if it is smaller, that means you don't have to get across two elements
        val_sft = val << cur_mod;//lsb, bwt in smaller position
        bwt[bucket_bwt_offset] |= val_sft;


    }else{
        // if it is bigger, that means u have a situation, need to store it across two elements
        val_sft = val << cur_mod;
        next_val_sft = val >> (8 - cur_mod);

        bwt[bucket_bwt_offset] |= val_sft;
        bwt[bucket_bwt_offset + 1] |= next_val_sft;

    }

}


/**
    With a loop boundary of FM_BP_RANGE^STEP_SIZE + 1, it checks which
    entry of i(x) to increment, by comparing to the permutation: struct BP_t

    @param  i       current i when looping across the BWT
    @param  bwt     compressed BWT
    @param  *bp_prmtn    the precomputed permutation for i(x)
*/
inline uint8_t  bwt_prmut_cmp(uint64_t i, BP_t * bwt, BP_t * bp_prmtn){
    uint8_t rtn_val = 0;

    for (uint8_t p = 0; p < FM_I_NUM + FM_STEP; p++){
#if FM_STEP == 1
        if (bwt[i].bp[0] == bp_prmtn[p].bp[0]){
            rtn_val = p;
            break;
        }
#elif FM_STEP == 2
        if (bwt[i].bp[0] == bp_prmtn[p].bp[0] && bwt[i].bp[1] == bp_prmtn[p].bp[1]){
            rtn_val = p;
            break;
        }
#elif FM_STEP == 3
        if (bwt[i].bp[0] == bp_prmtn[p].bp[0] && bwt[i].bp[1] == bp_prmtn[p].bp[1] && bwt[i].bp[2] == bp_prmtn[p].bp[2]){
            rtn_val = p;
            break;
        }
#endif
    }
    return rtn_val;
}

inline BP_t * gen_bp_prmut(BP_t* bwt, uint64_t* end_char_pos){

    // generate all the permutation
    char bp_end = '$';
    char bp_set[4] = { 'A', 'C', 'G', 'T'};

    BP_t * bp_prmtn = new BP_t[FM_I_NUM + FM_STEP];

#if FM_STEP == 1
    bp_prmtn[0].bp[0] = bp_end;
    for (uint8_t i = 1; i < FM_BP_RANGE+1; i++){
         bp_prmtn[i].bp[0] = bp_set[i-1];
    } 


#elif FM_STEP == 2
    for (uint8_t i = 0; i < FM_BP_RANGE; i++){
        for (uint8_t j = 0; j < FM_BP_RANGE; j++){
            bp_prmtn[i*FM_BP_RANGE + j + FM_STEP].bp[0] = bp_set[j];
            bp_prmtn[i*FM_BP_RANGE + j + FM_STEP].bp[1] = bp_set[i];
        }
    }

    for (uint8_t s = 0; s < FM_STEP; s++ ){
        bp_prmtn[s].bp[0] = bwt[end_char_pos[s]].bp[0];
        bp_prmtn[s].bp[1] = bwt[end_char_pos[s]].bp[1];       
    }


#elif FM_STEP == 3
    for (uint8_t i = 0; i < FM_BP_RANGE; i++){
        for (uint8_t j = 0; j < FM_BP_RANGE; j++){
            for (uint8_t k = 0; k < FM_BP_RANGE; k++){
                bp_prmtn[i*FM_BP_RANGE * FM_BP_RANGE + j * FM_BP_RANGE + FM_STEP].bp[0] = bp_set[j];
                bp_prmtn[i*FM_BP_RANGE * FM_BP_RANGE + j * FM_BP_RANGE + FM_STEP].bp[1] = bp_set[i];
                bp_prmtn[i*FM_BP_RANGE * FM_BP_RANGE + j * FM_BP_RANGE + FM_STEP].bp[2] = bp_set[2];
        }
    }

    for (uint8_t s = 0; s < FM_STEP; s++ ){
        bp_prmtn[s].bp[0] = bwt[end_char_pos[s]].bp[0];
        bp_prmtn[s].bp[1] = bwt[end_char_pos[s]].bp[1];  
        bp_prmtn[s].bp[2] = bwt[end_char_pos[s]].bp[2];      
    }
    

#endif
    return bp_prmtn;
}

void fmtToIdx(FILE *FM_fp, FILE * FM_meta_fp, FILE * SA_ref_fp, char *fmt, uint64_t fmt_len, uint32_t bucket_size, bool wpad){

    int n_threads = omp_get_max_threads();

    // suffix array
    int64_t *sai = NULL;

    // BWT in char and in binary
    BP_t * bwt = NULL;
    uint8_t * bwtM = NULL;

    uint32_t bucket_bwt_len;
    uint32_t n_buckets;

    // bool to indicate if 32bit is enough to hold the BWT of reference
    bool c32;

    // A bucket
    Index_t <u_int64_t> *bucket_indx_c64 = 0;
    Index_t <u_int32_t> *bucket_indx_c32 = 0;

    // Marker for $ pos
    uint64_t end_char_pos[FM_STEP] = {0};

#if FM_STEP == 1
    uint8_t endChar[1] = {0};
#elif FM_STEP == 2 
    uint8_t endChar[FM_STEP*FM_STEP] = {0};
#elif FM_STEP == 3
    uint8_t endChar[FM_STEP*FM_STEP*FM_STEP] = {0};
#endif

    uint32_t pad_size = 0;

    sai = new int64_t [fmt_len];

    // compute suffix array
    if (!sai) {
        fprintf(stderr, "error: unable to allocate memory for suffix array!\n");
        exit(1);
    }

    divsufsort64((uint8_t*)fmt, sai, (int64_t)fmt_len);

    // write SA/ref to a file
    for (int i = 0; i< fmt_len; i++){
        uint32_t sa_val = (uint32_t )sai[i];
        writeFile(SA_ref_fp, &sa_val, sizeof(uint32_t));
    }

    // compute BWT
    bwt = new BP_t [fmt_len];
    bwtM = new uint8_t [fmt_len];


    if (!bwt || !bwtM) {
        fprintf(stderr, "error: unable to allocate memory for BWT!\n");
        exit(1);
    }


#pragma omp parallel for num_threads(n_threads)
    for (uint64_t i = 0; i < fmt_len; i++) {
        uint8_t val = 0;

        //bwt[i] = sai[i] > 0 ? fmt[sai[i]-1] : fmt[fmt_len-1+sai[i]];
        // Obtain the BWT depending on the step size
        for (uint8_t s = 0; s < FM_STEP; s++){
            if (s == 0){
                // Go to the location pointed by SA, get the nucleotide right before
                bwt[i].bp[s] = sai[i] > 0 ? fmt[sai[i]-1-s] : fmt[fmt_len-1-s];
            }
            else if (s == 1){
                bwt[i].bp[s] = sai[i] > 1 ? fmt[sai[i]-1-s] : sai[i] == 1 ? fmt[fmt_len-1]: fmt[fmt_len-2];
            }
            else if (s == 2){
                bwt[i].bp[s] = sai[i] > 2 ? fmt[sai[i]-1-s] : sai[i] == 2 ? fmt[fmt_len-1]: sai[i] == 1 ? fmt[fmt_len-2]: fmt[fmt_len-3];
            }
        }

        // Compression, assign BWT nucleotide with 0 -> 3, random value for $
        for (uint8_t s = 0; s < FM_STEP; s++){
            if (s == 0){
                switch(bwt[i].bp[s]) {
                    case 'A': val = 0;  break;
                    case 'C': val = 1;  break;
                    case 'G': val = 2;  break;
                    case 'T': val = 3;  break;
                    case '$':
                        val += 16;
                        end_char_pos[s] = i;
                        break;
                    default:
                        break;
                }
            // for > 1 step size
            }else{
                switch(bwt[i].bp[s]) {
                    case 'A': val += 0;  break;
                    case 'C': val += 1 * (uint8_t) pow(FM_BP_RANGE, s);  break;
                    case 'G': val += 2 * (uint8_t) pow(FM_BP_RANGE, s);  break;
                    case 'T': val += 3 * (uint8_t) pow(FM_BP_RANGE, s);  break;
                    case '$':
                        val += 16;
                        end_char_pos[s] = i;
                        break;
                    default:
                        break;
                }
            }
        bwtM[i] = val;
        }        
    }

    //record all the char in end char
    for (uint8_t s = 0; s < FM_STEP; s++){
        for (uint8_t si = 0; si < FM_STEP; si++){
            
            uint8_t val = 0;
            switch(bwt[end_char_pos[s]].bp[si]) {
                case 'A': val = 0;  break;
                case 'C': val = 1;  break;
                case 'G': val = 2;  break;
                case 'T': val = 3;  break;
                case '$': val = 4;  break;
                default:
                    break;
            }
            endChar[s * FM_STEP + si] = val;
        }
    }


    // compute i(x) 
    uint64_t cnt[FM_I_NUM + FM_STEP] = {0};

    // generate all the permutation of i(x) 
    BP_t* bp_prmtn = gen_bp_prmut(bwt, end_char_pos);

    for (uint64_t i = 0; i < fmt_len; i++) {
        uint8_t val = 0;

        // get the entry number for i(x) that I need to increment
        val = bwt_prmut_cmp(i, bwt, bp_prmtn);
        cnt[val]++;
    }
    
    //free up some of the memory
    delete bwt;

    // sum counters, need to do from the end
    for (int i = FM_I_NUM + FM_STEP - 1; i > 0; i--)
        cnt[i] = sumArray(cnt, i);
    cnt[0] = 0;

    // remove the first $ counters, because it will always be 0
    for (int i = 0; i < FM_I_NUM; i++){
        cnt[i] = cnt[i+FM_STEP];
    }
    cnt[FM_I_NUM] = fmt_len;

    //DEBUG START:
    //for (int i = 0; i < FM_I_NUM+1; i++){
    //    cout<<"cnt"<<cnt[i]<<"\n";
    //}
    //DEBUG END:
    
    // build FM-index structure
    // 1. calculate the number of bits needed for C Matrix
    uint32_t bitCcnt = bitCounting(fmt_len);
    if (bitCcnt == 32)
        c32 = true;
    else if (bitCcnt == 64)
        c32 = false;

    printf("Bit required to store i-array: %d \n",bitCcnt);

    // 2. calculate the number of buckets needed and the bwt_length in each bucket
    int64_t bucket_size_bit = bucket_size * 8;
    bucket_bwt_len = (bucket_size_bit - FM_I_NUM * bitCcnt) / FM_BP_BIT;

    if (bucket_bwt_len < 0){
        fprintf(stderr, "error: Bucket size is less than i-array!\n");
        exit(1);
    }

    // calculate the closest power of 2
    if (wpad == true) {
        uint32_t bucket_bwt_len_power = (uint32_t) floor(log(float(bucket_bwt_len))/log(2.0));
        bucket_bwt_len = 1 << bucket_bwt_len_power;
    }
    else{
        // pop counter on FPGA can only do one for every 64 entries
        bucket_bwt_len = (bucket_bwt_len / 64 ) * 64;
    }

    printf("\t The length of BWT within a bucket: %d \n", bucket_bwt_len );

    n_buckets = CEIL(fmt_len, bucket_bwt_len);

    uint32_t tmp_bwt_size_byte = (FM_I_NUM * bitCcnt + bucket_bwt_len * FM_BP_BIT) / 8;
    pad_size = (tmp_bwt_size_byte >= BUCKET_SIZE) ? (0) : (BUCKET_SIZE - tmp_bwt_size_byte);

    printf("\t The padding size within a bucket: %d \n", pad_size );

    
    //3. build the index

    //calclate exactly how big the bwt array should be for type uint8_t
    uint32_t bwt_length_uint8 = (bucket_bwt_len * FM_BP_BIT) / 8;
    if ((bucket_bwt_len * FM_BP_BIT) % 8 > 0) bwt_length_uint8++;

    printf("\t Number of buckets %d: \n", n_buckets);
    printf("\t Byte needed for each BTW in a bucket: %d\n", bwt_length_uint8);


    if (c32) {
        // declare and initalize a bucket object
        bucket_indx_c32 = new Index_t<uint32_t>[n_buckets];
        memset(bucket_indx_c32, 0, n_buckets * sizeof(bucket_indx_c32));

        uint32_t cnt_tmp[FM_I_NUM] = {0};
        for (uint32_t i = 0; i < fmt_len; i++) {
            uint32_t bucket_i = i / bucket_bwt_len;

            if (i % (bucket_bwt_len) == 0) {
                // this should prevent the seg fault when accessing
                bucket_indx_c32[bucket_i].initalize(bwt_length_uint8);

                // assign the value for i(x) within a bucket
                for (int j = 0; j < FM_I_NUM; j++){
                    bucket_indx_c32[bucket_i].count[j] = cnt_tmp[j] + cnt[j];
                }

            }

            if (bwtM[i] < 16) {
                // 16 means $ sign
                cnt_tmp[bwtM[i]]++;
            }else{
                // entry with $ sign will always be 0 
                bwtM[i] = 0;
            }

            bucket_setVal<uint32_t>(bucket_indx_c32[bucket_i].bwt, bwtM[i], bucket_bwt_len, i);
        }
    }

    //4. write the index to a file
    for (int i = 0; i< FM_I_NUM + 1; i++){
        writeFile(FM_fp, & cnt[i], sizeof(uint64_t));
    }

    
    for (int i = 0; i < n_buckets; i++){
        // writing the c-matrix
        writeFile(FM_fp, &bucket_indx_c32[i].count, (FM_BP_RANGE * bitCcnt) / 8);
        // writing the bwt
        writeFile(FM_fp, bucket_indx_c32[i].bwt, bwt_length_uint8);

        //write pad
        if (pad_size != 0){
            for (int j = 0; j < pad_size; j++){
                char zero = 0;
                writeFile(FM_fp, &zero, sizeof(char));
            }
        }
    }

    write_meta(FM_meta_fp, 
                c32, fmt_len, 
                bwt, bucket_bwt_len, 
                wpad, pad_size, 
                end_char_pos,
                endChar);

    delete [] bp_prmtn;
    delete [] bwtM;
    delete [] sai;
    c32 ? delete [] bucket_indx_c32: delete [] bucket_indx_c64;
}







