#ifndef FM_SA_ALIGN_HPP
#define FM_SA_ALIGN_HPP

#include <omp.h>
#include <vector>

#include "reads.h"
#include "readFM.h"



/**
    Perform exact match alignment

    @param  *index              (pointer) index
    @param  *idx                (pointer) the compressed FM-index in consecutive memory space
    @param  *cnt                (pointer) i(x)
    @param  *sym                (pointer) read symbols
    @param  read_len            read length. 
                                This requires padding but forgos the use of divider.
    @param  *read_low           (pointer) To store up the final low
    @param  *read_high          (pointer) To store up the final high
    @param  bucket_bwt_len      Length of the BWT within a bucket
    @param  end_char_bucket     bucket index that contains $
    @param  end_char_bucketi    the offset withing a bucket for $
    @param  backward            bool for backward search

*/
template <class index_t, class cnt_t>
void exactSearch(index_t * index,  
                cnt_t * cnt, 
                char *sym, 
                uint8_t read_len, 
                bool *is_align, 
                uint32_t * read_low, 
                uint32_t * read_high, 
                uint32_t bucket_bwt_len, 
                uint32_t end_char_bucket, 
                uint32_t end_char_bucketi,
                bool backward) {

    uint32_t low = 0;
    uint32_t high = cnt[4];


    // get read symbol
    uint8_t val;
    for (uint8_t i = read_len; i > 0; i--) {

        uint8_t sym_idx = backward == 1? i - 1: read_len - i;
        // read from the end
        if (backward == 1){
            switch (sym[sym_idx]) {
                case 'A' : val = 0; break;
                case 'C' : val = 1; break;
                case 'G' : val = 2; break;
                case 'T' : val = 3; break;
                case 'a' : val = 0; break;
                case 'c' : val = 1; break;
                case 'g' : val = 2; break;
                case 't' : val = 3; break;
                default : val = 0;
            }
        }else{
            // reverse complement, A <-> T, G <-> C
             switch (sym[sym_idx]) {
                case 'A' : val = 3; break;
                case 'C' : val = 2; break;
                case 'G' : val = 1; break;
                case 'T' : val = 0; break;
                case 'a' : val = 3; break;
                case 'c' : val = 2; break;
                case 'g' : val = 1; break;
                case 't' : val = 0; break;
                default : val = 0;
            }           
        }

        //#TODO: Not sure if this is correct
        //uint32_t low_tmp = low == 0 ? 0 : low - 1;
        bool same_bucket = high - low < bucket_bwt_len ? true : false;


        uint32_t low_addr = low / bucket_bwt_len;
        uint32_t low_idx = low % bucket_bwt_len;

        uint32_t low_count = index[low_addr].count[val];
        uint32_t low_bwtcount = getOcc(val, index[low_addr].bwt, 0, low_idx);

        //DEBUG START   
        //cout<< "lowCCount "<<low_count<<" low_bwtcount "<<low_bwtcount<<" low_idx "<<low_idx<<"\n";
        //DEBUG END  

        //handle endchar, as it will add 1 for 'A', because in the index $ and A are both represented as 00
        if (end_char_bucket == low_addr &&  end_char_bucketi < low_idx && val == 0) {
            low_bwtcount--;
        }



        low = low_count + low_bwtcount;

        uint32_t high_addr = high / bucket_bwt_len;
        uint32_t high_idx = high % bucket_bwt_len;

        uint32_t high_count = index[high_addr].count[val];
        uint32_t high_bwtcount = getOcc(val, index[high_addr].bwt, 0, high_idx);

        //handle endchar, as it will add 1 for 'A', because in the index $ and A are both represented as 00
        if (end_char_bucket == high_addr &&  end_char_bucketi < high_idx && val == 0) {
            high_bwtcount--;
        }

        high = high_count + high_bwtcount;

        //DEBUG START        
        //cout<< "highCCount "<<high_count<<" high_bwtcount "<<high_bwtcount<<" high_idx "<<high_idx<<"\n\n";
        //cout<<"symVal: "<<sym[i - 1]<<" lowNew: "<<(long) low<<" highNew: "<<(long)high<<"\n";
        //DEBUG END   

        if (low >= high)
            return;

    }

    // store hit
    *is_align = true;
    * read_low = low;
    * read_high = high;

}

/**
    Perform exact searching for both backward and reverse complement

    @param  &reads              vector that contains the short read read_t
    @param  *idx                the compressed FM-index in consecutive memory space
    @param  *cnt                i(x)
    @param  bitCcnt             32bit or 64bit for i(x)
    @param  *wpad               bool to make the length of BWT within a bucket a multiple of 2. 
                                This requires padding but forgos the use of divider.
    @param  bucket_pad_size     The padding size, if wpad == true
    @param  end_char_bucket     bucket index that contains $
    @param  end_char_bucketi    the offset withing a bucket for $
*/
template <class index_t, class cnt_t>
void exactAlign(std::vector<read_t> &reads, 
                index_t * idx, cnt_t * cnt, 
                uint32_t bitCcnt, uint32_t bucket_pad_size,
                uint32_t end_char_bucket, uint32_t end_char_bucketi){

    uint32_t bucket_bwt_len = (BUCKET_SIZE * 8 - FM_BP_RANGE * bitCcnt - bucket_pad_size * 8) / FM_BP_BIT ;


    if (bitCcnt == 32){

     int n_threads = omp_get_max_threads();    
#pragma omp parallel for num_threads(n_threads)
        for (uint32_t i = 0; i < reads.size(); i++) {

            // only peform exact match when there is no N char
            if (reads[i].has_N == false) {
                // exact align read
                //DEBUG START   
                //cout<<"\nBackward: "<<"\n";
                //DEBUG END
                exactSearch<index32_t, uint32_t>(idx, cnt, 
                                                reads[i].seq, reads[i].seq_len, 
                                                &reads[i].is_b_align, &reads[i].b_low, &reads[i].b_high, 
                                                bucket_bwt_len, end_char_bucket, end_char_bucketi,
                                                true);
                //DEBUG START 
                //cout<<"\nForward : Reverse Complement "<<"\n";
                //DEBUG END
                exactSearch<index32_t, uint32_t>(idx, cnt, 
                                                reads[i].seq, reads[i].seq_len, 
                                                &reads[i].is_f_align, &reads[i].f_low, &reads[i].f_high, 
                                                bucket_bwt_len, end_char_bucket, end_char_bucketi,
                                                false);

            }
        }

    }else{
        cout<<"64-bit is to be done later";
    }


}

#endif
