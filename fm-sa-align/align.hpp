#ifndef FM_SA_ALIGN_HPP
#define FM_SA_ALIGN_HPP

#include <omp.h>
#include <vector>

#include "reads.h"
#include "readFM.h"


template <class index_t, class cnt_t>
void exactSearch(index_t * index,  cnt_t * cnt, 
                Read_t & read,
                uint32_t bucket_bwt_len, uint32_t endCharBucket, uint32_t endCharBucketi,
                bool backward) {

    char *sym = read.seq;
    uint16_t read_len = read.seq_len;


    uint8_t val;
    uint16_t cur_entry_cnt;
    for (uint16_t i = 0; i < read_len - TRIM_READ_LEN ; i = i + SEED_LEN){

        uint32_t low = 0;
        uint32_t high = cnt[4];

        uint16_t start_pos = backward == 1? i + TRIM_READ_LEN - 1: read_len - i - TRIM_READ_LEN;
        uint16_t sym_idx = start_pos;

        for (uint8_t j = 0; j < TRIM_READ_LEN; j++){
        
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

            bool same_bucket = high - low < bucket_bwt_len ? true : false;

            uint32_t low_addr = low / bucket_bwt_len;
            uint32_t low_idx = low % bucket_bwt_len;

            uint32_t low_count = index[low_addr].count[val];
            uint32_t low_bwtcount = getOcc(val, index[low_addr].bwt, 0, low_idx);

            cout<< "lowCCount "<<low_count<<" low_bwtcount "<<low_bwtcount<<" low_idx "<<low_idx<<"\n";


            //handle endchar, as it will add 1 for 'A', because in the index $ and A are both represented as 00
            if (endCharBucket == low_addr &&  endCharBucketi < low_idx && val == 0) {
                low_bwtcount--;
            }


            low = low_count + low_bwtcount;

            uint32_t high_addr = high / bucket_bwt_len;
            uint32_t high_idx = high % bucket_bwt_len;

            uint32_t high_count = index[high_addr].count[val];
            uint32_t high_bwtcount = getOcc(val, index[high_addr].bwt, 0, high_idx);

            //handle endchar, as it will add 1 for 'A', because in the index $ and A are both represented as 00
            if (endCharBucket == high_addr &&  endCharBucketi < high_idx && val == 0) {
                high_bwtcount--;
            }

            high = high_count + high_bwtcount;
                    cout<< "highCCount "<<high_count<<" high_bwtcount "<<high_bwtcount<<" high_idx "<<high_idx<<"\n\n";
                                                                            
            cout<<(int)sym_idx<<"symVal: "<<sym[sym_idx]<<" lowNew: "<<(long) low<<" highNew: "<<(long)high<<"\n";

            sym_idx = backward == 1? sym_idx - 1: sym_idx + 1;

            if (low >= high){
                break;
            }
        }
        if (low < high){
            read.setSeedRslt(cur_entry_cnt, start_pos, backward, low, high);
        }

        cur_entry_cnt++;
    }

    /*exit(1);



    // store hit
    *is_align = true;
    * read_low = low;
    * read_high = high;*/

}



template <class index_t, class cnt_t>
void exactAlign(std::vector<Read_t> &reads, index_t * idx, cnt_t * cnt, uint32_t bitCcnt, uint32_t bucket_pad_size,
                uint32_t endCharBucket, uint32_t endCharBucketi){

    uint32_t bucket_bwt_len = (BUCKET_SIZE * 8 - FM_BP_RANGE * bitCcnt - bucket_pad_size * 8) / FM_BP_BIT ;

    //uint32_t bucket_bwt_len = (BUCKET_SIZE * 8 - FM_BP_RANGE * bitCcnt) / FM_BP_BIT;

    int n_threads = omp_get_max_threads();

    if (bitCcnt == 32){

        //#pragma omp parallel for num_threads(n_threads)
        for (uint32_t i = 0; i < reads.size(); i++) {

            // only peform exact match when there is no N char
            if (reads[i].has_N == false) {
                // exact align read
                cout<<"\nBackward: "<<"\n";
                exactSearch<index32_t, uint32_t>(idx, cnt, 
                                                reads[i],
                                                bucket_bwt_len, endCharBucket, endCharBucketi,
                                                true);

                cout<<"\nForward : Reverse Complement "<<"\n";
                exactSearch<index32_t, uint32_t>(idx, cnt, 
                                                reads[i],
                                                bucket_bwt_len, endCharBucket, endCharBucketi,
                                                false);
                reads[i].callDestructor = true;

            }
        }

    }else{
        cout<<"64-bit is to be done later, after LSA";
    }


}





#endif
