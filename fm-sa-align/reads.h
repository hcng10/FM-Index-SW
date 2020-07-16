//
// Created by hcng on 9/20/18.
//

#ifndef FM_SA_ALIGN_READS_H
#define FM_SA_ALIGN_READS_H


#include <vector>
#include "file_op.h"
#include <sys/time.h>


struct sub_read_t{
    uint32_t start_pos;
    
    uint32_t low;
    uint32_t high;

    bool is_b_align;
    bool is_f_align;
    
};

class Read_t {
    private:
        uint16_t sub_read_cnt;
        sub_read_t * sub_reads;
        uint16_t * seed_rank;
        uint16_t seed_aligned_cnt;

    public:
        // trick the destructor
        bool callDestructor;

        char seq[MAX_READ_LEN+1];
        uint32_t seq_len;

        // symbol in binary
        uint8_t pck_sym[CEIL(MAX_READ_LEN, FM_BP_RANGE)];

        std::string at_data;
        std::string q_score;

        char plus_data;

        uint32_t b_low;//needed to be comment
        uint32_t b_high;

        uint32_t f_low;
        uint32_t f_high;//needed to be comment

        bool is_b_align;
        bool is_f_align;
        bool has_N;

        Read_t(){
            sub_read_cnt = 0;
            sub_reads = NULL;
            seed_rank = NULL;
            seed_aligned_cnt = 0;

            callDestructor = false;
        }

        
        void setSeqLen(uint32_t seq_len){
            uint16_t tmp = (seq_len / SEED_LEN) - 1;
            this->sub_read_cnt = (tmp * SEED_LEN + TRIM_READ_LEN <= seq_len ? tmp + 1 : tmp) * 2;

            this->seq_len = seq_len;
            sub_reads = new sub_read_t[sub_read_cnt];
            seed_rank = new uint16_t[sub_read_cnt];

            memset(sub_reads, 0, sizeof(sub_read_t) * sub_read_cnt);
            memset(seed_rank, 0, sizeof(uint16_t) * sub_read_cnt);

        }

        void setSeedRslt(uint16_t cur_entry_cnt, uint32_t start_pos, bool backward, uint32_t low, uint32_t high){
            uint16_t entry_num = backward ? cur_entry_cnt: (this->sub_read_cnt /2) + cur_entry_cnt;
            
            sub_reads[entry_num].start_pos = start_pos;
            if (backward == 1){
                sub_reads[entry_num].low = low;
                sub_reads[entry_num].high = high;
                sub_reads[entry_num].is_b_align = true;

                sub_reads[entry_num].start_pos = start_pos;

            }
            else{
                sub_reads[entry_num].low = low;
                sub_reads[entry_num].high = high;
                sub_reads[entry_num].is_f_align = true;

                sub_reads[entry_num].start_pos = start_pos;
            }

            if (seed_aligned_cnt == 0){
                seed_rank[0] = entry_num;
            }
            else{
                bool start_sft = false;
                uint16_t chk_seed_num = 0;
                for (uint16_t i = 0; i <seed_aligned_cnt; i++){

                    if (start_sft){
                        uint16_t tmp_seed_num = seed_rank[i];
                        seed_rank[i] = chk_seed_num;
                        chk_seed_num = tmp_seed_num;
                    }
                    else if (sub_reads[seed_rank[i]].high - sub_reads[seed_rank[i]].low >
                            high - low){
                        // insert here
                        chk_seed_num = seed_rank[i];
                        seed_rank[i] = entry_num;
                        start_sft = true;
                    
                    }
                }

                if (start_sft == false){
                    seed_rank[seed_aligned_cnt] = entry_num;
                }
            }

            seed_aligned_cnt++;

            for (int k = 0; k< sub_read_cnt; k++){
                cout<<seed_rank[k]<<" ";
            }
        }

        // we don't use destructor here because we didn't do deep copy
        ~Read_t(){
            if (callDestructor && sub_reads != NULL && seed_rank != NULL){
                delete[] sub_reads;
                delete[] seed_rank;
            }
        }


};




void loadReads(FILE *fp, std::vector<Read_t> &reads, char *buffer, uint64_t size_r, bool r_ctrl, uint64_t *bytes);
uint64_t writeReads(FILE *fp, std::vector<Read_t> &reads, char *buffer);


#endif //FM_SA_ALIGN_READS_H
