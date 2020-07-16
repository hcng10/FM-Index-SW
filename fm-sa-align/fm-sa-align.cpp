#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <thread>
#include <bitset>


#include "file_op.h"
#include "readFM.h"
#include "reads.h"
#include "align.hpp"


using namespace std;

int main(int argc, char *argv[]) {

    string ext = "ral";

    // variables for file
    FILE *N_fp = NULL;
    FILE *FM_fp = NULL;
    FILE * FM_meta_fp = NULL;

    FILE *rev_N_fp = NULL;
    FILE *rev_FM_fp = NULL;
    FILE * rev_FM_meta_fp = NULL;

    FILE * FM_SA = NULL;
    FILE * rev_FM_SA = NULL;

    FILE *in_fp = NULL;
    FILE *out_fp = NULL;

    // variables from meta-file
    uint64_t fmt_len;
    bool c32;
    uint32_t bucket_bwt_len;
    uint32_t bucket_pad_size;
    uint64_t end_char_pos[FM_STEP];
    uint8_t end_char[FM_STEP*FM_STEP];


    uint32_t end_char_bucket;
    uint32_t end_char_bucketi;

    // variables for bwt/ FM-index
    bool wpad;
    uint32_t n_buckets;
    index32_t * idx32;


    uint32_t cnt32[FM_I_NUM + 1] = {0};
    uint32_t bitCcnt;

    uint32_t * sai;

    std::thread align_thread;

    // program usage
    if (argc != 3) {
        printf("usage: %s <index basename> <Reads file>\n", argv[0]);
        exit(1);
    }


    string s1 = string(argv[1]) + ".1." + ext;
    string s2 = string(argv[1]) + ".2." + ext;
    string s3 = string(argv[1]) + ".3." + ext;

    string s4 = string(argv[1]) + ".4." + ext;
    string s5 = string(argv[1]) + ".5." + ext;
    string s6 = string(argv[1]) + ".6." + ext;

    string s7 = string(argv[1]) + ".7." + ext;
    string s8 = string(argv[1]) + ".8." + ext;

    string s9 = string(argv[1]) + ".9." + ext;

    string r1 = "result.fq";


    printf("loading index ... ");
    printf("Reading meta data ... ");fflush(stdout);

    openFile(&FM_fp, s1, "r");
    openFile(&FM_meta_fp, s2, "r");

    uint64_t f_size = fileSizeBytes(FM_fp);
    read_meta(FM_meta_fp,
                &c32, &fmt_len, 
                &bucket_bwt_len,
                &wpad, &bucket_pad_size, 
                end_char_pos, end_char);
                

#if FM_STEP == 1
    end_char_bucket = (uint32_t) (end_char_pos[0] / bucket_bwt_len);
    end_char_bucketi = (uint32_t) (end_char_pos[0] % bucket_bwt_len);
#endif

    n_buckets = CEIL(fmt_len, bucket_bwt_len);

    if (c32){
        //read the i-table (note it is always stored at 64bit)

        for (int i = 0; i < FM_BP_RANGE + 1; i++){
            uint64_t tmp_cnt;
            readFile(FM_fp, & tmp_cnt, sizeof(uint64_t));
            cnt32[i] = (uint32_t) tmp_cnt;

            //DEBUG START:
            //cout<<"cnt32 "<<i<<" "<< cnt32[i]<<"\n";
            //DEBUG END:
        }

        idx32 = new index32_t [n_buckets];
        // read the index one-by-one
        for (int i = 0; i < n_buckets; i++){
            readFile(FM_fp, & idx32[i], FM_BP_RANGE * 4 + BWT32_LEN_BYTE);
        }

        //readFile(FM_fp, idx32, f_size - sizeof(uint64_t) * (FM_BP_RANGE + 1));
        //cerr<<"f_size - sizeof(uint64_t) * (FM_BP_RANGE + 1) "<< sizeof(index32_t) <<" "<< fseek(FM_fp, 0, SEEK_CUR)<<"\n";

        //cout<<"n_buckets "<<n_buckets<<" fmt_len "<<fmt_len<<" bucket_bwt_len "<<bucket_bwt_len<<"\n\n";

        //cout<<"f_size: "<<f_size<<"Data structure size: "<<n_buckets*sizeof(index32_t);
        bitCcnt = 32;

        for (int c=0;c<4;c++){
            cout<<"The c index"<<(int)c<<" index: "<<(int) idx32[0].count[c]<<"\n";
        }

        //for (int b=0;b<(BUCKET_SIZE - 32 * FM_BP_RANGE / 8);b++)
            //cout<<"The final "<<(int)b<<" val "<<std::bitset<8>(idx32[0].bwt[b])<<"\n";

    }

    fclose(FM_fp);
    fclose(FM_meta_fp);

    printf("OK!\n");

   
    // read SA
    printf("Reading SA ... \n");fflush(stdout);

    openFile(&FM_SA, s7, "r");
    
    sai = new uint32_t[fmt_len];
    cerr<<fmt_len;
    if (fread(sai, sizeof(uint32_t), fmt_len, FM_SA) != fmt_len) {
        fprintf(stderr, "error: unable to read SA file!\n");
        exit(1);
    }

    /*for (int i= 0; i<fmt_len; i++){
        cout<<"SA "<<sai[i]<<"\n";
    }*/

    fclose(FM_SA);
    printf("OK!\n");


    printf("Reading reads ... \n");fflush(stdout);

    // allocate I/O buffers
    char * in_buff = new char [BUFF_SIZE + 512];
    char * out_buff = new char [BUFF_SIZE];

    std::vector<Read_t> reads1, reads2;


    openFile(&out_fp, r1, "w+");
    // read first batch
    openFile(&in_fp, argv[2], "r");
    f_size = fileSizeBytes(in_fp);

    uint64_t bytes_r = 0;
    uint64_t size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

    loadReads(in_fp, reads1, in_buff, size_r, true, &bytes_r);



    for (uint32_t i = 0; i < 5; i++) {
        printf("%u: %s\n", i, reads1[i].seq);
    }

    uint64_t aligned_cnt1 = 0;
    uint64_t aligned_cnt2 = 0;
    uint32_t cnt = 0;
    for (int i = 0; ; i++) {

        bool r_ctrl = bytes_r < f_size ? true : false;
        size_r = bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

        // read to r1, process r0
        if (!(i%2)) {
            align_thread = std::thread(loadReads, in_fp, std::ref(reads2), in_buff, size_r, r_ctrl, &bytes_r);

            if (reads1.size() > 0) {
                cnt += reads1.size();

                exactAlign <index32_t, uint32_t > (reads1, idx32, cnt32, bitCcnt, bucket_pad_size, end_char_bucket, end_char_bucketi);
                aligned_cnt1 = aligned_cnt1 + writeReads(out_fp, reads1, out_buff);
                

            }
            else
                break;
        }else{
            align_thread = std::thread(loadReads, in_fp, std::ref(reads1), in_buff, size_r, r_ctrl, &bytes_r);

            if (reads2.size() > 0) {
                cnt += reads2.size();

                exactAlign <index32_t, uint32_t > (reads2, idx32, cnt32, bitCcnt, bucket_pad_size, end_char_bucket, end_char_bucketi);
                aligned_cnt2 = aligned_cnt2 + writeReads(out_fp, reads2, out_buff);

            }
            else
                break;

        }

        align_thread.join();
        printf("processed %u reads\n", cnt);
        printf("aligned %lu reads\n", aligned_cnt1 + aligned_cnt2);

    }

    align_thread.join();

    delete [] idx32;
    delete [] in_buff;
    delete [] out_buff;
    delete [] sai;

    fclose(in_fp);
    fclose(out_fp);

}
