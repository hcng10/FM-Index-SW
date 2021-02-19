#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <thread>
#include <bitset>
#include <omp.h>

#include "file_op.h"
#include "readFM.h"
#include "reads.h"
#include "format.h"
#include "align.hpp"
#include "writeSAM.h"


using namespace std;

int main(int argc, char *argv[]) {
    
    // stands for reconfigurable alignment
    string ext = "ral";

    // file descriptor for the input files
    FILE * N_fp = NULL;
    FILE * FM_fp = NULL;
    FILE * FM_meta_fp = NULL;

    FILE * rev_N_fp = NULL;
    FILE * rev_FM_fp = NULL;
    FILE * rev_FM_meta_fp = NULL;

    FILE * SA_fp = NULL;
    FILE * rev_SA_fp = NULL;

    FILE * chr_fp = NULL;

    FILE * in_fp = NULL;
    FILE * out_fp = NULL;
    FILE * sam_fp = NULL;

    // fmt_len is the length of the reference without the N nucleotide
    uint64_t fmt_len;

    // bool to indicate if 32bit is enough to hold the BWT of reference
    bool c32;

    uint32_t bucket_bwt_len;
    uint32_t n_buckets;
    bool wpad;
    uint32_t bucket_pad_size;

    // Marker for $ pos
    uint64_t end_char_pos[FM_STEP];
    uint8_t end_char[FM_STEP*FM_STEP];

    uint32_t end_char_bucket;
    uint32_t end_char_bucketi;

    // variables for bwt/FM-index
    // bitwidth for i(x)
    uint32_t bitCcnt;
    // i(x)
    uint32_t cnt32[FM_I_NUM + 1] = {0};
    // bucket
    index32_t * idx32;

    // suffix array
    uint32_t * sai;

    // number of chromosomes
    vector<chr_t> chrs;
    uint16_t chrs_num;
    
    // N char info
    vector<nchar_cluster_t> nchar_clusters;
    uint64_t N_cluster;

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
    string s10 = string(argv[1]) + ".10." + ext;

    string r1 = "result.fq";
    string rs1 = "aligned.sam";


    printf("Reading meta data about the index ... \n");fflush(stdout);

    openFile(&FM_fp, s1, "r");
    openFile(&FM_meta_fp, s2, "r");

    uint64_t f_size = fileSizeBytes(FM_fp);
    read_meta(FM_meta_fp,
                &c32, 
                &fmt_len, 
                &bucket_bwt_len,
                &wpad, 
                &bucket_pad_size, 
                end_char_pos, 
                end_char,
                &N_cluster,
                &chrs_num);

#if FM_STEP == 1
    end_char_bucket = (uint32_t) (end_char_pos[0] / bucket_bwt_len);
    end_char_bucketi = (uint32_t) (end_char_pos[0] % bucket_bwt_len);
#endif

    n_buckets = CEIL(fmt_len, bucket_bwt_len);

    if (c32){
        bitCcnt = 32;

        // read i(x)
        for (int i = 0; i < FM_BP_RANGE + 1; i++){
            uint64_t tmp_cnt;
            readFile(FM_fp, & tmp_cnt, sizeof(uint64_t));
            cnt32[i] = (uint32_t) tmp_cnt;

            //cout<<"cnt32 "<<i<<" "<< cnt32[i]<<"\n";
        }

        idx32 = new index32_t [n_buckets];
        // read the index one-by-one
        for (int i = 0; i < n_buckets; i++){
            readFile(FM_fp, & idx32[i], FM_BP_RANGE * 4 + BWT32_LEN_BYTE);
        }
    }
    fclose(FM_fp);
    fclose(FM_meta_fp);
    printf("FINISH ---> Reading meta data!\n\n");


   
    // read SA
    printf("Reading SA ... \n");fflush(stdout);

    openFile(&SA_fp, s7, "r");
    sai = new uint32_t[fmt_len];

    if (fread(sai, sizeof(uint32_t), fmt_len, SA_fp) != fmt_len) {
        fprintf(stderr, "error: unable to read SA file!\n");
        exit(1);
    }
    fclose(SA_fp);
    printf("FINISH ---> Reading SA\n\n");



    // read N char info
    // printf("Reading N char info in reference ... \n");fflush(stdout);

    // openFile(&N_fp, s3, "r");
    // readNinfo(N_fp, N_cluster, nchar_clusters);

    // fclose(N_fp);
    // printf("FINISH ---> Reading N char info\n\n");



    // read the name of chromosome
    printf("Reading the name of chromosome ... \n");fflush(stdout);

    openFile(&chr_fp, s9, "r");
    readChrName(chr_fp, chrs_num, chrs);

    fclose(chr_fp);
    printf("FINISH ---> Reading name of chromosome\n\n");


    printf("Reading reads ... \n");fflush(stdout);

    // allocate I/O buffers
    char * in_buff = new char [BUFF_SIZE + 512];
    char * out_buff = new char [BUFF_SIZE];
    char * sam_buff = new char [BUFF_SIZE];

    std::vector<read_t> reads1, reads2;


    openFile(&out_fp, r1, "w+");
    openFile(&sam_fp, rs1, "w+");
    // read first batch
    openFile(&in_fp, argv[2], "r");
    f_size = fileSizeBytes(in_fp);

    uint64_t bytes_r = 0;
    uint64_t size_r =  bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

    uint64_t n_Nchar = 0;
    loadReads(in_fp, reads1, in_buff, size_r, true, &bytes_r, &n_Nchar);

    // print the first 5 reads
    for (uint32_t i = 0; i < 5; i++) {
        printf("%u: %s\n", i, reads1[i].seq);
    }

    printf("FINISH ---> Reading reads\n\n");fflush(stdout);

    // some debugging variables
    uint64_t aligned_cnt1 = 0;
    uint64_t aligned_cnt2 = 0;
    uint32_t cnt = 0;

    printf("Aligning reads ... \n ");

    // ping-pong processing
    for (int i = 0; ; i++) {

        bool r_ctrl = bytes_r < f_size ? true : false;
        size_r = bytes_r + BUFF_SIZE <= f_size ? BUFF_SIZE : f_size - bytes_r;

        // read to r1, process r0
        if (!(i%2)) {
            align_thread = std::thread(loadReads, in_fp, std::ref(reads2), in_buff, size_r, r_ctrl, &bytes_r, &n_Nchar);

            if (reads1.size() > 0) {
                cnt += reads1.size();

                exactAlign <index32_t, uint32_t > (reads1, 
                                                    idx32, cnt32, 
                                                    bitCcnt, bucket_pad_size, 
                                                    end_char_bucket, end_char_bucketi);

                convertSAM(sam_fp, reads1, sai, sam_buff, chrs);

                aligned_cnt1 = aligned_cnt1 + writeReads(out_fp, reads1, out_buff);

            }
            else
                break;
        }else{
            align_thread = std::thread(loadReads, in_fp, std::ref(reads1), in_buff, size_r, r_ctrl, &bytes_r, &n_Nchar);

            if (reads2.size() > 0) {
                cnt += reads2.size();

                exactAlign <index32_t, uint32_t > (reads2, 
                                                    idx32, cnt32, 
                                                    bitCcnt, bucket_pad_size, 
                                                    end_char_bucket, end_char_bucketi);

                aligned_cnt2 = aligned_cnt2 + writeReads(out_fp, reads2, out_buff);

            }
            else
                break;

        }

        align_thread.join();

    }

    align_thread.join();
    

    printf("FINISH ---> Checking N\n\n");
    printf("processed %u reads\n", cnt);
    printf("%u reads have N char\n", n_Nchar);

    printf("aligned %lu reads\n", aligned_cnt1 + aligned_cnt2);

    delete [] idx32;
    delete [] in_buff;
    delete [] out_buff;
    delete [] sam_buff;
    delete [] sai;

    fclose(in_fp);
    fclose(out_fp);
    fclose(sam_fp);

}
