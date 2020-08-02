#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdint.h>
#include <iostream>
#include <bitset>



#include "../fm-sa-def.h"
#include "file_op.h"
#include "format.h"
#include "buildFM.h"



using namespace std;

void reverse_N_fp(FILE * N_fp, FILE * rev_N_fp, uint64_t fmt_len, uint64_t ref_len, uint64_t N_cluster);
void reverse(char *fmt, uint32_t len);

int main(int argc, char *argv[]) {

    // stands for reconfigurable alignment
    string ext = "ral";

    // the input refernce fasta file
    FILE * fasta_fp = NULL;

    // temp storage for the format, in char 
    char * fmt = NULL;
    // convert the nucleotide to 2 bits, pack them
    uint8_t * pck_fmt = NULL;

    // file descriptor for the output files
    FILE * N_fp = NULL;
    FILE * FM_fp = NULL;
    FILE * FM_meta_fp = NULL;

    FILE * rev_N_fp = NULL;
    FILE * rev_FM_fp = NULL;
    FILE * rev_FM_meta_fp = NULL;

    FILE * SA_fp = NULL;
    FILE * rev_SA_fp = NULL;

    FILE * ref_bin_fp = NULL;

    // fmt_len is the length of the reference without the N nucleotide
    uint64_t fmt_len = 0;
    uint64_t ref_len = 0;
    uint64_t N_cluster = 0;


    // program usage
    if (argc != 4) {
        printf("usage: %s <FASTA file> <Output basename> <with padding: true or false>\n", argv[0]);
        exit(1);
    }

    //**
    // The files format
    // Forward reference
    // 1: FM-index in bucket format
    // 2: Meta file, e.g. endCharPos, fmt_len...
    // 3: Info for N nucleotide: pos, length...
    // Reverse reference
    // Same for 4-6
    //
    // 7: Suffix Array
    // 8: Reverse Suffix Array
    // 9: The reference in binary
    //**

    string s_in = string(argv[1]);
    string s1 = string(argv[2]) + ".1." + ext;
    string s2 = string(argv[2]) + ".2." + ext;
    string s3 = string(argv[2]) + ".3." + ext;

    string s4 = string(argv[2]) + ".4." + ext;
    string s5 = string(argv[2]) + ".5." + ext;
    string s6 = string(argv[2]) + ".6." + ext;

    string s7 = string(argv[2]) + ".7." + ext;
    string s8 = string(argv[2]) + ".8." + ext;
    string s9 = string(argv[2]) + ".8." + ext;

    // extra parameter that enables padding in a bucket
    bool wpad = strcmp(argv[3], "true") == 0 ? true: false;

    // 1. construct reference sequence in binary format
    printf("Getting reference sequence ... \n"); fflush(stdout);

    openFile(&fasta_fp, s_in, "r");
    uint64_t fa_len = fileSizeBytes(fasta_fp);

    // new char array to hold all the nucleotide, 1 byte 1 nucleotide
    fmt = new char [fa_len];
    if (!fmt){
        fprintf(stderr, "error: unable to allocate memory for reference sequence!\n");
        exit(1);
    }

    openFile(&N_fp, s3, "w+");
    faToFmt(fasta_fp, N_fp , fmt, FM_BP_BIT, fmt_len, ref_len, N_cluster);

    fmt[fmt_len++] = '$';
    //fmt[fmt_len] = '\0';
    fclose(fasta_fp);

    printf("\tOriginal Reference + '$' Length: %ld\n", ref_len+1);
    printf("\tConstructed Reference + '$' Length: %ld\n", fmt_len);
    printf("\tNumber of N nucleotide clusters: %ld\n", N_cluster);
        
    printf("FINISH ---> Getting reference\n\n");fflush(stdout);
    
    // 2. Store up the sequence in binary format
    //      Only need to store the forward part as it is used in Smith-Waterman only
    printf("Converting reference into binary ... \n"); fflush(stdout);
    
    pck_fmt = new uint8_t[CEIL(fmt_len-1, FM_BP_RANGE)];
    memset(pck_fmt, 0, CEIL(fmt_len-1, FM_BP_RANGE) * sizeof(uint8_t));
    packSymbols(fmt, pck_fmt, fmt_len-1);

    openFile(&ref_bin_fp, s9, "w+");
    writeFile(ref_bin_fp, pck_fmt, sizeof(uint8_t) * CEIL(fmt_len-1, FM_BP_RANGE));
    fclose(ref_bin_fp);

    delete [] pck_fmt;


    printf("FINISH ---> Converting reference into binary\n\n");fflush(stdout);

    // 3. construct the FM-index
    printf("constructing FM-index of reference sequence ... \n"); fflush(stdout);

    openFile(&FM_fp, s1, "w+");
    openFile(&FM_meta_fp, s2, "w+");
    openFile(&SA_fp, s7, "w+");


    fmtToIdx(FM_fp, FM_meta_fp, SA_fp, fmt, fmt_len, BUCKET_SIZE, wpad);

    fclose(FM_fp);
    fclose(FM_meta_fp);

    printf("FINISH: --> Constructing FM-index\n\n");fflush(stdout);


    // 4. construct the FM-index for the reverse of the reference
    printf("reversing reference sequence ... "); fflush(stdout);

    reverse(fmt, fmt_len-1);
    printf("OK!\n");


    printf("constructing FM-index of reversed reference sequence ... \n"); fflush(stdout);
    openFile(&rev_FM_fp, s4, "w+");
    openFile(&rev_FM_meta_fp, s5, "w+");
    openFile(&rev_SA_fp, s8, "w+");

    fmtToIdx(rev_FM_fp, rev_FM_meta_fp, rev_SA_fp, fmt, fmt_len, BUCKET_SIZE, wpad);

    fclose(rev_FM_fp);
    fclose(rev_FM_meta_fp);


    openFile(&rev_N_fp, s6, "w+");
    reverse_N_fp(N_fp, rev_N_fp, fmt_len, ref_len, N_cluster);
    fclose(N_fp);
    fclose(rev_N_fp);

    printf("FINISH: --> Constructing the reversed FM-index\n\n");fflush(stdout);

    delete[] fmt;

}

void reverse_N_fp(FILE * N_fp, FILE * rev_N_fp, uint64_t fmt_len, uint64_t ref_len, uint64_t N_cluster){


    u_int64_t ref_cnt;
    u_int64_t fmt_cnt;
    u_int32_t un_cnt;


    int pos = sizeof(u_int64_t) + sizeof(u_int64_t) + sizeof(un_cnt);
    for (int i = N_cluster - 1; i >= 0; i--){
        fseek(N_fp , pos * i , SEEK_SET);

        fread(&ref_cnt, sizeof(u_int64_t), 1,  N_fp);
        fread(&fmt_cnt, sizeof(u_int64_t), 1, N_fp);
        fread(&un_cnt, sizeof(u_int32_t), 1, N_fp);

        ref_cnt = ref_len - ref_cnt - 1;
        fmt_cnt = fmt_len - fmt_cnt -1;

        writeNinfo(rev_N_fp, ref_cnt, fmt_cnt, un_cnt);

    }

}



// reverse reference sequence
void reverse(char *fmt, uint32_t len)
{
    char *p1 = &fmt[0];
    char *p2 = &fmt[len-1];
    char tmp;

    while (p2 > p1) {
        tmp = *p1;
        *p1++ = *p2;
        *p2-- = tmp;
    }
}
