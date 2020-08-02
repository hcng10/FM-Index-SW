//
// Created by hcng on 6/29/18.
//

#include "format.h"

#include "iostream"

void faToFmt(FILE *fasta_fp, FILE * N_fp, char * fmt, 
                    uint8_t bit_for_bp, uint64_t & fmt_cnt,
                    uint64_t & ref_cnt, uint64_t & N_cluster){

    // all the counting starts from zero
    char line[LINE_SIZE];
    uint32_t un_cnt = 0;

    // read FASTA lines
    while (fgets(line, LINE_SIZE, fasta_fp) != NULL) {
        // not header line
        if (line[0] != '>') {

            int i = 0;

            // parse line
            while (line[i] > 32) {

                char sym = toupper(line[i]);

                if ((sym != 'A' && sym != 'C' && sym != 'G' && sym != 'T') && bit_for_bp == 2) {
                    un_cnt++;

                } else {
                    // find an actual nucleotide, write the previously 
                    // found N nucleotide into file
                    fmt[fmt_cnt] = sym;
                    if (un_cnt > 0) {
                        writeNinfo(N_fp, ref_cnt - un_cnt, fmt_cnt, un_cnt);
                        N_cluster++;
                        un_cnt = 0;
                    }

                    fmt_cnt++;
                }
                i++;
                ref_cnt++;
            }


        }
    }
    //return fmt_cnt;

}

inline void setVal(uint8_t *pck, uint64_t idx, uint8_t idx_bp_range, uint8_t val){
  uint8_t tmp = val << ((idx_bp_range * FM_BP_BIT) % (sizeof(uint8_t)*8));
  pck[idx] |= tmp;
}


void packSymbols(char *sym, uint8_t *pck, uint8_t len){

    int n_threads = omp_get_max_threads();
    #pragma omp parallel for num_threads(n_threads)

    for (uint64_t i = 0; i < CEIL(len, FM_BP_RANGE); i++){
        for (uint8_t j = 0; j < FM_BP_RANGE; j++){
            switch(sym[i * FM_BP_RANGE + j]) {
                case 'A': setVal(pck, i, j, 0); break;
                case 'C': setVal(pck, i, j, 1); break;
                case 'G': setVal(pck, i, j, 2); break;
                case 'T': setVal(pck, i, j, 3); break;
                default : setVal(pck, i, j, 0);           
            }
        }
    }
}