//
// Created by hcng on 6/29/18.
//

#include "format.h"



void faToFmt(FILE *fasta_fp, 
            FILE * N_fp, 
            char * fmt, 
            uint8_t bit_for_bp, 
            uint64_t & fmt_cnt,
            uint64_t & ref_cnt, 
            std::vector<nchar_cluster_t> &nchar_clusters,
            std::vector<chr_t> &chrs){

    // all the counting starts from zero
    char line[LINE_SIZE];
    uint32_t un_cnt = 0;
    uint64_t cum_un_cnt = 0;

    bool chr_valid = false;
    chr_t tmp;
    nchar_cluster_t nchar_cluster;

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
                        nchar_cluster.ref_cnt = ref_cnt - un_cnt;
                        nchar_cluster.fmt_cnt = fmt_cnt;
                        nchar_cluster.un_cnt = un_cnt;

                        nchar_cluster.cum_un_cnt = cum_un_cnt + un_cnt;

                        nchar_clusters.push_back(nchar_cluster);
                        cum_un_cnt += un_cnt;

                        writeNinfo(N_fp, ref_cnt - un_cnt, fmt_cnt, un_cnt, cum_un_cnt);

                        un_cnt = 0;
                    }

                    fmt_cnt++;
                }
                i++;
                ref_cnt++;
            }
        }
        // get the chromosome name
        else{
            if (chr_valid){
                chrs.push_back(tmp);
            }
            tmp.name.assign(line, strlen(line));
            tmp.begin = ref_cnt;
            chr_valid = true;
        }
        tmp.end = ref_cnt;
    }

    if (chr_valid){
        chrs.push_back(tmp);
    }
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


void writeChrName(FILE * fp, std::vector<chr_t> &chrs){

    for (uint16_t i = 0; i < chrs.size(); i++){
        fwrite(&(chrs[i].begin), sizeof(uint64_t), 1, fp);
        fwrite(&(chrs[i].end), sizeof(uint64_t), 1, fp);

        uint64_t str_length = (uint64_t) chrs[i].name.length();
        fwrite(&str_length, sizeof(uint64_t), 1, fp);
        fwrite(chrs[i].name.c_str(), sizeof(char), str_length, fp);
    }
    
    chrs.clear();

}