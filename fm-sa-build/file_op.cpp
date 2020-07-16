//
// Created by hcng on 6/28/18.
//

#include "file_op.h"


void openFile(FILE **fp, string f_name, const char *mode) {

    *fp = fopen(f_name.c_str(), mode);
    if (!(*fp)) {
        fprintf(stderr, "error: unable to open file '%s'!\n", f_name.c_str());
        exit(1);
    }
}

// get file size in bytes
uint64_t fileSizeBytes(FILE *fp)
{
    struct stat st;
    uint64_t len;
    int fd;

    if ((fd = fileno(fp)) == -1) {
        printf("error: unable to get file size!\n");
        exit(1);
    }

    if(fstat(fd, &st)) {
        printf("error: unable to get file size!\n");
        exit(1);
    }

    len = st.st_size;

    return len;
}

void writeFile(FILE *fp, void *data, uint64_t n_bytes)
{
    if (fwrite(data, n_bytes, 1, fp) != 1) {
        fprintf(stderr, "error: unable to write file!");
        exit(1);
    }
}

uint64_t writeNinfo(FILE *fp, u_int64_t ref_cnt, u_int64_t fmt_cnt, u_int32_t un_cnt){

    //printf("%ld, %ld, %d  %ld  %ld  %ld\n", ref_cnt, fmt_cnt, un_cnt, sizeof(u_int64_t), sizeof(u_int64_t), sizeof(u_int32_t));
    
    //Recording info: 
    // 1. the starting point of N chars in the original reference
    // 2. the point where N char got cut off and replaced with actual nucleotie
    // 3. Numner of N chars
    fwrite(&ref_cnt, sizeof(char), sizeof(u_int64_t), fp);
    fwrite(&fmt_cnt, sizeof(char), sizeof(u_int64_t), fp);
    fwrite(&un_cnt, sizeof(char), sizeof(u_int32_t), fp);

}