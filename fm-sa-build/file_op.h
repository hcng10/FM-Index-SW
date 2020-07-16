//
// Created by hcng on 6/27/18.
//

#ifndef FM_SA_FILE_OP_H
#define FM_SA_FILE_OP_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <bitset>

using namespace std;

// open file
void openFile(FILE **fp, string f_name, const char *mode);

// get file size in bytes
uint64_t fileSizeBytes(FILE *fp);
void writeFile(FILE *fp, void *data, uint64_t n_bytes);
uint64_t writeNinfo(FILE *fp, u_int64_t ref_cnt, u_int64_t fmt_cnt, u_int32_t un_cnt);

#endif //FM_SA_H
