//
// Created by hcng on 9/19/18.
//

#ifndef FM_SA_ALIGN_FILE_OP_H
#define FM_SA_ALIGN_FILE_OP_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>

#include <iostream>

#include "../fm-sa-def.h"


using namespace std;

// open file
void openFile(FILE **fp, string f_name, const char *mode);
void readFile(FILE *fp, void *a, uint64_t n_bytes);
uint64_t fileSizeBytes(FILE *fp);
void writeFile(FILE *fp, void *a, uint64_t n_bytes);

#endif //FM_SA_ALIGN_FILE_OP_H
