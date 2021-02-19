//
// Created by hcng on 9/20/18.
//

#include "reads.h"



// get remaining entry symbols if '@' is a quality score
uint64_t getOneMoreEntry(FILE *fp, char *buffer, uint64_t len)
{
    long int f_pos = ftell(fp);
    int cnt = 0;
    int c;

    // test if entry
    while ((c = fgetc(fp)) != EOF) {
        cnt += 1;
        if (c == '\n') {
            c = fgetc(fp);

            // It can be some random stuff?
            if (c != '@' && c != EOF)
                cnt = 0;
            break;
        }
    }

    // reset file pointer
    fseek(fp, f_pos, SEEK_SET);

    // read remaining positions
    for (int i = 0; i < cnt; i++)
        buffer[len++] = fgetc(fp);

    return len;
}



void loadReads(FILE *fp, std::vector<Read_t> &reads, char *buffer, uint64_t size_r, bool r_ctrl, uint64_t *bytes){

    struct timeval  tv1, tv2;
    gettimeofday(&tv1, NULL);

    reads.clear();

    if (r_ctrl == true) {

        // read a bunch of things from file first
        if (fread(buffer, size_r, 1, fp) != 1) {
            fprintf(stderr, "error: unable to read file!\n");
            exit(1);
        }

        // Then read until next record starts or eof
        int c;
        while ((c = fgetc(fp)) != EOF) {
            if (c == '@') {
                ungetc(c, fp);
                size_r = getOneMoreEntry(fp, buffer, size_r);
                break;
            }
            buffer[size_r++] = (char)c;
        }

        // update bytes read
        *bytes += size_r;


        // parse the buffer
        Read_t tmp;

        //tmp.is_b_align = false;
        //tmp.is_f_align = false;

        //tmp.b_low = 0;
        //tmp.b_high = 0;
        //tmp.f_low = 0;
        //tmp.f_high = 0;

        uint64_t i = 0;
        uint64_t s_position;

        int s_len;

        while (i < size_r) {
            tmp.has_N = false;

            //// Step 1: get meta data line, which is <@r337>
            s_position = i;
            s_len = 0;

            while (buffer[i++] != '\n') {
                s_len += 1;
            }

            tmp.at_data.assign(&buffer[s_position], s_len);


            ////Step 2: get the sequence
            //bool unknown = false;
            s_position = i;
            s_len = 0;

            while (buffer[i] != '\n') {
                if (buffer[i++] == 'N')
                    tmp.has_N = true;
                s_len += 1;
            }

            i += 1;

            //check if the length of the sequence is longer than the allocated size (not likely to happen)
            s_len = s_len > MAX_READ_LEN ? MAX_READ_LEN : s_len;

            memcpy(tmp.seq, &buffer[s_position], s_len * sizeof(char));
            tmp.seq[s_len] = '\0';
            tmp.setSeqLen(s_len);
            //tmp.seq_len = s_len;

            ////Step 3: get the line with '+' marker
            tmp.plus_data = buffer[i];
            while (buffer[i++] != '\n');

            ////Step 4: get the quality score
            s_position = i;
            s_len = 0;

            while (buffer[i++] != '\n') {
                s_len += 1;
            }
            tmp.q_score.assign(&buffer[s_position], s_len);

            //if (unknown == false)
            // will call the implicit copy constructor
            reads.push_back(tmp);
            
        }
    }
    gettimeofday(&tv2, NULL);
    printf("OK Load Reads Time [%.2f s]\n", (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
                                            (double) (tv2.tv_sec - tv1.tv_sec));
}



uint64_t writeReads(FILE *fp, std::vector<Read_t> &reads, char *buffer)
{
    uint64_t bytes = 0;
    uint64_t aligned_cnt = 0;

    // write hits to buffer
    for (uint32_t i = 0; i < reads.size(); i++) {
        if (reads[i].is_b_align == false || reads[i].is_f_align == false ) {

            // write buffer before overflow
            if ((bytes + 512) > BUFF_SIZE) {
                writeFile(fp, buffer, bytes);
                bytes = 0;
            }

            // write meta data
            memcpy(&buffer[bytes], reads[i].at_data.c_str(), reads[i].at_data.size());
            bytes += reads[i].at_data.size();
            buffer[bytes++] = '\n';

            // write seqeunce data
            memcpy(&buffer[bytes], reads[i].seq, reads[i].seq_len);
            bytes += reads[i].seq_len;
            buffer[bytes++] = '\n';

            // write strand
            buffer[bytes++] = '+';
            buffer[bytes++] = '\n';

            // write quality scores
            memcpy(&buffer[bytes], reads[i].q_score.c_str(), reads[i].q_score.size());
            bytes += reads[i].q_score.size();
            buffer[bytes++] = '\n';

        }else{
            aligned_cnt = aligned_cnt + 1;
        }
    }

    if (bytes > 0) {
        writeFile(fp, buffer, bytes);
    }

    return aligned_cnt;
}

