#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "utils.h"
#include "kseq.h"
#include "khash.h"
#include "argparser.c"

KSEQ_INIT(gzFile, gzread)

int scan_file(char *infile, long n_stop, struct results *results){

    int err = 0;
    if(n_stop > 0 && results->n_seq >= n_stop){
        return err;
    }
    
    results->filename = realloc(results->filename, strlen(results->filename) + strlen(infile) + 3);
    if( ! results->filename){
        fprintf(stderr, "Cannot allocate memory for filename string\n");
        exit(1);
    }
    strcat(results->filename, infile);
    strcat(results->filename, ", ");

    FILE *stream = NULL;
    if(strcmp(infile, "-") == 0){
        stream = stdin;
    } else {
        if( access(infile, R_OK) == -1 ){
            fprintf(stderr, "File '%s' not found or not readable\n", infile);
            return 1;
        }
        stream = fopen(infile, "r");
    }
    gzFile fh = gzdopen(fileno(stream), "r");
    
    kseq_t *seq;
    int l;

    seq = kseq_init(fh);
    while ((l = kseq_read(seq)) >= 0) {
        
        results->n_seq += 1;
        int seq_len = strlen(seq->seq.s);
        if(results->n_seq % 1000000 == 0){
            fprintf(stderr, "\r%ld reads processed", results->n_seq);
            fflush(stderr);
        }
        count_nt(seq->seq.s, seq_len, results->nt_counter);
        
        if(seq->qual.l){
            results->sum_read_quality += mean_quality(seq->qual.s);
        }
       
        h_length_update(results->h_length, seq_len);

        if(n_stop > 0 && results->n_seq >= n_stop){
            break;
        }
    } 
    
    kseq_destroy(seq);
   
    const char *err_msg = gzerror(fh, &err);
    if(strlen(err_msg) > 0){
        fprintf(stderr, "\rError: %s\n", err_msg);
        return err;
    }
    fclose(stream);
    gzclose(fh);

    return err;
}

int main(int argc, char *argv[])
{
    time_t t0 = time(NULL);

    struct args args = argparser(argc, &argv);

    FILE *fout;
    if(strcmp(args.outfile, "-") == 0){
        fout = stdout;
    } else {
        fout = fopen(args.outfile, "w");
        if(access(args.outfile, W_OK) == -1){
            fprintf(stderr, "Unable to write to '%s'\n", args.outfile);
            return 1;
        }
    }

    int nfiles = 0;
    int err;
    struct results results;
    while(args.infile[nfiles] != NULL){
        if(nfiles == 0 || args.pull_files == 0){
            results = init_results(args.nx_ints);
        }
        err = scan_file(args.infile[nfiles], args.n_stop, &results);
        if(err != 0){
            return err;
        }
        if(results.n_seq > 1000000){
            fprintf(stderr, "\n");
        }
        if(args.pull_files == 0){
            flush_results(fout, results);
        }
        nfiles++;
    }
    
    if(args.pull_files == 1){
        flush_results(fout, results);
    }
    fclose(fout);

    for(int i = 0; i <= nfiles; i++){
        free(args.infile[i]);
    }
    free(args.infile);

    time_t t1= time(NULL);
    fprintf(stderr, "# Proc time %s\n", format_seconds(t1-t0));
    return err;
}
