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
    
    results->filename = strdup(infile); // malloc(strlen(infile) + 1);

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
    
    khash_t(int2long) *h_length = kh_init(int2long); // create a hashtable

    kseq_t *seq;
    int l;
    float sum_read_quality = 0;

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
            sum_read_quality += mean_quality(seq->qual.s);
        }

        results->n_bases += seq_len;
       
        h_length_update(h_length, seq_len);

        if(n_stop > 0 && results->n_seq >= n_stop){
            break;
        }
    } 
    
    kseq_destroy(seq);
   
    int err;
    const char *err_msg = gzerror(fh, &err);
    if(strlen(err_msg) > 0){
        fprintf(stderr, "\rError: %s\n", err_msg);
        return err;
    }
    fclose(stream);
    gzclose(fh);

    struct int_count *len_counts;
    len_counts = sort_int_table(h_length); 


    int n_lens = kh_size(h_length);
    
    results->median_length = median_length(&len_counts, n_lens);

    int nnx = 0;
    while(results->nx[nnx].N != -1){
        results->nx[nnx].length = get_nx(&len_counts, n_lens, results->nx[nnx].N);
        nnx++;
    }

    kh_destroy_int2long(h_length);
    results->mean_read_quality = (float) sum_read_quality / results->n_seq;
    free(len_counts);
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
    while(args.infile[nfiles] != NULL){
        struct results results = init_results(args.nx_ints);
        err = scan_file(args.infile[nfiles], args.n_stop, &results);
        if(err != 0){
            return err;
        }
        if(results.n_seq > 1000000){
            fprintf(stderr, "\n");
        }
        print_results(fout, results);
        free(results.nx);        
        free(results.filename);
        nfiles++;
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
    /*
    char xseq[] = "ABCAACACACACABABABAACACBABABABABBABABA";
    
    int kmer_size = 1;
    int ret; 
    khiter_t k;
    khash_t(khStrLong) *h = kh_init(khStrLong); // create a hashtable

    for(int i = 0; i < 10000000; i++){
        int len = strlen(xseq);
        for(int j = 0; j < len - kmer_size + 1; j++){
            char kmer[kmer_size+1];  
            memcpy(kmer, xseq+j, kmer_size);
            kmer[kmer_size] = '\0';
            k = kh_get(khStrLong, h, kmer);
            if(k == kh_end(h)){
                k = kh_put(khStrLong, h, strdup(kmer), &ret);
                kh_value(h, k) = 0;
            }
            kh_value(h, k) += 1;
        }
    }
    for (k = kh_begin(h); k != kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            const char *key = kh_key(h,k);
            long tval = kh_value(h, k);
            printf("key=%s  val=%ld\n", key, tval);
            // free((char*)kh_key(h, k));
        }
    }
    kh_destroy(khStrLong, h);
    return 0;
    */    
    /*
    k = kh_put(khStrLong, h, 1, &ret); // add the key

    for(int i = 0; i < 10; i++){
        int key = 1;
        k = kh_get(int2long, h, key); 
        //if (k == kh_end(h)){
        //    // If a key is absent node == end of table
        //    k = kh_put(int2long, h, key, &ret); // add the key
        //    kh_value(h, k) = 1; // set the value of the key
        //} else {
            // If key found, increment value
            kh_value(h, k) += 1;
        //}
    }
    printf("Buckets %d\n", kh_size(h));

    for (k = kh_begin(h); k != kh_end(h); ++k) {
        if (kh_exist(h, k)) {
            const int key = kh_key(h,k);
            long tval = kh_value(h, k);
            printf("key=%d  val=%ld\n", key, tval);
        }
    }

    kh_destroy(int2long, h);              // deallocate the hash table
    return 0;
    */

