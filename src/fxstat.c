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

int main(int argc, char *argv[])
{
    time_t t0 = time(NULL);

    struct args args = argparser(argc, &argv);

    int err;
    FILE *stream = NULL;
    if(strcmp(args.infile, "-") == 0){
        stream = stdin;
    } else {
        if( access(args.infile, R_OK) == -1 ){
            fprintf(stderr, "File '%s' not found or not readable\n", args.infile);
            return 1;
        }
        stream = fopen(args.infile, "r");
    }
    gzFile fh = gzdopen(fileno(stream), "r");

    khash_t(int2long) *h_length = kh_init(int2long); // create a hashtable

    long nt_counter[256];
    for(int i = 0; i < 256; i++){
        nt_counter[i] = 0;
    }
    
    kseq_t *seq;
    int l;
    long n_seq = 0;
    long n_bases = 0;
    float sum_read_quality = 0;
    char spacer = '\0';

    seq = kseq_init(fh);
    while ((l = kseq_read(seq)) >= 0) {
        n_seq += 1;
        int seq_len = strlen(seq->seq.s);
        if(n_seq % 1000000 == 0){
            fprintf(stderr, "\r%ld reads processed", n_seq);
            fflush(stderr);
            spacer = '\n';
        }
        count_nt(seq->seq.s, seq_len, nt_counter);
        
        if(seq->qual.l){
            sum_read_quality += mean_quality(seq->qual.s);
        }

        n_bases += seq_len;
       
        h_length_update(h_length, seq_len);

        if(args.n_stop > 0 && n_seq >= args.n_stop){
            break;
        }
    } 
    kseq_destroy(seq);
    const char *err_msg = gzerror(fh, &err);
    if(strlen(err_msg) > 0){
        fprintf(stderr, "\rError: %s\n", err_msg);
    }
    fclose(stream);
    gzclose(fh);
    
    struct int_count *len_counts;
    len_counts = sort_int_table(h_length); 

    /*Collect stats*/

    int n_lens = kh_size(h_length);
   
    time_t t1= time(NULL);
    
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
    fprintf(stderr, "%c", spacer);
    fprintf(fout, "n_seq              %ld\n", n_seq);
    fprintf(fout, "n_bases            %ld\n", n_bases);
    fprintf(fout, "mean_length        %.2f\n", (float) n_bases / n_seq);
    fprintf(fout, "mean_read_quality  %.2f\n", sum_read_quality / n_seq);
    int nnx = 0;
    while(1){
        int nx = args.nx_ints[nnx];
        if(nx == -1){
            break;
        }
        fprintf(fout, "N%-5d             %d\n", nx, get_nx(&len_counts, n_lens, nx));
        nnx++;
    }
    fprintf(fout, "A                  %.2f%%\n", (float) 100 * (nt_counter['A'] + nt_counter['a']) / n_bases);
    fprintf(fout, "T                  %.2f%%\n", (float) 100 * (nt_counter['T'] + nt_counter['t']) / n_bases);
    fprintf(fout, "C                  %.2f%%\n", (float) 100 * (nt_counter['C'] + nt_counter['c']) / n_bases);
    fprintf(fout, "G                  %.2f%%\n", (float) 100 * (nt_counter['G'] + nt_counter['g']) / n_bases);
    fprintf(fout, "N                  %.2f%%\n", (float) 100 * (nt_counter['N'] + nt_counter['n']) / n_bases);
    fprintf(fout, "CG                 %.2f%%\n", (float) 100 * (nt_counter['G'] + nt_counter['g'] + nt_counter['C'] + nt_counter['c']) / n_bases);
    fclose(fout);

    fprintf(stderr, "# Proc time %s\n", format_seconds(t1-t0));

    kh_destroy_int2long(h_length);
    free(len_counts);
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

