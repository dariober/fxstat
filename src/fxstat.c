#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <zlib.h>
#include "utils.h"

int main(int argc, char *argv[])
{
    time_t t0 = time(NULL);
    
    // Default arguments 
    char *version = "0.1.0";
    int nx_cutoff = 50;
    char *infile = "-";
    char *outfile = "-";
    // int seed = 0;

    int index;
    int opt;
    opterr = 0;

    while ((opt = getopt (argc, argv, "N:hVo:")) != -1)
        switch (opt)
        {
            case 'N':
                nx_cutoff = atoi(optarg);
                break;
            case 'o':
                outfile = optarg;
                break;
            case 'V':
                printf("%s\n", version);
                return 0;
            case 'h':
                printf("Usage: fxstat [OPTION] ... [FILE]\n");
                printf("Collect sequence statistics from fastq/fasta file\n");
                printf("\n");
                printf("  -N  value for Nx statistics [%d]\n", nx_cutoff);
                printf("  -o  output file [%s]\n", outfile);
                //printf("  -s  seed for sampling. 0 for random seed [%d]\n", seed);
                //printf("  -b  number of bases for sample statistics to hold in memory [%ld]\n", max_bases_in_reservoir);
                printf("  -V  print version\n");
                printf("\n");
                printf("\
With no FILE, or when FILE is -, read standard input.\n\
Input may be gzip'd compressed. Base qualities must\n\
be in phred scale. There is no check for that!\n");
                printf("\n");
                printf("Version %s\n", version);
                return 0;
            case '?':
                if (optopt == 'N' || optopt == 's' || optopt == 'o'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else {
                    fprintf (stderr, "Unknown option `-%c'\n", optopt);
                }
            default:
                return 1;
        }

    for(index = optind; index < argc; index++){
        infile = argv[index];
    }

    if(nx_cutoff < 0 || nx_cutoff > 100){
        fprintf(stderr, "Argument to -N must be >= 0 and <= 100. Got %d\n", nx_cutoff);
        return 1;
    }
    
    int err;
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

    struct int_count *len_histogram = NULL;

    struct actg_count nt_counter = {.A = 0, .C = 0, .G = 0, .T = 0, .N = 0, .warns = 0, .MAX_WARN = 20};

    struct sequence_record *rec; 
    rec = malloc(sizeof(struct sequence_record));

    long n_seq = 0;
    long n_bases = 0;
    float sum_read_quality = 0;
    char *fasta_sequence_name;
    fasta_sequence_name = calloc(10, 1);
    char record_type = 'u'; 

    while(next_record(rec, fh, &record_type, &fasta_sequence_name) == 0) {
        n_seq += 1;
        if(n_seq % 1000000 == 0){
            fprintf(stderr, "%ld reads processed\n", n_seq);
        }
        count_nt(rec->sequence, &nt_counter);
        
        if(record_type == 'q'){
            sum_read_quality += mean_quality(rec->quality);
        }

        n_bases += rec->len;
        increment_int_count(rec->len, &len_histogram);
        
        free(rec->name);
        free(rec->sequence);
        free(rec->quality);
    } 
    const char *err_msg = gzerror(fh, &err);
    if(strlen(err_msg) > 0){
        fprintf(stderr, "Error: %s\n", err_msg);
    }
    fclose(stream);
    gzclose(fh);
    free(rec); 
    free(fasta_sequence_name);
    
    HASH_SORT(len_histogram, compare);

    /*Collect stats*/

    long nx_mass = 0;
    int nx = -1;
    float nx_target = n_bases * ((float) nx_cutoff / 100);
    int max_len = len_histogram->id;
   
    struct int_count *s;
    for(s = len_histogram; s != NULL; s = s->hh.next) {
        nx_mass += (s->id) * (s->count); 
        if(nx_mass >= nx_target && nx < 0){
            nx = s->id;
        }
    }

    delete_hash_table(&len_histogram);

    time_t t1= time(NULL);
    
    FILE *fout;
    if(strcmp(outfile, "-") == 0){
        fout = stdout;
    } else {
        fout = fopen(outfile, "w");
        if(access(outfile, W_OK) == -1){
            fprintf(stderr, "Unable to write to '%s'\n", outfile);
            return 1;
        }
    }

    fprintf(fout, "n_seq\t%ld\n", n_seq);
    fprintf(fout, "n_bases\t%ld\n", n_bases);
    fprintf(fout, "max_length\t%d\n", max_len);
    fprintf(fout, "mean_length\t%.2f\n", (float) n_bases / n_seq);
    fprintf(fout, "N%d\t%d\n", nx_cutoff, nx);
    fprintf(fout, "A\t%.2f%%\n", (float) 100 * nt_counter.A / n_bases);
    fprintf(fout, "T\t%.2f%%\n", (float) 100 * nt_counter.T / n_bases);
    fprintf(fout, "C\t%.2f%%\n", (float) 100 * nt_counter.C / n_bases);
    fprintf(fout, "G\t%.2f%%\n", (float) 100 * nt_counter.G / n_bases);
    fprintf(fout, "N\t%.2f%%\n", (float) 100 * nt_counter.N / n_bases);
    fprintf(fout, "CG\t%.2f%%\n", (float) 100 * (nt_counter.G + nt_counter.C) / n_bases);
    fprintf(fout, "mean_read_quality\t%.2f\n", sum_read_quality / n_seq);
    fclose(fout);

    fprintf(stderr, "# Proc time %s\n", format_seconds(t1-t0));
    return err;
}
        /* Reservoir sampling 
         * We start by filling up the reservoir by counting the number of
         * bases in it, not the number of sequences. We assume that in
         * terms of length distribution the first reads are representative
         * of the whole population.
         *
         * Once the reservoir is full, we use the number of sequences in it
         * to drive the sampling.
         *
         * */
        /*
        if(bases_in_reservoir < max_bases_in_reservoir){
            char **temp = realloc(reservoir, (reads_in_reservoir + 1) * sizeof(char*)); 
            if(temp == NULL){
                fprintf(stderr, "Failed to allocate memory\n");
                exit(1);
            }
            reservoir = temp;
            reservoir[reads_in_reservoir] = malloc((rec.len + 1) * sizeof(char));
            strcpy(reservoir[reads_in_reservoir], rec.sequence);
            
            reads_in_reservoir += 1;
            bases_in_reservoir += rec.len;
        } else {
            int j = rand() % n_seq; // Random int in range [0, n_seq]
            if(j <= reads_in_reservoir){
                reservoir[j] = strdup(rec.sequence);
            }
        }
        */
