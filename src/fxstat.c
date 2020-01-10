#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "utils.h"

int main(int argc, char *argv[])
{
    time_t t0 = time(NULL);
    
    // Default arguments 
    char *version = "0.1.0";
    int nx_cutoff = 50;
    char *infile = "-";
    int seed = 0;

    int index;
    int opt;
    opterr = 0;

    while ((opt = getopt (argc, argv, "N:hV")) != -1)
        switch (opt)
        {
            case 'N':
                nx_cutoff = atoi(optarg);
                break;
            case 'V':
                printf("%s\n", version);
                return 0;
            case 'h':
                printf("Usage: fxstat [OPTION] ... [FILE]\n");
                printf("Collect sequence statistics from fastq file\n");
                printf("\n");
                printf("  -N  value for Nx statistics [%d]\n", nx_cutoff);
                //printf("  -s  seed for sampling. 0 for random seed [%d]\n", seed);
                //printf("  -b  number of bases for sample statistics to hold in memory [%ld]\n", max_bases_in_reservoir);
                printf("  -V  print version\n");
                printf("\n");
                printf("With no FILE, or when FILE is -, read standard input\n");
                printf("Base qualities must be in phred scale. There is no check for that!\n");
                printf("\n");
                printf("Version %s\n", version);
                return 0;
            case '?':
                if (optopt == 'N' || optopt == 's'){
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

    FILE *fh = NULL;
    if(strcmp(infile, "-") == 0){
        fh = stdin;
    } else {
        fh = fopen(infile, "r");
        if( access(infile, R_OK) == -1 ){
            fprintf(stderr, "File '%s' not found or not readable\n", infile);
            return 1;
        }
    }
    
    if(seed == 0){
        srand(time(0) + getpid());
    } else {
        srand(seed);
    }

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
    fclose(fh);
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
    
    printf("n_seq\t%ld\n", n_seq);
    printf("n_bases\t%ld\n", n_bases);
    printf("max_length\t%d\n", max_len);
    printf("mean_length\t%.2f\n", (float) n_bases / n_seq);
    printf("N%d\t%d\n", nx_cutoff, nx);
    printf("A\t%.2f%%\n", (float) 100 * nt_counter.A / n_bases);
    printf("T\t%.2f%%\n", (float) 100 * nt_counter.T / n_bases);
    printf("C\t%.2f%%\n", (float) 100 * nt_counter.C / n_bases);
    printf("G\t%.2f%%\n", (float) 100 * nt_counter.G / n_bases);
    printf("N\t%.2f%%\n", (float) 100 * nt_counter.N / n_bases);
    printf("CG\t%.2f%%\n", (float) 100 * (nt_counter.G + nt_counter.C) / n_bases);
    printf("mean_read_quality\t%.2f\n", sum_read_quality / n_seq);

    // printf("# Sampled %d reads (%ld bases)\n", reads_in_reservoir, bases_in_reservoir);
    fprintf(stderr, "# Proc time %s\n", format_seconds(t1-t0));

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
