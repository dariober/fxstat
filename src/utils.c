#include <stdio.h>
#include <ctype.h> 
#include <stdlib.h>
#include <string.h>
#include "utils.h"

void increment_int_count(int id, struct int_count **counter) {
    struct int_count *s;
    HASH_FIND_INT(*counter, &id, s);  /* id already in the hash? */
    if(s == NULL) {
        s = (struct int_count *)malloc(sizeof(*s));
        s->id = id;
        HASH_ADD_INT(*counter, id, s);  /* id: name of key field */
        s->count = 0;
    }
    s->count = (s->count)+1;
}

long get_int_count(int id, struct int_count **counter) {
    struct int_count *s;
    HASH_FIND_INT(*counter, &id, s);  /* s: output pointer */
    if(s == NULL){
        return 0; 
    } else {
        return s->count;
    }
}

void delete_hash_table(struct int_count **counter) {
    struct int_count *current, *tmp;

    HASH_ITER(hh, *counter, current, tmp) {
        HASH_DEL(*counter, current);  /* delete it (users advances to next) */
        free(current);             /* free it */
    }
}

int compare(const void *a, const void *b)
{
     int int_a = * ( (int*) a );
     int int_b = * ( (int*) b );
     return(int_b - int_a);
}

void count_nt(char *seq, struct actg_count *counter)
{
    int seq_len = strlen(seq);
    for(int i = 0; i < seq_len; i++){
        char nuc = toupper(seq[i]);
        if(nuc == 'A'){
            counter->A += 1;
        } 
        else if(nuc == 'C') {
            counter->C += 1;
        }
        else if(nuc == 'G') {
            counter->G += 1;
        }
        else if(nuc == 'T') {
            counter->T += 1;
        }
        else if(nuc == 'N') {
            counter->N += 1;
        } else {
            counter->warns += 1;
            if(counter->warns <= counter->MAX_WARN){
                fprintf(stderr, "Warning: Unexpected nucleotide: %c\n", seq[i]);
            } 
            if(counter->warns == counter->MAX_WARN){
                fprintf(stderr, "Additional warning will not be shown\n");
            }
            counter->N += 1;
        }
    }
}

float mean_quality(char *sequence, int offset){
    int sum = 0;
    int len = strlen(sequence);
    for(int i = 0; i < len; i++){
        sum += (sequence[i] - offset);
    }
    return (float)sum/len;
}

char* format_seconds(int seconds){
	int h, m, s;
	h = (seconds/3600); 
	m = (seconds -(3600*h))/60;
	s = (seconds -(3600*h)-(m*60));
	static char fmt[100];
    snprintf(fmt, sizeof(fmt), "%02d:%02d:%02d", h, m, s);
    return fmt;
}

int is_fasta_comment(char *x){
    while(isspace(*x)){
        x++; 
    }
    if(x[0] == ';'){
        return 1; 
    } else {
        return 0;
    }
}

int is_blank(char *x){
    int len = strlen(x);
    for(int i = 0; i < len; i++){
        if(! isspace(x[i])){
            return 0; 
        }
    }
    return 1;
}

int next_record(struct sequence_record *rec, FILE *fh, char *record_type, char **sequence_name){
    int eof = 0;
    int add_fasta_name = 1;
    char *line = NULL;
    size_t len = 0;
    int mseq = 10000;
    rec->sequence = (char*)malloc(mseq);
    rec->quality = (char*)malloc(mseq);

    rec->len = 0;

    while(1){
        if(getline(&line, &len, fh) == -1){
            if(rec->len == 0){
                free(rec->sequence);
                free(rec->quality);
                eof = -1;
            }
            break;
        }

        line[strcspn(line, "\n\r")] = '\0';
        
        if(is_blank(line)){
            continue;
        }

        if(*record_type == 'u'){
            // Inspect first non-blank line to tell if record is fasta or fastq
            if(line[0] == '@'){
                *record_type = 'q'; 
            } else {
                *record_type = 'a';
            }
        }
        if(*record_type == 'q'){
            // FASTQ READER
            if(line[0] != '@'){
                fprintf(stderr, "Invalid fastq sequence name: '%s'\n", line);
                exit(1); 
            } 
            rec->name = strdup(line);

            getline(&line, &len, fh);    // Sequence
            line[strcspn(line, "\n\r")] = '\0';
            rec->len = strlen(line);
            if((rec->len + 1) > mseq){
                rec->sequence = (char*) realloc(rec->sequence, rec->len + 1);
            }
            strcpy(rec->sequence, line);

            getline(&line, &len, fh);    // Comment (ignored)

            getline(&line, &len, fh);    // Quality
            line[strcspn(line, "\n\r")] = '\0';
            if((rec->len + 1) > mseq){
                rec->quality = (char*) realloc(rec->quality, rec->len + 1);
            }
            strcpy(rec->quality, line);

            break;
        } else {
            // FASTA READER
            if(add_fasta_name && *sequence_name[0] != 0 && *sequence_name[0] != '1'){
                rec->name = strdup(*sequence_name); 
                add_fasta_name = 0;
            }

            if(is_fasta_comment(line)){
                continue;
            }

            if(line[0] == '>'){
                if(strlen(*sequence_name) == 0){
                    *sequence_name[0] = '1'; // Signal that we have read the first sequence header in the file
                    rec->name = strdup(line);
                } else {
                    *sequence_name = (char*) realloc(*sequence_name, strlen(line)+1);
                    if(!*sequence_name){
                        fprintf(stderr, "Failed to allocate memory for sequence header\n");
                        exit(1);
                    }
                    strcpy(*sequence_name, line);
                    break;
                }
            } else {
                rec->len += strlen(line);
                if((rec->len + 1) > mseq){
                    mseq = rec->len+1;
                    rec->sequence = (char*) realloc(rec->sequence, mseq);
                    if(!rec->sequence){
                        fprintf(stderr, "Failed to allocate memory for sequence of length %d\n", mseq-1);
                        exit(1);
                    }
                }
                sprintf(&rec->sequence[rec->len - strlen(line)], "%s", line);
            }
        } 
    }
    free(line); 
    return eof;
}

