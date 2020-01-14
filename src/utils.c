#include <stdio.h>
#include <ctype.h> 
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include <zlib.h>

// Conversion from Phred score to probability
static const float PHRED_PROB[] = {
    1.0,
    0.7943282347242815,
    0.6309573444801932,
    0.5011872336272722,
    0.3981071705534972,
    0.31622776601683794,
    0.251188643150958,
    0.19952623149688797,
    0.15848931924611134,
    0.12589254117941673,
    0.1,
    0.07943282347242814,
    0.06309573444801933,
    0.05011872336272722,
    0.039810717055349734,
    0.03162277660168379,
    0.025118864315095794,
    0.0199526231496888,
    0.015848931924611134,
    0.012589254117941675,
    0.01,
    0.007943282347242814,
    0.00630957344480193,
    0.005011872336272725,
    0.003981071705534973,
    0.0031622776601683794,
    0.0025118864315095794,
    0.001995262314968879,
    0.001584893192461114,
    0.0012589254117941675,
    0.001,
    0.0007943282347242813,
    0.000630957344480193,
    0.0005011872336272725,
    0.00039810717055349735,
    0.00031622776601683794,
    0.00025118864315095795,
    0.00019952623149688788,
    0.00015848931924611142,
    0.00012589254117941674,
    0.0001,
    7.943282347242822e-05,
    6.309573444801929e-05,
    5.011872336272725e-05,
    3.9810717055349695e-05,
    3.1622776601683795e-05,
    2.5118864315095822e-05,
    1.9952623149688786e-05,
    1.584893192461114e-05,
    1.2589254117941661e-05,
    1e-05,
    7.943282347242822e-06,
    6.30957344480193e-06,
    5.011872336272725e-06,
    3.981071705534969e-06,
    3.162277660168379e-06,
    2.5118864315095823e-06,
    1.9952623149688787e-06,
    1.584893192461114e-06,
    1.2589254117941661e-06,
    1e-06,
    7.943282347242822e-07,
    6.30957344480193e-07,
    5.011872336272725e-07,
    3.981071705534969e-07,
    3.162277660168379e-07,
    2.5118864315095823e-07,
    1.9952623149688787e-07,
    1.584893192461114e-07,
    1.2589254117941662e-07,
    1e-07,
    7.943282347242822e-08,
    6.30957344480193e-08,
    5.011872336272725e-08,
    3.981071705534969e-08,
    3.162277660168379e-08,
    2.511886431509582e-08,
    1.9952623149688786e-08,
    1.5848931924611143e-08,
    1.2589254117941661e-08,
    1e-08,
    7.943282347242822e-09,
    6.309573444801943e-09,
    5.011872336272715e-09,
    3.981071705534969e-09,
    3.1622776601683795e-09,
    2.511886431509582e-09,
    1.9952623149688828e-09,
    1.584893192461111e-09,
    1.2589254117941663e-09,
    1e-09,
    7.943282347242822e-10,
    6.309573444801942e-10,
    5.011872336272714e-10,
    3.9810717055349694e-10,
    3.1622776601683795e-10,
    2.511886431509582e-10,
    1.9952623149688828e-10,
    1.584893192461111e-10,
    1.2589254117941662e-10,
    1e-10,
    7.943282347242822e-11,
    6.309573444801942e-11,
    5.011872336272715e-11,
    3.9810717055349695e-11,
    3.1622776601683794e-11,
    2.5118864315095823e-11,
    1.9952623149688828e-11,
    1.5848931924611107e-11,
    1.2589254117941662e-11,
    1e-11,
    7.943282347242821e-12,
    6.309573444801943e-12,
    5.011872336272715e-12,
    3.9810717055349695e-12,
    3.1622776601683794e-12,
    2.5118864315095823e-12,
    1.9952623149688827e-12,
    1.584893192461111e-12,
    1.258925411794166e-12,
    1e-12,
    7.943282347242822e-13,
    6.309573444801942e-13,
    5.011872336272715e-13,
    3.981071705534969e-13,
    3.162277660168379e-13,
    2.511886431509582e-13,
    1.9952623149688827e-13};

char * zgetline(gzFile fh, char **line){
    int buffer_size = 1000;
    int line_length = buffer_size;

    *line = (char*) realloc(*line, buffer_size);
    if(*line == NULL){
        fprintf(stderr, "Failed to allocate memory for line of length %d\n", buffer_size);
        exit(1);
    }
    int pos = 0;
    char buffer[buffer_size];
    char *eof;
    while( (eof = gzgets(fh, buffer, buffer_size)) ){
        int len = strlen(buffer);
        if((pos + len + 1) >= line_length){
            line_length *= 2;
            *line = (char*) realloc(*line, line_length);
            if(*line == NULL){
                fprintf(stderr, "Failed to allocate memory for line of length %d\n", line_length);
                exit(1);
            }
        }
        memcpy((*line)+pos, buffer, len);
        pos += len;
        if(buffer[len-1] == '\n'){
            break;
        } 
    }
    (*line)[pos] = '\0';
    return eof;
}

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
    int i = 0;
    while(seq[i] != '\0'){
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
        i++;
    }
}

/* Faster implementation of log10 function from
 * http://www.machinedlearnings.com/2011/06/fast-approximate-logarithm-exponential.html
 * */
float fastlog10 (float x) {
    union { float f; uint32_t i; } vx = { x };
    union { uint32_t i; float f; } mx = { (vx.i & 0x007FFFFF) | (0x7e << 23) };
    float y = vx.i;
    y *= 1.0 / (1 << 23);
 
    float log_2 =  y - 124.22544637f - 1.498030302f * mx.f - 1.72587999f / (0.3520887068f + mx.f);
    // Convert log2 to log10. The constant 0.301 comes from log10(n) / log2(x)
    return 0.301029995f * log_2;
}

/* Qualities are in phred scale (log scale) so we cannot directly average
 * them.  We need to convert to probabilities, average, convert back to phred.
 * See also https://gigabaseorgigabyte.wordpress.com/2017/06/26/averaging-basecall-quality-scores-the-right-way/
 * */
float mean_quality(char *quality){
    float sum = 0;
    int len = strlen(quality);
    for(int i = 0; i < len; i++){
        sum += PHRED_PROB[quality[i] - 33];
    }
    float avg = sum / len;
    float phred = -10 * fastlog10(avg);
    return phred;
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

int next_record(struct sequence_record *rec, gzFile fh, char *record_type, char **sequence_name){
    int eof = 0;
    int add_fasta_name = 1;
    char *line = NULL;
    int mseq = 10000;
    rec->sequence = (char*)malloc(mseq);
    rec->quality = (char*)malloc(mseq);

    rec->len = 0;

    while(1){
        if(zgetline(fh, &line) == 0){
            if(rec->len == 0){
                free(rec->sequence);
                free(rec->quality);
                eof = -1;
            }
            break;
        }

        line[strcspn(line, "\n\r")] = '\0'; /* Strip trailing newline */
        
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
            
            char *p;
            
            p = zgetline(fh, &line);    // Sequence
            line[strcspn(line, "\n\r")] = '\0';
            rec->len = strlen(line);
            if((rec->len + 1) > mseq){
                rec->sequence = (char*) realloc(rec->sequence, rec->len + 1);
            }
            strcpy(rec->sequence, line);

            p= zgetline(fh, &line);    // Comment (ignored)

            p= zgetline(fh, &line);    // Quality
            if(p == 0){
                fprintf(stderr, "\nWarning: Fastq file may be truncated!\n\n");
            }
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

