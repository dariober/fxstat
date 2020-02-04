#include <stdio.h>
#include <ctype.h> 
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include "khash.h"

#ifndef UTILS
#define UTILS

KHASH_MAP_INIT_INT(int2long, long)  // instantiate structs and methods
                                     // with int key and long value

/*struct to hold Nx statistics */
struct Nx {
    int N;
    int length;
};

/* struct to hold results that will be finally presented to the user
 * */
struct results {
    char *filename;
    long n_seq;
    long n_bases; 
    float median_length;
    float mean_read_quality;
    struct Nx *nx;
    long nt_counter[256];
};

struct results init_results(int nx_ints[]){
    struct results results;
    results.n_seq = 0;
    results.n_bases = 0;
    results.mean_read_quality = 0.0;
    for(int i = 0; i < 256; i++){
        results.nt_counter[i] = 0;
    }

    int nnx = 0;
    while(nx_ints[nnx] != -1){
        nnx++;
    }
    results.nx = malloc((nnx + 1) * sizeof(struct results));
    for(int i = 0; i < nnx; i++){
        results.nx[i].N = nx_ints[i];
    }
    results.nx[nnx].N = -1; // Signal end of N cutoffs
    
    return results;
}

void print_results(FILE *fout, struct results results){
    fprintf(fout, "[%s]\n", results.filename);
    fprintf(fout, "n_seq              %ld\n", results.n_seq);
    fprintf(fout, "n_bases            %ld\n", results.n_bases);
    fprintf(fout, "mean_length        %.2f\n", (float) results.n_bases / results.n_seq);
    fprintf(fout, "median_length      %.2f\n", results.median_length);
    fprintf(fout, "mean_read_quality  %.2f\n", results.mean_read_quality);

    int nnx = 0;
    while((results.nx[nnx].N) != -1){
        char cmt[100];
        if(results.nx[nnx].N == 0){
            strncpy(cmt, "  # MAX length", sizeof(cmt));
        } else if (results.nx[nnx].N == 100){
            strncpy(cmt, "  # MIN length > 0", sizeof(cmt));
        } else {
            strncpy(cmt, "", sizeof(cmt));
        }
        fprintf(fout, "N%-5d             %d%s\n", results.nx[nnx].N, results.nx[nnx].length, cmt);
        nnx++;
    }
    fprintf(fout, "A                  %.2f%%\n", 
            (float) 100 * (results.nt_counter['A'] + results.nt_counter['a']) / results.n_bases);
    fprintf(fout, "T                  %.2f%%\n", 
            (float) 100 * (results.nt_counter['T'] + results.nt_counter['t']) / results.n_bases);
    fprintf(fout, "C                  %.2f%%\n", 
            (float) 100 * (results.nt_counter['C'] + results.nt_counter['c']) / results.n_bases);
    fprintf(fout, "G                  %.2f%%\n", 
            (float) 100 * (results.nt_counter['G'] + results.nt_counter['g']) / results.n_bases);
    fprintf(fout, "N                  %.2f%%\n", 
            (float) 100 * (results.nt_counter['N'] + results.nt_counter['n']) / results.n_bases);
    fprintf(fout, "GC                 %.2f%%\n", 
            (float) 100 * (results.nt_counter['G'] + results.nt_counter['g'] + results.nt_counter['C'] + results.nt_counter['c']) / results.n_bases);
}

struct int_count {
    int id;
    long count;
};

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

/* Get the Nx value (e.g. the N50) of the sequence length table len_counts We
 * assume the number of entries in len_counts is not huge so passing through
 * them multiple times is not too expensive.
 * */
int get_nx(struct int_count **sorted_len_counts, int n_lens, int nx){
    
    struct int_count s;
    long tot_bases = 0;
    for(int i = 0; i < n_lens; i++) {
        s = (*sorted_len_counts)[i];
        tot_bases += s.id * s.count;
        if(i > 0 && (*sorted_len_counts)[i-1].id > (*sorted_len_counts)[i].id){
            // A sanity check
            fprintf(stderr, "List of counts is not sorted by length. Length %d before %d\n", 
                    (*sorted_len_counts)[i-1].id, (*sorted_len_counts)[i].id);
            exit(1);
        }
    }

    float nx_target_mass = tot_bases * ((float) nx / 100);
    int nx_value = -1;
    long nx_mass = 0;
    for(int i = n_lens-1; i >= 0; i--) {
        s = (*sorted_len_counts)[i];
        nx_mass += s.id * s.count; 
        if(nx_mass >= nx_target_mass && nx_value < 0){
            nx_value = s.id;
            break;
        }
    }
    return nx_value;
}

float median_length(struct int_count **sorted_len_counts, int n_lens){
    // Get total number of reads
    struct int_count s;
    long tot_reads = 0;
    for(int i = 0; i < n_lens; i++) {
        s = (*sorted_len_counts)[i];
        tot_reads += s.count;
        if(i > 0 && (*sorted_len_counts)[i-1].id > (*sorted_len_counts)[i].id){
            // A sanity check
            fprintf(stderr, "List of counts is not sorted by length. Length %d before %d\n", 
                    (*sorted_len_counts)[i-1].id, (*sorted_len_counts)[i].id);
            exit(1);
        }
    }

    long idxlow, idxhigh;
    if(tot_reads % 2 == 0){
        idxlow = tot_reads / 2;
        idxhigh = idxlow + 1;
    } else {
        idxlow = (tot_reads+1) / 2;
        idxhigh = idxlow;
    }
    int low = -1;
    int high = -1;
    long sum = 0;
    for(int i = 0; i < n_lens; i++) {
        s = (*sorted_len_counts)[i];
        sum += s.count;
        if(sum >= idxlow && low == -1){
            low = s.id;
        }
        if(sum >= idxhigh && high == -1){
            high = s.id;
            break;
        }
    }
    return (low + high) / 2.0; 
}

int compare(const void *p, const void *q)  
{ 
    int l = ((struct int_count *)p)->id; 
    int r = ((struct int_count *)q)->id;  
    return (l - r); 
} 

static inline void count_nt(char *seq, int len, long *counter)
{
    for(int i = 0; i < len; i++){
        counter[(int)seq[i]] += 1;
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
static inline float mean_quality(char *quality){
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

static inline void h_length_update(void *_hash, int key){
    
    khash_t(int2long) *hash = (khash_t(int2long)*)_hash;
    
    int ret;
    khiter_t k = kh_get(int2long, hash, key);
    if(k == kh_end(hash)){
        k = kh_put(int2long, hash, key, &ret); // add the key
        kh_value(hash, k) = 1; // set the value of the key
    } else {
        // If key found, increment value
        kh_value(hash, k) += 1;
    }
}

struct int_count * sort_int_table(void *_hash){
    
    khash_t(int2long) *hash = (khash_t(int2long)*)_hash;

    int n = kh_size(hash);
    struct int_count *len_counts = malloc(n * sizeof(struct int_count));
    
    int i = 0;
    struct int_count s;
    khiter_t k;
    for (k = kh_begin(hash); k != kh_end(hash); ++k) {
        if (kh_exist(hash, k)) {
            s.id = kh_key(hash, k);
            s.count = kh_value(hash, k);
            len_counts[i] = s;
            i++;
        }
    }
    qsort(len_counts, n, sizeof(struct int_count), compare);
    return len_counts;
}

#endif
