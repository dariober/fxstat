#ifndef UTILS
#define UTILS
#include "uthash.h"

struct int_count {
    int id;         /* key */
    long count;
    UT_hash_handle hh;  /* makes this structure hashable */
};

struct sequence_record {
    char *name;
    char *sequence;
    char *quality;
    int len;
};

struct actg_count {
    long A;
    long C;
    long G;
    long T;
    long N;
    int warns;
    int MAX_WARN;
};

void increment_int_count(int id, struct int_count **counter);

long get_int_count(int id, struct int_count **counter);

void delete_hash_table(struct int_count **counter);

int compare(const void *a, const void *b);

void count_nt(char *seq, struct actg_count *counter);

float mean_quality(char *sequence, int offset);

char* format_seconds(int seconds);

int is_fasta_comment(char *x);

// int is_blank(char *x);

int next_record(struct sequence_record *rec, FILE *fh, char *record_type, char **sequence_name);

int next_fasta_record(struct sequence_record *rec, FILE *fh, char **sequence_name);

#endif
