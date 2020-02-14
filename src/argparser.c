#include <string.h>

# define nxstr "0,25,50,75,100"

int N_FILES = 0;

struct args {
    char *version;
    int nx_ints[101]; 
    char **infile;
    char *outfile;
    long n_stop;
    int pull_files;
};

// Default arguments
// We need the args object to be global becouse we pass args.infile to nftw
// function
struct args args ={
    .version = "0.1.0",
    .outfile = "-",
    .n_stop = -1,
    .pull_files = 0
};

struct file {
    char *userpath;
    char *fullpath;
};

void free_file(struct file f){
    free(f.userpath);
    free(f.fullpath);
}

void copy_file(struct file *src, struct file *dst){
    dst->fullpath = strdup(src->fullpath);
    dst->userpath = strdup(src->userpath);
}

/* Check fpath is in list file_list. n: Number files in file_list */
int file_found(char ***file_list, const char *fpath, int n){
    char x[PATH_MAX];
    char *qry = realpath(fpath, x);
    if(qry == NULL){
        // This could mean we are reading from stdin
        return 0;
    }
    for(int i = 0; i < n; i++){
        char y[PATH_MAX];
        char *sbj = realpath((*file_list)[i], y);
        if(strcmp(qry, sbj) == 0) {
            return 1;
        }
    }
    return 0;
}

int fullpath_cmp(const void *v1, const void *v2)
{
    const struct file *p1 = (struct file *)v1;
    const struct file *p2 = (struct file *)v2;
    return strcmp(p1->fullpath, p2->fullpath);
}

int userpath_cmp(const void *v1, const void *v2)
{
    const struct file *p1 = (struct file *)v1;
    const struct file *p2 = (struct file *)v2;
    return strcmp(p1->userpath, p2->userpath);
}

/* Remove duplicate files from array of file sorted by fullpath. Duplicates are
 * filepaths pointing the same file. n: number of elements in file_list */
void remove_duplicate_files(struct file **files, int *n){
    if(n == 0){
        return;
    }
    // Check if the next element in full_path list is different from the
    // previous one. If so, add it to the dedup list
    struct file *dedups = malloc((*n) * sizeof(struct file));
    copy_file(&(*files)[0], &(dedups[0]));
    int ndedup = 1;
    for(int i = 1; i < (*n); i++){
        int cmp = strcmp(dedups[ndedup-1].fullpath, (*files)[i].fullpath);
        if(cmp == 0){
            // Skip duplicate
        } else if(cmp < 0) {
            copy_file(&(*files)[i], &(dedups[ndedup]));
            ndedup++;
        } else {
            fprintf(stderr, "Files not sorted\n");
            exit(1);
        }
    }
    // Copy the dedup list into the original file_list and update the number 
    // of files in it.
    for(int i = 0; i < *n; i++){
        free_file((*files)[i]);
    }
    (*files) = realloc((*files), ndedup * sizeof(struct file));
    for(int i = 0; i < ndedup; i++){
        copy_file(&(dedups[i]), &(*files)[i]);
        free_file(dedups[i]);
    }
    (*n) = ndedup;
    free(dedups);
}

/* Add filename fpath to the array of pointers file_list. n is a counter of
 * files in file_list */
void add_file(char ***file_list, const char *fpath, int *n){
    if((*n) == 0){
        (*file_list) = malloc(sizeof(char *));
    } else {
        (*file_list) = realloc((*file_list), ((*n)+1) * sizeof(char *));
    }
    if(strcmp(fpath, "-") != 0 && access(fpath, R_OK) == -1){
        fprintf(stderr, "Invalid input file: '%s'\n", fpath);
        exit(1);
    }
    (*file_list)[(*n)] = strdup(fpath);
    (*n)++;
}

int isDirectory(const char *path) {
   struct stat statbuf;
   if (stat(path, &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}

/* From https://stackoverflow.com/questions/744766/how-to-compare-ends-of-strings-in-c */
int EndsWith(const char *str, const char *suffix)
{
    if (!str || !suffix)
        return 0;
    size_t lenstr = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix >  lenstr)
        return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}

/* Test if filename fpath is a fastq (1), fasta (2) or neither (0)
 * */
int is_sequence_filename(const char *fpath){
    
    char *f = strdup(fpath);
    for(int i = 0; i < strlen(f); i++){
        f[i] = tolower(f[i]);
    }
    if(EndsWith(f, ".gz")){
        f[strlen(f) - 3] = 0;
    }

    int type;
    if(EndsWith(f, ".fastq") || EndsWith(f, ".fq")){
        type = 1;
    }
    else if(EndsWith(f, ".fasta") || EndsWith(f, ".fa")) {
        type = 2;
    }
    else {
        type = 0;
    }
    free(f);
    return type;
}

int find_fastx_files(const char *fpath, const struct stat *sb, int tflag, struct FTW *ftwbuf) {
    if( (tflag == FTW_F || tflag == FTW_SL) && is_sequence_filename(fpath) ){
        add_file(&(args.infile), fpath, &N_FILES);
    }
    return 0;  /* To tell nftw() to continue */
}

void collect_files(char *fpath){
    if(strcmp(fpath, "-") == 0) {
        add_file(&(args.infile), fpath, &N_FILES);
    } else if(isDirectory(fpath) == 0) {
        add_file(&(args.infile), fpath, &N_FILES);
    } else {
        int x = nftw(fpath, find_fastx_files, 15, FTW_PHYS);
        if(x != 0){
            fprintf(stderr, "Cannot read file or directory '%s'\n", fpath);
            exit(x);
        }
    }
}

void dedup_files(char ***file_list, int *n_files)
{
    struct file *files = malloc(*n_files * sizeof(struct file));
    for(int i = 0; i < *n_files; i++){
        files[i].userpath = strdup((*file_list)[i]);
        char x[PATH_MAX];
        char *f = realpath((*file_list)[i], x);
        if(f == NULL){
            // This cold be we are reading from stdin
            files[i].fullpath = strdup((*file_list)[i]);
        } else {
            files[i].fullpath = strdup(f);
        }
        free((*file_list)[i]);
    }
    free((*file_list));

    qsort(files, *n_files, sizeof(struct file), fullpath_cmp);
    remove_duplicate_files(&files, n_files);
    qsort(files, *n_files, sizeof(struct file), userpath_cmp);
    
    // Copy the collected files to infile paramater:
    args.infile = malloc((*n_files + 1) * sizeof(char *));
    for(int i = 0; i < *n_files; i++){
        args.infile[i] = strdup(files[i].userpath);
    }
    args.infile[*n_files] = 0; // Signal end of file list

    for(int i = 0; i < *n_files; i++){
        free_file(files[i]);
    }
    free(files);
}

void tokenize_nx(char *str, int ints[]){
    char *token = strtok(str, ",");
    int nnx = 0;
    while (token != NULL) {
        int tok = atoi(token);
        int is_dupl = 0;
        for(int i = 0; i < nnx; i++){
            // Check for duplicates and skip them.
            if(tok == ints[i]){
                is_dupl = 1;
                break;
            }
        }
        if(! is_dupl){
            if(tok < 0 || tok > 100){
                fprintf(stderr, "Argument to -N must be >= 0 and <= 100\n");
                exit(1);
            }
            if(nnx > 101){
                fprintf(stderr, "Too many arguments to -s\n");
                exit(1);
            }
            ints[nnx] = tok;
            nnx++;
        }
        token = strtok(NULL, ",");
    }
    ints[nnx] = -1; // signal end of user's arguments
}

struct args argparser(int argc, char **argv[]){

    char nx_ints[1000] = nxstr;
    tokenize_nx(nx_ints, args.nx_ints); 

    int opt;
    opterr = 0;

    while ((opt = getopt (argc, *argv, "s:N:hVo:p")) != -1)
        switch (opt)
        {
            case 's':
                args.n_stop = atol(optarg);
                break;
            case 'N':
                strcpy(nx_ints, optarg);
                tokenize_nx(nx_ints, args.nx_ints); 
                break;
            case 'o':
                args.outfile = optarg;
                break;
            case 'p':
                args.pull_files = 1;
                break;
            case 'V':
                printf("%s\n", args.version);
                exit(0);
            case 'h':
                printf("Usage: fxstat [OPTION]... [FILE|DIR]...\n");
                printf("Collect sequence statistics from fastq/fasta files.\n");
                printf("\n");
                printf("  -p  Pull all FILEs in a single summary\n");
                printf("  -s  Stop after this many sequences [%ld]\n", args.n_stop);
                printf("  -N  Comma-separated thresholds for Nx statistics [%s]\n", nxstr);
                printf("  -o  output file [%s]\n", args.outfile);
                printf("  -V  print version\n");
                printf("\n");
                printf("\
With no FILE, or when FILE is -, read standard input.\n\
Input may be gzip'd compressed. Base qualities must\n\
be in phred scale. There is no check for that!\n\
\n\
If input is a directory, recursively find all sequence\n\
files based on filename extension.\n\
");
                printf("\n");
                printf("Version %s\n", args.version);
                exit(0);
            case '?':
                if (optopt == 'N' || optopt == 's' || optopt == 'o'){
                    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
                } else {
                    fprintf (stderr, "Unknown option `-%c'\n", optopt);
                }
            default:
                exit(1);
        }
    
    int num_pos_args = argc - optind;
    
    if(num_pos_args == 0){
        args.infile = malloc(2 * sizeof(char *));
        args.infile[0] = strdup("-"); // Default for input
        args.infile[1] = 0; // Signal end of file list 
    } else {
        for(int i = 0; i < num_pos_args; ++i) {
            collect_files((*argv)[i + optind]);
        }
        dedup_files(&(args.infile), &N_FILES);
    }
    return args; 
}

