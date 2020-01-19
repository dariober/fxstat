#include <string.h>

struct args {
    char *version;
    int nx_ints[101]; 
    char *infile;
    char *outfile;
    int n_stop;
};

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
    
    // Default arguments 
    struct args args ={
        .version = "0.1.0", 
        .infile = "-",
        .outfile = "-", 
        .n_stop = -1
    };
    char nx_ints[1000] = "0,25,50,75,90";

    int index;
    int opt;
    opterr = 0;

    while ((opt = getopt (argc, *argv, "s:N:hVo:")) != -1)
        switch (opt)
        {
            case 's':
                args.n_stop = atoi(optarg);
                break;
            case 'N':
                strcpy(nx_ints, optarg);
                break;
            case 'o':
                args.outfile = optarg;
                break;
            case 'V':
                printf("%s\n", args.version);
                exit(0);
            case 'h':
                printf("Usage: fxstat [OPTION] ... [FILE]\n");
                printf("Collect sequence statistics from fastq/fasta file\n");
                printf("\n");
                printf("  -s  Stop after this many sequences [%d]\n", args.n_stop);
                printf("  -N  Comma-separated thresholds for Nx statistics [%s]\n", nx_ints);
                printf("  -o  output file [%s]\n", args.outfile);
                printf("  -V  print version\n");
                printf("\n");
                printf("\
With no FILE, or when FILE is -, read standard input.\n\
Input may be gzip'd compressed. Base qualities must\n\
be in phred scale. There is no check for that!\n");
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

    for(index = optind; index < argc; index++){
        args.infile = (*argv)[index];
    }
    tokenize_nx(nx_ints, args.nx_ints); 
    return args; 
}

