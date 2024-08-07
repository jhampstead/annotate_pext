#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>
#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/hts.h>

#define MAX_LINE_LENGTH 1024
#define MAX_TISSUES 50 // GTEx isoform expression matrix = 43; alter this value if using a different isoform expression matrix
#define MAX_CONSEQUENCES 45 // Current VEP consequences: 41
#define MAX_CONSEQUENCE_LENGTH 40

// Structure to hold information about each transcript
typedef struct {
    char *id;
    double *expression_levels;
} Transcript;

// Structure to hold information about each gene
typedef struct {
    char *id;
    Transcript *transcripts;
    int num_transcripts;
} Gene;

// Structure to hold isoform expression matrix
typedef struct {
    Gene *genes;
    int num_genes;
    char **tissue_names;
    int num_tissues;
} ExpressionMatrix;

// Structure to hold information about each annotation
typedef struct {
    char *gene_id;
    char *transcript_id;
    char *consequence;
    char *loftee;
    bool calculate_pext;
} Annotation;

// Structure to hold information about each variant
typedef struct {
    int variant_id;
    Annotation *annotations;
    int num_annotations;
    int num_allocations;
} Variant;

// Structure to hold unique annotations over which to calculate PEXT
typedef struct {
    char *gene_id;
    char *consequence;
    char *loftee;
} UniqueAnnotation;

// Structure to hold PEXT annotation
typedef struct {
    double *mean_expression_proportion;
    int num_unique_annotations;
    UniqueAnnotation *unique_annotations;
} PextAnnotation;

int index_of(char *s, char delimiter) {
    char *occurrence = strchr(s, delimiter);
    return (int)(occurrence - s);
}

// Strip gene and transcript versions from Ensembl identifiers
char *strip_version(char *id) {
    char *underscore = strchr(id, '.');
    if (underscore != NULL) {
        *underscore = '\0';  // Truncate the string at the first .
    }
    return id;
}

ExpressionMatrix* read_expression_matrix(const char *filename, char *tissue_group, char **tissue_groups, int num_tissue_groups) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        perror("Could not open isoform expression matrix");
        return NULL;
    }

    ExpressionMatrix *matrix = malloc(sizeof(ExpressionMatrix));
    matrix->num_genes = 0;
    matrix->num_tissues = 0;
    matrix->tissue_names = NULL;
    matrix->genes = NULL;

    char line[MAX_LINE_LENGTH];
    bool header = true;

    while (fgets(line, sizeof(line), file)) {
        line[strcspn(line, "\n")] = '\0';

        if (header) {
            // The first line contains tissue names
            char *token = strtok(line, "\t");
            // First two columns contain gene_id, transcript_id
            token = strtok(NULL, "\t");
            token = strtok(NULL, "\t");

            // Count tissues
            while (token) {
                // Match tissues to tissue_groups
                if (tissue_group) {
                    for (int i = 0; i < num_tissue_groups; i++) {
                        if (strncmp(tissue_groups[i], token, strlen(tissue_groups[i])) == 0) {
                            matrix->num_tissues++;
                            matrix->tissue_names = realloc(matrix->tissue_names, matrix->num_tissues * sizeof(char *));
                            matrix->tissue_names[matrix->num_tissues - 1] = strdup(token);
                            break;  // Exit the for loop once matched
                        }
                    }
                }
                if (!tissue_group) {
                    matrix->num_tissues++;
                    matrix->tissue_names = realloc(matrix->tissue_names, matrix->num_tissues * sizeof(char *));
                    matrix->tissue_names[matrix->num_tissues - 1] = strdup(token);
                }
                token = strtok(NULL, "\t");
            }
            header = false;
            continue;
        }

        char *token = strtok(line, "\t");
        char *gene_id = strdup(strip_version(token));
        token = strtok(NULL, "\t");
        char *transcript_id = strdup(strip_version(token));

        // Store gene in matrix
        if (matrix->num_genes == 0 || strcmp(gene_id, matrix->genes[matrix->num_genes-1].id) != 0) {
            matrix->genes = realloc(matrix->genes, ++(matrix->num_genes) * sizeof(Gene));
            matrix->genes[matrix->num_genes-1] = (Gene) { gene_id, NULL, 0 };
        } else {
            free(gene_id);
        }

        // Prepare gene for storing transcripts
        Gene *gene = &(matrix->genes[matrix->num_genes-1]);
        gene->num_transcripts++;
        gene->transcripts = realloc(gene->transcripts, gene->num_transcripts * sizeof(Transcript));

        // Create transcript
        Transcript transcript = {transcript_id, NULL};
        transcript.expression_levels = malloc(matrix->num_tissues * sizeof(double));
        for (int i = 0; i < matrix->num_tissues; i++) {
            token = strtok(NULL, "\t");
            transcript.expression_levels[i] = token ? atof(token) : 0.0;
        }

        gene->transcripts[gene->num_transcripts-1] = transcript;

    }

    fclose(file);
    return matrix;
}

void free_expression_matrix(ExpressionMatrix *matrix) {
    if (matrix) {
        for (int i = 0; i < matrix->num_genes; i++) {
            for (int j = 0; j < matrix->genes[i].num_transcripts; j++) {
                free(matrix->genes[i].transcripts[j].id);
                free(matrix->genes[i].transcripts[j].expression_levels);
            }

            free(matrix->genes[i].id);
            free(matrix->genes[i].transcripts);
        }
        free(matrix->genes);

        for (int i = 0; i < matrix->num_tissues; i++) {
            free(matrix->tissue_names[i]);
        }
        free(matrix->tissue_names);

        free(matrix);
    }
}

void parse_most_severe_csq(char **ordered_consequences, int num_ordered_consequences, UniqueAnnotation *unique_annotations, int num_unique_annotations, Variant *variant) {
    UniqueAnnotation most_severe_annotation;
    int min = -1;
    for (int i = 0; i < num_unique_annotations; ++i) {
        for (int j = 0; j < num_ordered_consequences; ++j) {
            if (strcmp(unique_annotations[i].consequence, ordered_consequences[j]) == 0) {
                // Set if current does not already exist
                if (min == -1) {
                    most_severe_annotation = unique_annotations[i];
                    min = j;
                    continue;
                }

                if (j >= min) continue;

                if (strcmp(unique_annotations[i].loftee, "HC") == 0) {
                    if (strcmp(most_severe_annotation.loftee, "LC") == 0 || strcmp(most_severe_annotation.loftee, "") == 0) {
                        most_severe_annotation = unique_annotations[i];
                        min = j;
                        continue;
                    } else if (strcmp(most_severe_annotation.loftee, "HC") != 0) {
                        fprintf(stderr, "Unknown LOFTEE annotation encountered: %s\n", unique_annotations[i].loftee);
                        exit(1);
                    }
                } else if (strcmp(unique_annotations[i].loftee, "LC") == 0 || strcmp(unique_annotations[i].loftee, "") == 0) {
                    if (strcmp(most_severe_annotation.loftee, "HC") == 0) continue;
                } else {
                    fprintf(stderr, "Unknown LOFTEE annotation encountered: %s\n", unique_annotations[i].loftee);
                    exit(1);
                }

                // Annotations match between current and target annotation
                if (j < min) {
                    most_severe_annotation = unique_annotations[i];
                    min = j;
                }
            }
        }
    }

    for (int i = 0; i < variant->num_annotations; i++) {
        if (variant->annotations[i].calculate_pext && strcmp(variant->annotations[i].consequence, most_severe_annotation.consequence) != 0 && strcmp(variant->annotations[i].loftee, most_severe_annotation.loftee) != 0) {
            variant->annotations[i].calculate_pext = false;
        }
    }
}

void parse_csq(const char *csq_field, int gene_idx, int transcript_idx, int consequence_idx, int biotype_idx, int cds_pos_idx, int loftee_idx, Variant *variant) {
    char *csq = strdup(csq_field);
    char *csq_copy = csq;
    int clen = 0;
    int maxlen = strlen(csq);
    while (maxlen > 0 && *csq != '\0') {
        clen = 0;
        while (csq[clen] && csq[clen] != ',' && csq[clen] != '\0') clen++;
        csq[clen] = '\0';

        int flen = 0;
        char *fields = strdup(csq);
        char *fields_copy = fields;
        Annotation annotation = {NULL, NULL, NULL, NULL, true};
        for (int i = 0; *fields && *fields != '\0'; i++) {
            flen = 0;
            while (fields[flen] && fields[flen] != '|' && fields[flen] != '\0') flen++;
            fields[flen] = '\0';

            if (i == gene_idx) annotation.gene_id = strdup(strip_version(fields));
            else if (i == transcript_idx) annotation.transcript_id = strdup(strip_version(fields));
            else if (i == consequence_idx) {
                char *end = strchr(fields, '&'); // Remove secondary/modifier consequences
                if (end != NULL) *end = '\0';
                annotation.consequence = strdup(fields);
            } else if (i == loftee_idx) annotation.loftee = strdup(fields);
            else if (i == biotype_idx && strcmp(fields, "protein_coding") != 0) annotation.calculate_pext = false; // Filter to protein_coding biotypes
            else if (i == cds_pos_idx && strlen(fields) == 0) annotation.calculate_pext = false; // Filter to variants within CDS

            fields += flen + 1;
        }
        free(fields_copy);

        // Store variant annotations
        if (annotation.gene_id && annotation.transcript_id && annotation.consequence) {
            if (variant->num_annotations >= variant->num_allocations) {
                variant->num_allocations += 5;
                variant->annotations = realloc(variant->annotations, variant->num_allocations * sizeof(Annotation));
            }
            variant->annotations[variant->num_annotations++] = annotation;
        }

        csq += clen + 1;
        maxlen -= (clen + 1);
    }

    free(csq_copy);
}

void free_variant(Variant *variant) {
    if(variant) {
        for (int i = 0; i < variant->num_annotations; i++) {
            free(variant->annotations[i].gene_id);
            free(variant->annotations[i].transcript_id);
            free(variant->annotations[i].consequence);
            free(variant->annotations[i].loftee);
        }
        free(variant->annotations);
    } 
}

void get_unique_gene_set(int num_elements, char **list, int *len, char **unique_set[]) {
    char **ids = malloc(sizeof(char *));
    int nlen = 0, mlen = 1;

    for (int i = 0; i < num_elements; i++) {
        char *id = list[i];

        bool within_set = false;
        for (int j = 0; j < nlen; j++) {
            if (strcmp(id, ids[j]) == 0) {
                within_set = true;
                break;
            }
        }

        if (within_set) continue;

        if (nlen >= mlen) {
            mlen++;
            ids = realloc(ids, mlen * sizeof(char *));
        }

        ids[nlen++] = strdup(id);
    }

    *len = nlen;  // Set the output length to the number of unique IDs
    *unique_set = ids;
}

void get_unique_annotation_set(char *gene, Variant *variant, int *len, UniqueAnnotation *unique_annotations[]) {
    UniqueAnnotation *set = malloc(sizeof(UniqueAnnotation));
    int nlen = 0;

    for (int i = 0; i < variant->num_annotations; i++) {
        if (strcmp(variant->annotations[i].gene_id, gene) != 0) continue;

        char *gene_id = variant->annotations[i].gene_id;
        char *csq = variant->annotations[i].consequence;
        char *lof = variant->annotations[i].loftee;

        bool within_set = false;
        for (int j = 0; j < nlen; j++) {
            if (strcmp(csq, set[j].consequence) == 0 && strcmp(lof, set[j].loftee) == 0) {
                within_set = true;
                break;
            } 
        }

        if (within_set) continue;

        set = realloc(set, (nlen+1) * sizeof(UniqueAnnotation));

        set[nlen].gene_id = strdup(gene_id);
        set[nlen].consequence = strdup(csq);
        set[nlen++].loftee = strdup(lof);
    }

    *len = nlen;
    *unique_annotations = set;
}

int get_matrix_gene_subset(ExpressionMatrix *matrix, int ngenes, char **unique_gene_ids, ExpressionMatrix *out) {
    ExpressionMatrix subset = {NULL, 0, matrix->tissue_names, matrix->num_tissues}; 
    subset.genes = malloc(ngenes * sizeof(Gene));

    int i;
    for (i = 0; i < ngenes; i++) {
        for (int j = 0; j < matrix->num_genes; j++) {
            if (strcmp(matrix->genes[j].id, unique_gene_ids[i]) != 0) continue;

            subset.genes[subset.num_genes++] = matrix->genes[j];
            break;
        } 

        if (subset.num_genes <= i) {
            subset.genes[subset.num_genes++] = (Gene) {"", NULL, 0};
        }
    }

    *out = subset;

    return 0;
} 


// Compute annotation-level PEXT scores
    // total_expression = sum of expression_level for all transcripts per gene in a tissue
    // base_expression = sum of expression_level for all transcripts per gene touching a variant position with a given consequence in a tissue
    // PEXT (mean_expression_proportion) = mean of base_expression/total_expression for each consequence across all included tissues
PextAnnotation *calculate_pext(char **ordered_consequences, int num_ordered_consequences, ExpressionMatrix *matrix, Variant *variant) {
    // Get unique genes
    int num_genes = 0;
    char **gene_ids = malloc(variant->num_annotations * sizeof(char *));
    for (int i = 0; i < variant->num_annotations; i++, num_genes++) {
        gene_ids[i] = variant->annotations[i].gene_id;
    }

    int num_unique_genes;
    char **unique_gene_ids;
    get_unique_gene_set(num_genes, gene_ids, &num_unique_genes, &unique_gene_ids);

    ExpressionMatrix matrix_subset;
    get_matrix_gene_subset(matrix, num_unique_genes, unique_gene_ids, &matrix_subset);

    // Initialise PextAnnotation struct
    PextAnnotation *pext = malloc(sizeof(PextAnnotation));
    if (!pext) {
        fprintf(stderr, "Memory allocation failed: PextAnnotation\n");
        exit(EXIT_FAILURE);
    }

    pext->mean_expression_proportion = malloc(0);
    pext->num_unique_annotations = 0;
    pext->unique_annotations = malloc(0);

    // Loop over unique genes
    for (int i = 0; i < num_unique_genes; i++) {
        // Calculate total expression for all transcripts for all transcripts per gene in a tissue
        double *total_expression = calloc(matrix_subset.num_tissues, sizeof(double));

        for (int j = 0; j < matrix_subset.num_tissues; j++) {
            for (int l = 0; l < matrix_subset.genes[i].num_transcripts; l++) {
                total_expression[j] += matrix_subset.genes[i].transcripts[l].expression_levels[j];
            }
        }

        // Get unique consequences per gene
        int num_unique_annotations;
        UniqueAnnotation *unique_annotations;
        get_unique_annotation_set(unique_gene_ids[i], variant, &num_unique_annotations, &unique_annotations);

        // Update bool calculate_pext to most severe consequence if option is specified
        if (num_ordered_consequences > 0) {
            parse_most_severe_csq(ordered_consequences, num_ordered_consequences, unique_annotations, num_unique_annotations, variant);
        }

        // Initalise base_expression
        double **base_expression = malloc(matrix_subset.num_tissues * sizeof(double *));
        for (int j = 0; j < matrix_subset.num_tissues; j++) {
            base_expression[j] = calloc(num_unique_annotations, sizeof(double));
        }

        // Loop over unique consequences and calculate base_expression
        for (int l = 0; l < matrix_subset.genes[i].num_transcripts; l++) {
            for (int m = 0; m < variant->num_annotations; m++) {
                if (!variant->annotations[m].calculate_pext) continue;
                if (strcmp(variant->annotations[m].transcript_id, matrix_subset.genes[i].transcripts[l].id) == 0) {
                    for (int n = 0; n < num_unique_annotations; n++) {
                        if (strcmp(variant->annotations[m].consequence, unique_annotations[n].consequence) == 0 && strcmp(variant->annotations[m].loftee, unique_annotations[n].loftee) == 0) {
                            for (int j = 0; j < matrix_subset.num_tissues; j++) {
                                base_expression[j][n] += matrix_subset.genes[i].transcripts[l].expression_levels[j];
                            }

                            break;
                        }
                    }

                    break;
                }
            }
        }

        // Calculate normalised base_expression by dividing by total_expression per tissue
        for (int j = 0; j < matrix_subset.num_tissues; j++) {
            for (int k = 0; k < num_unique_annotations; k++) {
                base_expression[j][k] = total_expression[j] > 0 ? base_expression[j][k] / total_expression[j] : 0.0;
            }
        }

        double *mean_expression_proportion = calloc(num_unique_annotations, sizeof(double));

        // Calculate mean normalised base_expression across tissues
        for (int j = 0; j < num_unique_annotations; j++) {
            double sum = 0.0;
            for (int k = 0; k < matrix_subset.num_tissues; k++) {
                sum += base_expression[k][j];
            }
            mean_expression_proportion[j] = sum / matrix_subset.num_tissues;
        }

        // Store PextAnnotation
        int num_pext_scores = pext->num_unique_annotations + num_unique_annotations;
        pext->mean_expression_proportion = realloc(pext->mean_expression_proportion, num_pext_scores * sizeof(double));
        pext->unique_annotations = realloc(pext->unique_annotations, num_pext_scores * sizeof(UniqueAnnotation));
        for (int i = 0; i < num_unique_annotations; i++) {
            int index = pext->num_unique_annotations + i;
            pext->mean_expression_proportion[index] = mean_expression_proportion[i];
            pext->unique_annotations[index] = unique_annotations[i];
        }
        pext->num_unique_annotations = num_pext_scores;

        // Free memory 
        free(total_expression);
        for (int j = 0; j < matrix_subset.num_tissues; j++) {
            free(base_expression[j]);
        }
        free(base_expression);

    }

    free(matrix_subset.genes);
    free(gene_ids);
    for (int i = 0; i < num_unique_genes; i++) {
        free(unique_gene_ids[i]);
    }
    free(unique_gene_ids);

    return pext;
}

void parse_info_tag(bcf_hdr_t *hdr, bcf1_t *rec, Variant *variant, PextAnnotation *pext) {
    kstring_t str = {0, 0, NULL};

    for (int i = 0; i < variant->num_annotations; i++) {
        Annotation a = variant->annotations[i];

        bool is_matched = false;
        for (int j = 0; j < pext->num_unique_annotations; j++) {
            UniqueAnnotation u = pext->unique_annotations[j];

            if (strcmp(a.gene_id, u.gene_id) != 0 || strcmp(a.loftee, u.loftee) != 0 || strcmp(a.consequence, u.consequence) != 0) continue;
            is_matched = true;

            if (str.l) kputc(',', &str);
            if (variant->annotations[i].calculate_pext) {
                ksprintf(&str, "%.3f", pext->mean_expression_proportion[j]);
            } else {
                ksprintf(&str, ".");
            }

            break;
        }

        if (!is_matched) {
            if (str.l) kputc(',', &str);
            ksprintf(&str, ".");
        }
    }

    bcf_update_info_string(hdr, rec, "PEXT", str.s);
    free(str.s);
}

int main(int argc, char *argv[]) {
    // Parse option arguments
    char *most_severe = NULL;  // String for --most-severe
    char *tissue_group = NULL;  // String for --tissue-group
    char *tissue_groups[MAX_TISSUES];  // Array to hold individual tissue groups
    int num_tissue_groups = 0;  // Counter for the number of tissue groups

    struct option long_options[] = {
        {"most-severe", required_argument, 0, 'm'},
        {"tissue-group", required_argument, 0, 't'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "mt:", long_options, &option_index)) != -1) {
        switch (c) {
            case 'm':
                most_severe = optarg;
                break;
            case 't':
                tissue_group = optarg;
                break;
            case '?':
                // getopt_long already printed an error message.
                break;
            default:
                abort();
        }
    }

    argc -= optind;
    argv += optind;

    // Parse non-option arguments
    if (argc != 3) {
        fprintf(stderr, "Usage: %s [--most-severe <ordered_consequences.txt>] [--tissue-group GROUP] <isoform_expression_matrix.tsv> <input.vcf> <output.vcf>\n", argv[-optind]);
        return EXIT_FAILURE;
    }

    // Print parsed options
    char *ordered_consequences[MAX_CONSEQUENCES];
    int num_ordered_consequences = 0;
    if (most_severe) {
        printf("Most severe: Only calculating PEXT scores across transcripts with the most severe consequence\n");

        FILE *file = fopen(most_severe, "r");
        if (file == NULL) {
            perror("Error opening file of ordered consequences");
            EXIT_FAILURE;
        }

        char line[MAX_CONSEQUENCE_LENGTH];  // Adjust buffer size as needed
        while (fgets(line, sizeof(line), file)) {
            line[strcspn(line, "\n")] = '\0';
            ordered_consequences[num_ordered_consequences] = strdup(line);
            num_ordered_consequences++;

            if (num_ordered_consequences >= MAX_CONSEQUENCES) {
                fprintf(stderr, "Exceeded maximum number of consequences: set MAX_CONSEQUENCES parameter\n");
                break;
            }
        }
        fclose(file);
    }

    if (tissue_group) {
        // Split tissue_group into individual tissue groups
        char *token = strtok(tissue_group, ",");
        while (token != NULL && num_tissue_groups < MAX_TISSUES) {
            tissue_groups[num_tissue_groups++] = token;
            token = strtok(NULL, ",");
        }

        printf("Tissue groups:\n");
        for (int i = 0; i < num_tissue_groups; i++) {
            printf("  %s\n", tissue_groups[i]);
        }
    }

    // Access the remaining non-option arguments
    const char *isoform_expression_matrix_file = argv[0];
    const char *input_vcf_file = argv[1];
    const char *output_vcf_file = argv[2];

    // Print input
    printf("Isoform expression matrix: %s\n", isoform_expression_matrix_file);
    printf("Input VCF: %s\n", input_vcf_file);
    printf("Output VCF: %s\n", output_vcf_file);

    // Read input data into memory
    ExpressionMatrix *matrix = read_expression_matrix(isoform_expression_matrix_file, tissue_group, tissue_groups, num_tissue_groups);
    if (!matrix) {
        fprintf(stderr, "Error opening isoform expression matrix: %s\n", argv[1]);
        return EXIT_FAILURE;
    }

    htsFile *input_vcf = bcf_open(input_vcf_file, "r");
    if (!input_vcf) {
        fprintf(stderr, "Error opening input VCF file: %s\n", argv[2]);
        return EXIT_FAILURE;
    }

    bcf_hdr_t *hdr = bcf_hdr_read(input_vcf);
    if (!hdr) {
        fprintf(stderr, "Error reading input VCF header: %s\n", argv[2]);
        bcf_close(input_vcf);
        return EXIT_FAILURE;
    }

    // Append a new INFO field to the VCF header for PEXT annotations
    const char *pext_info_tag = "##INFO=<ID=PEXT,Number=1,Type=Float,Description=\"Consequence-aware proportion expressed across transcripts, a measure of the transcriptional output affected by a variant\">";

    htsFile *output_vcf = bcf_open(output_vcf_file, "w");
    if (!output_vcf) {
        fprintf(stderr, "Failed to open output VCF file for writing: %s\n", argv[3]);
        return EXIT_FAILURE;
    }

    if (bcf_hdr_append(hdr, pext_info_tag) != 0) {
        fprintf(stderr, "Error appending PEXT INFO field to VCF header\n");
        bcf_hdr_destroy(hdr);
        bcf_close(input_vcf);
        return EXIT_FAILURE;
    }

    if (bcf_hdr_write(output_vcf, hdr) < 0) {
        fprintf(stderr, "Failed to write updated VCF header\n");
        bcf_hdr_destroy(hdr);
        bcf_close(output_vcf);
        return EXIT_FAILURE;
    }

    // Calculate PEXT scores per VCF line
    int gene_idx = -1, transcript_idx = -1, consequence_idx = -1, biotype_idx = -1, cds_pos_idx = -1, loftee_idx = -1;

    // Find the CSQ INFO field
    int csq_id = bcf_hdr_id2int(hdr, BCF_DT_ID, "CSQ");
    if (csq_id < 0) {
        fprintf(stderr, "CSQ INFO field not found in header\n");
        bcf_hdr_destroy(hdr);
        bcf_close(input_vcf);
        return EXIT_FAILURE;
    }

    const char *csq_desc = bcf_hdr_int2id(hdr, BCF_DT_ID, csq_id);

    if (csq_desc) {
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, NULL, csq_desc, NULL);
        int ret = bcf_hrec_find_key(hrec, "Description");
        char *csq_format = strstr(hrec->vals[ret], "Format: ");

        if (csq_format) {
            csq_format += 8; // Skip "Format: "
            char *csq_copy = strdup(csq_format);
            int num_fields = 0;
            char *field = strtok(csq_copy, "|");

            while (field != NULL) {
                if (strcmp(field, "Gene") == 0) gene_idx = num_fields;
                else if (strcmp(field, "Feature") == 0) transcript_idx = num_fields;
                else if (strcmp(field, "Consequence") == 0) consequence_idx = num_fields;
                else if (strcmp(field, "BIOTYPE") == 0) biotype_idx = num_fields;
                else if (strcmp(field, "CDS_position") == 0) cds_pos_idx = num_fields;
                else if (strcmp(field, "LoF") == 0) loftee_idx = num_fields;
                num_fields++;
                field = strtok(NULL, "|");
            }
            free(csq_copy);
        }
    }

    bcf1_t *rec = bcf_init();
    Variant variant;
    variant.annotations = NULL;
    variant.num_annotations = 0;
    variant.num_allocations = 0;

    while (bcf_read(input_vcf, hdr, rec) == 0) {
        bcf_unpack(rec, BCF_UN_INFO);
        char *csq = NULL;
        int n_csq_values = 0;
        bcf_get_info_string(hdr, rec, "CSQ", &csq, &n_csq_values);

        if (n_csq_values > 0) {
            parse_csq(csq, gene_idx, transcript_idx, consequence_idx, biotype_idx, cds_pos_idx, loftee_idx, &variant);

            if (variant.num_annotations == 0) continue;

            PextAnnotation *pext = calculate_pext(ordered_consequences, num_ordered_consequences, matrix, &variant);
            if (pext) {
                parse_info_tag(hdr, rec, &variant, pext);
                bcf_write(output_vcf, hdr, rec);
                bcf_clear(rec);

                // Free memory after tag is written
                free(pext->mean_expression_proportion);
                for (int i = 0; i < pext->num_unique_annotations; i++) {
                    free(pext->unique_annotations[i].gene_id);
                    free(pext->unique_annotations[i].consequence);
                    free(pext->unique_annotations[i].loftee);
                }
                free(pext->unique_annotations);
                free(pext);
            }

            free_variant(&variant);
            variant.annotations = NULL;
            variant.num_annotations = 0;
            variant.num_allocations = 0;
        }

        free(csq);
    }

    free_expression_matrix(matrix);

    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    bcf_close(input_vcf);
    bcf_close(output_vcf);

    return EXIT_SUCCESS;
}

