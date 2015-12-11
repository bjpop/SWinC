/* This module will be the main module in which
 * the SWinC aligner (Smith-Waterman alignmer) will
 * be implemented. 
 * In addition to the vanilla Smith-Waterman 
 * (bottom-up Dynamic Programming with
 * gap-scoring scheme), the algorithm will be modified
 * to suit its use in the context of DNA nucleotide 
 * sequence alignment. In particular, thermodynamic
 * consideration would be incorporated into the
 * scoring scheme. 
 *
 * The algorithm is described in the following paper:
 * Kaderali, L, "Primer Design for Multiplexed Genotyping",
 * Methods In Molecular Biology, Humana Press, Totowa,
 * New Jersy, pg. 269-285
 */

#include <stdio.h>
#include <string.h>
#include <ctype.h>

/* define maximum primer + heel length acceptable */
#define MAX_SEQ_LEN 60

/* produce record of nearest neighbour thermodynamic
 * data here. It should be in "hash table" like 
 * structure
 */


// XXX should we use another symbol for mismatch? maybe X?

/* Record score and decision at every sw_matrix entry
 * the decision are 
 *   'M' for match
 *   'm' for mismatch
 *   'D' for deletion
 *   'I' for insertion
 *   '\0' for termimate alignment */
typedef struct {
    float score;
    char decision;
} SW_entry;

float swalign(char *ref, char *query);
float fill_matrix(int ref_len, int query_len);
float score(int row, int col);
SW_entry score_mm(int row, int col);
SW_entry score_insert(int row, int col);
SW_entry score_delete(int row, int col);
float penalise_gap(int gap_len);
SW_entry max(SW_entry list[], int list_len);
void print_sw_matrix(int nrow, int ncol);

void print_interaction_matrix(int nrow, int ncol);
void align_pool(int pool_size);
int get_primers(char *filename);

void rev_complement(char *seq, char *new);

// XXX These should either be const or preferably #define constants
// use UPPER_CASE_NAMES for constants
float match_score = 2.0;
float mismatch_penalty = -1.0;
float gap_open = -1.0;
float gap_extension = 0.0;

// XXX Why are these declared static?
static SW_entry null_entry = {0.0, '\0'};
/* initialise sw_matrix as (MAX_SEQ_LEN+1) x (MAX_SEQ_LEN+1)
 * all entries will be initialised with SW_entry = {0.0,'\0'} */
static SW_entry sw_matrix[MAX_SEQ_LEN +1][MAX_SEQ_LEN +1];


// XXX I think you mean "global" instead of "external",
// XXX "external" means something else in C
/* external variable recording the reference and query sequence
 * visible to all the alignment/scoring related functions */
char *g_ref;
char *g_query;


/******************************************************************************/
#define MAX_POOL_SIZE 100
float interaction_matrix[MAX_POOL_SIZE][MAX_POOL_SIZE];
char pool[MAX_POOL_SIZE][MAX_SEQ_LEN];
/* main() */
int main(int argc, char **argv){
    // XXX check the value of argc before indexing into argv
    char *ref = argv[1];
    char *query = argv[2];
    int ref_len = strlen(ref);
    int query_len = strlen(query);
    float align_score;
    printf("ref=%s\nquery=%s\nref_len=%i query_len=%i\n", ref, query, ref_len, query_len);
    align_score = swalign(ref, query);
    // XXX put single space on either side of plus sign: query_len + 1
    print_sw_matrix(query_len +1, ref_len +1);
    printf("Alignment score = %.2f\n\n", align_score);

    // XXX Do all your command line argument processing in one place.
    // XXX Consider using getopts 
    char *filename = argv[3];
    int num_primer;
    num_primer = get_primers(filename);
    printf("num_primer = %i\n", num_primer);

    align_pool(num_primer);
    print_interaction_matrix(num_primer, num_primer);
    return 0;
}

// XXX use of global pool variable should probably be avoided
int get_primers(char *filename){
    FILE *file_handle = fopen(filename, "r");
    // XXX test if the fopen succeeded or failed.
    int i = 0;
    // XXX Should this be MAX_SEQ_LEN + 1 to accommodate the terminating null byte?
    char buff[MAX_SEQ_LEN];
    char *temp;
    while ((temp = fgets(buff, MAX_SEQ_LEN, file_handle)) != NULL){
        strcpy(pool[i], temp);
        i++;
    }
    return i;
}

/* align_pool:
 * given an array of strings, pool, return a matrix of size
 * |pool|x|pool| with matrix[i][j] representing the alignment score
 * between string_i and string_j (a complete graph).
 * The matrix will be symmetric and the diagonal element will be 0.0, i.e.
 * by default we don't want information about self-alignment.*/
void align_pool(int pool_size){
    int i, j;
    // XXX Put braces around all blocks, including those that only
    // XXX contain one statement.
    // XXX Put a single space on either side of all binary operators
    // XXX including the ? 
    for (i = 0; i < pool_size; i++)
        for (j = 0; j < pool_size; j++)
            interaction_matrix[i][j] = (i>=j)? 0.0 : swalign(pool[i], pool[j]);
}

/* swalign: produce the best alignment score for the 2 input string */
float swalign(char *ref, char *query){
    // XXX avoid globals here if possible
    g_ref = ref;
    g_query = query;
    float best_score;
    int ref_len = strlen(ref);
    int query_len = strlen(query);

    best_score = fill_matrix(ref_len, query_len);
    return best_score;
}

/* fill_matrix:
 * fill the matrix with the the best scores and decisions that lead
 * to those scores. It return the best scored found in the matrix */
float fill_matrix(int ref_len, int query_len){
    register int row, col;
    float best_score = 0.0;
    float new_score;
    for (row = 1; row < query_len +1; row++)
        for (col = 1; col < ref_len +1; col++){
            // XXX Pass the matrix as an argument instead of using a global
            new_score = score(row, col);
            best_score = (new_score > best_score)? new_score : best_score;
        }
    return best_score;
}


/* score:
 * for the specifed row and col in sw_matrix
 * assign the best  */
float score(int row, int col){
    SW_entry best_choice;
    // XXX I think you can omit the array size (4) if you define its elements statically?
    SW_entry choices[4] = 
                       {null_entry,
                        score_mm(row, col),
                        score_insert(row, col),
                        score_delete(row, col)};
    best_choice = max(choices, 4);
    sw_matrix[row][col] = best_choice;
    return best_choice.score;
}

// XXX comment on what this is doing
// XXX maybe use more meaningful name than score_mm, perhaps
// XXX score_match_mismatch?
SW_entry score_mm(int row, int col){
    float prefix_score = sw_matrix[row-1][col-1].score;
    SW_entry mm_entry;
    if (g_ref[col -1] == g_query[row -1]){
        mm_entry.score = prefix_score + match_score;
        mm_entry.decision = 'M';
    } else{
        mm_entry.score = prefix_score + mismatch_penalty;
        mm_entry.decision = 'm';
    }; // XXX why is there a semicolon here?
    return mm_entry;
}

SW_entry score_insert(int row, int col){
    int gap_len;
    SW_entry insert_entry = {0.0, 'I'};
    float new_score;
    // XXX I don't understand the need for a loop here
    // XXX can this be computed in constant time?
    for (gap_len = col; gap_len > 0; gap_len--){
        new_score = sw_matrix[row][col-gap_len].score + penalise_gap(gap_len);
        if (new_score > insert_entry.score)
            insert_entry.score = new_score;
    }; // XXX why is there a semicolon here?
    return insert_entry;
}

SW_entry score_delete(int row, int col){
    int gap_len;
    SW_entry delete_entry = {0.0, 'D'};
    float new_score;
    // XXX I don't understand the need for a loop here
    // XXX can this be computed in constant time?
    for (gap_len = row; gap_len > 0; gap_len--){
        new_score = sw_matrix[row-gap_len][col].score + penalise_gap(gap_len);
        if (new_score > delete_entry.score)
            delete_entry.score = new_score;
    }; // XXX why is there a semicolon here?
    return delete_entry;
}

float penalise_gap(int gap_len){
    // XXX probably could avoid the variable assignment and return
    // XXX the result directly?
    float penalty = gap_open + gap_extension * gap_len;
    return penalty;
}




/***** Utilities *********/

/* max:
 * return the entry in list[] with maximal score */
SW_entry max(SW_entry list[], int list_len){
    SW_entry maximum = list[0]; // assume list_len > 0
    register int i;
    for (i = 1; i < list_len; i++)
        if (list[i].score > maximum.score)
            maximum = list[i];
    return maximum;
}


void print_sw_matrix(int nrow, int ncol){
    int row, col;
    SW_entry entry;
    for (row = 0; row < nrow; row++){
        for (col = 0; col < ncol; col++){
            entry = sw_matrix[row][col];
            printf("%.2f%c ", entry.score, entry.decision);
        };
        putchar('\n');
    }
}

void print_interaction_matrix(int nrow, int ncol){
    register int row, col;
    for (row = 0; row < nrow; row++){
        for (col = 0; col < ncol; col++){
            printf("%.2f ", interaction_matrix[row][col]);
        }
        putchar('\n');
    }
}

char complement(char base){
    base = toupper(base);
    switch (base){
        case 'A':
            return 'T';
            break;
        case 'T':
            return 'A';
            break;
        case 'G':
            return 'C';
            break;
        case 'C':
            return 'G';
            break;
        default:
            return '\0';
            break;
    }
}


// XXX new must have sufficient memory to hold the data
// XXX maybe it should be malloced in this function instead?
void rev_complement(char *seq, char *new){
    int len = strlen(seq);
    int i;
    for (i = 0; i < len; i++)
        new[i] = complement(seq[len-i-1]);
    new[len] = '\0';
}

