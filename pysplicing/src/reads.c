
#include "splicing.h"
#include "splicing_memory.h"

void splicing_reads_destroy(splicing_reads_t *reads) {
  splicing_vector_int_destroy(&reads->pos);
  splicing_strvector_destroy(&reads->cigar);
}

void splicing_replicate_reads_destroy(splicing_replicate_reads_t *reads) {
  int i, n = splicing_vector_ptr_size(&reads->reads);
  for (i = 0; i < n; i++) {
    splicing_reads_t *it = VECTOR(reads->reads)[i];
    if (!it) break;
    splicing_reads_destroy(it);
    splicing_Free(it);
  }
  splicing_vector_ptr_destroy(&reads->reads);
}

const splicing_vector_int_t *
splicing_replicate_reads_pos(const splicing_replicate_reads_t *reads,
			     int rep_num) {
  splicing_reads_t *reads1 = VECTOR(reads->reads)[rep_num];
  return &reads1->pos;
}

const char **
splicing_replicate_reads_cigar(const splicing_replicate_reads_t *reads,
			       int rep_num) {

  splicing_reads_t *reads1 = VECTOR(reads->reads)[rep_num];
  return (const char **) reads1->cigar.table;
}
