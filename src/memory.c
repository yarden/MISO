/* -*- mode: C -*-  */

#include "splicing_memory.h"

/**
 * \function splicing_free
 * Deallocate memory that was allocated by splicing functions
 * 
 * Some splicing functions return a pointer vector (splicing_vector_ptr_t)
 * containing pointers to other splicing or other data types. These data
 * types are dynamically allocated and have to be deallocated
 * manually, if the user does not need them any more. This can be done
 * by calling splicing_free on them.
 * 
 * </para><para>
 * Here is a complete example on how to use \c splicing_free properly.
 * <programlisting>
 * <![CDATA[#include <splicing.h>
 * 
 * int main(void)
 * {
 *    splicing_t graph;
 *    splicing_vector_ptr_t seps;
 *    long int i;
 * 
 *    splicing_famous(&graph, "tutte");
 *    splicing_vector_ptr_init(&seps, 0);
 *    splicing_minimum_size_separators(&graph, &seps);
 *    
 *    for (i=0; i<splicing_vector_ptr_size(&seps); i++) {
 *      splicing_vector_t *v=VECTOR(seps)[i];
 *      splicing_vector_print(v);
 *      splicing_vector_destroy(v);
 *      splicing_free(v);
 *    }
 *
 *    splicing_vector_ptr_destroy(&seps);
 *    splicing_destroy(&graph);
 *    return 0;
 * }]]>
 * </programlisting>
 * 
 * 
 * \param p Pointer to the piece of memory to be deallocated.
 * \return Error code, currently always zero, meaning success.
 * 
 * Time complexity: platform dependent, ideally it should be O(1).
 */

int splicing_free(void *p) {
  splicing_Free(p);
  return 0;
}
