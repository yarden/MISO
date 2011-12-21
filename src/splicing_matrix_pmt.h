/* -*- mode: C -*-  */

typedef struct TYPE(splicing_matrix) {
  TYPE(splicing_vector) data;
  long int nrow, ncol;
} TYPE(splicing_matrix);

/*---------------*/
/* Allocation    */
/*---------------*/

int FUNCTION(splicing_matrix,init)(TYPE(splicing_matrix) *m, 
				 long int nrow, long int ncol);
int FUNCTION(splicing_matrix,copy)(TYPE(splicing_matrix) *to, 
				 const TYPE(splicing_matrix) *from);
void FUNCTION(splicing_matrix,destroy)(TYPE(splicing_matrix) *m);
int FUNCTION(splicing_matrix,reserve)(TYPE(splicing_matrix) *m,
				      long int nrow, long int ncol);

/*--------------------*/
/* Accessing elements */
/*--------------------*/

/* MATRIX */
BASE FUNCTION(splicing_matrix,e)(const TYPE(splicing_matrix) *m, 
			       long int row, long int col);
BASE* FUNCTION(splicing_matrix,e_ptr)(const TYPE(splicing_matrix) *m,
				    long int row, long int col);
void FUNCTION(splicing_matrix,set)(TYPE(splicing_matrix)* m, long int row, long int col,
				 BASE value);

/*------------------------------*/
/* Initializing matrix elements */
/*------------------------------*/

void FUNCTION(splicing_matrix,null)(TYPE(splicing_matrix) *m);
void FUNCTION(splicing_matrix,fill)(TYPE(splicing_matrix) *m, BASE e);

/*------------------*/
/* Copying matrices */
/*------------------*/

void FUNCTION(splicing_matrix,copy_to)(const TYPE(splicing_matrix) *m, BASE *to);
int FUNCTION(splicing_matrix,update)(TYPE(splicing_matrix) *to, 
				   const TYPE(splicing_matrix) *from);
int FUNCTION(splicing_matrix,rbind)(TYPE(splicing_matrix) *to,
				  const TYPE(splicing_matrix) *from);
int FUNCTION(splicing_matrix,cbind)(TYPE(splicing_matrix) *to,
				  const TYPE(splicing_matrix) *from);
int FUNCTION(splicing_matrix,swap)(TYPE(splicing_matrix) *m1, TYPE(splicing_matrix) *m2);

/*--------------------------*/
/* Copying rows and columns */
/*--------------------------*/

int FUNCTION(splicing_matrix,get_row)(const TYPE(splicing_matrix) *m, 
				    TYPE(splicing_vector) *res, long int index);
int FUNCTION(splicing_matrix,get_col)(const TYPE(splicing_matrix) *m, 
				    TYPE(splicing_vector) *res, long int index);
int FUNCTION(splicing_matrix,set_row)(TYPE(splicing_matrix) *m,
				     const TYPE(splicing_vector) *v, long int index);
int FUNCTION(splicing_matrix,set_col)(TYPE(splicing_matrix) *m,
				    const TYPE(splicing_vector) *v, long int index);
int FUNCTION(splicing_matrix,select_rows)(const TYPE(splicing_matrix) *m,
					TYPE(splicing_matrix) *res, 
					const splicing_vector_t *rows);
int FUNCTION(splicing_matrix,select_cols)(const TYPE(splicing_matrix) *m,
					TYPE(splicing_matrix) *res, 
					const splicing_vector_t *cols);
int FUNCTION(splicing_matrix,select_rows_cols)(const TYPE(splicing_matrix) *m,
					TYPE(splicing_matrix) *res, 
					const splicing_vector_t *rows,
					const splicing_vector_t *cols);

/*-----------------------------*/
/* Exchanging rows and columns */
/*-----------------------------*/

int FUNCTION(splicing_matrix,swap_rows)(TYPE(splicing_matrix) *m, 
				      long int i, long int j);
int FUNCTION(splicing_matrix,swap_cols)(TYPE(splicing_matrix) *m, 
				      long int i, long int j);
int FUNCTION(splicing_matrix,swap_rowcol)(TYPE(splicing_matrix) *m,
				       long int i, long int j);
int FUNCTION(splicing_matrix,transpose)(TYPE(splicing_matrix) *m);

/*-----------------------------*/
/* Matrix operations           */
/*-----------------------------*/

int FUNCTION(splicing_matrix,add)(TYPE(splicing_matrix) *m1, 
				const TYPE(splicing_matrix) *m2);
int FUNCTION(splicing_matrix,sub)(TYPE(splicing_matrix) *m1, 
				const TYPE(splicing_matrix) *m2);
int FUNCTION(splicing_matrix,mul_elements)(TYPE(splicing_matrix) *m1, 
					 const TYPE(splicing_matrix) *m2);
int FUNCTION(splicing_matrix,div_elements)(TYPE(splicing_matrix) *m1, 
					 const TYPE(splicing_matrix) *m2);
void FUNCTION(splicing_matrix,scale)(TYPE(splicing_matrix) *m, BASE by);
void FUNCTION(splicing_matrix,add_constant)(TYPE(splicing_matrix) *m, BASE plus);

/*-----------------------------*/
/* Finding minimum and maximum */
/*-----------------------------*/

double FUNCTION(splicing_matrix,min)(const TYPE(splicing_matrix) *m);
double FUNCTION(splicing_matrix,max)(const TYPE(splicing_matrix) *m);
int FUNCTION(splicing_matrix,which_min)(const TYPE(splicing_matrix) *m,
				      long int *i, long int *j);
int FUNCTION(splicing_matrix,which_max)(const TYPE(splicing_matrix) *m,
				      long int *i, long int *j);
int FUNCTION(splicing_matrix,minmax)(const TYPE(splicing_matrix) *m,
				   BASE *min, BASE *max);
int FUNCTION(splicing_matrix,which_minmax)(const TYPE(splicing_matrix) *m,
					 long int *imin, long int *jmin,
					 long int *imax, long int *jmax);

/*-------------------*/
/* Ordering columns  */
/*-------------------*/

int FUNCTION(splicing_matrix,order_cols)(const TYPE(splicing_matrix) *m, 
					 splicing_vector_int_t *order);
int FUNCTION(splicing_matrix,binorder_cols)(const TYPE(splicing_matrix) *m, 
					    splicing_vector_int_t *order);
int FUNCTION(splicing_matrix,sort_cols)(TYPE(splicing_matrix) *m);

/*-------------------*/
/* Matrix properties */
/*-------------------*/

int FUNCTION(splicing_matrix,isnull)(const TYPE(splicing_matrix) *m);
int FUNCTION(splicing_matrix,empty)(const TYPE(splicing_matrix) *m);
long int FUNCTION(splicing_matrix,size)(const TYPE(splicing_matrix) *m);
long int FUNCTION(splicing_matrix,nrow)(const TYPE(splicing_matrix) *m);
long int FUNCTION(splicing_matrix,ncol)(const TYPE(splicing_matrix) *m);
int FUNCTION(splicing_matrix,is_symmetric)(const TYPE(splicing_matrix) *m);
BASE FUNCTION(splicing_matrix,sum)(const TYPE(splicing_matrix) *m);
BASE FUNCTION(splicing_matrix,prod)(const TYPE(splicing_matrix) *m);
int FUNCTION(splicing_matrix,rowsum)(const TYPE(splicing_matrix) *m,
				   TYPE(splicing_vector) *res);
int FUNCTION(splicing_matrix,colsum)(const TYPE(splicing_matrix) *m,
				   TYPE(splicing_vector) *res);
int FUNCTION(splicing_matrix,is_equal)(const TYPE(splicing_matrix) *m1, 
					       const TYPE(splicing_matrix) *m2);
BASE FUNCTION(splicing_matrix,maxdifference)(const TYPE(splicing_matrix) *m1,
						    const TYPE(splicing_matrix) *m2);

/*------------------------*/
/* Searching for elements */
/*------------------------*/

int FUNCTION(splicing_matrix,contains)(const TYPE(splicing_matrix) *m,
					       BASE e);
int FUNCTION(splicing_matrix,search)(const TYPE(splicing_matrix) *m,
					     long int from, BASE what, 
					     long int *pos, 
					     long int *row, long int *col);

/*------------------------*/
/* Resizing operations    */
/*------------------------*/

int FUNCTION(splicing_matrix,resize)(TYPE(splicing_matrix) *m, 
				   long int nrow, long int ncol);
int FUNCTION(splicing_matrix,add_cols)(TYPE(splicing_matrix) *m, long int n);
int FUNCTION(splicing_matrix,add_rows)(TYPE(splicing_matrix) *m, long int n);
int FUNCTION(splicing_matrix,remove_col)(TYPE(splicing_matrix) *m, long int col);
int FUNCTION(splicing_matrix,remove_row)(TYPE(splicing_matrix) *m, long int row);

/*------------------------*/
/* Print as text          */
/*------------------------*/

int FUNCTION(splicing_matrix,print)(const TYPE(splicing_matrix) *m);
int FUNCTION(splicing_matrix,fprint)(const TYPE(splicing_matrix) *m,
				   FILE *file);

#ifdef BASE_COMPLEX

int splicing_matrix_complex_real(const splicing_matrix_complex_t *v, 
			       splicing_matrix_t *real);
int splicing_matrix_complex_imag(const splicing_matrix_complex_t *v, 
			       splicing_matrix_t *imag);
int splicing_matrix_complex_realimag(const splicing_matrix_complex_t *v, 
				   splicing_matrix_t *real, 
				   splicing_matrix_t *imag);
int splicing_matrix_complex_create(splicing_matrix_complex_t *v,
				 const splicing_matrix_t *real,
				 const splicing_matrix_t *imag);
int splicing_matrix_complex_create_polar(splicing_matrix_complex_t *v,
				       const splicing_matrix_t *r,
				       const splicing_matrix_t *theta);

#endif

/* ----------------------------------------------------------------------------*/
/* For internal use only, may be removed, rewritten ... */
/* ----------------------------------------------------------------------------*/

int FUNCTION(splicing_matrix,permdelete_rows)(TYPE(splicing_matrix) *m, 
					    long int *index, long int nremove);
int FUNCTION(splicing_matrix,delete_rows_neg)(TYPE(splicing_matrix) *m, 
					    const splicing_vector_t *neg, 
					    long int nremove);

