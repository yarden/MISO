/* -*- mode: C -*-  */

#ifndef SPLICING_ERROR_H
#define SPLICING_ERROR_H

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

/* This file contains the splicing error handling.
 * Most bits are taken literally from the GSL library (with the GSL_
 * prefix renamed to SPLICING_), as I couldn't find a better way to do
 * them. */

/**
 * \section errorhandlingbasics Error handling basics
 * 
 * <para>\a splicing functions can run into various problems preventing them 
 * from normal operation. The user might have supplied invalid arguments,
 * e.g. a non-square matrix when a square-matrix was expected, or the program 
 * has run out of memory while some more memory allocation is required, etc.
 * </para>
 * 
 * <para>By default \a splicing aborts the program when it runs into an 
 * error. While this behavior might be good enough for smaller programs, 
 * it is without doubt avoidable in larger projects. Please read further
 * if your project requires more sophisticated error handling. You can 
 * safely skip the rest of this chapter otherwise.
 * </para>
 */

/**
 * \section errorhandlers Error handlers
 *
 * <para>
 * If \a splicing runs into an error - an invalid argument was supplied
 * to a function, or we've ran out of memory - the control is
 * transferred to the \emb error handler \eme function.
 * </para><para>
 * The default error handler is \ref splicing_error_handler_abort which
 * prints an error message and aborts the program.
 * </para>
 * <para>
 * The \ref splicing_set_error_handler() function can be used to set a new
 * error handler function of type \ref splicing_error_handler_t; see the
 * documentation of this type for details.
 * </para>
 * <para>
 * There are two other predefined error handler functions,
 * \ref splicing_error_handler_ignore and \ref splicing_error_handler_printignore.
 * These deallocate the temporarily allocated memory (more about this
 * later) and return with the error code. The latter also prints an
 * error message. If you use these error handlers you need to take
 * care about possible errors yourself by checking the return value of
 * (almost) every non-void \a splicing function.
 * </para><para>
 * Independently of the error handler installed, all functions in the
 * library do their best to leave their arguments
 * \em semantically unchanged if an error
 * happens. By semantically we mean that the implementation of an
 * object supplied as an argument might change, but its
 * \quote meaning \endquote in most cases does not. The rare occasions
 * when this rule is violated are documented in this manual.
 * </para>
 */

/**
 * \section errorcodes Error codes
 * 
 * <para>Every \a splicing function which can fail return a
 * single integer error code. Some functions are very simple and
 * cannot run into any error, these may return other types, or
 * \type void as well. The error codes are defined by the
 * \ref splicing_error_type_t enumeration.
 * </para>
 */

/**
 * \section writing_error_handlers Writing error handlers
 *
 * <para>
 * The contents of the rest of this chapter might be useful only
 * for those who want to create an interface to \a splicing from another
 * language. Most readers can safely skip to the next chapter.
 * </para>
 * 
 * <para>
 * You can write and install error handlers simply by defining a
 * function of type \ref splicing_error_handler_t and calling
 * \ref splicing_set_error_handler(). This feature is useful for interface
 * writers, as \a splicing will have the chance to
 * signal errors the appropriate way, eg. the R interface defines an
 * error handler which calls the <function>error()</function>
 * function, as required by R, while the Python interface has an error
 * handler which raises an exception according to the Python way.
 * </para>
 * <para> 
 * If you want to write an error handler, your error handler should
 * call \ref SPLICING_FINALLY_FREE() to deallocate all temporary memory to
 * prevent memory leaks.
 * </para>
 */

/**
 * \section error_handling_internals Error handling internals
 *
 * <para>
 * If an error happens, the functions in the library call the
 * \ref SPLICING_ERROR macro with a textual description of the error and an
 * \a splicing error code. This macro calls (through the \ref
 * splicing_error() function) the installed error handler. Another useful
 * macro is \ref SPLICING_CHECK(). This checks the return value of its
 * argument, which is normally a function call, and calls \ref
 * SPLICING_ERROR if it is not \c SPLICING_SUCCESS. 
 * </para>
 */

/** 
 * \section deallocating_memory Deallocating memory
 *
 * <para>
 * If a function runs into an error (and the program is not aborted)
 * the error handler should deallocate all temporary memory. This is
 * done by storing the address and the destroy function of all temporary
 * objects in a stack. The \ref SPLICING_FINALLY function declares an object as
 * temporary by placing its address in the stack. If an \a splicing function returns
 * with success it calls \ref SPLICING_FINALLY_CLEAN() with the
 * number of objects to remove from the stack. If an error happens
 * however, the error handler should call \ref SPLICING_FINALLY_FREE() to
 * deallocate each object added to the stack. This means that the
 * temporary objects allocated in the calling function (and etc.) will
 * be freed as well.
 * </para>
 */

/**
 * \section writing_functions_error_handling Writing \a splicing functions with
 * proper error handling
 *
 * <para>
 * There are some simple rules to keep in order to have functions
 * behaving well in erroneous situations. First, check the arguments
 * of the functions and call \ref SPLICING_ERROR if they are invalid. Second,
 * call \ref SPLICING_FINALLY on each dynamically allocated object and call
 * \ref SPLICING_FINALLY_CLEAN() with the proper argument before returning. Third, use
 * \ref SPLICING_CHECK on all \a splicing function calls which can generate errors.
 * </para>
 * <para>
 * The size of the stack used for this bookkeeping is fixed, and
 * small. If you want to allocate several objects, write a destroy
 * function which can deallocate all of these. See the
 * <filename>adjlist.c</filename> file in the
 * \a splicing source for an example.
 * </para>
 * <para> 
 * For some functions these mechanisms are simply not flexible
 * enough. These functions should define their own error handlers and
 * restore the error handler before they return.
 * </para>
 */

/**
 * \section error_handling_threads Error handling and threads
 *
 * <para>
 * It is likely that the \a splicing error handling
 * method is \em not thread-safe, mainly because of
 * the static global stack which is used to store the address of the
 * temporarily allocated objects. This issue might be addressed in a
 * later version of \a splicing.
 * </para>
 */

/**
 * \typedef splicing_error_handler_t
 * \brief Type of error handler functions.
 * 
 * This is the type of the error handler functions.
 * \param reason Textual description of the error.
 * \param file The source file in which the error is noticed.
 * \param line The number of the line in the source file which triggered
 *   the error
 * \param splicing_errno The \a splicing error code.
 */

typedef void splicing_error_handler_t (const char * reason, const char * file,
				     int line, int splicing_errno);

/**
 * \var splicing_error_handler_abort
 * \brief Abort program in case of error.
 *
 * The default error handler, prints an error message and aborts the
 * program. 
 */

extern splicing_error_handler_t splicing_error_handler_abort;

/**
 * \var splicing_error_handler_ignore
 * \brief Ignore errors.
 *
 * This error handler frees the temporarily allocated memory and returns
 * with the error code. 
 */

extern splicing_error_handler_t splicing_error_handler_ignore;

/**
 * \var splicing_error_handler_printignore
 * \brief Print and ignore errors.
 * 
 * Frees temporarily allocated memory, prints an error message to the
 * standard error and returns with the error code. 
 */

extern splicing_error_handler_t splicing_error_handler_printignore;

/**
 * \function splicing_set_error_handler
 * \brief Set a new error handler.
 *
 * Installs a new error handler. If called with 0, it installs the
 * default error handler (which is currently
 * \ref splicing_error_handler_abort). 
 * \param new_handler The error handler function to install.
 * \return The old error handler function. This should be saved and
 *   restored if \p new_handler is not needed any
 *   more.
 */

splicing_error_handler_t*
splicing_set_error_handler(splicing_error_handler_t* new_handler);

/**
 * \typedef splicing_error_type_t
 * \brief Error code type.
 * These are the possible values returned by \a splicing functions.
 * Note that these are interesting only if you defined an error handler
 * with \ref splicing_set_error_handler(). Otherwise the program is aborted 
 * and the function causing the error never returns.
 * 
 * \enumval SPLICING_SUCCESS The function successfully completed its task.
 * \enumval SPLICING_FAILURE Something went wrong. You'll almost never
 *    meet this error as normally more specific error codes are used.
 * \enumval SPLICING_ENOMEM There wasn't enough memory to allocate
 *    on the heap. 
 * \enumval SPLICING_PARSEERROR A parse error was found in a file.
 * \enumval SPLICING_EINVAL A parameter's value is invalid. Eg. negative
 *    number was specified as the number of vertices.
 * \enumval SPLICING_EXISTS A graph/vertex/edge attribute is already
 *    installed with the given name.
 * \enumval SPLICING_EINVEVECTOR Invalid vector of vertex ids. A vertex id
 *    is either negative or bigger than the number of vertices minus one.
 * \enumval SPLICING_EINVVID Invalid vertex id, negative or too big.
 * \enumval SPLICING_NONSQUARE A non-square matrix was received while a
 *    square matrix was expected.
 * \enumval SPLICING_EINVMODE Invalid mode parameter.
 * \enumval SPLICING_EFILE A file operation failed. Eg. a file doesn't exist,
 *   or the user has no rights to open it.
 * \enumval SPLICING_UNIMPLEMENTED Attempted to call an unimplemented or
 *   disabled (at compile-time) function.
 * \enumval SPLICING_DIVERGED A numeric algorithm failed to converge.
 * \enumval SPLICING_ARPACK_PROD Matrix-vector product failed.
 * \enumval SPLICING_ARPACK_NPOS N must be positive.
 * \enumval SPLICING_ARPACK_NEVNPOS NEV must be positive.
 * \enumval SPLICING_ARPACK_NCVSMALL NCV must be bigger.
 * \enumval SPLICING_ARPACK_NONPOSI Maximum number of iterations should be positive.
 * \enumval SPLICING_ARPACK_WHICHINV Invalid WHICH parameter.
 * \enumval SPLICING_ARPACK_BMATINV Invalid BMAT parameter.
 * \enumval SPLICING_ARPACK_WORKLSMALL WORKL is too small.
 * \enumval SPLICING_ARPACK_TRIDERR LAPACK error in tridiagonal eigenvalue calculation.
 * \enumval SPLICING_ARPACK_ZEROSTART Starting vector is zero.
 * \enumval SPLICING_ARPACK_MODEINV MODE is invalid.
 * \enumval SPLICING_ARPACK_MODEBMAT MODE and BMAT are not compatible.
 * \enumval SPLICING_ARPACK_ISHIFT ISHIFT must be 0 or 1.
 * \enumval SPLICING_ARPACK_NEVBE NEV and WHICH='BE' are incompatible.
 * \enumval SPLICING_ARPACK_NOFACT Could not build an Arnoldi factorization.
 * \enumval SPLICING_ARPACK_FAILED No eigenvalues to sufficient accuracy.
 * \enumval SPLICING_ARPACK_HOWMNY HOWMNY is invalid.
 * \enumval SPLICING_ARPACK_HOWMNYS HOWMNY='S' is not implemented.
 * \enumval SPLICING_ARPACK_EVDIFF Different number of converged Ritz values.
 * \enumval SPLICING_ARPACK_SHUR Error from calculation of a real Schur form.
 * \enumval SPLICING_ARPACK_LAPACK LAPACK (dtrevc) error for calculating eigenvectors.
 * \enumval SPLICING_ARPACK_UNKNOWN Unknown ARPACK error.
 * \enumval SPLICING_ENEGLOOP Negative loop detected while calculating shortest paths.
 * \enumval SPLICING_EINTERNAL Internal error, likely a bug in splicing.
 * \enumval SPLICING_EDIVZERO Big integer division by zero.
 * \enumval IGARPH_GLP_EBOUND GLPK error (GLP_EBOUND).
 * \enumval IGARPH_GLP_EROOT GLPK error (GLP_EROOT).
 * \enumval IGARPH_GLP_ENOPFS GLPK error (GLP_ENOPFS).
 * \enumval IGARPH_GLP_ENODFS GLPK error (GLP_ENODFS).
 * \enumval IGARPH_GLP_EFAIL GLPK error (GLP_EFAIL).
 * \enumval IGARPH_GLP_EMIPGAP GLPK error (GLP_EMIPGAP).
 * \enumval IGARPH_GLP_ETMLIM GLPK error (GLP_ETMLIM).
 * \enumval IGARPH_GLP_ESTOP GLPK error (GLP_ESTOP).
 * \enumval SPLICING_EATTRIBUTES Attribute handler error. The user is not 
 *   expected to find this; it is signalled if some splicing function is
 *   not using the attribute handler interface properly.
 * \enumval SPLICING_EATTRCOMBINE Unimplemented attribute combination 
 *   method for the given attribute type.
 * \enumval SPLICING_ELAPACK A LAPACK call resulted an error.
 */

typedef enum {
  SPLICING_SUCCESS       = 0,
  SPLICING_FAILURE       = 1,
  SPLICING_ENOMEM        = 2,
  SPLICING_PARSEERROR    = 3,
  SPLICING_EINVAL        = 4,
  SPLICING_EXISTS        = 5,
  SPLICING_EINVEVECTOR   = 6,
  SPLICING_EINVVID       = 7,
  SPLICING_NONSQUARE     = 8,
  SPLICING_EINVMODE      = 9,
  SPLICING_EFILE         = 10,
  SPLICING_UNIMPLEMENTED = 12,
  SPLICING_INTERRUPTED   = 13,
  SPLICING_DIVERGED      = 14,
  SPLICING_ARPACK_PROD      = 15,
  SPLICING_ARPACK_NPOS      = 16,
  SPLICING_ARPACK_NEVNPOS   = 17,
  SPLICING_ARPACK_NCVSMALL  = 18,
  SPLICING_ARPACK_NONPOSI   = 19,
  SPLICING_ARPACK_WHICHINV  = 20,
  SPLICING_ARPACK_BMATINV   = 21,
  SPLICING_ARPACK_WORKLSMALL= 22,
  SPLICING_ARPACK_TRIDERR   = 23,
  SPLICING_ARPACK_ZEROSTART = 24,
  SPLICING_ARPACK_MODEINV   = 25,
  SPLICING_ARPACK_MODEBMAT  = 26,
  SPLICING_ARPACK_ISHIFT    = 27,
  SPLICING_ARPACK_NEVBE     = 28,
  SPLICING_ARPACK_NOFACT    = 29,
  SPLICING_ARPACK_FAILED    = 30,
  SPLICING_ARPACK_HOWMNY    = 31,
  SPLICING_ARPACK_HOWMNYS   = 32,
  SPLICING_ARPACK_EVDIFF    = 33,
  SPLICING_ARPACK_SHUR      = 34,
  SPLICING_ARPACK_LAPACK    = 35,
  SPLICING_ARPACK_UNKNOWN   = 36,
  SPLICING_ENEGLOOP         = 37,
  SPLICING_EINTERNAL        = 38,
  SPLICING_ARPACK_MAXIT     = 39,
  SPLICING_ARPACK_NOSHIFT   = 40,
  SPLICING_ARPACK_REORDER   = 41,
  SPLICING_EDIVZERO         = 42,
  SPLICING_GLP_EBOUND       = 43,
  SPLICING_GLP_EROOT        = 44,
  SPLICING_GLP_ENOPFS       = 45,
  SPLICING_GLP_ENODFS       = 46,
  SPLICING_GLP_EFAIL        = 47, 
  SPLICING_GLP_EMIPGAP      = 48,
  SPLICING_GLP_ETMLIM       = 49,
  SPLICING_GLP_ESTOP        = 50,
  SPLICING_EATTRIBUTES      = 51,
  SPLICING_EATTRCOMBINE     = 52,
  SPLICING_ELAPACK          = 53
} splicing_error_type_t;

/**
 * \define SPLICING_ERROR
 * \brief Trigger an error.
 * 
 * \a splicing functions usually use this macro when they notice an error.
 * It calls
 * \ref splicing_error() with the proper parameters and if that returns 
 * the macro returns the "calling" function as well, with the error
 * code. If for some (suspicious) reason you want to call the error
 * handler without returning from the current function, call
 * \ref splicing_error() directly.
 * \param reason Textual description of the error. This should be
 *   something more descriptive than the text associated with the error
 *   code. Eg. if the error code is \c SPLICING_EINVAL,
 *   its associated text (see  \ref splicing_strerror()) is "Invalid
 *   value" and this string should explain which parameter was invalid
 *   and maybe why. 
 * \param splicing_errno The \a splicing error code.
 */

#define SPLICING_ERROR(reason,splicing_errno) \
       do { \
       splicing_error (reason, __FILE__, __LINE__, splicing_errno) ; \
       return splicing_errno ; \
       } while (0)

/**
 * \function splicing_error
 * \brief Trigger an error.
 *
 * \a splicing functions usually call this function (most often via the 
 * \ref SPLICING_ERROR macro) if they notice an error.
 * It calls the currently installed error handler function with the
 * supplied arguments. 
 *
 * \param reason Textual description of the error.
 * \param file The source file in which the error was noticed.
 * \param line The number of line in the source file which triggered the
 *   error.
 * \param splicing_errno The \a splicing error code.
 * \return the error code (if it returns)
 */

int splicing_error(const char *reason, const char *file, int line,
		 int splicing_errno);

/**
 * \function splicing_strerror
 * \brief Textual description of an error.
 * 
 * This is a simple utility function, it gives a short general textual
 * description for an \a splicing error code.
 * 
 * \param splicing_errno The \a splicing error code.
 * \return pointer to the textual description of the error code.
 */

const char* splicing_strerror(const int splicing_errno);

#define SPLICING_ERROR_SELECT_2(a,b)       ((a) != SPLICING_SUCCESS ? (a) : ((b) != SPLICING_SUCCESS ? (b) : SPLICING_SUCCESS))
#define SPLICING_ERROR_SELECT_3(a,b,c)     ((a) != SPLICING_SUCCESS ? (a) : SPLICING_ERROR_SELECT_2(b,c))
#define SPLICING_ERROR_SELECT_4(a,b,c,d)   ((a) != SPLICING_SUCCESS ? (a) : SPLICING_ERROR_SELECT_3(b,c,d))
#define SPLICING_ERROR_SELECT_5(a,b,c,d,e) ((a) != SPLICING_SUCCESS ? (a) : SPLICING_ERROR_SELECT_4(b,c,d,e))

/* Now comes the more convenient error handling macro arsenal.
 * Ideas taken from exception.{h,c} by Laurent Deniau see
 * http://cern.ch/Laurent.Deniau/html/oopc/oopc.html#Exceptions for more 
 * information. We don't use the exception handling code though.  */

struct splicing_i_protectedPtr {
  int all;
  void *ptr;
  void (*func)(void*);
};

typedef void splicing_finally_func_t (void*);

void SPLICING_FINALLY_REAL(void (*func)(void*), void* ptr);

/**
 * \function SPLICING_FINALLY_CLEAN
 * \brief Signal clean deallocation of objects.
 * 
 * Removes the specified number of objects from the stack of
 * temporarily allocated objects. Most often this is called just
 * before returning from a function.
 * \param num The number of objects to remove from the bookkeeping
 *   stack. 
 */

void SPLICING_FINALLY_CLEAN(int num); 

/**
 * \function SPLICING_FINALLY_FREE
 * \brief Deallocate all registered objects.
 *
 * Calls the destroy function for all objects in the stack of
 * temporarily allocated objects. This is usually called only from an
 * error handler. It is \em not appropriate to use it
 * instead of destroying each unneeded object of a function, as it
 * destroys the temporary objects of the caller function (and so on)
 * as well.
 */

void SPLICING_FINALLY_FREE(void);

/**
 * \function SPLICING_FINALLY_STACK_SIZE
 * \brief Returns the number of registered objects.
 *
 * Returns the number of objects in the stack of temporarily allocated
 * objects. This function is handy if you write an own splicing routine and
 * you want to make sure it handles errors properly. A properly written
 * splicing routine should not leave pointers to temporarily allocated objects
 * in the finally stack, because otherwise an \ref SPLICING_FINALLY_FREE call
 * in another splicing function would result in freeing these objects as well
 * (and this is really hard to debug, since the error will be not in that
 * function that shows erroneous behaviour). Therefore, it is advised to
 * write your own test cases and examine \ref SPLICING_FINALLY_STACK_SIZE
 * before and after your test cases - the numbers should be equal.
 */
int SPLICING_FINALLY_STACK_SIZE(void);

/**
 * \define SPLICING_FINALLY_STACK_EMPTY
 * \brief Returns true if there are no registered objects, false otherwise.
 *
 * This is just a shorthand notation for checking that
 * \ref SPLICING_FINALLY_STACK_SIZE is zero.
 */
#define SPLICING_FINALLY_STACK_EMPTY (SPLICING_FINALLY_STACK_SIZE() == 0)

/**
 * \define SPLICING_FINALLY
 * \brief Register an object for deallocation.
 * \param func The address of the function which is normally called to
 *   destroy the object.
 * \param ptr Pointer to the object itself.
 * 
 * This macro places the address of an object, together with the
 * address of its destructor in a stack. This stack is used if an
 * error happens to deallocate temporarily allocated objects to
 * prevent memory leaks.
 */

#define SPLICING_FINALLY(func,ptr) \
  SPLICING_FINALLY_REAL((splicing_finally_func_t*)(func), (ptr))

#if (defined(__GNUC__) && GCC_VERSION_MAJOR >= 3)
#  define SPLICING_UNLIKELY(a) __builtin_expect((a), 0)
#  define SPLICING_LIKELY(a)   __builtin_expect((a), 1)
#else
#  define SPLICING_UNLIKELY(a) a
#  define SPLICING_LIKELY(a)   a
#endif

/**
 * \define SPLICING_CHECK
 * \brief Check the return value of a function call.
 *
 * \param a An expression, usually a function call.
 * 
 * Executes the expression and checks its value. If this is not
 * \c SPLICING_SUCCESS, it calls \ref SPLICING_ERROR with
 * the value as the error code. Here is an example usage:
 * \verbatim SPLICING_CHECK(vector_push_back(&amp;v, 100)); \endverbatim
 * 
 * </para><para>There is only one reason to use this macro when writing 
 * \a splicing functions. If the user installs an error handler which
 * returns to the auxiliary calling code (like \ref
 * splicing_error_handler_ignore and \ref
 * splicing_error_handler_printignore), and the \a splicing function
 * signalling the error is called from another \a splicing function 
 * then we need to make sure that the error is propagated back to 
 * the auxiliary (ie. non-splicing) calling function. This is achieved
 * by using <function>SPLICING_CHECK</function> on every \a splicing
 * call which can return an error code.
 */

#define SPLICING_CHECK(a) do { \
                 int splicing_i_ret=(a); \
                 if (SPLICING_UNLIKELY(splicing_i_ret != 0)) {\
                     SPLICING_ERROR("", splicing_i_ret); \
                 } } while (0)


typedef splicing_error_handler_t splicing_warning_handler_t;

splicing_warning_handler_t*
splicing_set_warning_handler(splicing_warning_handler_t* new_handler);

extern splicing_warning_handler_t splicing_warning_handler_ignore;
extern splicing_warning_handler_t splicing_warning_handler_print;

int splicing_warning(const char *reason, const char *file, int line,
		   int splicing_errno);

#define SPLICING_WARNING(reason) \
       do { \
         splicing_warning(reason, __FILE__, __LINE__, -1); \
       } while (0)

__END_DECLS

#endif
