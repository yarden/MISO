/* -*- mode: C -*-  */

#define CONCAT2x(a,b) a ## _ ## b 
#define CONCAT2(a,b) CONCAT2x(a,b)
#define CONCAT3x(a,b,c) a ## _ ## b ## _ ## c
#define CONCAT3(a,b,c) CONCAT3x(a,b,c)
#define CONCAT4x(a,b,c,d) a ## _ ## b ## _ ## c ## _ ## d
#define CONCAT4(a,b,c,d) CONCAT4x(a,b,c,d)

#if defined(BASE_DOUBLE)
#define BASE double
#define SHORT
#define OUT_FORMAT "%g"
#define ZERO 0.0
#define ONE 1.0
#define MULTIPLICITY 1

#elif defined(BASE_LONG)
#define BASE long
#define SHORT long
#define OUT_FORMAT "%ld"
#define ZERO 0L
#define ONE 1L
#define MULTIPLICITY 1

#elif defined(BASE_CHAR)
#define BASE char
#define SHORT char
#define OUT_FORMAT "%d"
#define ZERO 0
#define ONE 1
#define MULTIPLICITY 1

#elif defined(BASE_INT)
#define BASE int
#define SHORT int
#define OUT_FORMAT "%d"
#define ZERO 0
#define ONE 1
#define MULTIPLICITY 1

#elif defined(BASE_PTR)
#define BASE void*
#define SHORT ptr
#define ZERO 0
#define MULTIPLICITY 1

#else
#error unknown BASE_ directive
#endif

#if defined(BASE_DOUBLE)
#  define FUNCTION(dir,name) CONCAT2(dir,name)
#  define TYPE(dir) CONCAT2(dir,t)
#elif defined(BASE_BOOL)
   /* Special case because stdbool.h defines bool as a macro to _Bool which would
    * screw things up */
#  define FUNCTION(a,c) CONCAT3x(a,bool,c)
#  define TYPE(dir) CONCAT3x(dir,bool,t)
#else
#  define FUNCTION(a,c) CONCAT3(a,SHORT,c)
#  define TYPE(dir) CONCAT3(dir,SHORT,t)
#endif

#if defined(HEAP_TYPE_MIN)
#define HEAPMORE <
#define HEAPMOREEQ <=
#define HEAPLESS >
#define HEAPLESSEQ >=
#undef FUNCTION
#undef TYPE
#if defined(BASE_DOUBLE)
#define FUNCTION(dir,name) CONCAT3(dir,min,name)
#define TYPE(dir) CONCAT3(dir,min,t)
#else
#define FUNCTION(a,c) CONCAT4(a,min,SHORT,c)
#define TYPE(dir) CONCAT4(dir,min,SHORT,t)
#endif
#endif

#if defined(HEAP_TYPE_MAX)
#define HEAPMORE >
#define HEAPMOREEQ >=
#define HEAPLESS <
#define HEAPLESSEQ <=
#endif

