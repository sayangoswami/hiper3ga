//
// Created by sayan on 9/13/19.
//

#ifndef HIPER3GA_HPRELUDE_H
#define HIPER3GA_HPRELUDE_H

#define GASNET_PAR 1
#define GASNETT_THREAD_SAFE 1

#define LOG_LEVEL 3

/// uncomment the following line to use gasnet's malloc services
#define USE_DEBUG_MALLOC

/// uncomment the following line to enable sanity checks
#define SANITY_CHECKS

#include <gasnet.h>
#include <gasnet_tools.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <stddef.h>
#include <fcntl.h>
#include <time.h>
#include <errno.h>
#include <unistd.h>
#include <inttypes.h>
#include <malloc.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <stdlib.h>
#include <execinfo.h>
#include <math.h>
#include <signal.h>

/** Logging routines and macros */

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"

static char timeStr[9] = {0};

static char *time_str()
{
    time_t rawtime;
    struct tm * timeinfo;
    time ( &rawtime );
    timeinfo = localtime ( &rawtime );
    sprintf(timeStr, "%02d:%02d:%02d", timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);
    return timeStr;
}

static void print_trace() {
    void *array[16];
    int size;
    char **strings;
    size_t i;

    size = backtrace (array, 10);
    strings = backtrace_symbols (array, size);

    printf ("Stack frames (%d).\n", size);

    for (i = 0; i < size; i++)
        printf ("%s\n", strings[i]);

    free (strings);
}

#define LOG_ERR(fmt, ...) do {                                      \
      print_trace();                                                \
      fprintf(stderr, "[%s]" RED "ERROR: " fmt RESET " at %s:%i\n", \
            time_str(), ##__VA_ARGS__, __FILE__, __LINE__);         \
      fflush(stderr);                                                \
      abort();                                                        \
} while(0)

#define LOG_INFO(fmt, ...) do {                             \
      fprintf(stdout, "[%s]" GRN "INFO: " fmt RESET "\n",   \
            time_str(), ##__VA_ARGS__);                     \
      fflush(stdout);                                 \
} while(0)

#define LOG_DEBUG(lvl, fmt, ...) do {                         \
    if (lvl <= LOG_LEVEL) {                                   \
        fprintf(stdout, "[%s]" CYN "INFO: " fmt RESET "\n",   \
            time_str(), ##__VA_ARGS__);                       \
        fflush(stdout);                                       \
    }                                                         \
} while(0)

#define LOG_WARN(fmt, ...) do {                                      \
      fprintf(stdout, "[%s]" YEL "WARN: " fmt RESET " at %s:%i\n",   \
            time_str(), ##__VA_ARGS__, __FILE__, __LINE__);          \
      fflush(stdout);                                          \
} while(0)

#define GASNET_SAFE(fncall) do {                                \
      int err;                                                  \
      if ((err = fncall) != GASNET_OK) {                        \
         fprintf(stderr, "[%s]" RED "ERROR calling %s" RESET    \
               " at %s:%i\n"                                    \
               "Cause: %s (%s)\n",                              \
               time_str(), #fncall, __FILE__, __LINE__,         \
               gasnet_ErrorName(err),                           \
               gasnet_ErrorDesc(err)                            \
         );                                                     \
         fflush(stderr);                                        \
         gasnet_exit(err);                                      \
      }                                                         \
} while(0)

#define SEND_PTR(ptr) \
    (((uintptr_t)(ptr))>>32u), \
    ((uintptr_t)(ptr))
#define RECV_PTR(a0, a1) \
    ((void*)(((uint64_t)((unsigned)(a0)))<<32u | (uint64_t)((unsigned)(a1))))
#define PtrArg(i) \
    gasnet_handlerarg_t a##i##0, gasnet_handlerarg_t a##i##1
#define PtrParam(i) \
    a##i##0, a##i##1
#define Ptr(i) \
    RECV_PTR(a##i##0, a##i##1)
#define GASNET_CALL(fn, ...) \
    GASNET_SAFE(fn(__VA_ARGS__))

/** memory allocation/deallcoation utils */

#define KiB <<10u
#define MiB <<20u
#define GiB <<30u

//size_t available_memory = 4 GiB;
//
//static inline void *
//safe_malloc (size_t size, const char *file, unsigned line)
//{
//    if (size > available_memory) {
//        fprintf(stderr, "[%s]" MAG "ERROR: Cannot allocate %ld bytes (available %zd bytes)." RESET " at %s:%i\n",
//                    time_str(), size, available_memory, file, line);
//        fflush (stderr);
//        abort();
//    }
//    void *mem = calloc(1, size);
//    if (mem == NULL) {
//        fprintf (stderr, "FATAL ERROR at %s:%u\n", file, line);
//        fprintf (stderr, "OUT OF MEMORY (calloc returned NULL)\n");
//        fflush (stderr);
//        abort ();
//    }
//    available_memory -= size;
//    return mem;
//}
//
//static inline void*
//safe_realloc(void *ptr, size_t size, const char *file, unsigned line)
//{
//    if (size > available_memory) {
//        fprintf(stderr, "[%s]" MAG "ERROR: Cannot allocate %ld bytes (available %zd bytes)." RESET " at %s:%i\n",
//                time_str(), size, available_memory, file, line);
//        fflush (stderr);
//        abort();
//    }
//    void *mem = realloc(ptr, size);
//    if (mem == NULL)
//    {
//        fprintf (stderr, "FATAL ERROR at %s:%u\n", file, line);
//        fprintf (stderr, "OUT OF MEMORY (realloc returned NULL)\n");
//        fflush (stderr);
//        abort ();
//    }
//    available_memory -= size;
//    return mem;
//}
//
//static inline void
//safe_free(void *ptr)

#if GASNET_DEBUGMALLOC && defined USE_DEBUG_MALLOC
#   define hmalloc(type, count) \
        ((count)*sizeof(type) <= 2 GiB)? \
            ((type*)gasnett_debug_calloc((count), sizeof(type))) : \
                (fprintf(stderr, "[%s]" MAG "ERROR: Cannot allocate %zd bytes." RESET " at %s:%i\n", \
                    time_str(), (size_t)(count)*sizeof(type), __FILE__, __LINE__), raise(SIGABRT), NULL)
#else
#   define hmalloc(type, count) ((type*)safe_malloc((count)*sizeof(type), __FILE__, __LINE__))
#endif

#if GASNET_DEBUGMALLOC && defined USE_DEBUG_MALLOC
#   define hrealloc(ptr, type, count) \
        ((count)*sizeof(type) <= 2 GiB)? \
            ((type*)gasnett_debug_realloc((ptr), (count)*sizeof(type))) : \
                (fprintf(stderr, "[%s]" MAG "ERROR: Cannot allocate %zd bytes." RESET " at %s:%i\n", \
                    time_str(), (size_t)(count)*sizeof(type), __FILE__, __LINE__), raise(SIGABRT), NULL)
#else
#   define hrealloc(ptr, type, count) ((type*)safe_realloc((ptr), (count)*sizeof(type), __FILE__, __LINE__))
#endif

#if GASNET_DEBUGMALLOC && defined USE_DEBUG_MALLOC
#   define hfree(x) do {gasnett_debug_free(x); (x) = NULL;} while(0)
#else
#   define hfree(x) do {free((x)); (x) = NULL;} while(0)
#endif

#define anonymous_barrier (gasnet_barrier_notify(0,GASNET_BARRIERFLAG_ANONYMOUS), \
                            gasnet_barrier_wait(0,GASNET_BARRIERFLAG_ANONYMOUS))

/*
#   define hmalloc(type, count) ((type*)safe_malloc((count)*sizeof(type), __FILE__, __LINE__))
#   define hrealloc(ptr, type, count) ((type*)safe_realloc((ptr), (count)*sizeof(type), __FILE__, __LINE__))
#   define hfree(x) do {free((x)); (x) = NULL;} while(0)
*/

//#define HPAIR_NANO_SLEEP_TIME 25
//#define cpu_relax() asm volatile("pause\n": : :"memory")

#define streq(s1, s2) (!(strcmp(s1, s2)))

#endif //HIPER3GA_HPRELUDE_H
