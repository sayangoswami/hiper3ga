//
// Created by sayan on 2/16/19.
//

#ifndef HIPER3GA_0_1_HFILE_H
#define HIPER3GA_0_1_HFILE_H

#include <sys/mman.h>
#include "hiper3ga.h"

#define MAX_FNAME_LEN 512
typedef struct {
    char fname[MAX_FNAME_LEN];
    FILE *fp;
    struct stat st;
    void *mmaddr;
    size_t mmlen;
    bool read, write;
} hfile_t;

#define hfile_init_input_stream(hfile, fmt, ...) do {                     \
    sprintf((hfile).fname, fmt, ##__VA_ARGS__);                           \
    (hfile).fp = fopen((hfile).fname, "r");                               \
    if (!(hfile).fp)                                                      \
        LOG_ERR("%s while opening %s.", strerror(errno), (hfile).fname);  \
    stat((hfile).fname, &(hfile).st);                                     \
    (hfile).read = true, (hfile).write = false;                           \
    (hfile).mmaddr = NULL, (hfile).mmlen = 0;                             \
} while(0)

#define hfile_init_output_stream(hfile, fmt, ...) do {                    \
    sprintf((hfile).fname, fmt, ##__VA_ARGS__);                           \
    (hfile).fp = fopen((hfile).fname, "w");                               \
    if (!(hfile).fp)                                                      \
        LOG_ERR("%s while opening %s.", strerror(errno), (hfile).fname);  \
    stat((hfile).fname, &(hfile).st);                                     \
    (hfile).read = false, (hfile).write = true;                           \
    (hfile).mmaddr = NULL, (hfile).mmlen = 0;                             \
} while(0)

#define hfile_close_stream(hfile) do {                        \
    if((hfile).fp) {fclose((hfile).fp); (hfile).fp = NULL;}   \
    if((hfile).mmaddr) hfile_unmap(hfile);                    \
    (hfile).read = false, (hfile).write = false;              \
} while(0)

#define hfile_read(hfile, buf, type, count) fread((buf), sizeof(type), (count), (hfile).fp)

#define hfile_write(hfile, buf, type, count) (fwrite((buf), sizeof(type), (count), (hfile).fp) == count)

#define hfile_forwrd(hfile, nbytes) (fseek((hfile).fp, (nbytes), SEEK_CUR))

#define hfile_reset(hfile, nbytes) (fseek((hfile).fp, (nbytes), SEEK_SET) == 0)

#define hfile_bytes(hfile) ((hfile).st.st_size)

#define hfile_unmap(hfile) do {                                     \
    if (munmap((hfile).mmaddr, (hfile).mmlen))                      \
        LOG_ERR("Could not unmap because %s.", strerror(errno));    \
    else {(hfile).mmaddr = NULL, (hfile).mmlen = 0;}                \
} while(0)

static inline void* hfile_map(hfile_t *hfile, off64_t offset, size_t *length_p) {
    if (offset >= hfile_bytes(*hfile))
        LOG_ERR("offset is past end of file");
    if (hfile->mmaddr) hfile_unmap(*hfile);
    int fd = fileno(hfile->fp);
    size_t pa_offset = offset & ~(sysconf(_SC_PAGE_SIZE) - 1);
    if (offset + *length_p > hfile_bytes(*hfile))
        *length_p = hfile_bytes(*hfile) - offset;
    hfile->mmlen = *length_p + offset - pa_offset;
    hfile->mmaddr = mmap(NULL, hfile->mmlen, PROT_READ, MAP_PRIVATE, fd, pa_offset);
    if (hfile->mmaddr == MAP_FAILED)
        LOG_ERR("Could not map because %s.", strerror(errno));
    return hfile->mmaddr + offset - pa_offset;
}

#endif //HIPER3GA_0_1_HFILE_H
