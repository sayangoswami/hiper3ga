//
// Created by sayan on 9/22/19.
//

#ifndef HIPER3GA_UTILS_H
#define HIPER3GA_UTILS_H

static inline int nearest_power2(int v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

static inline unsigned int hashint32(unsigned int x) {
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    return x;
}

#include <dirent.h>
static inline bool dir_exisis(const char *path) {
    DIR* dir = opendir(path);
    if (dir) {
        /* Directory exists. */
        closedir(dir);
        return true;
    } else if (ENOENT == errno) {
        /* Directory does not exist. */
        return false;
    } else {
        /* opendir() failed for some other reason. */
        LOG_ERR("Could not verify the presence or absence of dir %s because %s (%d).", path, strerror(errno), errno);
    }
}

#define UTILS_INIT(name, type_t, __sort_lt)                 \
    /**
     * Returns the number of elements in arr that are less than or equal to key
     * @param arr - sorted (non-decreasing) array
     * @param n - number of elements in array
     * @param key - element to query for
     * @return The largest index where key can be inserted in arr without violating the non-decreasing property of arr
     */\
    static inline size_t upper_bound_##name(type_t *arr, const size_t n, type_t key) {      \
        size_t start = 0, end = n, i;                       \
        while(start < end) {                                \
            i = (start + end) / 2;                          \
            if (__sort_lt(key, arr[i])) end = i;            \
            else start = i + 1;                             \
        }                                                   \
        return start;                                       \
    }                                                       \
    /**
     * Returns the number of elements in arr that are strictly less than key
     * @param arr - sorted (non-decreasing) array
     * @param n - number of elements in array
     * @param key - element to query for
     * @return The smallest index where key can be inserted in arr without violating the non-decreasing property of arr
     */\
    static inline size_t lower_bound_##name(type_t *arr, const size_t n, type_t key) {      \
        size_t start = 0, end = n, i;                       \
        while(start < end) {                                \
            i = (start + end) / 2;                          \
            if (__sort_lt(arr[i], key)) start = i + 1;      \
            else end = i;                                   \
        }                                                   \
        return start;                                       \
    }                                                       \
    /**
     * merges 2 sorted partitions into a single partition
     * @param beg_A
     * @param end_A
     * @param beg_B
     * @param end_B
     * @param result
     */\
    static inline void two_way_merge_##name(                                                    \
            type_t *beg_A, type_t *end_A, type_t *beg_B, type_t *end_B, type_t *result) {       \
        type_t *A = beg_A, *B = beg_B, *C = result;         \
        while (A < end_A && B < end_B) {                    \
            if (__sort_lt(*A, *B)) *C = *A, A++;            \
            else *C = *B, B++;                              \
            C++;                                            \
        }                                                   \
        if (A != end_A)                                     \
            memcpy(C, A, sizeof(type_t) * (end_A - A));     \
        else if (B != end_B)                                \
            memcpy(C, B, sizeof(type_t) * (end_B - B));     \
    }                                                       \

#endif //HIPER3GA_UTILS_H
