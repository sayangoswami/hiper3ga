//
// Created by sayan on 9/13/19.
//

#include "hclasses.h"
static const char *TAGSTR = "CONF";

#include "ketopt.h"
#define MAX_PATH_LEN 512

#include "utils.h"

struct _hconf_t {
    uint32_t tag;
    int k, w, t, iot;
    size_t max_nbases;
    char input[MAX_PATH_LEN],
         output[MAX_PATH_LEN],
         tempdir[MAX_PATH_LEN];
};
bool sync_requests = false;
bool bi_dir = true;
int min_span = 2000;
int min_match = 100;
int min_dp = 3;
float min_iden = .05f;

int max_hang = 1000;
int min_ovlp = 2000;
float int_frac = .8f;

int gap_fuzz = 1000;
int n_rounds = 2;
int bub_dist = 50000;
int max_ext = 4;
float min_ovlp_drop_ratio = .5f;
float max_ovlp_drop_ratio = .7f;
float final_ovlp_drop_ratio = .8f;

static inline size_t hmsize2bytes(char *hsize) {
    char unit = hsize[strlen(hsize)-1];
    if (unit >= 48 && unit <= 57)
        return (size_t) strtol(hsize, (char **)NULL, 10);
    else {
        hsize[strlen(hsize)-1] = '\0';
        size_t size = strtol(hsize, (char **)NULL, 10);;
        switch (unit) {
            case 'k':
            case 'K':
                size = size KiB;
                break;
            case 'm':
            case 'M':
                size = size MiB;
                break;
            case 'g':
            case 'G':
                size = size GiB;
                break;
            default:
                LOG_ERR("Unknown unit %c.", unit);
        }
        return size;
    }
}

static void parse_args(hconf_t *self, int argc, char **argv, int node_id) {
    const char* phase = argv[1];
    argc--, argv++;
    self->tempdir[0] = '\0';
    self->k = 15;
    self->w = 5;
    self->t = gasnett_cpu_count()? gasnett_cpu_count() : 1;
    self->iot = -1;
    self->max_nbases = 64ULL MiB;
    static ko_longopt_t longopts[] = {
            { "tempdir", ko_required_argument, 301 },
            { "block-nbase", ko_required_argument, 302 },
            { "io-threads", ko_required_argument, 308 },
            { "sync-requests", ko_no_argument, 310},
            { NULL, 0, 0 }
    };
    ketopt_t opt = KETOPT_INIT;
    int i, c;
    while ((c = ketopt(&opt, argc, argv, 1, "kwt:", longopts)) >= 0) {
        if (c == 'k') self->k = (int) strtol(opt.arg, (char **)NULL, 10);
        else if (c == 'w') self->w = (int) strtol(opt.arg, (char **)NULL, 10);
        else if (c == 't') self->t = (int) strtol(opt.arg, (char **)NULL, 10);
        else if (c == 308) self->iot = (int) strtol(opt.arg, (char **)NULL, 10);
        else if (c == 310) sync_requests = true;
        else if (c == 301) strcpy(self->tempdir, opt.arg);
        else if (c == 302) {
            self->max_nbases = hmsize2bytes(opt.arg);
            self->max_nbases = (self->max_nbases > 1u GiB)? (1u GiB) : (self->max_nbases);
        }
        else if (c == '?') LOG_ERR("Unknown option -%c", opt.opt? opt.opt : ':');
        else if (c == ':') LOG_ERR("missing arg: -%c", opt.opt? opt.opt : ':');
    }
    i = opt.ind;
    if (i != argc-2) LOG_ERR("Usage: Hiper3ga [options] <input> <output>");
    strcpy(self->input, argv[i++]);
    strcpy(self->output, argv[i]);
    if (self->tempdir[0] == '\0')
        strcpy(self->tempdir, "/tmp/hiper3ga");
    if (streq(phase, "align")) {
        if (node_id == 0) {
            if (!dir_exisis(self->tempdir)) {
                if (mkdir(self->tempdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
                    LOG_ERR("Could not create temp directory %s because %s (%d)", self->tempdir, strerror(errno),
                            errno);
            } else
                LOG_DEBUG(2, "Temp directory %s exists.", self->tempdir);
        }
    }
    anonymous_barrier;
    sprintf(self->tempdir + strlen(self->tempdir), "/%d", node_id);
}

hconf_t *hconf_create(int argc, char **argv, hctx_t *ctx) {
    hconf_t *self = hmalloc(hconf_t, 1);
    self->tag = *((uint32_t *)TAGSTR);
    parse_args(self, argc, argv, hctx_mynode(ctx));

    if (streq(argv[1], "align")) {
        if (!dir_exisis(self->tempdir)) {
            if (mkdir(self->tempdir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) != 0)
                LOG_ERR("Could not create temp sub-directory %s because %s (%d)", self->tempdir, strerror(errno),
                        errno);
        } else
            LOG_DEBUG(2, "Temp sub-directory %s exists.", self->tempdir);
    }

    return self;
}

void hconf_destroy(hconf_t **self_p) {
    if (self_p && *self_p) {
        hconf_t *self = *self_p;
        assert(hconf_is(self));
        hfree(self);
        *self_p = NULL;
    }
}

bool hconf_is(void *self) {
    assert(self);
    return ((hconf_t*)self)->tag == *((uint32_t *)TAGSTR);
}

char* hconf_input(hconf_t *self) {
    assert(hconf_is(self));
    return self->input;
}

char* hconf_output(hconf_t *self) {
    assert(hconf_is(self));
    return self->output;
}

size_t hconf_max_bases_per_batch(hconf_t *self) {
    assert(hconf_is(self));
    return self->max_nbases;
}

char* hconf_tempdir(hconf_t *self) {
    assert(hconf_is(self));
    return self->tempdir;
}

int hconf_nthreads(hconf_t *self) {
    assert(hconf_is(self));
    return self->t;
}

int hconf_niothreads(hconf_t *self) {
    assert(hconf_is(self));
    return self->iot;
}

int hconf_window_len(hconf_t *self) {
    assert(hconf_is(self));
    return self->w;
}

int hconf_kmer_len(hconf_t *self) {
    assert(hconf_is(self));
    return self->k;
}

bool hconf_sync_requests() {
    return sync_requests;
}
