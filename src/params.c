#include "params.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static const char *BC_NAMES[BC_COUNT]   = { "free", "pinned", "clamped", "periodic" };
static const char *DIS_NAMES[DIS_COUNT] = { "none", "light", "heavy", "random" };

const char *bc_name(Boundary bc)       { return (bc >= 0 && bc < BC_COUNT) ? BC_NAMES[bc] : "?"; }
const char *disorder_name(Disorder d)  { return (d >= 0 && d < DIS_COUNT) ? DIS_NAMES[d] : "?"; }

void params_default(Params *p)
{
    *p = (Params){
        .cols = 64, .rows = 64,
        .mass = 1.0, .k = 400.0, .kdiag = -1.0, .diagonal = 1,
        .damping = 0.4, .dt = 5e-4, .speed = 1.0,
        .bc = BC_FREE, .prestrain = 0.0,
        .disorder = DIS_NONE, .defect_light = 0.2, .defect_heavy = 5.0,
        .disorder_width = 0.5,
        .absorb = 0, .absorb_width = 10, .absorb_gamma = -1.0,
        .fput_strain = 0.15,
        .nt = 256,
        .seed = 0,
        .win_w = 0, .win_h = 0,
    };
}

static void usage(const char *argv0)
{
    printf(
"usage: %s [options]\n"
"\n"
"Lattice (lattice units: spacing a = 1)\n"
"  -c, --cols N            columns                         (64)\n"
"  -r, --rows N            rows                            (64)\n"
"  -m, --mass M            node mass, kg                   (1.0)\n"
"  -k, --k K               cardinal spring constant        (400)\n"
"      --kdiag K           diagonal spring constant        (= k)\n"
"      --no-diag           start with diagonal springs off\n"
"  -b, --bc NAME           free | pinned | clamped | periodic   (free)\n"
"      --prestrain E       L0 = (1-E)*a; acts under clamped/periodic (0)\n"
"\n"
"Dynamics\n"
"  -g, --damping G         uniform drag, kg/s              (0.4)\n"
"      --dt DT             integration step, s             (5e-4)\n"
"  -s, --speed S           sim seconds per wall second     (1.0)\n"
"\n"
"Disorder and boundaries\n"
"      --disorder NAME     none | light | heavy | random   (none)\n"
"      --defect-light F    light defect mass ratio         (0.2)\n"
"      --defect-heavy F    heavy defect mass ratio         (5.0)\n"
"      --disorder-width W  random masses m*(1+W*U[-1,1])   (0.5)\n"
"      --absorb            start with the absorbing layer on\n"
"      --absorb-width N    absorbing layer width, nodes    (10)\n"
"      --absorb-gamma G    peak layer drag, kg/s           (sqrt(k*m))  \n"
"\n"
"Analysis\n"
"      --fput-strain S     peak strain of the F/Shift-F mode (0.15)\n"
"      --nt N              time samples for omega-q, pow2  (256)\n"
"      --seed N            RNG seed                        (time)\n"
"\n"
"Window\n"
"      --width W           initial window width, points    (auto)\n"
"      --height H          initial window height, points   (auto)\n"
"  -h, --help\n",
        argv0);
}

static int parse_int(const char *s, int *out)
{
    char *end; long v = strtol(s, &end, 10);
    if (*s == '\0' || *end != '\0') return -1;
    *out = (int)v; return 0;
}

static int parse_dbl(const char *s, double *out)
{
    char *end; double v = strtod(s, &end);
    if (*s == '\0' || *end != '\0' || !isfinite(v)) return -1;
    *out = v; return 0;
}

enum {
    O_KDIAG = 256, O_NODIAG, O_PRESTRAIN, O_DT, O_DISORDER, O_DLIGHT, O_DHEAVY,
    O_DWIDTH, O_ABSORB, O_AWIDTH, O_AGAMMA, O_FPUT, O_NT, O_SEED, O_WIDTH, O_HEIGHT
};

int params_parse(Params *p, int argc, char **argv)
{
    static const struct option opts[] = {
        { "cols",           required_argument, 0, 'c' },
        { "rows",           required_argument, 0, 'r' },
        { "mass",           required_argument, 0, 'm' },
        { "k",              required_argument, 0, 'k' },
        { "kdiag",          required_argument, 0, O_KDIAG },
        { "no-diag",        no_argument,       0, O_NODIAG },
        { "bc",             required_argument, 0, 'b' },
        { "prestrain",      required_argument, 0, O_PRESTRAIN },
        { "damping",        required_argument, 0, 'g' },
        { "dt",             required_argument, 0, O_DT },
        { "speed",          required_argument, 0, 's' },
        { "disorder",       required_argument, 0, O_DISORDER },
        { "defect-light",   required_argument, 0, O_DLIGHT },
        { "defect-heavy",   required_argument, 0, O_DHEAVY },
        { "disorder-width", required_argument, 0, O_DWIDTH },
        { "absorb",         no_argument,       0, O_ABSORB },
        { "absorb-width",   required_argument, 0, O_AWIDTH },
        { "absorb-gamma",   required_argument, 0, O_AGAMMA },
        { "fput-strain",    required_argument, 0, O_FPUT },
        { "nt",             required_argument, 0, O_NT },
        { "seed",           required_argument, 0, O_SEED },
        { "width",          required_argument, 0, O_WIDTH },
        { "height",         required_argument, 0, O_HEIGHT },
        { "help",           no_argument,       0, 'h' },
        { 0, 0, 0, 0 }
    };

    int seed_given = 0, c, bad = 0;
    while ((c = getopt_long(argc, argv, "c:r:m:k:b:g:s:h", opts, NULL)) != -1) {
        switch (c) {
        case 'c':  bad |= parse_int(optarg, &p->cols); break;
        case 'r':  bad |= parse_int(optarg, &p->rows); break;
        case 'm':  bad |= parse_dbl(optarg, &p->mass); break;
        case 'k':  bad |= parse_dbl(optarg, &p->k); break;
        case O_KDIAG:     bad |= parse_dbl(optarg, &p->kdiag); break;
        case O_NODIAG:    p->diagonal = 0; break;
        case 'b': {
            int found = 0;
            for (int i = 0; i < BC_COUNT; i++)
                if (strcmp(optarg, BC_NAMES[i]) == 0) { p->bc = (Boundary)i; found = 1; }
            if (!found) { fprintf(stderr, "unknown --bc '%s'\n", optarg); bad = -1; }
            break; }
        case O_PRESTRAIN: bad |= parse_dbl(optarg, &p->prestrain); break;
        case 'g':  bad |= parse_dbl(optarg, &p->damping); break;
        case O_DT:        bad |= parse_dbl(optarg, &p->dt); break;
        case 's':  bad |= parse_dbl(optarg, &p->speed); break;
        case O_DISORDER: {
            int found = 0;
            for (int i = 0; i < DIS_COUNT; i++)
                if (strcmp(optarg, DIS_NAMES[i]) == 0) { p->disorder = (Disorder)i; found = 1; }
            if (!found) { fprintf(stderr, "unknown --disorder '%s'\n", optarg); bad = -1; }
            break; }
        case O_DLIGHT:    bad |= parse_dbl(optarg, &p->defect_light); break;
        case O_DHEAVY:    bad |= parse_dbl(optarg, &p->defect_heavy); break;
        case O_DWIDTH:    bad |= parse_dbl(optarg, &p->disorder_width); break;
        case O_ABSORB:    p->absorb = 1; break;
        case O_AWIDTH:    bad |= parse_int(optarg, &p->absorb_width); break;
        case O_AGAMMA:    bad |= parse_dbl(optarg, &p->absorb_gamma); break;
        case O_FPUT:      bad |= parse_dbl(optarg, &p->fput_strain); break;
        case O_NT:        bad |= parse_int(optarg, &p->nt); break;
        case O_SEED: { int s = 0; bad |= parse_int(optarg, &s); p->seed = (unsigned)s; seed_given = 1; break; }
        case O_WIDTH:     bad |= parse_int(optarg, &p->win_w); break;
        case O_HEIGHT:    bad |= parse_int(optarg, &p->win_h); break;
        case 'h':  usage(argv[0]); return 1;
        default:   usage(argv[0]); return -1;
        }
        if (bad) { fprintf(stderr, "bad value for option near '%s'\n", optarg ? optarg : ""); return -1; }
    }
    if (optind < argc) { fprintf(stderr, "unexpected argument '%s'\n", argv[optind]); return -1; }

    if (p->kdiag < 0) p->kdiag = p->k;
    if (p->absorb_gamma < 0) p->absorb_gamma = sqrt(p->k * p->mass);
    if (!seed_given) p->seed = (unsigned)time(NULL);

    const char *err = NULL;
    if (p->cols < 4 || p->rows < 4 || p->cols > 1024 || p->rows > 1024) err = "cols/rows must be in [4, 1024]";
    else if (p->bc == BC_PERIODIC && (p->cols % 2 || p->rows % 2)) err = "periodic BC needs even cols and rows";
    else if (p->mass <= 0 || p->k <= 0 || p->kdiag < 0) err = "mass and k must be positive";
    else if (p->damping < 0 || p->absorb_gamma < 0) err = "damping must be >= 0";
    else if (p->dt <= 0 || p->dt > 0.1) err = "dt must be in (0, 0.1]";
    else if (p->speed <= 0 || p->speed > 64) err = "speed must be in (0, 64]";
    else if (p->prestrain < 0 || p->prestrain >= 0.5) err = "prestrain must be in [0, 0.5)";
    else if (p->defect_light <= 0 || p->defect_heavy <= 0) err = "defect mass ratios must be positive";
    else if (p->disorder_width < 0 || p->disorder_width >= 1) err = "disorder width must be in [0, 1)";
    else if (p->absorb_width < 1) err = "absorb width must be >= 1";
    else if (p->nt < 16 || (p->nt & (p->nt - 1))) err = "nt must be a power of two >= 16";
    if (err) { fprintf(stderr, "%s\n", err); return -1; }
    return 0;
}
