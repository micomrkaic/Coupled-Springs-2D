/*
 * params.h — run-time parameters and command-line parsing
 *
 * Units: lattice units throughout.  The equilibrium spacing a = 1, mass is
 * in kg, time in s, spring constants in N/m with lengths measured in units
 * of a.  Pixels appear only in the renderer.
 */
#ifndef PARAMS_H
#define PARAMS_H

#define PI 3.14159265358979323846

typedef enum { BC_FREE, BC_PINNED, BC_CLAMPED, BC_PERIODIC, BC_COUNT } Boundary;
typedef enum { DIS_NONE, DIS_LIGHT, DIS_HEAVY, DIS_RANDOM, DIS_COUNT } Disorder;

typedef struct {
    int      cols, rows;
    double   mass;            /* nominal node mass                          */
    double   k;               /* cardinal spring constant                   */
    double   kdiag;           /* diagonal spring constant                   */
    int      diagonal;        /* diagonal springs on/off                    */
    int      rod_card;        /* cardinal bonds are rigid rods              */
    int      rod_diag;        /* diagonal bonds are rigid rods              */
    double   damping;         /* uniform linear drag gamma (kg/s)           */
    double   dt;              /* integration step (s)                       */
    double   speed;           /* simulated seconds per wall-clock second    */
    Boundary bc;
    double   prestrain;       /* L0 = (1 - prestrain)*|r|; clamped/periodic */
    Disorder disorder;
    double   defect_light;    /* mass ratio of the light defect             */
    double   defect_heavy;    /* mass ratio of the heavy defect             */
    double   disorder_width;  /* random masses: m*(1 + w*U[-1,1])           */
    int      absorb;          /* graded absorbing layer on/off              */
    int      absorb_width;    /* layer thickness in nodes                   */
    double   absorb_gamma;    /* peak extra drag in the layer (kg/s)        */
    double   fput_strain;     /* peak strain of the F-key mode excitation   */
    int      nt;              /* time samples for the omega-q spectrum      */
    unsigned seed;
    int      win_w, win_h;    /* initial window size in points (0 = auto)   */
} Params;

void        params_default(Params *p);
/* 0 = run, 1 = exit successfully (--help), -1 = error */
int         params_parse(Params *p, int argc, char **argv);
const char *bc_name(Boundary bc);
const char *disorder_name(Disorder d);

#endif
