#include "modes.h"

#include <math.h>

int line_w0(Boundary bc) { return bc == BC_CLAMPED ? 1 : 0; }

int line_nw(Boundary bc, int N)
{
    switch (bc) {
    case BC_CLAMPED:  return N - 1;          /* w = 1 .. N-2 */
    case BC_PERIODIC: return N / 2 + 1;      /* w = 0 .. N/2 */
    default:          return N;              /* w = 0 .. N-1 */
    }
}

double line_q(Boundary bc, int N, int w)
{
    switch (bc) {
    case BC_CLAMPED:  return PI * w / (N - 1);
    case BC_PERIODIC: return 2.0 * PI * w / N;
    default:          return PI * w / N;
    }
}

int line_has_sin(Boundary bc, int N, int w)
{
    return bc == BC_PERIODIC && w > 0 && 2 * w < N;
}

double line_phi(Boundary bc, int N, int w, int s, int j)
{
    switch (bc) {
    case BC_CLAMPED:
        if (j <= 0 || j >= N - 1) return 0.0;
        return sqrt(2.0 / (N - 1)) * sin(PI * w * j / (N - 1));
    case BC_PERIODIC: {
        double th = 2.0 * PI * w * j / N;
        if (w == 0 || 2 * w == N) return cos(th) / sqrt((double)N);
        return sqrt(2.0 / N) * (s ? sin(th) : cos(th));
    }
    default:
        if (w == 0) return 1.0 / sqrt((double)N);
        return sqrt(2.0 / N) * cos(PI * w * (j + 0.5) / N);
    }
}

void dynmat(const Lattice *L, double qx, double qy, double D[3])
{
    double xx = 0, yy = 0, xy = 0, c;

    c = 1.0 - cos(qx);                 /* bond (1,0): n = x-hat */
    xx += L->k * c;  yy += L->tc * c;
    c = 1.0 - cos(qy);                 /* bond (0,1): n = y-hat */
    yy += L->k * c;  xx += L->tc * c;

    if (L->diag) {
        double s = 0.5 * (L->kd + L->td), d = 0.5 * (L->kd - L->td);
        c = 1.0 - cos(qx + qy);        /* bond (1,1)  */
        xx += s * c;  yy += s * c;  xy += d * c;
        c = 1.0 - cos(qx - qy);        /* bond (1,-1) */
        xx += s * c;  yy += s * c;  xy -= d * c;
    }
    double f = 2.0 / L->m;
    D[0] = f * xx;  D[1] = f * yy;  D[2] = f * xy;
}

void branches(const Lattice *L, double qx, double qy, double w[2], double e[2][2])
{
    double D[3];
    dynmat(L, qx, qy, D);
    double a = D[0], b = D[1], c = D[2];
    double mean = 0.5 * (a + b), half = 0.5 * (a - b);
    double rad = sqrt(half * half + c * c);
    double l1 = mean + rad, l2 = mean - rad;

    double qn = sqrt(qx * qx + qy * qy);
    double qhx = qn > 1e-12 ? qx / qn : 1.0, qhy = qn > 1e-12 ? qy / qn : 0.0;

    double e1x, e1y;
    if (rad <= 1e-12 * (fabs(a) + fabs(b) + 1e-30)) {   /* degenerate */
        e1x = qhx; e1y = qhy;
    } else {
        double th = 0.5 * atan2(2.0 * c, a - b);
        e1x = cos(th); e1y = sin(th);
    }
    double e2x = -e1y, e2y = e1x;

    /* label: quasi-L has the larger overlap with q-hat */
    double o1 = fabs(e1x * qhx + e1y * qhy), o2 = fabs(e2x * qhx + e2y * qhy);
    double lx, ly, lL, lT;
    if (o1 >= o2) { lx = e1x; ly = e1y; lL = l1; lT = l2; }
    else          { lx = e2x; ly = e2y; lL = l2; lT = l1; }
    if (lx * qhx + ly * qhy < 0) { lx = -lx; ly = -ly; }

    e[0][0] = lx;  e[0][1] = ly;
    e[1][0] = -ly; e[1][1] = lx;
    w[0] = sqrt(lL > 0 ? lL : 0);
    w[1] = sqrt(lT > 0 ? lT : 0);
}

double band_top(const Lattice *L)
{
    double top = 0;
    const int n = 64;
    for (int i = 0; i <= n; i++)
        for (int j = -n; j <= n; j++) {
            double w[2], e[2][2];
            branches(L, PI * i / n, PI * j / n, w, e);
            if (w[0] > top) top = w[0];
            if (w[1] > top) top = w[1];
        }
    return top;
}
