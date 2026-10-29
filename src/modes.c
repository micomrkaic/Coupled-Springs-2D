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

static int rod_active(double qx, double qy, double rx, double ry)
{
    return 1.0 - cos(qx * rx + qy * ry) > 1e-9;
}

void branches(const Lattice *L, double qx, double qy, double w[2], double e[2][2])
{
    double D[3];
    dynmat(L, qx, qy, D);

    /* rod constraints: collect active normals */
    double nx = 0, ny = 0;
    int nact = 0, spans = 0;
    const double s2 = sqrt(0.5);
    double bonds[4][4] = { { 1, 0, 1, 0 }, { 0, 1, 0, 1 }, { 1, 1, s2, s2 }, { 1, -1, s2, -s2 } };
    for (int k = 0; k < 4; k++) {
        if (k < 2 ? !L->rod_card : !(L->rod_diag && L->diag)) continue;
        if (!rod_active(qx, qy, bonds[k][0], bonds[k][1])) continue;
        if (nact == 0) { nx = bonds[k][2]; ny = bonds[k][3]; }
        else if (fabs(nx * bonds[k][3] - ny * bonds[k][2]) > 1e-9) spans = 1;
        nact++;
    }
    if (nact) {
        double qn = sqrt(qx * qx + qy * qy);
        double qhx = qn > 1e-12 ? qx / qn : 1.0, qhy = qn > 1e-12 ? qy / qn : 0.0;
        double tx = -ny, ty = nx;                       /* allowed direction */
        double lx, ly;
        int t_is_L = fabs(tx * qhx + ty * qhy) >= fabs(nx * qhx + ny * qhy);
        if (t_is_L) { lx = tx; ly = ty; } else { lx = nx; ly = ny; }
        if (lx * qhx + ly * qhy < 0) { lx = -lx; ly = -ly; }
        e[0][0] = lx;  e[0][1] = ly;
        e[1][0] = -ly; e[1][1] = lx;
        double wt = INFINITY;
        if (!spans) {
            double l = tx * tx * D[0] + ty * ty * D[1] + 2 * tx * ty * D[2];
            wt = sqrt(l > 0 ? l : 0);
        }
        w[0] = t_is_L ? wt : INFINITY;
        w[1] = t_is_L ? INFINITY : wt;
        return;
    }
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
            if (isfinite(w[0]) && w[0] > top) top = w[0];
            if (isfinite(w[1]) && w[1] > top) top = w[1];
        }
    return top;
}
