/*
 * spring_grid.c — 2D NxM spring-mass lattice simulation
 *
 * Window layout:
 *   [ simulation | spectrum | help panel ]
 *
 * Physics options:
 *   - Cardinal + optional diagonal springs between neighbours
 *   - Optional fixed-edge boundary: edge masses pinned to equilibrium
 *     by springs of the same stiffness k
 *
 * Spectrum panel:
 *   Real-time 2D DFT of the transverse (y) velocity field.
 *   Displays |ṽ(qx,qy)|² on a log scale — the phonon spectral weight.
 *   For small displacements this traces the dispersion surface
 *       ω²(q) = (4k/m)[sin²(qx·a/2) + sin²(qy·a/2)]
 *   plus additional branches when diagonal springs are on.
 *
 * Build:
 *   gcc -O2 -std=c17 spring_grid.c -o spring_grid \
 *       $(sdl2-config --cflags --libs) -lfftw3 -lm
 *
 * Controls:  see the help panel, or the table in README.md
 */

#include <SDL.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ───────────────────────── tuneable parameters ─────────────────────────── */
#define COLS            64
#define ROWS            64
#define MASS            1.0
#define SPRING_K        400.0
#define DAMPING         0.4
#define DT              0.0005
#define STEPS_PER_FRAME 16
#define SPACING         16
#define MASS_RADIUS     2
#define WINDOW_TITLE    "Spring-Mass Grid"

/* panel widths */
#define SPEC_W          320
#define HELP_W          220
/* ─────────────────────────────────────────────────────────────────────────── */

#define NMASSES  (ROWS * COLS)
#define MAX_SPRINGS (NMASSES * 8)

typedef struct { double x,y, vx,vy, ax,ay, ox,oy; } Mass;
typedef struct { int a,b; double L0; } Spring;

static Mass   masses[NMASSES];
static Spring springs[MAX_SPRINGS];
static int    nsprings = 0;

static int    win_w, win_h, sim_w;
static int    diagonal    = 1;
static int    fixed_edges = 0;
static double damping     = DAMPING;

static int    drag_idx = -1;
static double drag_x, drag_y;

static fftw_complex *fft_in  = NULL;
static fftw_complex *fft_out = NULL;
static fftw_plan     fft_plan;

static double spec[ROWS][COLS];

/* ─────────────────────────────────────────────────────────────── helpers ── */

static inline int idx(int c, int r){ return r*COLS+c; }

static void build_springs(void)
{
    nsprings = 0;
    double diag_L0 = SPACING * sqrt(2.0);

    for (int r = 0; r < ROWS; r++){
        for (int c = 0; c < COLS; c++){
            int a = idx(c,r);
            if (c+1 < COLS)                       springs[nsprings++]=(Spring){a,idx(c+1,r),  SPACING};
            if (r+1 < ROWS)                       springs[nsprings++]=(Spring){a,idx(c,r+1),  SPACING};
            if (diagonal && c+1<COLS && r+1<ROWS) springs[nsprings++]=(Spring){a,idx(c+1,r+1),diag_L0};
            if (diagonal && c-1>=0  && r+1<ROWS)  springs[nsprings++]=(Spring){a,idx(c-1,r+1),diag_L0};
        }
    }

    if (fixed_edges){
        for (int c = 0; c < COLS; c++){
            springs[nsprings++]=(Spring){idx(c,0),      -(idx(c,0)+1),      0.0};
            springs[nsprings++]=(Spring){idx(c,ROWS-1), -(idx(c,ROWS-1)+1), 0.0};
        }
        for (int r = 1; r < ROWS-1; r++){
            springs[nsprings++]=(Spring){idx(0,r),      -(idx(0,r)+1),      0.0};
            springs[nsprings++]=(Spring){idx(COLS-1,r), -(idx(COLS-1,r)+1), 0.0};
        }
    }
}

static void reset(void)
{
    int mx = (sim_w  - (COLS-1)*SPACING) / 2;
    int my = (win_h  - (ROWS-1)*SPACING) / 2;
    for (int r=0;r<ROWS;r++) for (int c=0;c<COLS;c++){
        Mass *m = &masses[idx(c,r)];
        m->ox = mx + c*SPACING;
        m->oy = my + r*SPACING;
        m->x=m->ox; m->y=m->oy;
        m->vx=m->vy=m->ax=m->ay=0.0;
    }
    build_springs();
    memset(spec, 0, sizeof spec);
}

/* ─────────────────────────────────────────── force integration ──────── */

static void compute_forces(void)
{
    for (int i=0;i<NMASSES;i++) masses[i].ax=masses[i].ay=0.0;

    for (int s=0;s<nsprings;s++){
        int ai = springs[s].a;
        int bi = springs[s].b;

        if (bi < 0){
            int real_b = (-bi) - 1;
            double dx = masses[real_b].ox - masses[ai].x;
            double dy = masses[real_b].oy - masses[ai].y;
            double dist = sqrt(dx*dx+dy*dy);
            if (dist < 1e-9) continue;
            double f = SPRING_K * dist;
            masses[ai].ax += f*dx/(dist*MASS);
            masses[ai].ay += f*dy/(dist*MASS);
            continue;
        }

        double dx = masses[bi].x - masses[ai].x;
        double dy = masses[bi].y - masses[ai].y;
        double dist = sqrt(dx*dx+dy*dy);
        if (dist < 1e-9) continue;
        double f = SPRING_K * (dist - springs[s].L0) / dist;
        double fx = f*dx, fy = f*dy;
        masses[ai].ax += fx/MASS;
        masses[ai].ay += fy/MASS;
        masses[bi].ax -= fx/MASS;
        masses[bi].ay -= fy/MASS;
    }

    double dc = damping/MASS;
    for (int i=0;i<NMASSES;i++){
        if (i==drag_idx) continue;
        masses[i].ax -= dc*masses[i].vx;
        masses[i].ay -= dc*masses[i].vy;
    }
}

static void step(void)
{
    for (int i=0;i<NMASSES;i++){
        if (i==drag_idx) continue;
        masses[i].vx += 0.5*masses[i].ax*DT;
        masses[i].vy += 0.5*masses[i].ay*DT;
        masses[i].x  += masses[i].vx*DT;
        masses[i].y  += masses[i].vy*DT;
    }
    compute_forces();
    for (int i=0;i<NMASSES;i++){
        if (i==drag_idx) continue;
        masses[i].vx += 0.5*masses[i].ax*DT;
        masses[i].vy += 0.5*masses[i].ay*DT;
    }
    if (drag_idx>=0){
        masses[drag_idx].x=drag_x; masses[drag_idx].y=drag_y;
        masses[drag_idx].vx=masses[drag_idx].vy=0.0;
    }
}

/* ────────────────────────────────────────── spectrum ─────────── */

static void update_spectrum(void)
{
    for (int r=0;r<ROWS;r++) for (int c=0;c<COLS;c++){
        fft_in[r*COLS+c][0] = masses[idx(c,r)].vy;
        fft_in[r*COLS+c][1] = 0.0;
    }
    fftw_execute(fft_plan);

    double alpha = 0.15;
    for (int r=0;r<ROWS;r++) for (int c=0;c<COLS;c++){
        double re = fft_out[r*COLS+c][0];
        double im = fft_out[r*COLS+c][1];
        spec[r][c] = (1.0-alpha)*spec[r][c] + alpha*(re*re+im*im);
    }
}

/* ─────────────────────────────────────────────────── bitmap font ─────────── */

static const uint8_t FONT5x7[][5] = {
    {0x00,0x00,0x00,0x00,0x00},{0x00,0x00,0x5F,0x00,0x00},{0x00,0x07,0x00,0x07,0x00},
    {0x14,0x7F,0x14,0x7F,0x14},{0x24,0x2A,0x7F,0x2A,0x12},{0x23,0x13,0x08,0x64,0x62},
    {0x36,0x49,0x55,0x22,0x50},{0x00,0x05,0x03,0x00,0x00},{0x00,0x1C,0x22,0x41,0x00},
    {0x00,0x41,0x22,0x1C,0x00},{0x08,0x2A,0x1C,0x2A,0x08},{0x08,0x08,0x3E,0x08,0x08},
    {0x00,0x50,0x30,0x00,0x00},{0x08,0x08,0x08,0x08,0x08},{0x00,0x60,0x60,0x00,0x00},
    {0x20,0x10,0x08,0x04,0x02},{0x3E,0x51,0x49,0x45,0x3E},{0x00,0x42,0x7F,0x40,0x00},
    {0x42,0x61,0x51,0x49,0x46},{0x21,0x41,0x45,0x4B,0x31},{0x18,0x14,0x12,0x7F,0x10},
    {0x27,0x45,0x45,0x45,0x39},{0x3C,0x4A,0x49,0x49,0x30},{0x01,0x71,0x09,0x05,0x03},
    {0x36,0x49,0x49,0x49,0x36},{0x06,0x49,0x49,0x29,0x1E},{0x00,0x36,0x36,0x00,0x00},
    {0x00,0x56,0x36,0x00,0x00},{0x00,0x08,0x14,0x22,0x41},{0x14,0x14,0x14,0x14,0x14},
    {0x41,0x22,0x14,0x08,0x00},{0x02,0x01,0x51,0x09,0x06},{0x32,0x49,0x79,0x41,0x3E},
    {0x7E,0x11,0x11,0x11,0x7E},{0x7F,0x49,0x49,0x49,0x36},{0x3E,0x41,0x41,0x41,0x22},
    {0x7F,0x41,0x41,0x22,0x1C},{0x7F,0x49,0x49,0x49,0x41},{0x7F,0x09,0x09,0x01,0x01},
    {0x3E,0x41,0x41,0x51,0x32},{0x7F,0x08,0x08,0x08,0x7F},{0x00,0x41,0x7F,0x41,0x00},
    {0x20,0x40,0x41,0x3F,0x01},{0x7F,0x08,0x14,0x22,0x41},{0x7F,0x40,0x40,0x40,0x40},
    {0x7F,0x02,0x04,0x02,0x7F},{0x7F,0x04,0x08,0x10,0x7F},{0x3E,0x41,0x41,0x41,0x3E},
    {0x7F,0x09,0x09,0x09,0x06},{0x3E,0x41,0x51,0x21,0x5E},{0x7F,0x09,0x19,0x29,0x46},
    {0x46,0x49,0x49,0x49,0x31},{0x01,0x01,0x7F,0x01,0x01},{0x3F,0x40,0x40,0x40,0x3F},
    {0x1F,0x20,0x40,0x20,0x1F},{0x3F,0x40,0x38,0x40,0x3F},{0x63,0x14,0x08,0x14,0x63},
    {0x03,0x04,0x78,0x04,0x03},{0x61,0x51,0x49,0x45,0x43},{0x00,0x00,0x7F,0x41,0x41},
    {0x02,0x04,0x08,0x10,0x20},{0x41,0x41,0x7F,0x00,0x00},{0x04,0x02,0x01,0x02,0x04},
    {0x40,0x40,0x40,0x40,0x40},{0x00,0x01,0x02,0x04,0x00},{0x20,0x54,0x54,0x54,0x78},
    {0x7F,0x48,0x44,0x44,0x38},{0x38,0x44,0x44,0x44,0x20},{0x38,0x44,0x44,0x48,0x7F},
    {0x38,0x54,0x54,0x54,0x18},{0x08,0x7E,0x09,0x01,0x02},{0x08,0x14,0x54,0x54,0x3C},
    {0x7F,0x08,0x04,0x04,0x78},{0x00,0x44,0x7D,0x40,0x00},{0x20,0x40,0x44,0x3D,0x00},
    {0x00,0x7F,0x10,0x28,0x44},{0x00,0x41,0x7F,0x40,0x00},{0x7C,0x04,0x18,0x04,0x78},
    {0x7C,0x08,0x04,0x04,0x78},{0x38,0x44,0x44,0x44,0x38},{0x7C,0x14,0x14,0x14,0x08},
    {0x08,0x14,0x14,0x18,0x7C},{0x7C,0x08,0x04,0x04,0x08},{0x48,0x54,0x54,0x54,0x20},
    {0x04,0x3F,0x44,0x40,0x20},{0x3C,0x40,0x40,0x20,0x7C},{0x1C,0x20,0x40,0x20,0x1C},
    {0x3C,0x40,0x30,0x40,0x3C},{0x44,0x28,0x10,0x28,0x44},{0x0C,0x50,0x50,0x50,0x3C},
    {0x44,0x64,0x54,0x4C,0x44},{0x00,0x08,0x36,0x41,0x00},{0x00,0x00,0x7F,0x00,0x00},
    {0x00,0x41,0x36,0x08,0x00},{0x08,0x08,0x2A,0x1C,0x08},
};

static void draw_char(SDL_Renderer *ren, int x, int y, char ch, int scale)
{
    if (ch < 32 || ch > 126) ch = '?';
    const uint8_t *glyph = FONT5x7[(unsigned char)ch - 32];
    for (int col=0; col<5; col++){
        uint8_t bits = glyph[col];
        for (int row=0; row<7; row++){
            if (bits & (1<<row)){
                SDL_Rect r={x+col*scale, y+row*scale, scale, scale};
                SDL_RenderFillRect(ren,&r);
            }
        }
    }
}

static void draw_string(SDL_Renderer *ren, int x, int y, const char *s, int scale)
{
    for(;*s;s++,x+=6*scale) draw_char(ren,x,y,*s,scale);
}

/* ─────────────────────────────────────────── rendering helpers ───────── */

static void draw_circle_outline(SDL_Renderer *ren, int cx, int cy, int r)
{
    int x=r,y=0,err=0;
    while(x>=y){
        SDL_RenderDrawPoint(ren,cx+x,cy+y); SDL_RenderDrawPoint(ren,cx+y,cy+x);
        SDL_RenderDrawPoint(ren,cx-y,cy+x); SDL_RenderDrawPoint(ren,cx-x,cy+y);
        SDL_RenderDrawPoint(ren,cx-x,cy-y); SDL_RenderDrawPoint(ren,cx-y,cy-x);
        SDL_RenderDrawPoint(ren,cx+y,cy-x); SDL_RenderDrawPoint(ren,cx+x,cy-y);
        if(err<=0){y++;err+=2*y+1;} if(err>0){x--;err-=2*x+1;}
    }
}

static void fill_circle(SDL_Renderer *ren, int cx, int cy, int r)
{
    for(int dy=-r;dy<=r;dy++){
        int dx=(int)sqrt((double)(r*r-dy*dy));
        SDL_RenderDrawLine(ren,cx-dx,cy+dy,cx+dx,cy+dy);
    }
}

static void speed_colour(SDL_Renderer *ren, int i)
{
    double spd=sqrt(masses[i].vx*masses[i].vx+masses[i].vy*masses[i].vy);
    double t=spd/600.0; if(t>1.0)t=1.0;
    SDL_SetRenderDrawColor(ren,(Uint8)(50+205*t),(Uint8)(130-80*t),(Uint8)(220-200*t),255);
}

static void strain_colour(SDL_Renderer *ren, int s)
{
    if(springs[s].b<0){SDL_SetRenderDrawColor(ren,80,80,200,120);return;}
    Mass *a=&masses[springs[s].a],*b=&masses[springs[s].b];
    double dx=b->x-a->x,dy=b->y-a->y,dist=sqrt(dx*dx+dy*dy);
    double strain=(dist-springs[s].L0)/(springs[s].L0+1e-9);
    double t=fabs(strain)*8.0; if(t>1.0)t=1.0;
    if(strain>0){Uint8 v=(Uint8)(255-200*t);SDL_SetRenderDrawColor(ren,255,v,v,80+(Uint8)(100*t));}
    else        {Uint8 v=(Uint8)(255-200*t);SDL_SetRenderDrawColor(ren,v,v,255,80+(Uint8)(100*t));}
}

/* ─────────────────────────────────────────── spectrum panel ───────── */

static void render_spectrum(SDL_Renderer *ren, int ox, int oy, int pw, int ph)
{
    SDL_SetRenderDrawColor(ren,8,8,14,255);
    SDL_Rect bg={ox,oy,pw,ph}; SDL_RenderFillRect(ren,&bg);

    double mx=1e-30;
    for(int r=0;r<ROWS;r++) for(int c=0;c<COLS;c++) if(spec[r][c]>mx) mx=spec[r][c];

    int cw=pw/COLS, ch=ph/ROWS;
    if(cw<1)cw=1;
    if(ch<1)ch=1;

    for(int r=0;r<ROWS;r++) for(int c=0;c<COLS;c++){
        int sc=(c+COLS/2)%COLS, sr=(r+ROWS/2)%ROWS;
        double v=spec[sr][sc];
        double t=(v>1e-30)?log(v/mx)/10.0+1.0:0.0;
        if(t<0)t=0;
        Uint8 R,G,B;
        if(t<0.5){double s=t*2.0;R=(Uint8)(120*s);G=0;B=(Uint8)(180*s);}
        else     {double s=(t-0.5)*2.0;R=(Uint8)(120+135*s);G=(Uint8)(220*s);B=(Uint8)(180+75*s);}
        SDL_SetRenderDrawColor(ren,R,G,B,255);
        SDL_Rect cell={ox+c*cw,oy+r*ch,cw,ch};
        SDL_RenderFillRect(ren,&cell);
    }

    /* axis label */
    SDL_SetRenderDrawColor(ren,130,130,150,255);
    draw_string(ren,ox+pw/2-12,oy+ph+4,"qx",1);
    draw_string(ren,ox+2,oy+2,"qy",1);

    SDL_SetRenderDrawColor(ren,60,60,80,255);
    SDL_RenderDrawRect(ren,&bg);

    SDL_SetRenderDrawColor(ren,180,180,200,255);
    draw_string(ren,ox+4,oy-13,"|v(q)|^2  log scale",1);
}

/* ─────────────────────────────────────────── help panel ───────── */

static void render_help(SDL_Renderer *ren, int ox, int oy, int pw, int ph)
{
    SDL_SetRenderDrawColor(ren,12,12,18,255);
    SDL_Rect bg={ox,oy,pw,ph}; SDL_RenderFillRect(ren,&bg);
    SDL_SetRenderDrawColor(ren,40,40,60,255);
    SDL_RenderDrawRect(ren,&bg);

    const int xs=ox+8, lh=11;
    int y=oy+8;

    SDL_SetRenderDrawColor(ren,200,200,255,255);
    draw_string(ren,xs,y,"CONTROLS",1); y+=lh+4;
    SDL_SetRenderDrawColor(ren,50,50,70,255);
    SDL_RenderDrawLine(ren,ox+4,y,ox+pw-4,y); y+=6;

    static const struct { const char *key, *desc; } help[] = {
        {"L-drag", "pull a mass"},
        {"R-click","random impulse"},
        {"Space",  "pluck centre"},
        {"D",      "diag springs"},
        {"E",      "fix edges"},
        {"G",      "chaos impulse"},
        {"R",      "reset"},
        {"+/-",    "damping"},
        {"Q/Esc",  "quit"},
    };
    for(int i=0;i<(int)(sizeof help/sizeof*help);i++){
        SDL_SetRenderDrawColor(ren,100,210,120,255);
        draw_string(ren,xs,y,help[i].key,1);
        SDL_SetRenderDrawColor(ren,150,150,170,255);
        draw_string(ren,xs+40,y,help[i].desc,1);
        y+=lh;
    }

    y+=6;
    SDL_SetRenderDrawColor(ren,50,50,70,255);
    SDL_RenderDrawLine(ren,ox+4,y,ox+pw-4,y); y+=8;

    SDL_SetRenderDrawColor(ren,200,200,255,255);
    draw_string(ren,xs,y,"STATE",1); y+=lh+2;

    char buf[40];
    SDL_SetRenderDrawColor(ren,diagonal?80:180, diagonal?200:80, 100,255);
    snprintf(buf,sizeof buf,"Diagonal: %s", diagonal?"ON ":"OFF");
    draw_string(ren,xs,y,buf,1); y+=lh;

    SDL_SetRenderDrawColor(ren,fixed_edges?80:180, fixed_edges?200:80, 100,255);
    snprintf(buf,sizeof buf,"FixEdges: %s", fixed_edges?"ON ":"OFF");
    draw_string(ren,xs,y,buf,1); y+=lh;

    SDL_SetRenderDrawColor(ren,150,150,170,255);
    snprintf(buf,sizeof buf,"Damping:  %.2f", damping);
    draw_string(ren,xs,y,buf,1); y+=lh;
    snprintf(buf,sizeof buf,"Grid: %dx%d", COLS, ROWS);
    draw_string(ren,xs,y,buf,1); y+=lh;
    snprintf(buf,sizeof buf,"k=%.0f  m=%.1f", SPRING_K, MASS);
    draw_string(ren,xs,y,buf,1);
}

/* ─────────────────────────────────────────── main render ──────── */

static void render(SDL_Renderer *ren)
{
    SDL_SetRenderDrawColor(ren,15,15,20,255);
    SDL_RenderClear(ren);

    SDL_Rect sim_clip={0,0,sim_w,win_h};
    SDL_RenderSetClipRect(ren,&sim_clip);

    for(int s=0;s<nsprings;s++){
        if(springs[s].b<0) continue;
        strain_colour(ren,s);
        SDL_RenderDrawLine(ren,
            (int)masses[springs[s].a].x,(int)masses[springs[s].a].y,
            (int)masses[springs[s].b].x,(int)masses[springs[s].b].y);
    }

    if(fixed_edges){
        for(int i=0;i<NMASSES;i++){
            int r=i/COLS,c=i%COLS;
            if(r==0||r==ROWS-1||c==0||c==COLS-1){
                SDL_SetRenderDrawColor(ren,50,50,140,160);
                SDL_RenderDrawLine(ren,(int)masses[i].x,(int)masses[i].y,
                                       (int)masses[i].ox,(int)masses[i].oy);
                SDL_SetRenderDrawColor(ren,30,30,100,200);
                fill_circle(ren,(int)masses[i].ox,(int)masses[i].oy,3);
            }
        }
    }

    for(int i=0;i<NMASSES;i++){
        int cx=(int)masses[i].x,cy=(int)masses[i].y;
        speed_colour(ren,i); fill_circle(ren,cx,cy,MASS_RADIUS);
        if(i==drag_idx) SDL_SetRenderDrawColor(ren,255,255,100,255);
        else            SDL_SetRenderDrawColor(ren,200,200,220,200);
        draw_circle_outline(ren,cx,cy,MASS_RADIUS);
    }

    SDL_RenderSetClipRect(ren,NULL);

    SDL_SetRenderDrawColor(ren,40,40,60,255);
    SDL_RenderDrawLine(ren,sim_w,0,sim_w,win_h);

    render_spectrum(ren, sim_w+4, 20, SPEC_W-8, win_h-40);

    SDL_RenderDrawLine(ren,sim_w+SPEC_W,0,sim_w+SPEC_W,win_h);

    render_help(ren, sim_w+SPEC_W+4, 4, HELP_W-8, win_h-8);

    SDL_RenderPresent(ren);
}

/* ─────────────────────────────────────────────── nearest mass ──── */

static int nearest_mass(int mx, int my, double max_dist)
{
    int best=-1; double best2=max_dist*max_dist;
    for(int i=0;i<NMASSES;i++){
        double dx=masses[i].x-mx,dy=masses[i].y-my,d2=dx*dx+dy*dy;
        if(d2<best2){best2=d2;best=i;}
    }
    return best;
}

/* ─────────────────────────────────────────────────── main ──────── */

int main(void)
{
    srand((unsigned)time(NULL));

    if(SDL_Init(SDL_INIT_VIDEO)){
        fprintf(stderr,"SDL_Init: %s\n",SDL_GetError()); return 1;
    }

    sim_w = (COLS+1)*SPACING;
    win_w = sim_w + SPEC_W + HELP_W;
    win_h = (ROWS+1)*SPACING;

    SDL_Window *win = SDL_CreateWindow(WINDOW_TITLE,
        SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,
        win_w,win_h,SDL_WINDOW_SHOWN);
    if(!win){fprintf(stderr,"SDL_CreateWindow: %s\n",SDL_GetError());return 1;}

    SDL_Renderer *ren = SDL_CreateRenderer(win,-1,
        SDL_RENDERER_ACCELERATED|SDL_RENDERER_PRESENTVSYNC);
    if(!ren){fprintf(stderr,"SDL_CreateRenderer: %s\n",SDL_GetError());return 1;}

    SDL_SetRenderDrawBlendMode(ren,SDL_BLENDMODE_BLEND);

    fft_in  = fftw_alloc_complex(NMASSES);
    fft_out = fftw_alloc_complex(NMASSES);
    fft_plan = fftw_plan_dft_2d(ROWS,COLS,fft_in,fft_out,FFTW_FORWARD,FFTW_MEASURE);

    reset();
    masses[idx(COLS/2,ROWS/2)].y -= SPACING*2.5;

    int running=1; SDL_Event ev;
    while(running){
        while(SDL_PollEvent(&ev)){
            switch(ev.type){
            case SDL_QUIT: running=0; break;
            case SDL_KEYDOWN:
                switch(ev.key.keysym.sym){
                case SDLK_ESCAPE: case SDLK_q: running=0; break;
                case SDLK_r:  reset(); drag_idx=-1; break;
                case SDLK_SPACE:{
                    int ci=idx(COLS/2,ROWS/2);
                    masses[ci].y-=SPACING*2.0; masses[ci].vx=masses[ci].vy=0.0;
                    break;}
                case SDLK_d: diagonal=!diagonal;    build_springs(); break;
                case SDLK_e: fixed_edges=!fixed_edges; build_springs(); break;
                case SDLK_g:
                    for(int i=0;i<NMASSES;i++){
                        masses[i].vx+=((rand()%400)-200);
                        masses[i].vy+=((rand()%400)-200);
                    }
                    break;
                case SDLK_PLUS: case SDLK_EQUALS: damping*=1.25; break;
                case SDLK_MINUS: damping/=1.25; if(damping<0)damping=0; break;
                }
                break;
            case SDL_MOUSEBUTTONDOWN:
                if(ev.button.x < sim_w){
                    if(ev.button.button==SDL_BUTTON_LEFT){
                        drag_idx=nearest_mass(ev.button.x,ev.button.y,SPACING*0.9);
                        if(drag_idx>=0){drag_x=ev.button.x;drag_y=ev.button.y;}
                    } else if(ev.button.button==SDL_BUTTON_RIGHT){
                        for(int i=0;i<NMASSES;i++){
                            masses[i].vy+=((rand()%600)-300);
                            masses[i].vx+=((rand()%200)-100);
                        }
                    }
                }
                break;
            case SDL_MOUSEMOTION:
                if(drag_idx>=0){drag_x=ev.motion.x;drag_y=ev.motion.y;}
                break;
            case SDL_MOUSEBUTTONUP:
                if(ev.button.button==SDL_BUTTON_LEFT && drag_idx>=0){
                    masses[drag_idx].vx=ev.button.x-drag_x;
                    masses[drag_idx].vy=ev.button.y-drag_y;
                    drag_idx=-1;
                }
                break;
            }
        }

        for(int s=0;s<STEPS_PER_FRAME;s++) step();
        update_spectrum();
        render(ren);
    }

    fftw_destroy_plan(fft_plan);
    fftw_free(fft_in); fftw_free(fft_out);
    SDL_DestroyRenderer(ren); SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
