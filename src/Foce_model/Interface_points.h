#ifndef INTERFACE_POINTS_H
#define INTERFACE_POINTS_H

#include <stdlib.h>

typedef struct {
    coord x;// interface centroid
    coord n;// interface normal
    double Delta;// local grid size
    int tag;// droplet ID
    double alpha;// PLIC plane constant
    int level; // Grid level
} InterfacePoint;

#define MAX_INTERFACE_POINTS 5000

static InterfacePoint interface_points[MAX_INTERFACE_POINTS];
static int n_interface_points = 0;

//construct interface points, store and return the number of interface points
static int build_interface_points (scalar f, scalar droplet)
{
    n_interface_points = 0;

    foreach(serial) {
        if (f[] > 1e-6 && f[] < 1. - 1e-6) {

            if (n_interface_points >= MAX_INTERFACE_POINTS) {
                fprintf(stderr, "Warning: MAX_INTERFACE_POINTS exceeded\n");
                continue;
            }

            InterfacePoint *ip = &interface_points[n_interface_points++];

            coord n = interface_normal(point,f);
            double alpha = plane_alpha(f[],n);

            coord p;
            plane_area_center(n,alpha,&p);

            ip->x.x = x + Delta*p.x;
            ip->x.y = y + Delta*p.y;
            ip->x.z = 0.;
#if dimension == 3
            ip->x.z = z + Delta*p.z;
#endif
            ip->n = n;
            ip->Delta = Delta;
            ip->tag = (int) droplet[];
            ip->alpha = alpha;
            ip->level = level;
        }
    }

    return n_interface_points;
}

// Neighbhour search

typedef struct {
    int i, j;// indices for interface_points[]
    double d;// centroid-to-centroid distance
} CandidatePair;

#define MAX_CANDIDATE_PAIRS 20000
static CandidatePair candidate_pairs[MAX_CANDIDATE_PAIRS];
static int n_candidate_pairs = 0;

static void find_neighbour_pairs (double cutoff)
{
    n_candidate_pairs = 0;

    for (int a = 0; a < n_interface_points; a++) {
        for (int b = a + 1; b < n_interface_points; b++) {

            if (interface_points[a].tag == interface_points[b].tag)
                continue;// same droplet condition

            double dx = interface_points[a].x.x - interface_points[b].x.x;
            double dy = interface_points[a].x.y - interface_points[b].x.y;
            double d = sqrt(dx*dx + dy*dy);

            if (d <= cutoff && n_candidate_pairs < MAX_CANDIDATE_PAIRS) {
                candidate_pairs[n_candidate_pairs].i = a;
                candidate_pairs[n_candidate_pairs].j = b;
                candidate_pairs[n_candidate_pairs].d = d;
                n_candidate_pairs++;
            }
        }
    }
}

// reduce to min gap PER droplet pair

#define MAX_DROPLET_PAIRS 100

typedef struct {
    int tagA, tagB;
    double h_min;
    int point_i, point_j;
} DropletPairGap;

static DropletPairGap droplet_pair_gaps[MAX_DROPLET_PAIRS];
static int n_droplet_pairs = 0;

static int find_or_create_pair_gap (int tagA, int tagB)
{
    if (tagA > tagB) { int tmp = tagA; tagA = tagB; tagB = tmp; }

    for (int k = 0; k < n_droplet_pairs; k++)
        if (droplet_pair_gaps[k].tagA == tagA && droplet_pair_gaps[k].tagB == tagB)
            return k;

    if (n_droplet_pairs >= MAX_DROPLET_PAIRS) {
        fprintf(stderr, "Warning: MAX_DROPLET_PAIRS exceeded\n");
        return -1;
    }

    int k = n_droplet_pairs++;
    droplet_pair_gaps[k].tagA = tagA;
    droplet_pair_gaps[k].tagB = tagB;
    droplet_pair_gaps[k].h_min = 1e30;
    droplet_pair_gaps[k].point_i = -1;
    droplet_pair_gaps[k].point_j = -1;
    return k;
}

static DropletPairGap prev_pair_gaps[MAX_DROPLET_PAIRS];
static int n_prev_pairs = 0;

// remembers last step's gaps before compute_minimum_gaps() rebuilds them
static void snapshot_previous_gaps()
{
    for (int i = 0; i < n_droplet_pairs; i++)
        prev_pair_gaps[i] = droplet_pair_gaps[i];
    n_prev_pairs = n_droplet_pairs;
}

// returns last step's h_min for this tag pair, or -1 if not seen last step
static double lookup_prev_h_min (int tagA, int tagB)
{
    if (tagA > tagB) { int tmp = tagA; tagA = tagB; tagB = tmp; }

    for (int i = 0; i < n_prev_pairs; i++)
        if (prev_pair_gaps[i].tagA == tagA && prev_pair_gaps[i].tagB == tagB)
            return prev_pair_gaps[i].h_min;

    return -1.;
}

static void compute_minimum_gaps()
{
    snapshot_previous_gaps();

    n_droplet_pairs = 0;

    for (int c = 0; c < n_candidate_pairs; c++) {
        int a = candidate_pairs[c].i;
        int b = candidate_pairs[c].j;
        double d = candidate_pairs[c].d;

        int k = find_or_create_pair_gap(interface_points[a].tag,
                                         interface_points[b].tag);
        if (k < 0) continue;

        if (d < droplet_pair_gaps[k].h_min) {
            droplet_pair_gaps[k].h_min = d;
            droplet_pair_gaps[k].point_i = a;
            droplet_pair_gaps[k].point_j = b;
        }
    }
}

#endif

// Force Model Constants

double lub_epsilon = 10.; // proportionality constant for lubrication force
double lub_h0; 
double lub_hcutoff; // repulsive radius, set from grid spacing

double lub_attract_epsilon = 1700.; // strength of attractive force
double lub_damping = 1050.; // damping coefficient
double lub_attract_cutoff; // attraction reach, set from grid spacing
double lub_damping_cutoff; 


double lub_h0_cells = 1.5; 
double lub_hcutoff_cells   = 5.0; 
double lub_attract_cells   = 15.0; 
double lub_damping_cells   = 25.0; 


static void set_lubrication_lengthscales (double Delta_min)
{
    lub_h0 = lub_h0_cells      * Delta_min;
    lub_hcutoff = lub_hcutoff_cells * Delta_min;
    lub_attract_cutoff = lub_attract_cells * Delta_min;
    lub_damping_cutoff = lub_damping_cells * Delta_min;
}

static double lubrication_force_magnitude (double h)
{
    if (h >= lub_hcutoff || h <= 0.)
        return 0.;

    double hr  = h + lub_h0;
    double hcr = lub_hcutoff + lub_h0;

    return lub_epsilon * (1./(hr*hr*hr) - 1./(hcr*hcr*hcr));
}

static double attractive_force_magnitude (double h)
{
    if (h < lub_hcutoff || h >= lub_attract_cutoff)
        return 0.;

    double s = (h - lub_hcutoff) / (lub_attract_cutoff - lub_hcutoff); // 0..1 across the band
    double shape = sin(M_PI * s); // 0 at s=0 and s=1, peak of 1 at s=0.5

    return -lub_attract_epsilon * shape;
}

// total force: repulsive core + attractive
static double total_force_magnitude (double h)
{
    return lubrication_force_magnitude(h) + attractive_force_magnitude(h);
}
//storage of droplet pair 
typedef struct {
    coord xa, xb;// interface point positions
    coord dir;// unit vector from B to A
    double Fmag;// force magnitude
} LubForceEntry;

#define MAX_LUB_FORCES 100
static LubForceEntry lub_forces[MAX_LUB_FORCES];
static int n_lub_forces = 0;

//force list
static void prepare_lubrication_forces()
{
    n_lub_forces = 0;

    for (int k = 0; k < n_droplet_pairs; k++) {

        double h = droplet_pair_gaps[k].h_min;
        if (h <= 0. || h >= lub_damping_cutoff) continue;
        if (n_lub_forces >= MAX_LUB_FORCES) break;

        int ia = droplet_pair_gaps[k].point_i;
        int ib = droplet_pair_gaps[k].point_j;

        coord xa = interface_points[ia].x;
        coord xb = interface_points[ib].x;

        double dx = xa.x - xb.x;
        double dy = xa.y - xb.y;
        double dist = sqrt(dx*dx + dy*dy);
        if (dist < 1e-12) continue;

        // Damping 
        double h_prev = lookup_prev_h_min(droplet_pair_gaps[k].tagA, droplet_pair_gaps[k].tagB);
        double Fdamp = 0.;
        if (h_prev >= 0. && dt > 0.) {
            double dhdt = (h - h_prev) / dt;
            Fdamp = -lub_damping * dhdt;
        }

        LubForceEntry *lf = &lub_forces[n_lub_forces++];
        lf->xa = xa;
        lf->xb = xb;
        lf->dir.x = dx/dist;
        lf->dir.y = dy/dist;
        lf->Fmag = total_force_magnitude(h) + Fdamp;
    }
}

// point force to cell face force conversion
static void inject_face_force (face vector a, coord x0, coord dir, double Fmag, double sign)
{
    foreach_face(x) {
        double dx = x - x0.x;
        double dy = y - x0.y;
        double r2 = dx*dx + dy*dy;
        double w = exp(-r2 / (2.*sq(0.3*Delta)));
        a.x[] += sign * Fmag * dir.x * w;   
    }
    foreach_face(y) {
        double dx = x - x0.x;
        double dy = y - x0.y;
        double r2 = dx*dx + dy*dy;
        double w = exp(-r2 / (2.*sq(0.3*Delta)));
        a.y[] += sign * Fmag * dir.y * w;  
    }
}
static void apply_lubrication_forces (face vector a)
{    
    for (int k = 0; k < n_lub_forces; k++) {
        LubForceEntry *lf = &lub_forces[k];

        //equal and opposite forces on the two droplets
        inject_face_force (a, lf->xa, lf->dir, lf->Fmag, +1.);
        inject_face_force (a, lf->xb, lf->dir, lf->Fmag, -1.);
    }
}