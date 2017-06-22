// expanded smallppm (code is exactly the same as smallppm.cpp but with more comments)

#include <math.h>   // smallppm, Progressive Photon Mapping by T. Hachisuka
#include <stdlib.h> // originally smallpt, a path tracer by Kevin Beason, 2008
#include <stdio.h>  // Usage: ./smallppm 100000 && xv image.ppm 
#include <ctime> 
#include <vector>
#include "fstream"
#include "Bondingbox.h"
#include "Crash.h"

#define PI ((double)3.14159265358979) // ^^^^^^:number of photons emitted
#define ALPHA ((double)0.7) // the alpha parameter of PPM

// Halton 序列反向排序 
int primes[61] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79,
        83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181,
        191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283
};

inline int rev(const int i, const int p) {
    if (i == 0) return i; 
	else return p - i;
}

double hal(const int b, int j) {
    const int p = primes[b];
    double h = 0.0, f = 1.0 / (double) p, fct = f;
    while (j > 0) {
        h += rev(j % p, p) * fct;
        j /= p;
        fct *= f;
    }
    return h;
}

List *ListAdd(Crash_point *point, List *L){
	List *_L = new List;
	_L->id = point;
	_L->next = L;
	return _L;
}


unsigned int num_hash, pixel_index, num_photon;
double hash_s;
List **hash_grid;
List *hitpoints = NULL;
AABB hpbbox;

// spatial hash function  //空间hash  
inline unsigned int hash(const int ix, const int iy, const int iz) {
    return (unsigned int) ((ix * 73856093) ^ (iy * 19349663) ^ (iz * 83492791)) % num_hash;
}

void build_hash_grid(const int w, const int h) {  //hash 网格 
    // find the bounding box of all the measurement points
    hpbbox.reset();
    List *lst = hitpoints;
    while (lst != NULL) {
        Crash_point *hp = lst->id; 
        lst = lst->next; 
        hpbbox.fit(hp->pos); 
    }

    // heuristic for initial radius  启发式 
    vector3 ssize = hpbbox.max - hpbbox.min;
    double irad = ((ssize.x + ssize.y + ssize.z) / 3.0) / ((w + h) / 2.0) * 2.0;

    // determine hash table size
    // we now find the bounding box of all the measurement points inflated by the initial radius
    hpbbox.reset();
    lst = hitpoints;
    int vphoton = 0;
    while (lst != NULL) {
        Crash_point *hp = lst->id;
        lst = lst->next;
        hp->r2 = irad * irad;
        hp->n = 0;
        hp->flux = vector3();
        vphoton++;
        hpbbox.fit(hp->pos - irad);
        hpbbox.fit(hp->pos + irad);
    }

    // make each grid cell two times larger than the initial radius
    hash_s = 1.0 / (irad * 2.0);
    num_hash = vphoton;

    // build the hash table
    hash_grid = new List *[num_hash];
    for (unsigned int i = 0; i < num_hash; i++) hash_grid[i] = NULL;
    lst = hitpoints;
    while (lst != NULL) {
        Crash_point *hp = lst->id;
        lst = lst->next;
        vector3 BMin = ((hp->pos - irad) - hpbbox.min) * hash_s;
        vector3 BMax = ((hp->pos + irad) - hpbbox.min) * hash_s;
        for (int iz = abs(int(BMin.z)); iz <= abs(int(BMax.z)); iz++) {
            for (int iy = abs(int(BMin.y)); iy <= abs(int(BMax.y)); iy++) {
                for (int ix = abs(int(BMin.x)); ix <= abs(int(BMax.x)); ix++) {
                    int hv = hash(ix, iy, iz);
                    hash_grid[hv] = ListAdd(hp, hash_grid[hv]);
                }
            }
        }
    }
}

struct Ray {  //ray  
    vector3 o, d;

    Ray() {};

    Ray(vector3 o_, vector3 d_) : o(o_), d(d_) {}
};

enum Refl_t {  
    DIFF, SPEC, REFR  //漫反射 镜面反射 折射  
};  // material types, used in radiance()

struct Sphere {  //球类 
    double rad;
    vector3 p, c;
    Refl_t refl;

    Sphere(double r_, vector3 p_, vector3 c_, Refl_t re_) : rad(r_), p(p_), c(c_), refl(re_) {}  //radius  postion color reflect_type 

    inline double intersect(const Ray &r) const {  //球体求交 
        // ray-sphere intersection returns distance
        vector3 op = p - r.o;
        double t, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        if (det < 0) {
            return 1e20;
        } else {
            det = sqrt(det);
        }
        return (t = b - det) > 1e-4 ? t : ((t = b + det) > 1e-4 ? t : 1e20);
    }
};

Sphere sph[] = { // Scene: radius, position, color, material
        Sphere(1e5, vector3(1e5 + 1, 40.8, 81.6), vector3(.75, .25, .25), DIFF),//Left
        Sphere(1e5, vector3(-1e5 + 99, 40.8, 81.6), vector3(.25, .25, .75), DIFF),//Right
        Sphere(1e5, vector3(50, 40.8, 1e5), vector3(.75, .75, .75), DIFF),//Back
        Sphere(1e5, vector3(50, 40.8, -1e5 + 170), vector3(), DIFF),//Front
        Sphere(1e5, vector3(50, 1e5, 81.6), vector3(.75, .75, .75), DIFF),//Bottomm
        Sphere(1e5, vector3(50, -1e5 + 81.6, 81.6), vector3(.75, .75, .75), DIFF),//Top
        Sphere(16.5, vector3(27, 16.5, 47), vector3(1, 1, 1) * .999, SPEC),//Mirror
        Sphere(16.5, vector3(73, 16.5, 88), vector3(1, 1, 1) * .999, REFR),//Glass
        Sphere(8.5, vector3(50, 8.5, 60), vector3(1, 1, 1) * .999, DIFF)//Middle
};

struct Triangle {
    vector3 p0, p1, p2, e, c;  //color  emission 
    Refl_t refl;  //折射类型 

    Triangle(vector3 p0_, vector3 p1_, vector3 p2_, vector3 e_, vector3 c_, Refl_t refl_) :
            p0(p0_), p1(p1_), p2(p2_), e(e_), c(c_), refl(refl_) {}

    double intersect(const Ray &r) const { // return distance
        vector3 e1 = p0 - p1;
        vector3 e2 = p0 - p2;
        vector3 s = p0 - r.o;
        double t, eps = 1e-4, beta, gamma;
        double det1 = det(r.d, e1, e2);
        t = det(s, e1, e2) / det1;
        beta = det(r.d, s, e2) / det1;
        gamma = det(r.d, e1, s) / det1;
        if (beta > 1 || beta < 0 || gamma > 1 || gamma < 0 || (gamma + beta) > 1)
            return 0;
        return t > eps ? t : 0;
    }

    double det(vector3 a, vector3 b, vector3 c) const {
        double a1 = a.x * b.y * c.z;
        double a2 = b.x * c.y * a.z;
        double a3 = c.x * a.y * b.z;
        double b1 = c.x * b.y * a.z;
        double b2 = b.x * a.y * c.z;
        double b3 = a.x * c.y * b.z;
        return a1 + a2 + a3 - b1 - b2 - b3;
    }
};

std::vector<Triangle> tri;  //存三角面片 

void loadObj(const char* filename) {
    std::ifstream fin(filename);
    std::vector<vector3> vertices;
    std::string vof;
    double x, y, z;
    int a, b, c;
    while(!fin.eof()) {
        fin >> vof;
        if(vof == "v") {
            fin >> x >> y >> z;
            x=(x+4)*10, y=(y+4)*10, z=(z+4)*10;
            vertices.push_back(vector3(x,y,z));
        } else {
            fin >> a >> b >> c;
            a--, b--, c--;
            tri.push_back(Triangle(vertices[a], vertices[b], vertices[c], vector3(), vector3(.25, .25, .25), DIFF));
        }
    }
}

// tone mapping and gamma correction  色调映射与伽玛校正
int toInt(double x) {
    return int(pow(1 - exp(-x), 1 / 2.2) * 255 + .5);
}

double m, n;
// find the closet interection  找三角面片的 最近的交点 
inline bool intersect(const Ray &r, double &t, int &id) {
    double d, inf = 1e20;
    t = inf;
    for (int i = int(m + n - 1); i > n - 1; i--)
        if ((d = tri[i - int(n)].intersect(r)) && d < t) {
            t = d;
            id = i;
        }
    for (int i = 0; i < n; i++) {
        d = sph[i].intersect(r);
        if (d < t) {
            t = d;
            id = i;
        }
    }
    return t < inf;
}

// generate a photon ray from the point light source with QMC  生成光子射线 
void genp(Ray *pr, vector3 *f, int i) {
    *f = vector3(2500, 2500, 2500) * (PI * 4.0); // flux
    double p = 2. * PI * hal(0, i), t = 2. * acos(sqrt(1. - hal(1, i)));
    double st = sin(t);
    pr->d = vector3(cos(p) * st, cos(t), sin(p) * st);
    pr->o = vector3(50, 60, 85);
}

void trace(const Ray &r, int dpt, bool m, const vector3 &fl, const vector3 &adj, int i) {
    double t;
    int id;

    dpt++;
    if (!intersect(r, t, id) || (dpt >= 20))return;

    int d3 = dpt * 3;
    if(id < n) {
        const Sphere &obj = sph[id];
        vector3 x = r.o + r.d * t, n = (x - obj.p).normllize(), f = obj.c;
        vector3 nl = n.dot(r.d) < 0 ? n : n * -1;
        double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

        if (obj.refl == DIFF) {
            // Lambertian

            // use QMC to sample the next direction
            double r1 = 2. * PI * hal(d3 - 1, i), r2 = hal(d3 + 0, i);
            double r2s = sqrt(r2);
            vector3 w = nl, u = ((fabs(w.x) > .1 ? vector3(0, 1) : vector3(1)) % w).normllize();
            vector3 v = w % u, d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normllize();

            if (m) {
                // eye ray
                // store the measurment point
                Crash_point *hp = new Crash_point;
                hp->f = f.mul(adj);
                hp->pos = x;
                hp->nrm = n;
                hp->pix = pixel_index;
                hitpoints = ListAdd(hp, hitpoints);
            } else {
                // photon ray
                // find neighboring measurement points and accumulate flux via progressive density estimation
                vector3 hh = (x - hpbbox.min) * hash_s;
                int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
                // strictly speaking, we should use #pragma omp critical here.
                // it usually works without an artifact due to the fact that photons are
                // rarely accumulated to the same measurement points at the same time (especially with QMC).
                // it is also significantly faster.
                {
                    List *hp = hash_grid[hash(ix, iy, iz)];
                    while (hp != NULL) {
                        Crash_point *hitpoint = hp->id;
                        hp = hp->next;
                        vector3 v = hitpoint->pos - x;
                        // check normllizeals to be closer than 90 degree (avoids some edge brightning)
                        if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= hitpoint->r2)) {
                            // unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
                            double g = (hitpoint->n * ALPHA + ALPHA) / (hitpoint->n * ALPHA + 1.0);
                            hitpoint->r2 = hitpoint->r2 * g;
                            hitpoint->n++;
                            hitpoint->flux = (hitpoint->flux + hitpoint->f.mul(fl) * (1. / PI)) * g;
                        }
                    }
                }
                if (hal(d3 + 1, i) < p) trace(Ray(x, d), dpt, m, f.mul(fl) * (1. / p), adj, i);
            }

        } else if (obj.refl == SPEC) {
            // mirror
            trace(Ray(x, r.d - n * 2.0 * n.dot(r.d)), dpt, m, f.mul(fl), f.mul(adj), i);

        } else {
            // glass
            Ray lr(x, r.d - n * 2.0 * n.dot(r.d));
            bool into = (n.dot(nl) > 0.0);
            double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

            // total internal reflection
            if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) return trace(lr, dpt, m, fl, adj, i);

            vector3 td = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normllize();
            double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.dot(n));
            double Re = R0 + (1 - R0) * c * c * c * c * c, P = Re;
            Ray rr(x, td);
            vector3 fa = f.mul(adj);
            if (m) {
                // eye ray (trace both rays)
                trace(lr, dpt, m, fl, fa * Re, i);
                trace(rr, dpt, m, fl, fa * (1.0 - Re), i);
            } else {
                // photon ray (pick one via Russian roulette)
                (hal(d3 - 1, i) < P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
            }
        }

    } else {
        printf("%s\n", "hit tri");
        const Triangle &obj = tri[id - int(n)];        // the hit object
        vector3 x = r.o + r.d * t, n = (x - obj.p0).normllize(), f = obj.c;  // 三角形 的两条边 叉乘 需要改 
        vector3 nl = n.dot(r.d) < 0 ? n : n * -1;
        double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z;

        if (obj.refl == DIFF) {
            // Lambertian

            // use QMC to sample the next direction
            double r1 = 2. * PI * hal(d3 - 1, i), r2 = hal(d3 + 0, i);
            double r2s = sqrt(r2);
            vector3 w = nl, u = ((fabs(w.x) > .1 ? vector3(0, 1) : vector3(1)) % w).normllize();
            vector3 v = w % u, d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).normllize();

            if (m) {
                // eye ray
                // store the measurment point
                Crash_point *hp = new Crash_point;
                hp->f = f.mul(adj);
                hp->pos = x;
                hp->nrm = n;
                hp->pix = pixel_index;
                hitpoints = ListAdd(hp, hitpoints);
            } else {
                // photon ray
                // find neighboring measurement points and accumulate flux via progressive density estimation
                vector3 hh = (x - hpbbox.min) * hash_s;
                int ix = abs(int(hh.x)), iy = abs(int(hh.y)), iz = abs(int(hh.z));
                // strictly speaking, we should use #pragma omp critical here.
                // it usually works without an artifact due to the fact that photons are
                // rarely accumulated to the same measurement points at the same time (especially with QMC).
                // it is also significantly faster.
                {
                    List *hp = hash_grid[hash(ix, iy, iz)];
                    while (hp != NULL) {
                        Crash_point *hitpoint = hp->id;
                        hp = hp->next;
                        vector3 v = hitpoint->pos - x;
                        // check normllizeals to be closer than 90 degree (avoids some edge brightning)
                        if ((hitpoint->nrm.dot(n) > 1e-3) && (v.dot(v) <= hitpoint->r2)) {
                            // unlike N in the paper, hitpoint->n stores "N / ALPHA" to make it an integer value
                            double g = (hitpoint->n * ALPHA + ALPHA) / (hitpoint->n * ALPHA + 1.0);
                            hitpoint->r2 = hitpoint->r2 * g;
                            hitpoint->n++;
                            hitpoint->flux = (hitpoint->flux + hitpoint->f.mul(fl) * (1. / PI)) * g;
                        }
                    }
                }
                if (hal(d3 + 1, i) < p) trace(Ray(x, d), dpt, m, f.mul(fl) * (1. / p), adj, i);
            }

        } else if (obj.refl == SPEC) {
            // mirror
            trace(Ray(x, r.d - n * 2.0 * n.dot(r.d)), dpt, m, f.mul(fl), f.mul(adj), i);

        } else {
            // glass
            Ray lr(x, r.d - n * 2.0 * n.dot(r.d));
            bool into = (n.dot(nl) > 0.0);
            double nc = 1.0, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;

            // total internal reflection
            if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) return trace(lr, dpt, m, fl, adj, i);

            vector3 td = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).normllize();
            double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : td.dot(n));
            double Re = R0 + (1 - R0) * c * c * c * c * c, P = Re;
            Ray rr(x, td);
            vector3 fa = f.mul(adj);
            if (m) {
                // eye ray (trace both rays)
                trace(lr, dpt, m, fl, fa * Re, i);
                trace(rr, dpt, m, fl, fa * (1.0 - Re), i);
            } else {
                // photon ray (pick one via Russian roulette)
                (hal(d3 - 1, i) < P) ? trace(lr, dpt, m, fl, fa, i) : trace(rr, dpt, m, fl, fa, i);
            }
        }
    }
}

int main(int argc, char *argv[]) {
    // samps * 1000 photon paths will be traced
    int w = 300, h = 300, samps = (argc == 2) ? MAX(atoi(argv[1]) / 1000, 1) : 1000;

    // trace eye rays and store measurement points
    Ray cam(vector3(50, 48, 295.6), vector3(0, -0.042612, -1).normllize());
    vector3 cx = vector3(w * .5135 / h), cy = (cx % cam.d).normllize() * .5135, *c = new vector3[w * h], vw;
//    tri.push_back(Triangle(vector3(2, 40, 2), vector3(2, 20, 75), vector3(35, 20, 45), vector3(), vector3(), DIFF));
//    tri.push_back(Triangle(vector3(2, 30, 2), vector3(2, 50, 75), vector3(22, 20, 45), vector3(), vector3(), DIFF));
    loadObj("cup.obj");
    n = sizeof(sph) / sizeof(Sphere);
    m = tri.size();
//#pragma omp parallel for schedule(dynamic, 1)
    for (int y = 0; y < h; y++) {
        fprintf(stderr, "\rHitPointPass %5.2f%%", 100.0 * y / (h - 1));
        for (int x = 0; x < w; x++) {
            pixel_index = x + y * w;
            vector3 d = cx * ((x + 0.5) / w - 0.5) + cy * (-(y + 0.5) / h + 0.5) + cam.d;
            trace(Ray(cam.o + d * 140, d.normllize()), 0, true, vector3(), vector3(1, 1, 1), 0);
        }
    }
    //fprintf(stderr, "\n");  刚打的注释 

    // build the hash table over the measurement points
    build_hash_grid(w, h);

    // trace photon rays with multi-threading
    num_photon = samps;
    vw = vector3(1, 1, 1);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < num_photon; i++) {
        double p = 100. * (i + 1) / num_photon;
        fprintf(stderr, "\rPhotonPass %5.2f%%", p);
        int m = 1000 * i;
        Ray r;
        vector3 f;
        for (int j = 0; j < 1000; j++) {
            genp(&r, &f, m + j);
            trace(r, 0, 0 > 1, f, vw, m + j);
        }
    }

    // density estimation
    List *lst = hitpoints;
    while (lst != NULL) {
        Crash_point *hp = lst->id;
        lst = lst->next;
        int i = hp->pix;
        c[i] = c[i] + hp->flux * (1.0 / (PI * hp->r2 * num_photon * 1000.0));
    }

    // save the image after tone mapping and gamma correction
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"%m-%d-%H-%M-%S.ppm",now);
    FILE *f = fopen(buffer, "w");         // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++) {
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
    }
}
