#ifndef SPLASHBACK_H
#define SPLASHBACK_H

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <fstream>
#include "cosmology.h"
#include "kdtree2.h"
#include <vector>
#include <ctime>
#include <cstring>
#include <cmath>
#include <boost/multi_array.hpp>
typedef boost::multi_array<float,1> array1dfloat;
#include <boost/timer.hpp>

class splashback;

class splashback
{
    protected:
        cosmology *cosmo;
        bool verbose;
        float rmin;
        float rmax;
        float logrmin;
        float logrmax;
        float logrdiff;
        float mag_limit;
        float zmax, Dcomzmax;
        bool deproject;
        int rbins;
        int colored;
        float *ra, *dec, *zred, *lam, *datacen, *c_lra, *s_lra, *c_ldec, *s_ldec, *Dcom, *maglimz, *wt;
        double *Rarr;
        int *jackreg;
        int Ncen, Nfill;
        float zredmin, zredmax, Dcommin;
        double *count_num;
        double *denom_count_num;
        char outfile[1000];
        char outfile_denom[1006];
        boost::timer *timer;
        long int Niter;
        int Njack;
        bool lenses_alloc, lenses_finalized;
        array1dfloat xper;
        int dim;
        bool class_initialized;

        double gee;
        double cee;

    public:
        kdtree2::KDTree *tree;
        //splashback();
        splashback(float xrmin=0.5, float xrmax=15.0, int xrbins=15, char *xoutfile=(char *)"Debug.dat", float xmag_limit=25.0, float zmax=0.5, int Njack=25, bool deproject=false, int colored=0, bool verbose=false);
        ~splashback();
        int allocate_lens_memory(int xNcen);
        int process_lens(float xra, float xdec, float xzred, int xjackreg, float wt=1.0);
        int finalize_lenses();

        int process_source(float sra, float sdec, float smag, bool disable_magcheck, float color=1.E30);
        int finalize_results(bool writeok=false);
        int test_searchrecord();
        double deprojection_kernel(double x);
};

#endif
