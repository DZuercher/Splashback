#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include "splashback.h"
#include "spline.h"

splashback::splashback(float xrmin, float xrmax, int xrbins, char *xoutfile, float xmag_limit, float xzmax, int xNjack, bool xdeproject, int xcolored, bool xverbose)
{
    fprintf(stderr,"Using right code");
    class_initialized=true;
    rmin=xrmin;
    rmax=xrmax;
    rbins=xrbins;
    logrmin=log10(rmin);
    logrmax=log10(rmax);
    logrdiff=(logrmax-logrmin)/rbins;
    mag_limit=xmag_limit;
    zmax = xzmax;
    Njack=xNjack;
    sprintf(outfile, "%s", xoutfile);

    deproject=xdeproject;
    verbose=xverbose;
    colored=xcolored;

    count_num=NULL;
    denom_count_num=NULL;
    Rarr=NULL;
    count_num = (double *)calloc(rbins*Njack,sizeof(double));
    Rarr = (double *)calloc(rbins+1, sizeof(double));
    for (int i=0;i<=rbins;i++)
    {
        float rrmin = pow(10., logrmin + i*logrdiff);
        Rarr[i] = rrmin;
    }

    // ========================================================
    // Initialize cosmology
    // ========================================================
    cosmo=NULL;
    cosmo = new cosmology(0.27,0.0,-1.0,0.0,0.0476,0.7,2.726,0.8,0.96,log10(8.0),1.0);
    lenses_alloc=false;
    lenses_finalized=false;
    zredmin=1.E30;
    zredmax=-1.E30;
    gee = 4.2994e-9;
    cee = 3.E5;
    Dcomzmax = cosmo->Dcofz(double(zmax));

    if (verbose){
        fprintf(stderr, "Initialized cosmology, and delta sigma sum arrays\n");
        fprintf(stderr, "Options passed rmin:%f rmax:%f rbins:%d outfile:%s \n", rmin, rmax, rbins, outfile);
    }
}

splashback::~splashback(){
    if (class_initialized){
	free(Rarr);
	free(count_num);
	delete cosmo;
    /*    */
	class_initialized=false;
    }
    if (verbose){
        fprintf(stderr, "Freed memory related to cosmology and the dsigma arrays\n");
    }
    if (lenses_alloc){
        free(ra      );
        free(dec     );
        free(zred    );
        free(wt      );
        free(datacen );
        free(c_lra   );
        free(s_lra   );
        free(s_ldec  );
        free(c_ldec  );
        free(Dcom    );
        free(maglimz );
        free(jackreg );
        free(denom_count_num);
        if (verbose){
            fprintf(stderr, "Freed memory related to lenses as lenses were allocated\n");
        }
        lenses_alloc=false;
    }else{
        if (verbose){
            fprintf(stderr, "Did not free memory related to lenses as lenses were not finalized\n");
        }
    }
    if (lenses_finalized){
        delete tree;
        delete timer;
        if (verbose){
            fprintf(stderr, "Freed memory related to tree and timer as lenses were finalized\n");
        }
        lenses_finalized=false;
    }else{
        if (verbose){
            fprintf(stderr, "Did not free memory related to tree and timer as lenses were not finalized\n");
        }
    }
}

int splashback::allocate_lens_memory(int xNcen)
{
    Ncen = xNcen;
    ra = dec = zred = datacen = c_lra = s_lra = s_ldec = c_ldec = Dcom = maglimz = NULL;
    jackreg = NULL;
    ra      =(float *)  calloc(Ncen,sizeof(float));
    dec     =(float *)  calloc(Ncen,sizeof(float));
    zred    =(float *)  calloc(Ncen,sizeof(float));
    wt      =(float *)  calloc(Ncen,sizeof(float));
    datacen =(float *)  calloc(Ncen*3,sizeof(float));
    c_lra   =(float *)  calloc(Ncen,sizeof(float));
    s_lra   =(float *)  calloc(Ncen,sizeof(float));
    s_ldec  =(float *)  calloc(Ncen,sizeof(float));
    c_ldec  =(float *)  calloc(Ncen,sizeof(float));
    Dcom    =(float *)  calloc(Ncen,sizeof(float));
    maglimz =(float *)  calloc(Ncen,sizeof(float));
    jackreg =(int   *)  calloc(Ncen,sizeof(int));
    denom_count_num = (double *)calloc(Ncen*8, sizeof(double));
    Nfill=0;
    lenses_alloc=true;
    if (verbose){
        fprintf(stderr, "Initialized lenses Ncen=%d\n", Ncen);
    }
    return 0;
}

int splashback::process_lens(float xra, float xdec, float xzred, int xjackreg, float xwt)
{
    if (Nfill==Ncen){
        fprintf(stderr, "Lens arrays were not allocated correctly, please call allocate_lens_memory once again with the correct number of lenses\n");
        //exit(1);
        return 1;
    }
    
    xra = xra*M_PI/180.;
    xdec = xdec*M_PI/180.;
    ra     [Nfill] = xra; 
    dec    [Nfill] = xdec;
    zred   [Nfill] = xzred;
    wt     [Nfill] = xwt;
    c_lra  [Nfill] = cos(xra);
    s_lra  [Nfill] = sin(xra);
    c_ldec [Nfill] = cos(xdec);
    s_ldec [Nfill] = sin(xdec);
    datacen[3*Nfill] = c_ldec[Nfill]*c_lra[Nfill];
    datacen[3*Nfill+1] = c_ldec[Nfill]*s_lra[Nfill];
    datacen[3*Nfill+2] = s_ldec[Nfill];
    Dcom   [Nfill] = cosmo->Dcofz(double(xzred));
    maglimz[Nfill] = mag_limit - 5.0*log10(cosmo->Dlofz(zmax)/cosmo->Dlofz(xzred));
    denom_count_num[Nfill] = 0.0;
    //maglimz[Nfill] = mag_limit - 5.0*log10(Dcomzmax/Dcom[Nfill]*(1.+zmax)/(1.+zred[Nfill]));
    //maglimz[Nfill] = mag_limit - 5.0*log10(Dcom[Nfill]/Dcomzmax*(1.+zmax)/(1.+zred[Nfill]));
    jackreg[Nfill] = xjackreg;
    if (xzred<zredmin) zredmin = xzred;
    if (xzred>zredmax) zredmax = xzred;
    Nfill++;
    if (verbose){
        fprintf(stderr, "Filled %d-th lenses with ra:%f dec:%f\n", Nfill, xra*180./M_PI, xdec*180./M_PI);
    }

    return 0;
}

int splashback::finalize_lenses()
{
    Dcommin = cosmo->Dcofz(zredmin);
    if (Ncen>Nfill){
        Ncen=Nfill;
        ra      =(float *)  realloc(ra     , Ncen*sizeof(float));
        dec     =(float *)  realloc(dec    , Ncen*sizeof(float));
        zred    =(float *)  realloc(zred   , Ncen*sizeof(float));
        datacen =(float *)  realloc(datacen, Ncen*3*sizeof(float));
        c_lra   =(float *)  realloc(c_lra  , Ncen*sizeof(float));
        s_lra   =(float *)  realloc(s_lra  , Ncen*sizeof(float));
        s_ldec  =(float *)  realloc(s_ldec , Ncen*sizeof(float));
        c_ldec  =(float *)  realloc(c_ldec , Ncen*sizeof(float));
        Dcom    =(float *)  realloc(Dcom   , Ncen*sizeof(float));
        jackreg =(int   *)  realloc(jackreg, Ncen*sizeof(int));
        denom_count_num = (double *) realloc(denom_count_num, Ncen*sizeof(double));
    }

    // Ok now use kdtree2
    dim=3;
    xper.resize(boost::extents[dim]);
    for (int j=0; j<dim; j++)
      xper[j] = -1.0;

    // rearrange is true
    timer = NULL;
    tree = NULL;
    timer = new boost::timer;
    tree = new kdtree2::KDTree(datacen,xper,Ncen,dim);
    tree->sort_results = true;
    printf("Tree done %.1lf seconds\n", timer->elapsed());

    Niter=0;

    lenses_finalized=true;

    if (verbose){
        fprintf(stderr, "Finalized %d lenses, now ready to ingest source data\n", Ncen);
    }

    return 0;

}

int splashback::test_searchrecord(){
}

double splashback::deprojection_kernel(double R_by_a){
    if (R_by_a<1.0) return 1.0;
    double x = 1./sqrt(pow(R_by_a,2)-1.0);
    return 2.0/M_PI*(atan(x)-x);
}


int splashback::process_source(float sra, float sdec, float smag, bool disable_magcheck, float color)
{

   //For spline cut
   //------------------------------------------------------------------------------------
   std::string line, sub_1, sub_2;
   std::istringstream os;
   std::vector<double> redshifts, cuts;
   double foo;
   std::ifstream inputFile("/work/dominik.zuercher/Output/match_PS/spline_data.dat");

   if (inputFile.is_open()){
        while ( getline(inputFile, line) ){
            sub_1 = line.substr(0, line.find(' '));
            sub_2 = line.substr(line.find(' '));
            os.str(sub_1);
            os >> foo;
            redshifts.push_back(foo);
            os.str("");
            os.clear();
            os.str(sub_2);
            os >> foo;
            cuts.push_back(foo);
            os.str("");
            os.clear();
        }
        inputFile.close();
    }

    tk::spline spl;
    spl.set_points(redshifts, cuts);
   //------------------------------------------------------------------------------------


    if(std::isnan(sra) || std::isnan(sdec) || std::isnan(smag)){
        return 1;
    }

    if(!lenses_finalized){
        fprintf(stderr, "Lenses were not finalized, please call finalize_lenses() before processing sources\n");
        //exit(1);
        return 1;
    }
    Niter++;
    sra *= (M_PI/180.);
    sdec *= (M_PI/180.);

    double c_sdec=cos(sdec);
    double s_sdec=sin(sdec);
    double s_sra=sin(sra);
    double c_sra=cos(sra);

    double sx=c_sdec*c_sra;
    double sy=c_sdec*s_sra;
    double sz=s_sdec;

    // Ok now get results
    kdtree2::KDTreeResultVector res;
    std::vector<float> qv(3);
    qv[0]=sx;
    qv[1]=sy;
    qv[2]=sz;
    float dis2 = 0.0;
    if (!deproject){
        dis2 = pow(rmax/Dcommin, 2.);
    }else{
        dis2 = pow(10.0*rmax/Dcommin, 2.);
    }
    tree->r_nearest(qv, dis2, res);
    //0 0.034433 0.260000 435.623138 0.150020
    if (verbose)
        fprintf(stderr, "%d %f %f %f %f\n", res.size(), sqrt(dis2), Dcommin, zredmin);
    //FILE *dfile =fopen("debug.dat", "w");


    //Old cut : 0.65 + (z-0.1)* 0.8/0.23

    for(int k=0;k<res.size();k++){
        int xx = res.at(k).idx;
        if (smag>maglimz[xx] && (!disable_magcheck)) continue;
        if (colored==1 && color< spl(zred[xx])) continue;
        if (colored==-1 && color>= spl(zred[xx])) continue;

        float dis2=sqrt(res.at(k).dis);
        float rp=dis2*Dcom[xx];
        float logrp = log10(dis2*Dcom[xx]);

        if ((!deproject) && (logrp<=logrmin || logrp>=logrmax)) continue;
        //std::cout<<szbest<<" "<<zred[xx]<<std::endl;

        int rpbin = int((logrp-logrmin)/logrdiff);
        int jackbin = jackreg[xx];
        if (!deproject){
            double add=wt[xx];
            for(int jj=0; jj<Njack; jj++)
                count_num[Njack*rpbin+jj]+=add;
            count_num[Njack*rpbin+jackbin]-=add;
        }else{
            int rpbindeproject = rpbin<rbins ? rpbin : rbins;
            for(int ii=0; ii<rpbindeproject; ii++)
            {
                double add = deprojection_kernel(rp/Rarr[ii+1]) - deprojection_kernel(rp/Rarr[ii]);
                for(int jj=0; jj<Njack; jj++)
                    count_num[Njack*ii+jj]+=add;
                count_num[Njack*ii+jackbin]-=add;
            }
        }

    }

    for(int xx=0;xx<Ncen;xx++){
        if (smag>maglimz[xx] && (!disable_magcheck)) continue;
        if (colored==1 && color<0.8 + (zred[xx]-0.1) * 0.8/0.23) continue;
        if (colored==-1 && color>=0.8 + (zred[xx]-0.1) * 0.8/0.23) continue;
        denom_count_num[xx]++;
    }
       //exit(11);


    // Now go through the indexing
    //if(Niter%100000==0) printf("%ld done %.1lf seconds\n", Niter, timer->elapsed());

    return 0;
}

int splashback::finalize_results(bool writeok){
    fprintf(stderr, "%e %e %e\n", logrmin, logrmax, logrdiff);
    FILE *fout=fopen(outfile, "w");

    for (int i=0;i<rbins;i++)
    {
        float rr = logrmin + (i+0.5)*logrdiff;
        float rrmin = pow(10., logrmin + i*logrdiff);
        float rrmax = pow(10., logrmin + (i+1)*logrdiff);
        float area = M_PI*(rrmax*rrmax-rrmin*rrmin);
        for (int jj=0; jj<Njack; jj++)
            fprintf(fout, "%e %le %e %e \n", rr, count_num[Njack*i+jj], rrmin/2+rrmax/2, area);
    }
    if (writeok) fprintf(fout, "#OK \n");
    fclose(fout);

    sprintf(outfile_denom, "%s.denom", outfile);
    fout=fopen(outfile_denom, "w");

    for (int i=0;i<Ncen;i++)
    {
        fprintf(fout, "%le \n", denom_count_num[i]);
    }
    //if (writeok) fprintf(fout, "#OK \n");
    fclose(fout);
    return 0;
}
