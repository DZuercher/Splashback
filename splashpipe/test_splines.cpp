#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include "spline.h"

/*
#include "./alglib/src/ap.h"
//disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC) && !defined(AE_ALL_WARNINGS)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4611)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
#include "./alglib/src/alglibinternal.h"
#include "./alglib/src/alglibmisc.h"
//#include "./alglib/src/diffequations.h"
#include "./alglib/src/linalg.h"
#include "./alglib/src/optimization.h"
//#include "./alglib/src/solvers.h"
#include "./alglib/src/statistics.h"
#include "./alglib/src/dataanalysis.h"
//#include "./alglib/src/specialfunctions.h"
#include "./alglib/src/integration.h"
//#include "./alglib/src/fasttransforms.h"
#include "./alglib/src/interpolation.h"
*/


int main(){

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
            //std::cout<<x<<'\n';
            cuts.push_back(foo);
	    os.str("");
	    os.clear();
	}
	inputFile.close();
    }
    
    std::cout << "redshifts \n";
    for (std::vector<double>::const_iterator i = redshifts.begin(); i != redshifts.end(); ++i){
        std::cout << *i << ' ';
    }
    std::cout << "colors \n";
    for (std::vector<double>::const_iterator i = cuts.begin(); i != cuts.end(); ++i){
        std::cout << *i << ' ';
    }
    std::cout << '\n';
    tk::spline spl;
    spl.set_points(redshifts, cuts);
    
    printf("spline at %f is %f\n", 0.25, double(spl(0.25)));

    /*
    alglib::real_1d_array x;
    alglib::real_1d_array y;



    x.setcontent(redshifts.size(), &redshifts[0]);
    y.setcontent(cuts.size(), &cuts[0]);

    std::cout << int(x.length()) << '\n';

    for (int i = 0; i < int(y.length());++i){
        std::cout << y[i] << '\n';
	}

    alglib::spline1dinterpolant spl;
    alglib::spline1dbuildcubic(x, y, spl);
    printf("spline at %f is %f\n", 0.25, double(alglib::spline1dcalc(spl, 0.25)));
    */
 
    return 0;
}

