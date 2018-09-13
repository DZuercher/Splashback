#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "splashback.h"
#include <vector>
#include <ctime>
#include <cstring>
#include <cmath>

int main(int argc, char** argv)
{
    /*
    cosmology a = cosmology(0.30, 0.0, -1.0, 0.0, 0.0476, 0.7, 2.726, 0.8, 0.96, log10(8.0), 1.0);
    printf("%.7e\n",a.Dlofz(0.5));
    */
    splashback a = splashback(0.10, 10.0, 10, (char *)"debug.dat", 21.0, 0.33, 100, false, 0, true);
    a.allocate_lens_memory(10000);
    for (int i=0; i<100; i++){
        for (int j=0; j<100; j++){
            a.process_lens(i/100.*2*M_PI, j/100.*2*M_PI-M_PI, 0.5, i);
        }
    }
    a.finalize_lenses();

    for (int i=0; i<100; i++){
        a.process_source(0.6*2*M_PI, 0.1*M_PI, 21.0, 0);
    }

    a.finalize_results();
    a.finalize_results(true);
}

