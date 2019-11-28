#include "Acoustic2d.h"

// indent

int main()
{

    try
    {

        const double len = 1.0;
        const int nodes = 51;
        const double time_lim = 5.0, time_step = 0.01;
        const double dens = 1.0, el_rat = 1.0;
        const double wave_sp = 1.0;

        const Coord <int> pml_depth(10, 10);
        const Coord <double> pml_coef(1.0, 1.0);
/*
        Acoustic2d test1(TVD, PML, len, nodes,
                         time_lim, time_step,
                         dens, el_rat,
                         wave_sp,
                         pml_depth, pml_coef);
*/
        Acoustic2d test2(FD, PML, len, nodes,
                         time_lim, time_step,
                         dens, el_rat,
                         wave_sp,
                         pml_depth, pml_coef);

        Acoustic2d test3(FD, MUR, len, nodes,
                         time_lim, time_step,
                         dens, el_rat,
                         wave_sp,
                         pml_depth, pml_coef);
        //test1.Solver();
        test2.Solver();
        //test3.Solver();

    }
    catch(string& err)
    {

        cout << "Exception: " << err << "\n";

    }
    catch(const char* err)
    {

        cout << "Exception: " << err << "\n";

    }
    catch(...)
    {

        cout << "Unknown exception\n" << "\n";

    }

    return 0;

}
