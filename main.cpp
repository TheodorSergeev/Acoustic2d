#include "Acoustic2d.h"

// indent

int main()
{

    try
    {

        const double len = 1.0;
        const int nodes = 51;
        const double time_lim = 1.0, time_step = 0.005;
        const double dens = 1.0, el_rat = 1.0;
        const double wave_sp = 2.0;

        const Coord <int> pml_depth(10, 10);
        const Coord <double> pml_coef(100.0, 100.0);

        Acoustic2d test(FD, PML, len, nodes,
                        time_lim, time_step,
                        dens, el_rat,
                        wave_sp,
                        pml_depth, pml_coef);
        test.Solver();

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
