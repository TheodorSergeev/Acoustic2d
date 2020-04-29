#include "Acoustic2d.h"

// indent

int main()
{

    try
    {

        const int PML_DEPTH = 40;

        const double len = 1.8;
        const int nodes = 101;
        const double time_step = 0.005;
        const double time_lim = 1000 * time_step;//2 * 600 * time_step;
        const double dens = 1.0;
        const double wave_sp = 1.0;

        // time_step <= x_step / (sqrt(2) * wave_sp)

        const Coord <int> pml_depth(PML_DEPTH, PML_DEPTH);
        const double Rmax = 20.0 * dens * wave_sp;
        const Coord <double> pml_coef(Rmax, Rmax);

        Acoustic2d test4(FD, SPLID_PML, len, nodes + 2 * PML_DEPTH,
                         time_lim, time_step,
                         dens, wave_sp,
                         pml_depth, pml_coef);
        test4.Solver();

        Acoustic2d test5(TVD, SPLID_PML, len, nodes + 2 * PML_DEPTH,
                         time_lim, time_step,
                         dens, wave_sp,
                         pml_depth, pml_coef);
        test5.Solver();

        Acoustic2d test2(FD, PML, len, nodes + 2 * PML_DEPTH,
                         time_lim, time_step,
                         dens, wave_sp,
                         pml_depth, pml_coef);
        test2.Solver();

        Acoustic2d test1(TVD, PML, 1.8, nodes + 2 * PML_DEPTH,
                         time_lim, time_step,
                         dens, wave_sp,
                         pml_depth, pml_coef);

        test1.Solver();

        const Coord <int> pml_depth_mur(0, 0);
        const Coord <double> pml_coef_mur(0.0, 0.0);

        Acoustic2d test3(FD, MUR, 1.0, nodes,
                         time_lim, time_step,
                         dens, wave_sp,
                         pml_depth_mur, pml_coef_mur);
        test3.Solver();


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
