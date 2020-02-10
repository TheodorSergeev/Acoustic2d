    #include "Acoustic2d.h"

// indent

int main()
{

    try
    {

        const double len = 1.4;
        const int nodes = 101 + 40;
        const double time_step = 0.005;
        const double time_lim = 400 * time_step;
        const double dens = 1.0;
        const double wave_sp = 1.0;
        const double el_rat = dens * wave_sp * wave_sp;

        // time_step <= x_step / (sqrt(2) * wave_sp)

        const Coord <int> pml_depth(20, 20);
        // 0.8 / time_step / wave_sp
        const double Rmax = 0.8 * dens * wave_sp;
        const Coord <double> pml_coef(Rmax, Rmax);

        /*Acoustic2d test1(TVD, PML, len, nodes,
                         time_lim, time_step,
                         dens, el_rat,
                         wave_sp,
                         pml_depth, pml_coef);*/

        Acoustic2d test2(FD, PML, len, nodes,
                         time_lim, time_step,
                         dens, el_rat,
                         wave_sp,
                         pml_depth, pml_coef);

        const Coord <int> pml_depth_mur(0, 0);
        const Coord <double> pml_coef_mur(0.0, 0.0);

        Acoustic2d test3(FD, MUR, 1.0, 101,
                         time_lim, time_step,
                         dens, el_rat,
                         wave_sp,
                         pml_depth_mur, pml_coef_mur);

        //test1.Solver();
        test3.Solver();
        test2.Solver();

        /*for(int i = 0; i < 141; ++i) {
            printf("%d %d %lf\n", i, i, test2.sigma_x(i));
        }*/

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
