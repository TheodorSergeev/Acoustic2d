#include "Acoustic2d.h"

// indent

int main()
{

    try
    {

        const double len = 1.0;
        const int nodes = 51;
        const double time_lim = 1.0, time_step = 0.01;
        const double dens = 1.0, el_rat = 1.0, wave_sp = 1.0;

        Acoustic2d test(len, nodes,
                        time_lim, time_step,
                        dens, el_rat, wave_sp);
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
