#include "Acoustic2d.h"

#define EPS 1e-14

void Acoustic2d::CheckVal(double val, double left_lim, double right_lim)
{

    if(val < left_lim - EPS || val > right_lim + EPS)
        throw "Value val = " + to_string(val) + " but cannot be < " +
              to_string(left_lim) + " or > " + to_string(right_lim);

}


void Acoustic2d::CheckCoord(Coord<double>& r)
{

    CheckVal(r.x, 0.0, length.x);
    CheckVal(r.y, 0.0, length.y);

}


void Acoustic2d::Check()
{

    if(next_sol.size() != grid_size.x)
        throw "grid_size.x = " + to_string(grid_size.x) + " is not equal to next_sol.size() = " + to_string(next_sol.size());

    if(curr_sol.size() > 0)
    {

        if(curr_sol[0].size() != grid_size.y)
            throw "grid_size.y = " + to_string(grid_size.y) + " is not equal to curr_sol.size() = " + to_string(curr_sol.size());

    }

    if(length.x     <= 0) throw "length.x = "     + to_string(length.x)     + " cannot be < 0";
    if(length.y     <= 0) throw "length.y = "     + to_string(length.y)     + " cannot be < 0";
    if(time_lim     <= 0) throw "time_lim = "     + to_string(time_lim)     + " cannot be < 0";
    if(t_step       <= 0) throw "t_step = "       + to_string(t_step)       + " cannot be < 0";
    if(step.x       <= 0) throw "step.x = " + to_string(step.x) + " cannot be < 0";
    if(step.y       <= 0) throw "step.y = " + to_string(step.y) + " cannot be < 0";
    if(courant_num  <= 0) throw "courant_num = "  + to_string(courant_num)  + " cannot be < 0";

}

void Acoustic2d::Dump()
{

    cout << "length.x = "    << length.x     << "\n"
         << "length.y = "    << length.y     << "\n"
         << "time_lim = "    << time_lim     << "\n"
         << "t_step = "      << t_step       << "\n"
         << "x_step = "      << step.x << "\n"
         << "y_step = "      << step.y << "\n"
         << "courant_num = " << courant_num  << "\n"
         << "grid_size.x = " << grid_size.x  << "\n"
         << "grid_size.y = " << grid_size.y  << "\n";

    cout << "\n";

    /*cout << "p:\n";

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y; ++j)
            printf("%.3f ", curr_sol[i][j].p);
            //printf("(%.1f,%.1f,%.1f) ", curr_sol[i][j].u, curr_sol[i][j].v, curr_sol[i][j].p);

        cout << "\n";

    }*/

}


Acoustic2d::Acoustic2d(SchemeType scheme_type_, BoundCondType bound_cond_,
                       double len, int nodes,
                       double time_lim, double time_step,
                       double density, double elastic_ratio,
                       double w_sp,
                       const Coord <int>& pml_depth_,
                       const Coord <double>& pml_coef_
                       ):
    scheme_type(scheme_type_), bound_cond(bound_cond_),
    length(len, len), step(0.0, 0.0),  grid_size(nodes, nodes),
    time_lim(time_lim), t_step(time_step),
    next_sol_arr(), curr_sol_arr(),
    next_sol(next_sol_arr), curr_sol(curr_sol_arr),
    den(density), el_rat(elastic_ratio), wave_sp(w_sp),
    pml_depth(pml_depth_), pml_coef(pml_coef_)
{

    if(nodes < 2)
    {

        throw "Number of nodes = " + to_string(nodes) +  " cannot be < 2\n";

    }
    else
    {

        step = Coord <double> (len / (nodes - 1), len / (nodes - 1));

    }

    AcVars zero(0.0, 0.0, 0.0);
    vector <AcVars> line;
    line.assign(grid_size.y, zero);
    curr_sol_arr.assign(grid_size.x, line);
    next_sol_arr = curr_sol_arr;
    courant_num = t_step * wave_sp / step.x;

    printf("pml_depth = %d %d\n", pml_depth.x, pml_depth.y);
    printf("pml_coef = %lf %lf\n", pml_coef.x, pml_coef.y);

    Check();
    Dump();

}

Acoustic2d::Acoustic2d():
    length(0.0, 0.0), time_lim(0.0), t_step(0.0),
    step(0.0, 0.0), courant_num(0.0), grid_size(0, 0),
    next_sol_arr(), curr_sol_arr(),
    next_sol(next_sol_arr), curr_sol(curr_sol_arr),
    den(0.0), el_rat(0.0), pml_depth(0, 0), pml_coef(0.0, 0.0)
{

    Dump();

}


void Acoustic2d::record_curr_sol(string& prefix, double t)
{

    std::ofstream output;
    //string fname = "data/" + prefix + "_" + to_string(t) + ".txt";

    string fname = "data/test.csv." + to_string(int(t / t_step));
    output.open(fname, std::ifstream::out);

    if(!output.is_open() || !output.good())
        throw "Bad output filename " + fname + " in RecordCurrSol (check if the directory exists).\n";
    output << "x coord, y coord, z coord, scalar\n";

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y; ++j)
        {

            output << i * step.x << ", "
                   << j * step.y << ", "
                   << "0.0"                  << ", "
                   << curr_sol[i][j].p << "\n";
        }

    }

    output.close();

}

void Acoustic2d::record_energy(string& predix, double t)
{

    double energy = 0.0;

    for(unsigned int i = pml_depth.x; i < grid_size.x - pml_depth.x; ++i)
    {

        for(unsigned int j = pml_depth.y; j < grid_size.y - pml_depth.y; ++j)
        {

            energy += (curr_sol[i][j].u * curr_sol[i][j].u +
                       curr_sol[i][j].v * curr_sol[i][j].v);

            //energy = energy + curr_sol[i][j].p * 1.5 * (step.x * step.y * 1.0);  // K = 3/2 pV

        }

    }

    // E = mv^2 / 2
    // m = V * rho
    //energy *= (step.x * step.y * 1.0) * den;

    std::ofstream output;
    //string fname = "data/" + prefix + "_" + to_string(t) + ".txt";

    string fname = "energy/energy_" + SCHEME_TYPE_NAMES[scheme_type] +
                   "_" + BOUND_COND_NAMES[bound_cond] + ".csv";

    if(t < t_step)
        output.open(fname, std::ifstream::out | std::ifstream::trunc);
    else
        output.open(fname, std::ifstream::out | std::ifstream::app);

    if(!output.is_open() || !output.good())
        throw "Bad output filename " + fname + " in RecordCurrSol (check if the directory exists).\n";
    //output << "E, t\n";
    output << energy << " " << t << "\n";

    output.close();

}

void Acoustic2d::record_curr_sol_shifted_mesh(string& prefix, double t)
{

    std::ofstream output;
    //string fname = "data/" + prefix + "_" + to_string(t) + ".txt";

    string fname = "data/test.csv." + to_string(int(t / t_step));
    output.open(fname, std::ifstream::out);

    if(!output.is_open() || !output.good())
        throw "Bad output filename " + fname + " in RecordCurrSol (check if the directory exists).\n";
    output << "x coord, y coord, z coord, scalar\n";

    for(unsigned int i = 0; i < grid_size.x - 1; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y - 1; ++j)
        {

            output << (i + 0.5) * step.x << ", "
                   << (j + 0.5) * step.y << ", "
                   << "0.0"              << ", "
                   << curr_sol[i][j].p   << ", "
                   << sigma_x(i)         << ", "
                   << sigma_y(j)         << "\n";

        }

    }

    output.close();

}

void Acoustic2d::user_action(double t)
{

    string pref = "";

    switch(scheme_type)
    {

        case TVD:
            record_curr_sol(pref, t);
            break;

        case FD:
            record_curr_sol_shifted_mesh(pref, t);
            break;

    }

    record_energy(pref, t);

}

AcVars Acoustic2d::InitCond(Coord<double>& r)
{

    CheckCoord(r);
    AcVars val(0.0, 0.0, 0.0); // cos(pow(r.x * r.x + r.y * r.y, 0.5))
    double X0 = 0.0;

    switch(scheme_type)
    {

        case TVD:
            X0 = length.x * 0.5 /*+ 0.5 * step.x*/;
            break;

        case FD:
            X0 = length.x * 0.5 - 0.5 * step.x;
            break;

    }

    const double rad = pow((r.x - X0) * (r.x - X0) + (r.y - X0) * (r.y - X0), 0.5);
    const double amp = 1.0;

    val.p = 0.0;

    const double init_rad = step.x * 30;

    if(rad < init_rad)
        val.p += amp * cos(rad * M_PI / 2 / (init_rad));

    return val;

}

void Acoustic2d::Solver()
{

    Check();

    int t_steps_num = (int) round(time_lim / t_step);

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y; ++j)
        {

            Coord <double> pos(step.x * i, step.y * j);
            curr_sol[i][j] = InitCond(pos);

        }

    }

    user_action(0.0);

    for(int j = 1; j <= t_steps_num; ++j)
    {

        Iteration(j);
        user_action(j * t_step);

    }

}

void Acoustic2d::ReflBoundCond()
{

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        next_sol[i][0].u = 0.0;
        next_sol[i][0].v = 0.0;
        next_sol[i][grid_size.y - 1].u = 0.0;
        next_sol[i][grid_size.y - 1].v = 0.0;

    }

    for(unsigned int j = 0; j < grid_size.y; ++j)
    {

        next_sol[0][j].u = 0.0;
        next_sol[0][j].v = 0.0;
        next_sol[grid_size.x - 1][j].u = 0.0;
        next_sol[grid_size.x - 1][j].v = 0.0;

    }

}

void Acoustic2d::MurBoundCond()
{

    double coef_x = (wave_sp * t_step - step.x) / (wave_sp * t_step + step.x);
    double coef_y = (wave_sp * t_step - step.y) / (wave_sp * t_step + step.y);

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        next_sol[i][0].u = curr_sol[i][1].u - coef_x * (curr_sol[i][0].u - next_sol[i][1].u);
        next_sol[i][0].v = curr_sol[i][1].v - coef_y * (curr_sol[i][0].v - next_sol[i][1].v);

        next_sol[i][grid_size.y - 1].u = curr_sol[i][grid_size.y - 2].u +
                                         coef_x * (next_sol[i][grid_size.y - 2].u -
                                                   curr_sol[i][grid_size.y - 1].u);
        next_sol[i][grid_size.y - 1].v = curr_sol[i][grid_size.y - 2].v +
                                         coef_y * (next_sol[i][grid_size.y - 2].v -
                                                   curr_sol[i][grid_size.y - 1].v);

    }

    for(unsigned int j = 0; j < grid_size.y; ++j)
    {

        next_sol[0][j].u = curr_sol[1][j].u - coef_x * (curr_sol[0][j].u - next_sol[1][j].u);
        next_sol[0][j].v = curr_sol[1][j].v - coef_y * (curr_sol[0][j].v - next_sol[1][j].v);

        next_sol[grid_size.x - 1][j].u = curr_sol[grid_size.x - 2][j].u +
                                         coef_x * (next_sol[grid_size.x - 2][j].u -
                                                   curr_sol[grid_size.x - 1][j].u);
        next_sol[grid_size.x - 1][j].v = curr_sol[grid_size.x - 2][j].v +
                                         coef_y * (next_sol[grid_size.x - 2][j].v -
                                                   curr_sol[grid_size.x - 1][j].v);

    }

}

double Acoustic2d::R(int node_pos, int pml_depth, double R_max, double power) // left boundaries
{

    if(node_pos > pml_depth || node_pos < 0)
        throw "node_pos in R = " + to_string(node_pos);

    double R_val = R_max * pow((double) node_pos / pml_depth, power);

    return R_val;

}

double Acoustic2d::sigma_x(int node)
{

    double j = (double) node;
    double x_size = (double) grid_size.x;

    double penetration_depth = 0.0;

    if(j <= pml_depth.x)
        penetration_depth = pml_depth.x - j;
    else if(j > x_size - 1.0 - pml_depth.x)
        penetration_depth = j - (x_size - 1.0 - pml_depth.x);

    double sigma_val = R(penetration_depth,
                         pml_depth.x, pml_coef.x, 2.0);

    return sigma_val;

}

double Acoustic2d::sigma_y(int node)
{

    double j = (double) node;
    double y_size = (double) grid_size.y;

    double penetration_depth = 0.0;

    if(j <= pml_depth.y)
        penetration_depth = pml_depth.y - j;
    else if(j > y_size - 1.0 - pml_depth.y)
        penetration_depth = j - (y_size - 1.0 - pml_depth.y);

    double sigma_val = R(penetration_depth,
                         pml_depth.y, pml_coef.y, 2.0);

    return sigma_val;

}


void Acoustic2d::PmlBoundCond()
{
    for (unsigned int i = 1; i < grid_size.x - 1; ++i){

        double Rx = sigma_x(i);

        for (unsigned int j = 1; j < grid_size.y - 1; ++j){

            double Ry = sigma_y(j);

            next_sol[i][j].u = curr_sol[i][j].u * (1.0 - Rx * t_step * wave_sp) -
                               t_step / den / step.x * (curr_sol[i][j].p - curr_sol[i - 1][j].p);
            next_sol[i][j].v = curr_sol[i][j].v * (1.0 - Ry * t_step * wave_sp) -
                               t_step / den / step.y * (curr_sol[i][j].p - curr_sol[i][j - 1].p);
        }
    }

    MurBoundCond();

    for (unsigned int i = 0; i < grid_size.x; ++i) {

        double Rx = sigma_x(i);

        for (unsigned int j = 0; j < grid_size.y; ++j) {

            double Ry = sigma_y(j);

            next_sol[i][j].Q = curr_sol[i][j].Q + t_step * wave_sp * curr_sol[i][j].u;
            next_sol[i][j].R = curr_sol[i][j].R + t_step * wave_sp * curr_sol[i][j].v;
        }
    }

    for (unsigned int i = 0; i < grid_size.x - 1; ++i) {

        double Rx = sigma_x(i);

        for (unsigned int j = 0; j < grid_size.y - 1; ++j) {

            double Ry = sigma_y(j);

            next_sol[i][j].p =  curr_sol[i][j].p -
                                t_step * (den * wave_sp * wave_sp *
                                ((next_sol[i + 1][j].u - next_sol[i][j].u) / step.x +
                                 (next_sol[i][j + 1].v - next_sol[i][j].v) / step.y +
                                 (next_sol[i + 1][j].Q - next_sol[i][j].Q) / step.y * Rx +
                                 (next_sol[i][j + 1].R - next_sol[i][j].R) / step.x * Ry) +
                                  wave_sp * (Rx + Ry) * curr_sol[i][j].p);
        }
    }
}

void Acoustic2d::FD_Iteration(int t_step_num)
{

    double Rx = 0.0;
    double Ry = 0.0;

    for(unsigned int i = 1; i < grid_size.x - 1; ++i)
    {

        for(unsigned int j = 1; j < grid_size.y - 1; ++j)
        {

            next_sol[i][j].u = curr_sol[i][j].u * (1.0 - Rx * t_step * wave_sp) -
                               t_step / den / step.x * (curr_sol[i][j].p - curr_sol[i - 1][j].p);
            next_sol[i][j].v = curr_sol[i][j].v * (1.0 - Ry * t_step * wave_sp) -
                               t_step / den / step.y * (curr_sol[i][j].p - curr_sol[i][j - 1].p);

        }

    }

    MurBoundCond();

    for(unsigned int i = 0; i < grid_size.x - 1; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y - 1; ++j)
        {

            next_sol[i][j].p =  curr_sol[i][j].p -
                                t_step * (den * wave_sp * wave_sp *
                                ((next_sol[i + 1][j].u - next_sol[i][j].u) / step.x +
                                 (next_sol[i][j + 1].v - next_sol[i][j].v) / step.y +
                                 (next_sol[i + 1][j].Q - next_sol[i][j].Q) / step.y * Rx +
                                 (next_sol[i][j + 1].R - next_sol[i][j].R) / step.x * Ry) +
                                  wave_sp * (Rx + Ry) * curr_sol[i][j].p);

        }

    }

}

void Acoustic2d::Iteration(int t_step_num)
{

    /*if(scheme_type == TVD && bound_cond == PML)
    {

        TvdPmlBoundCond(t_step_num);
        //ReflBoundCond();
        swap(curr_sol, next_sol);
        return;

    }*/

    switch(scheme_type)
    {

        case TVD:
            TVD_Iteration(t_step_num);
            break;

        case FD:
            if(bound_cond != PML)
                FD_Iteration(t_step_num);
            break;

    }

    switch(bound_cond)
    {

        case REFL:
            ReflBoundCond();
            break;
        case MUR:
            //MurBoundCond();
            break;
        case PML:
            PmlBoundCond();
            break;
        case NONE:
            break;

    }

    swap(curr_sol, next_sol);

}

#define real double
inline real le_min(real a, real b) { return a > b ? b : a; }
inline real le_max(real a, real b) { return a > b ? a : b; }
inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

inline real tvd2(const double c, const double u_2, const double u_1, const double u, const double u1)
{

    #define TVD2_EPS 1e-6
    #define limiter limiter_superbee
    #define limiter_superbee(r) (le_max3(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r)))

    double r1 = (u - u_1);
    double r2 = (u1 - u) + TVD2_EPS;
    const double r = r1 / r2;

    r1 = (u_1 - u_2) + TVD2_EPS;
    r2 = (u   - u_1) + TVD2_EPS;
    const double r_1 = r1 / r2;

    const double f12  = u   + limiter(r)   / 2.0 * (1.0 - c) * (u1 - u);
    const double f_12 = u_1 + limiter(r_1) / 2.0 * (1.0 - c) * (u  - u_1);

    return u - c * (f12 - f_12);

}

void Acoustic2d::TvdPmlBoundCond(int t_step_num)
{

    // Step X

    const double cour_num_x = (t_step / 1.0) * wave_sp / step.x;

    for(int j = 0; j < grid_size.y; ++j)
    {

        double sig_y = sigma_y(j);

        RiemanInv nnu(curr_sol[2][j], PML_X, den, wave_sp, sig_y); // w2
        RiemanInv nu (curr_sol[1][j], PML_X, den, wave_sp, sig_y); // w1
        RiemanInv u  (curr_sol[0][j], PML_X, den, wave_sp, sig_y); // w
        RiemanInv pu (u);                                   // w_1 = w
        RiemanInv ppu(u);                                   // w_2 = w

        for(int i = 0; i < grid_size.x; ++i)
        {

            RiemanInv d(tvd2(0.0, ppu.w1, pu.w1, u.w1, nu.w1),
                        tvd2(0.0, ppu.w2, pu.w2, u.w2, nu.w2),
                        tvd2(0.0, ppu.w3, pu.w3, u.w3, nu.w3),
                        tvd2(cour_num_x, nnu.w4, nu.w4, u.w4, pu.w4),
                        tvd2(cour_num_x, ppu.w5, pu.w5, u.w5, nu.w5));

            next_sol[i][j].u = d.w3;
            next_sol[i][j].v = - sig_y * d.w2 + (d.w4 - d.w5) / den / wave_sp;
            next_sol[i][j].p = d.w4 + d.w5;
            next_sol[i][j].Q = d.w2;
            next_sol[i][j].R = d.w1;

            ppu = pu;
            pu  = u;
            u   = nu;
            nu  = nnu;

            if(i < grid_size.x - 3)
                nnu = RiemanInv(curr_sol[i + 3][j], PML_X, den, wave_sp, sig_y);
            else
                nnu = RiemanInv(curr_sol[grid_size.x - 1][j], PML_X, den, wave_sp, sig_y);

        }

    }
    swap(curr_sol, next_sol);

    // Step Y

    const double cour_num_y = (t_step / 1.0) * wave_sp / step.y;

    for(int i = 0; i < grid_size.x; ++i)
    {

        double sig_x = sigma_x(i);

        RiemanInv nnu(curr_sol[i][2], PML_Y, den, wave_sp, sig_x);
        RiemanInv nu (curr_sol[i][1], PML_Y, den, wave_sp, sig_x);
        RiemanInv u  (curr_sol[i][0], PML_Y, den, wave_sp, sig_x);
        RiemanInv pu (u);
        RiemanInv ppu(u);

        for(int j = 0; j < grid_size.y; ++j)
        {

            RiemanInv d(tvd2(0.0, ppu.w1, pu.w1, u.w1, nu.w1),
                        tvd2(0.0, ppu.w2, pu.w2, u.w2, nu.w2),
                        tvd2(0.0, ppu.w3, pu.w3, u.w3, nu.w3),
                        tvd2(cour_num_y, nnu.w4, nu.w4, u.w4, pu.w4),
                        tvd2(cour_num_y, ppu.w5, pu.w5, u.w5, nu.w5));

            next_sol[i][j].u = - sig_x * d.w1 + (d.w4 - d.w5) / den / wave_sp;
            next_sol[i][j].v = d.w3;
            next_sol[i][j].p = d.w4 + d.w5;
            next_sol[i][j].Q = d.w2;
            next_sol[i][j].R = d.w1;

            ppu = pu;
            pu  = u;
            u   = nu;
            nu  = nnu;

            if(j < grid_size.y - 3)
                nnu = RiemanInv(curr_sol[i][j + 3], PML_Y, den, wave_sp, sig_x);
            else
                nnu = RiemanInv(curr_sol[i][grid_size.y - 1], PML_Y, den, wave_sp, sig_x);

        }

    }

    swap(curr_sol, next_sol);

    // Step analytical

    for(int i = 0; i < grid_size.x; ++i)
    {

        for(int j = 0; j < grid_size.y; ++j)
        {

            if(i == 0 || j == 0 || i == grid_size.x - 1 || j == grid_size.y - 1)
            {

                // This gives us MUR ABC if pml_coef = 0.0

                next_sol[i][j].u = 0.0; // pml!
                next_sol[i][j].v = 0.0; // pml!
                next_sol[i][j].p = 0.0; // pml!
                next_sol[i][j].Q = 0.0;
                next_sol[i][j].R = 0.0;

            }
            else
            {

                double sig_x = sigma_x(i);
                double sig_y = sigma_y(j);

                sig_x += 1e-12;
                sig_y += 1e-12;

                RiemanInv w(curr_sol[i][j].u / sig_x + curr_sol[i][j].R,
                            curr_sol[i][j].v / sig_y + curr_sol[i][j].Q,
                            - curr_sol[i][j].u / sig_x,
                            - curr_sol[i][j].v / sig_y,
                            curr_sol[i][j].p);

                double t = t_step / 3.0;
                sig_x -= 1e-12;
                sig_y -= 1e-12;

                RiemanInv d(w.w1,
                            w.w2,
                            w.w3 * exp(- wave_sp * sig_x * t),
                            w.w4 * exp(- wave_sp * sig_y * t),
                            w.w5 * exp(- wave_sp * (sig_x + sig_y) * t));

                sig_x += 1e-12;
                sig_y += 1e-12;

                next_sol[i][j].u = - sig_x * d.w3;
                next_sol[i][j].v = - sig_y * d.w4;
                next_sol[i][j].p = d.w5;
                next_sol[i][j].Q = d.w2 + d.w4;
                next_sol[i][j].R = d.w1 + d.w3;

            }

        }

    }

}

void Acoustic2d::TVD_Iteration(int t_step_num)
{

    /*
     * We solve regular hyperbolic system of PDE (http://en.wikipedia.org/wiki/Hyperbolic_partial_differential_equation) in form:
     * du/dt + Ax * du/dx + Ay * du/dy = 0.
     *
     * During time integration we use dimension split method:
     * 1. Step:
     * Integrate system dv/dt + Ax * dv/dx = 0, get u = v^(n + 1).
     * 2. Step:
     * Integrate system du/dt + Ay * du/dy = 0, get on next time step u^(n + 1).
     */

    // curr = t
    // next = -
    TVD_Step_x(t_step_num);
    // curr = t
    // next = t+1/2
    swap(curr_sol, next_sol);
    // curr = t+1/2
    // next = t
    TVD_Step_y(t_step_num);
    // curr = t+1/2
    // next = t+1

}



void Acoustic2d::TVD_Step_x(int t_step_num)
{

    const double cour_num_x = (t_step / 2.0) * wave_sp / step.x;

    for(int j = 0; j < grid_size.y/* - 1*/; ++j)
    {

        RiemanInv nnu(curr_sol[2][j], X, den, wave_sp); // w2
        RiemanInv nu (curr_sol[1][j], X, den, wave_sp); // w1
        RiemanInv u  (curr_sol[0][j], X, den, wave_sp); // w
        RiemanInv pu (u); // w_1 = w - valid for the linear case
        RiemanInv ppu(u); // w_2 = w - valid for the linear case

        for(int i = 0; i < grid_size.x/* - 1*/; ++i)
        {

            RiemanInv d(tvd2(0.0, ppu.w1, pu.w1, u.w1, nu.w1),
                        tvd2(cour_num_x, nnu.w2, nu.w2, u.w2, pu.w2),
                        tvd2(cour_num_x, ppu.w3, pu.w3, u.w3, nu.w3));

            next_sol[i][j].u = (d.w2 - d.w3) / den / wave_sp;
            next_sol[i][j].v = d.w1;
            next_sol[i][j].p = d.w2 + d.w3;

            ppu = pu;
            pu  = u;
            u   = nu;
            nu  = nnu;

            if(i < grid_size.x - 3)
                nnu = RiemanInv(curr_sol[i + 3][j], X, den, wave_sp);
            else
                nnu = RiemanInv(curr_sol[grid_size.x - 1][j], X, den, wave_sp);

        }

    }

}

void Acoustic2d::TVD_Step_y(int t_step_num)
{

    const double cour_num_y = (t_step / 2.0) * wave_sp / step.y;

    for(int i = 0; i < grid_size.x/* - 1*/; ++i)
    {

        RiemanInv nnu(curr_sol[i][2], Y, den, wave_sp);
        RiemanInv nu (curr_sol[i][1], Y, den, wave_sp);
        RiemanInv u  (curr_sol[i][0], Y, den, wave_sp);
        RiemanInv pu (u); // Valid for the linear case
        RiemanInv ppu(u); // Valid for the linear case

        for(int j = 0; j < grid_size.y/* - 1*/; ++j)
        {

            RiemanInv d(tvd2(0.0, ppu.w1, pu.w1, u.w1, nu.w1),
                        tvd2(cour_num_y, nnu.w2, nu.w2, u.w2, pu.w2),
                        tvd2(cour_num_y, ppu.w3, pu.w3, u.w3, nu.w3));

            //inc_y(curr_sol_arr[i][j], d, den, p_wave_sp, s_wave_sp);
            next_sol[i][j].u = d.w1;
            next_sol[i][j].v = (d.w2 - d.w3) / den / wave_sp;
            next_sol[i][j].p = d.w2 + d.w3;

            ppu = pu;
            pu  = u;
            u   = nu;
            nu  = nnu;

            if(j < grid_size.y - 3)
                nnu = RiemanInv(curr_sol[i][j + 3], Y, den, wave_sp);
            else
                nnu = RiemanInv(curr_sol[i][grid_size.y - 1], Y, den, wave_sp);

        }

    }

}
