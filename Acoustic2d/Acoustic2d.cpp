#include "Acoustic2d.h"


void Acoustic2d::CheckVal(double val, double left_lim, double right_lim)
{

    if(val < left_lim || val > right_lim)
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


Acoustic2d::Acoustic2d(double len, int nodes,
                       double time_lim, double time_step,
                       double density, double elastic_ratio, double wave_speed):
    length(len, len), step(0, 0),  grid_size(nodes, nodes),
    time_lim(time_lim), t_step(time_step),
    next_sol_arr(), curr_sol_arr(),
    next_sol(next_sol_arr), curr_sol(curr_sol_arr),
    den(density), el_rat(elastic_ratio), wave_sp(wave_speed)
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

    Check();
    Dump();

}

Acoustic2d::Acoustic2d():
    length(0.0, 0.0), time_lim(0), t_step(0), step(0.0, 0.0), courant_num(0), grid_size(0, 0),
    next_sol_arr(), curr_sol_arr(),
    next_sol(next_sol_arr), curr_sol(curr_sol_arr),
    den(0.0), el_rat(0.0)
{

    Dump();

}


void Acoustic2d::RecordCurrSol(string& prefix, double t)
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

        output << (i + 0.5) * step.x << ", " << (j + 0.5) * step.y << ", 0.0, " << curr_sol[i][j].p << "\n";

    }

    output.close();

}

void Acoustic2d::user_action(double t)
{

    string pref = "";
    RecordCurrSol(pref, t);

}

#define X0 length.x / 2 + 0.5 * step.x

AcVars Acoustic2d::InitCond(Coord<double>& r)
{

    CheckCoord(r);
    AcVars val(0.0, 0.0, 0.0); // cos(pow(r.x * r.x + r.y * r.y, 0.5))

    double rad = pow((r.x - X0) * (r.x - X0) + (r.y - X0) * (r.y - X0), 0.5);

    if(rad < step.x * 10)
        val.p = 1.0 * cos(rad * M_PI / 2 / (step.x * 10));

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
        //PmlIteration(j); //tmp!!!!
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

    double coef = (wave_sp * t_step - step.x) / (wave_sp * t_step + step.x);

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        next_sol[i][0].u = curr_sol[i][1].u - coef * (curr_sol[i][0].u - next_sol[i][1].u);
        next_sol[i][0].v = curr_sol[i][1].v - coef * (curr_sol[i][0].v - next_sol[i][1].v);

        next_sol[i][grid_size.y - 1].u = curr_sol[i][grid_size.y - 2].u +
                                         coef * (next_sol[i][grid_size.y - 2].u - curr_sol[i][grid_size.y - 1].u);
        next_sol[i][grid_size.y - 1].v = curr_sol[i][grid_size.y - 2].v +
                                         coef * (next_sol[i][grid_size.y - 2].v - curr_sol[i][grid_size.y - 1].v);

    }

    for(unsigned int j = 0; j < grid_size.y; ++j)
    {

        next_sol[0][j].u = curr_sol[1][j].u - coef * (curr_sol[0][j].u - next_sol[1][j].u);
        next_sol[0][j].v = curr_sol[1][j].v - coef * (curr_sol[0][j].v - next_sol[1][j].v);

        next_sol[grid_size.x - 1][j].u = curr_sol[grid_size.x - 2][j].u +
                                         coef * (next_sol[grid_size.x - 2][j].u - curr_sol[grid_size.x - 1][j].u);
        next_sol[grid_size.x - 1][j].v = curr_sol[grid_size.x - 2][j].v +
                                         coef * (next_sol[grid_size.x - 2][j].v - curr_sol[grid_size.x - 1][j].v);

    }

}

double Acoustic2d::R(int node_pos, int pml_depth, double R_max, double power) // left boundaries
{

    if(node_pos >= pml_depth)
        return 0.0;

    return R_max * pow((double) node_pos / pml_depth, power);

}

void Acoustic2d::PmlBoundCond()
{

    const Coord <int> pml_depth(10, 10);
    const Coord <double> pml_coef(0.0, 0.0);
    double c = 1.0;

    ReflBoundCond();

    for(unsigned int i = 1; i < pml_depth.x; ++i)
    {

        for(unsigned int j = 1; j < pml_depth.y; ++j)
        {

            next_sol[i][j].u = curr_sol[i][j].u * (1.0 - R(i, pml_depth.x, pml_coef.x, 2) * t_step / den) -
                               t_step / den / step.x * (next_sol[i][j].p - next_sol[i - 1][j].p);
            next_sol[i][j].v = curr_sol[i][j].v * (1.0 - R(i, pml_depth.x, pml_coef.y, 2) * t_step / den) -
                               t_step / den / step.y * (next_sol[i][j].p - next_sol[i][j - 1].p);

            next_sol[i][j].p = curr_sol[i][j].p * (1.0 - R(i, pml_depth.x, pml_coef.x, 2)  * t_step / den) -
                               el_rat * t_step * ((curr_sol[i + 1][j].u - curr_sol[i][j].u) / step.x +
                               (curr_sol[i][j + 1].v - curr_sol[i][j].v) / step.y);

        }

    }

}

void Acoustic2d::PmlIteration(int t_step_num)
{
    const Coord <int> pml_depth(10, 10);
    const Coord <double> pml_coef(0.0, 20.0); //den / t_step * 0.8;
    double c = pow(el_rat / den, 0.5);

    for(unsigned int i = 0; i < grid_size.x; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y; ++j)
        {

            if(i == 0 || j == 0 || i == grid_size.x - 1 || j == grid_size.y - 1)
            {

                next_sol[i][j].u = 0.0; // pml!
                next_sol[i][j].v = 0.0; // pml!
                next_sol[i][j].Q = 0.0;
                next_sol[i][j].R = 0.0;

            }
            else
            {

                // 4 side - commented
                //double Rx = R(std::min(i, (grid_size.x - i)), pml_depth.x, pml_coef.x, 2);
                //double Ry = R(std::min(j, (grid_size.y - j)), pml_depth.y, pml_coef.y, 2);
                double Rx = R(i, pml_depth.x, pml_coef.x, 2);
                double Ry = R(j, pml_depth.y, pml_coef.y, 2);

                next_sol[i][j].u = curr_sol[i][j].u * (1.0 - Rx * t_step * wave_sp) -
                                   t_step / den * (curr_sol[i][j].p - curr_sol[i - 1][j].p) / step.x ;
                next_sol[i][j].v = curr_sol[i][j].v * (1.0 - Ry * t_step * wave_sp) -
                                   t_step / den * (curr_sol[i][j].p - curr_sol[i][j - 1].p) / step.y;

                next_sol[i][j].Q = curr_sol[i][j].Q + t_step * wave_sp * curr_sol[i][j].u;
                next_sol[i][j].R = curr_sol[i][j].R + t_step * wave_sp * curr_sol[i][j].v;

            }

        }

    }

    for(unsigned int i = 0; i < grid_size.x - 1; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y - 1; ++j)
        {

            //double Rx = R(std::min(i, (grid_size.x - i)), pml_depth.x, pml_coef.x, 2);
            //double Ry = R(std::min(j, (grid_size.y - j)), pml_depth.y, pml_coef.y, 2);
            double Rx = R(i, pml_depth.x, pml_coef.x, 2);
            double Ry = R(j, pml_depth.y, pml_coef.y, 2);

            next_sol[i][j].p =  curr_sol[i][j].p - t_step * (den * wave_sp * wave_sp *
                                                             ((next_sol[i + 1][j].u - next_sol[i][j].u) / step.x +
                                                              (next_sol[i][j + 1].v - next_sol[i][j].v) / step.y +
                                                              (next_sol[i + 1][j].Q - next_sol[i][j].Q) / step.y * Ry +
                                                              (next_sol[i][j + 1].R - next_sol[i][j].R) / step.x * Rx) +
                                                             wave_sp * (Rx + Ry) * curr_sol[i][j].p);

        }

    }

    swap(curr_sol, next_sol);

}


void Acoustic2d::Iteration(int t_step_num)
{
    const Coord <int> pml_depth(1, 1);


    for(unsigned int i = 0; i < grid_size.x - 1; ++i)
    {

        for(unsigned int j = 0; j < grid_size.y - 1; ++j)
        {

            //if(i != grid_size.x - 1 && j != grid_size.y - 1)
            {

                next_sol[i][j].p = curr_sol[i][j].p - el_rat * t_step * ((curr_sol[i + 1][j].u - curr_sol[i][j].u) / step.x +
                                                                        (curr_sol[i][j + 1].v - curr_sol[i][j].v) / step.y);

            }

        }

    }

    for(unsigned int i = 1; i < grid_size.x - 1; ++i)
    {

        for(unsigned int j = 1; j < grid_size.y - 1; ++j)
        {

            next_sol[i][j].u = curr_sol[i][j].u - t_step / den / step.x * (next_sol[i][j].p - next_sol[i - 1][j].p);
            next_sol[i][j].v = curr_sol[i][j].v - t_step / den / step.y * (next_sol[i][j].p - next_sol[i][j - 1].p);

        }

    }

    //ReflBoundCond();
    MurBoundCond();
    //PmlBoundCond();


    swap(curr_sol, next_sol);

}
