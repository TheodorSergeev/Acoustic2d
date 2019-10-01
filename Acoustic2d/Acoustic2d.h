#ifndef ACOUSTIC2D_H
#define ACOUSTIC2D_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::cout;
using std::to_string;

enum BoundCondType { REFL, MUR, PML };

template <class TYPE>
struct Coord
{

    TYPE x, y;

    Coord(TYPE x_, TYPE y_):
        x(x_), y(y_)
    {}

};

struct AcVars  //acoustic variables
{


    AcVars(double p_, double u_, double v_, double Q_ = 0.0, double R_ = 0.0):
        p(p_), u(u_), v(v_), Q(Q_), R(R_)
    {}

    double p, // sound pressure
           u, // particle velocity x
           v, // particle velocity y
           Q,
           R;

};

class Acoustic2d
{

protected:

    BoundCondType bound_cond; // bounary condition type

    double den;     // density
    double el_rat;  // volume elastic ratio
    double wave_sp; // wave speed

    vector < vector<AcVars> > next_sol_arr; // t+1
    vector < vector<AcVars> > curr_sol_arr; // t

    vector < vector<AcVars> >& next_sol;    // references to the arrays above
    vector < vector<AcVars> >& curr_sol;    // to avoid copying when swapping (time step)

    Coord <unsigned int> grid_size;  // grid size
    Coord <double>       length;     // x and y dimensions of the rectangular field
    Coord <double>       step;       // x and y spatial grid step sizes

    double time_lim;
    double t_step;

    double courant_num;

    void CheckCoord(Coord <double>& r);
    void CheckVal(double val, double left_lim, double right_lim); // check whether x is in [0, length]
    void Check();                 // check whether variables are defined correctly
    void Dump();                  // print techical info

    void Iteration(int t_step_num);

    void ReflBoundCond();
    void MurBoundCond();
    void PmlBoundCond();
    double R(int node_depth, int pml_depth, double R_max, double power);

public:

    virtual AcVars InitCond (Coord<double>& r);

    Acoustic2d();

    Acoustic2d(BoundCondType bound_cond_, double len, int nodes,
               double time_lim, double time_step,
               double density, double elastic_ratio, double wave_speed);

    void RecordCurrSol(string& prefix, double t);

    void Solver();

    virtual void user_action(double t);

};

#endif
