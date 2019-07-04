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

enum BoundCondType { NONE, DIRICHLET, NEUMANN, MUR };

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


    AcVars(double p_, double u_, double v_):
        p(p_), u(u_), v(v_)
    {}

    double p, // sound pressure
           u, // particle velocity x
           v; // particle velocity y

};

class Acoustic2d
{

protected:

    double den;     // density
    double el_rat;  // volume elastic ratio

    vector < vector<AcVars> > next_sol_arr; // t+1
    vector < vector<AcVars> > curr_sol_arr; // t

    vector < vector<AcVars> >& next_sol;    // references to arrays above
    vector < vector<AcVars> >& curr_sol;    // to avoid copying when swapping

    Coord <unsigned int> grid_size;  // grid size
    Coord <double>       length;     // x and y dimensiions of the rectangular field
    Coord <double>       step;       // x and y spatial grid sizes

    double time_lim;
    double t_step;

    double courant_num;

    void CheckCoord(Coord <double>& r);
    void CheckVal(double val, double left_lim, double right_lim); // check whether x is in [0, length]
    void Check();                 // check whether variables are defined correctly
    void Dump();                  // print techical info

public:

    virtual AcVars InitCond (Coord<double>& r);

    Acoustic2d();

    Acoustic2d(double len, int nodes,
               double time_lim, double time_step,
               double density, double elastic_ratio);

    void RecordCurrSol(string& prefix, double t);

    virtual void user_action(double t);

    void Solver();
    void Iteration(int t_step_num);

};

#endif
