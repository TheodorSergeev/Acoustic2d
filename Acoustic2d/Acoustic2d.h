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

enum BoundCondType { NONE, REFL, MUR, PML };
enum SchemeType    { FD, TVD };

template <class TYPE>
struct Coord
{

    TYPE x, y;

    Coord(TYPE x_, TYPE y_):
        x(x_), y(y_)
    {}

};

/*struct StressTensor
{

    double xx, xy, yy;

};*/

struct AcVars  //acoustic variables
{

    AcVars(double u_, double v_, double p_, double Q_ = 0.0, double R_ = 0.0):
        u(u_), v(v_), p(p_), Q(Q_), R(R_)
    {}

    double p, // sound pressure
           u, // particle velocity x
           v, // particle velocity y
           Q,
           R;

};


#define p_wave_sp wave_sp
#define s_wave_sp wave_sp

enum RiemanInvType {X, Y};

struct RiemanInv
{

    double w1, w2, w3;

    RiemanInv(const RiemanInv& cpy):
        w1(cpy.w1), w2(cpy.w2), w3(cpy.w3)
    {}

    RiemanInv(const AcVars& sol, RiemanInvType type, double den, double wave_sp)
    {


        if(type == X)
        {

            w1 = sol.v;
            w2 = + sol.u * den * wave_sp / 2 + sol.p / 2;
            w3 = - sol.u * den * wave_sp / 2 + sol.p / 2;

        }
        else if(type == Y)
        {

            w1 = sol.u;
            w2 = + sol.v * den * wave_sp / 2 + sol.p / 2;
            w3 = - sol.v * den * wave_sp / 2 + sol.p / 2;

        }

    }

    RiemanInv(double w1_, double w2_, double w3_):
        w1(w1_), w2(w2_), w3(w3_)
    {}

};


class Acoustic2d
{

protected:

    BoundCondType bound_cond; // bounary condition type
    SchemeType scheme_type;

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
    void FD_Iteration(int t_step_num);
    void TVD_Iteration(int t_step_num);
    void TVD_Step_x(int t_step_num);
    void TVD_Step_y(int t_step_num);

    void ReflBoundCond();
    void MurBoundCond();
    void PmlBoundCond();
    double R(int node_depth, int pml_depth, double R_max, double power);

public:

    virtual AcVars InitCond (Coord<double>& r);

    Acoustic2d();

    Acoustic2d(SchemeType scheme_type_, BoundCondType bound_cond_, 
               double len, int nodes,
               double time_lim, double time_step,
               double density, double elastic_ratio, 
               double wave_speed);

    void RecordCurrSol(string& prefix, double t);

    void Solver();

    virtual void user_action(double t);
};

#endif
