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

enum BoundCondType { NONE, REFL, MUR, PML, SPLID_PML };
enum SchemeType    { FD, TVD };

const string SCHEME_TYPE_NAMES[4] = {"FD", "TVD"};
const string BOUND_COND_NAMES [5] = {"NONE", "REFL", "MUR", "PML", "SPLID_PML"};

template <class TYPE>
struct Coord
{

    TYPE x, y;

    Coord(TYPE x_, TYPE y_):
        x(x_), y(y_)
    {}

    Coord(const Coord<TYPE>& cpy):
        x(cpy.x), y(cpy.y)
    {}

};

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

enum RiemanInvType {X, Y, PML_X, PML_Y, SPLID_PML_X, SPLID_PML_Y};

struct RiemanInv
{

    double w1, w2, w3, w4, w5;

    RiemanInv(const RiemanInv& cpy):
        w1(cpy.w1), w2(cpy.w2), w3(cpy.w3), w4(cpy.w4), w5(cpy.w5)
    {}

    RiemanInv(const AcVars& sol, RiemanInvType type, double den, double wave_sp, double pml_coef = 0.0)
    {

        double s = den * wave_sp / 2;
        double el_rat = den * wave_sp * wave_sp;

        if(type == X)
        {

            w1 = sol.v;
            w2 = + sol.u * s + sol.p / 2.0;
            w3 = - sol.u * s + sol.p / 2.0;

        }
        else if(type == Y)
        {

            w1 = sol.u;
            w2 = + sol.v * s + sol.p / 2.0;
            w3 = - sol.v * s + sol.p / 2.0;

        }
        else if(type == PML_X)
        {

            w1 = sol.R;
            w2 = sol.Q;
            w3 = sol.u;
            w4 = + sol.v * s + sol.p / 2.0 + s * pml_coef * sol.Q;
            w5 = - sol.v * s + sol.p / 2.0 - s * pml_coef * sol.Q;

        }
        else if(type == PML_Y)
        {

            w1 = sol.R;
            w2 = sol.Q;
            w3 = sol.v;
            w4 = + sol.u * s + sol.p / 2.0 + s * pml_coef * sol.R;
            w5 = - sol.u * s + sol.p / 2.0 - s * pml_coef * sol.R;

        }
        else if(type == SPLID_PML_X)
        {

            // p1 = sol.Q
            // p2 = sol.R

            w1 = + sol.R;
            w2 = sol.v;
            /*not required in calculations*/ w3 = 0.0;
            w4 = - sqrt(el_rat * den) * 0.5 * sol.u + 0.5 * (sol.Q + sol.R);
            w5 = + sqrt(el_rat * den) * 0.5 * sol.u + 0.5 * (sol.Q + sol.R);

        }
        else if(type == SPLID_PML_Y)
        {

            // p1 = sol.Q
            // p2 = sol.R

            w1 = - sol.Q;
            w2 = sol.u;
            /*not required in calculations*/ w3 = 0.0;
            w4 = - sqrt(el_rat * den) * 0.5 * sol.v + 0.5 * (sol.Q + sol.R);
            w5 = + sqrt(el_rat * den) * 0.5 * sol.v + 0.5 * (sol.Q + sol.R);

        }

    }

    RiemanInv(double w1_, double w2_, double w3_, double w4_ = 0.0, double w5_ = 0.0):
        w1(w1_), w2(w2_), w3(w3_), w4(w4_), w5(w5_)
    {}

};


class Acoustic2d
{

//protected:
public:
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
    void TvdPmlBoundCond(int t_step_num);
    void TvdSplidPmlBoundCond(int t_step_num);

    Coord <int> pml_depth;
    Coord <double> pml_coef;

    void ReflBoundCond();
    void MurBoundCond();
    void PmlBoundCond();
    void SplidPmlBoundCond();

    double R(int node_depth, int pml_depth, double R_max, double power);
    double sigma_x_shifted_mesh(double node);
    double sigma_y_shifted_mesh(double node);
    double sigma_x(double node);
    double sigma_y(double node);


    void record_curr_sol             (string& prefix, double t);
    void record_curr_sol_shifted_mesh(string& prefix, double t);
    void record_energy               (string& predix, double t);

public:

    virtual AcVars InitCond (Coord<double>& r);

    Acoustic2d();

    Acoustic2d(SchemeType scheme_type_, BoundCondType bound_cond_,
               double len, int nodes,
               double time_lim, double time_step,
               double density, double wave_speed,
               const Coord <int>& pml_depth_, const Coord <double>& pml_coef_
               );

    void Solver();

    virtual void user_action(double t);
};

#endif
