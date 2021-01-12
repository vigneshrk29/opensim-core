#include "Simbody.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int Nelt = 9; // # elements in a rotation matrix
constexpr int NR = 1;   // # outputs
/// typedef for CasAD int
typedef long long int casadi_int;
/// sparsity information
static const casadi_int sp_in1[] = { Nelt, 1, 0, Nelt, 0, 1, 2, 3, 4, 5, 6, 7, 8 };
static const casadi_int sp_in2[] = { Nelt, 1, 0, Nelt, 0, 1, 2, 3, 4, 5, 6, 7, 8 };
static const casadi_int sp_out[] = { NR, 1, 0, NR, 0 };
static std::vector<casadi_int> sp_jac(1);
/// nominal information
static const casadi_int sp_in1_null[] = { Nelt, 1, 0, 0 };
static const casadi_int sp_in2_null[] = { Nelt, 1, 0, 0 };
static const casadi_int sp_out_null[] = { NR, 1, 0, 0 };

#if defined _WIN32
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT
#endif // defined _WIN32

/* Function F, using templated type T
F(x,u) -> (tau)
*/

template<typename T>
T value(const double& e) { return e; }

template<typename T>
int F_generic(const T** arg, T** res) {

    // Read inputs
    std::vector<T> R_GS_vec(arg[0], arg[0] + Nelt);
    std::vector<T> R_GO_vec(arg[1], arg[1] + Nelt);

    // Construct rotation matrices from vectors.
    Mat33 R_GS_Mat33(
        R_GS_vec[0], R_GS_vec[1], R_GS_vec[2],
        R_GS_vec[3], R_GS_vec[4], R_GS_vec[5],
        R_GS_vec[6], R_GS_vec[7], R_GS_vec[8]);
    Rotation R_GS(R_GS_Mat33);

    Mat33 R_GO_Mat33(
        R_GO_vec[0], R_GO_vec[1], R_GO_vec[2],
        R_GO_vec[3], R_GO_vec[4], R_GO_vec[5],
        R_GO_vec[6], R_GO_vec[7], R_GO_vec[8]);
    Rotation R_GO(R_GO_Mat33);

    // Operations.
    auto R_SO = ~R_GS * R_GO;
    //auto Q_SO = R_SO.convertRotationToQuaternion();
    //// This part is inspired from convertQuaternionToAngleAxis().
    //// TODO: This is an attempt to make things work with AD.
    //// I removed the conditional statements and will instead generate
    //// two expression graphs. One for which the conditional statement
    //// canonicalizing the quaternion is active and one for which that
    //// statement is not active. Doing so should prevent the need for
    //// further conditional statements when calculating the angle from
    //// the angle-axis representation. Note that I also removed the part
    //// that assigns 0 when the value is very small, since it might not
    //// be needed for the tracking purpose.    
    //auto& ca2 = Q_SO[0];
    //auto& sa2v = Q_SO.getSubVec<3>(1);
    //auto sa2 = sa2v.norm();
    //// This should return an angle between [0, Pi].
    //auto angle = 2 * atan2(sa2, ca2);
    auto aa_SO = R_SO.convertRotationToAngleAxis();

    res[0][0] = value<T>(aa_SO[0]);


    return 0;
}

#ifdef __cplusplus
extern "C" {
#endif

    /*
    F(x,u,p) -> (Tau)
    */

    DLL_EXPORT casadi_int F(const double** arg, double** res) {
        F_generic<double>(arg, res);
        return 0;
    }

    DLL_EXPORT casadi_int F_n_in(void) { return n_in; }

    DLL_EXPORT casadi_int F_n_out(void) { return n_out; }

    DLL_EXPORT const casadi_int* F_sparsity_in(casadi_int i) {
        switch (i) {
        case 0: return sp_in1;
        case 1: return sp_in2;
        default: return 0;
        }
    }

    DLL_EXPORT const casadi_int* F_sparsity_out(casadi_int i) {
        switch (i) {
        case 0: return sp_out;
        default: return 0;
        }
    }

#ifdef __cplusplus
} /* extern "C" */
#endif
