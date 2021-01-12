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
constexpr int Nelt = 3; // # elements in a rotation matrix
constexpr int NR = 1;   // # outputs
/// typedef for CasAD int
typedef long long int casadi_int;
/// sparsity information
static const casadi_int sp_in1[] = { Nelt, 1, 0, Nelt, 0, 1, 2 };
static const casadi_int sp_in2[] = { Nelt, 1, 0, Nelt, 0, 1, 2 };
static const casadi_int sp_out[] = { NR, 1, 0, NR, 0 };
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
    std::vector<T> in1(arg[0], arg[0] + Nelt);
    std::vector<T> in2(arg[1], arg[1] + Nelt);

    // Construct rotation matrices from vectors.
    SimTK::Rotation R_GS(
        SimTK::BodyOrSpaceType::BodyRotationSequence,
        in1[0], SimTK::XAxis,
        in1[1], SimTK::YAxis,
        in1[2], SimTK::ZAxis);

    SimTK::Rotation R_GO(
        SimTK::BodyOrSpaceType::BodyRotationSequence,
        in2[0], SimTK::XAxis,
        in2[1], SimTK::YAxis,
        in2[2], SimTK::ZAxis);

    // Operations.
    auto R_SO = ~R_GS * R_GO;
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
