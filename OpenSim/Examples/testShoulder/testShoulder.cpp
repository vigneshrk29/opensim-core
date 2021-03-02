#include <algorithm>
#include <cassert>
#include <chrono>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/InverseDynamicsSolver.h>

using namespace SimTK;
using namespace OpenSim;

// Declare inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 17;     // degrees of freedom
constexpr int NX = 2 * ndof; // states
constexpr int NU = ndof;     // controls
constexpr int NR = ndof;     // residual forces

// OpenSim and Simbody use different indices for the states/controls when the
// kinematic chain has joints up and down the origin (e.g., lumbar joint/arms
// and legs with pelvis as origin).
// The two following functions allow getting the indices from one reference
// system to the other. These functions are inspired from
// createSystemYIndexMap() in Moco.
// getIndicesOSInSimbody() returns the indices of the OpenSim Qs in the Simbody
// reference system. Note that we only care about the order here so we divide
// by 2 because the states include both Qs and Qdots.
SimTK::Array_<int> getIndicesOSInSimbody(const Model& model) {
    auto s = model.getWorkingState();
    const auto svNames = model.getStateVariableNames();
    SimTK::Array_<int> idxOSInSimbody(s.getNU());
    s.updU() = 0;
    for (int iy = 0; iy < s.getNU(); ++iy) {
        s.updU()[iy] = SimTK::NaN;
        const auto svValues = model.getStateVariableValues(s);
        for (int isv = 0; isv < svNames.size(); ++isv) {
            if (SimTK::isNaN(svValues[isv])) {
                s.updU()[iy] = 0;
                idxOSInSimbody[iy] = isv / 2;
                break;
            }
        }
    }
    return idxOSInSimbody;
}
// getIndicesSimbodyInOS() returns the indices of the Simbody Qs in the OpenSim
// reference system.
SimTK::Array_<int> getIndicesSimbodyInOS(const Model& model) {
    auto idxOSInSimbody = getIndicesOSInSimbody(model);
    auto s = model.getWorkingState();
    SimTK::Array_<int> idxSimbodyInOS(s.getNU());
    for (int iy = 0; iy < s.getNU(); ++iy) {
        for (int iyy = 0; iyy < s.getNU(); ++iyy) {
            if (idxOSInSimbody[iyy] == iy) {
                idxSimbodyInOS[iy] = iyy;
                break;
            }
        }
    }
    return idxSimbodyInOS;
}
int F_generic(const double** arg, double** res) {
    Model* model = new Model("ThoracoscapularShoulderModel.osim");

    /// Initialize  system and state.
    SimTK::State* state;
    state = new State(model->initSystem());
    // Read inputs
    std::vector<double> x(arg[0], arg[0] + NX);
    std::vector<double> u(arg[1], arg[1] + NU);

    // States and controls
    double ua[NU];   /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    QsUs.setToZero();
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// Controls
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
    for (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

	const auto svValues = model->getStateVariableValues(*state);
    std::cout << svValues << std::endl;

    InverseDynamicsSolver* IDsolver;
    IDsolver = new InverseDynamicsSolver(*model);
    Vector residualMobilityForces = IDsolver->solve(*state);

    // Residual forces in OpenSim order
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < NU; ++i)
        res[0][i] = (residualMobilityForces[indicesSimbodyInOS[i]]);
    return 0;
}

int main() {
    double x[NX];
    double u[NU];
    double tau[NR];

    for (int i = 1; i < NX; i += 2) x[i] = 0;
    x[0] = 0.01;
    x[2] = 0.02;
    x[4] = 0.03;
    x[6] = 0.04;
    x[8] = 0.05;
    x[10] = 0.06;
    x[12] = 0.07;
    x[14] = 0.08;
    x[16] = 0.09;
    x[18] = 0.1;
    x[20] = 0.11;
    x[22] = 0.12;
    x[24] = 0.13;
    x[26] = 0.14;
    x[28] = 0.15;
    x[30] = 0.16;
    x[32] = 0.17;

    for (int i = 0; i < NU; ++i) u[i] = 0;

    const double* Recorder_arg[n_in] = {x, u};
    double* Recorder_res[n_out] = {tau};

    F_generic(Recorder_arg, Recorder_res);

    for (int i = 0; i < NR; ++i) std::cout << Recorder_res[0][i] << std::endl;

    return 0;
}
