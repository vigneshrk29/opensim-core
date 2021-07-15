/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/MultivariatePolynomialFunction.h>
#include <OpenSim/Common/PolynomialFunction.h>
#include <OpenSim/Common/MultiplierFunction.h>
#include "SimTKcommon/internal/recorder.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>

using namespace SimTK;
using namespace OpenSim;

/*  The function F describes the OpenSim model and, implicitly, the skeleton
    dynamics. F takes as inputs joint positions and velocities (states x),
    joint accelerations (controls u), and returns the joint torques as well as
    several variables for use in the optimal control problems. F is templatized
    using type T. F(x,u)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 10;        // # degrees of freedom
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NR = ndof;    // # residual torques + # joint origins

// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

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
    SimTK::Array_<int> idxOSInSimbody(s.getNQ());
    s.updQ() = 0;
    for (int iy = 0; iy < s.getNQ(); ++iy) {
        s.updQ()[iy] = SimTK::NaN;
        const auto svValues = model.getStateVariableValues(s);
        for (int isv = 0; isv < svNames.size(); ++isv) {
            if (SimTK::isNaN(svValues[isv])) {
                s.updQ()[iy] = 0;
                idxOSInSimbody[iy] = isv/2;
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
    SimTK::Array_<int> idxSimbodyInOS(s.getNQ());
	for (int iy = 0; iy < s.getNQ(); ++iy) {
		for (int iyy = 0; iyy < s.getNQ(); ++iyy) {
			if (idxOSInSimbody[iyy] == iy) {
				idxSimbodyInOS[iy] = iyy;
				break;
			}
		}
	}
    return idxSimbodyInOS;
}

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

    // OpenSim model: create components
    /// Model
    OpenSim::Model* model;
    /// Bodies
    OpenSim::Body* pelvis;
    OpenSim::Body* femur_r;
    OpenSim::Body* tibia_r;
    /// Joints
    OpenSim::CustomJoint* ground_pelvis;
    OpenSim::CustomJoint* hip_r;
    OpenSim::CustomJoint* knee_r;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
    pelvis = new OpenSim::Body("pelvis", 11.751210011095651, Vec3(-0.069729482228687481, 0, 0), Inertia(0.099778065737821386, 0.08453958682650041, 0.056197957258948036, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 9.2810312301269491, Vec3(0, -0.17281329846170712, 0), Inertia(0.13806543520006995, 0.036191910198076598, 0.14559252763442779, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 3.6993810916309018, Vec3(0, -0.20693800230674214, 0), Inertia(0.061783188539925718, 0.0062518702689210552, 0.062641288380758026, 0, 0, 0));
    /// Joints
    /// Ground-Pelvis transform
    SpatialTransform st_ground_pelvis;
    st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
    st_ground_pelvis[0].setFunction(new LinearFunction());
    st_ground_pelvis[0].setAxis(Vec3(0, 0, 1));
    st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_list", 1, 1));
    st_ground_pelvis[1].setFunction(new LinearFunction());
    st_ground_pelvis[1].setAxis(Vec3(1, 0, 0));
    st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_rotation", 1, 1));
    st_ground_pelvis[2].setFunction(new LinearFunction());
    st_ground_pelvis[2].setAxis(Vec3(0, 1, 0));
    st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
    st_ground_pelvis[3].setFunction(new LinearFunction());
    st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
    st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
    st_ground_pelvis[4].setFunction(new LinearFunction());
    st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));
    st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tz", 1, 1));
    st_ground_pelvis[5].setFunction(new LinearFunction());
    st_ground_pelvis[5].setAxis(Vec3(0, 0, 1));
    /// Hip_r transform
    SpatialTransform st_hip_r;
    st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
    st_hip_r[0].setFunction(new LinearFunction());
    st_hip_r[0].setAxis(Vec3(0, 0, 1));
    st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
    st_hip_r[1].setFunction(new LinearFunction());
    st_hip_r[1].setAxis(Vec3(1, 0, 0));
    st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
    st_hip_r[2].setFunction(new LinearFunction());
    st_hip_r[2].setAxis(Vec3(0, 1, 0));
    /// Knee_r transform
    SpatialTransform st_knee_r;
    st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[0].setFunction(new LinearFunction());
    st_knee_r[0].setAxis(Vec3(0,0,1));
    st_knee_r[1].setFunction(new Constant(0));
    st_knee_r[1].setAxis(Vec3(0,1,0));
    st_knee_r[2].setFunction(new Constant(0));
    st_knee_r[2].setAxis(Vec3(1,0,0));

    /*st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    osim_double_adouble knee_X_l_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
    osim_double_adouble knee_X_l_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
    OpenSim::SimmSpline* knee_X_l = new SimmSpline(12, knee_X_l_x, knee_X_l_y, "function_X");
    st_knee_r[3].setFunction(new MultiplierFunction(knee_X_l, 1.0165488144806301));
    st_knee_r[3].setAxis(Vec3(1,0,0));
    st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    osim_double_adouble knee_Y_l_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
    osim_double_adouble knee_Y_l_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
    OpenSim::SimmSpline* knee_Y_l = new SimmSpline(7, knee_Y_l_x, knee_Y_l_y, "function_Y");
    st_knee_r[4].setFunction(new MultiplierFunction(knee_Y_l, 1.0165488144806301));
    st_knee_r[4].setAxis(Vec3(0,1,0));*/
	
	// Example with MultivariatePolynomialFunction
	//st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	//osim_double_adouble coefficientsTr1[5] = { -0.0036, -0.0063, 0.0028, 0.0013, -0.0007 };
	//auto* tr1 = new MultivariatePolynomialFunction();
	//tr1->setDimension(1);
	//tr1->setOrder(4);
	//Vector coefficientsTr1_vec(5);
	//for (int i = 0; i < 5; ++i) coefficientsTr1_vec[i] = coefficientsTr1[i];
	//tr1->setCoefficients(coefficientsTr1_vec);
	//st_knee_r[3].setFunction(new MultiplierFunction(tr1, 1.0165488144806301));
	//st_knee_r[3].setAxis(Vec3(1, 0, 0));

	//// Example with PolynomialFunction
	//st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	//osim_double_adouble coefficientsTr1[5] = {-0.0007, 0.0013, 0.0028, -0.0063 -0.0036 };
	//auto* tr1 = new PolynomialFunction();
	//Vector coefficientsTr1_vec(5);
	//for (int i = 0; i < 5; ++i) coefficientsTr1_vec[i] = coefficientsTr1[i];
	//tr1->setCoefficients(coefficientsTr1_vec);
	//st_knee_r[3].setFunction(new MultiplierFunction(tr1, 1.0165488144806301));
	//st_knee_r[3].setAxis(Vec3(1, 0, 0));

	// Example with PolynomialFunction
	st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	osim_double_adouble coefficientsTr1[5] = { -0.0007, 0.0013, 0.0028, -0.0063 - 0.0036 };
	//auto* tr1 = new PolynomialFunction();
	Vector coefficientsTr1_vec(5);
	for (int i = 0; i < 5; ++i) coefficientsTr1_vec[i] = coefficientsTr1[i];
	//tr1->setCoefficients(coefficientsTr1_vec);
	st_knee_r[3].setFunction(new MultiplierFunction(new PolynomialFunction(coefficientsTr1_vec), 1.0165488144806301));
	st_knee_r[3].setAxis(Vec3(1, 0, 0));
	
	st_knee_r[4].setFunction(new MultiplierFunction(new Constant(0), 1));
	st_knee_r[4].setAxis(Vec3(0, 1, 0));
	
	st_knee_r[5].setFunction(new Constant(0));
    st_knee_r[5].setAxis(Vec3(0,0,1));

    /// Joint specifications
    ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.069729482228687481, -0.065192627656523949, 0.08235377321209908), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(0, 0, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
    /// Add bodies and joints to model
    model->addBody(pelvis);		    model->addJoint(ground_pelvis);
    model->addBody(femur_r);		model->addJoint(hip_r);
    model->addBody(tibia_r);		model->addJoint(knee_r);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);

    // States and controls
    T ua[NU]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    T ut[NU];
    for (int i = 0; i < NU; ++i) ut[i] = u[i];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
    for (int i = 0; i < ndof; ++i) ua[i] = ut[indicesOSInSimbody[i]];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    // Compute residual forces
    /// appliedMobilityForces (# mobilities)
    Vector appliedMobilityForces(ndof);
    appliedMobilityForces.setToZero();
    /// appliedBodyForces (# bodies + ground)
    Vector_<SpatialVec> appliedBodyForces;
    int nbodies = model->getBodySet().getSize() + 1;
    appliedBodyForces.resize(nbodies);
    appliedBodyForces.setToZero();
    /// Set gravity
    Vec3 gravity(0);
    gravity[1] = -9.80665;
    /// Add weights to appliedBodyForces
    for (int i = 0; i < model->getBodySet().getSize(); ++i) {
        model->getMatterSubsystem().addInStationForce(*state,
            model->getBodySet().get(i).getMobilizedBodyIndex(),
            model->getBodySet().get(i).getMassCenter(),
            model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
    }
    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
        *state, appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);
    // Residual forces in OpenSim order
    T res_os[ndof];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < ndof; ++i) res_os[i] =
            value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
    // Extract results
    int nc = 3;
    /// Residual forces
    /// We do want to extract the pro_sup torques (last two -> until NU)
    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];

    return 0;

}

/* In main(), the Recorder is used to save the expression graph of function F.
This expression graph is saved as a MATLAB function named foo.m. From this
function, a c-code can be generated via CasADi and then compiled as a dll. This
dll is then imported in MATLAB as an external function. With this workflow,
CasADi can use algorithmic differentiation to differentiate the function F.
*/
int main() {

    Recorder x[NX];
    Recorder u[NU];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 1;
    for (int i = 0; i < NU; ++i) u[i] <<= 1;

    const Recorder* Recorder_arg[n_in] = { x,u };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) {
        Recorder_res[0][i] >>= res[i];
        std::cout << res[i] << std::endl;
    }

    Recorder::stop_recording();

    return 0;

}
