#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
#include <OpenSim/Simulation/Model/MovingPathPoint.h>
#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
#include "SimTKcommon/internal/recorder.h"

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>
#include <chrono>

using namespace SimTK;
using namespace OpenSim;

// Declare inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 2;    // degrees of freedom
constexpr int NX = 2 * ndof;  // states
constexpr int NU = ndof;    // controls
constexpr int NR = ndof;      // number of returns: residual forces (ndof), maybe add 2 rotations, 4 velocities later. All in one output
							/// Value


// Helper function value
template<typename T>
T value(const Recorder& e) { return e; }
template<>
double value(const Recorder& e) { return e.getValue(); }

/* Function F, using templated type T
F(x,u) -> (Tau)
*/

template<typename T>
int F_generic(const T** arg, T** res) {
	/// bodies
	OpenSim::Model* model;
	OpenSim::Body* brach;
	OpenSim::Body* abrach;

	/// joints
	OpenSim::PinJoint* groundToBrach;
	OpenSim::PinJoint* brachToAbrach;

	/// contact elements
	/// none required for this model

	/// states
	SimTK::State* state;

	// OpenSim model
	/// Model
	model = new OpenSim::Model();
	/// Bodies - Definition
	osim_double_adouble brachMass = 1;
	osim_double_adouble abrachMass = 1;

	brach = new OpenSim::Body("Brachium", brachMass, Vec3(0), Inertia(1, 0.1, 0.1, 0, 0, 0));
	abrach = new OpenSim::Body("Antebrachium", abrachMass, Vec3(0), Inertia(1, 0.1, 0.1, 0, 0, 0));

	/// Joints - Transforms
	// none required, all defined as pin joints

	/// Joints - Definition
	osim_double_adouble brachLength = 0.5;
	osim_double_adouble abrachLength = 0.5;

	groundToBrach = new PinJoint("groundToBrach", model->getGround(), Vec3(0, 1.5, 0), Vec3(0), *brach, Vec3(0, -brachLength / 2, 0), Vec3(0));
	brachToAbrach = new PinJoint("brachToAbrach", *brach, Vec3(0, brachLength/2, 0), Vec3(0), *abrach, Vec3(0, -abrachLength / 2, 0), Vec3(0));

	/// bodies and joints
	model->addBody(brach);		    model->addJoint(groundToBrach);
	model->addBody(abrach);		model->addJoint(brachToAbrach);

	/// Contact elements
	/// none required


	//// Don't touch this section -- Begin
	/// Initialize  system and state.
	state = new State(model->initSystem());

	// Read inputs
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);

	// States and controls
	T ua[NU]; /// joint accelerations (Qdotdots), i.e. controls
	Vector QsUs(NX); /// joint positions (Qs) and velocities (Us), i.e. states

					 // Assign inputs to model variables
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];    // states
	for (int i = 0; i < NU; ++i) ua[i] = u[i];      // controls

	model->setStateVariableValues(*state, QsUs);
	model->realizeVelocity(*state);

	// Residual forces
	/// appliedMobilityForces (# mobilities)
	Vector appliedMobilityForces(ndof);
	appliedMobilityForces.setToZero();
	/// appliedBodyForces (# bodies + ground)
	Vector_<SpatialVec> appliedBodyForces;
	int nbodies = model->getBodySet().getSize() + 1; // including ground
	appliedBodyForces.resize(nbodies);
	appliedBodyForces.setToZero();
	/// Gravity
	Vec3 gravity(0);
	gravity[1] = -9.80665;    // Y direction is up
						   /// Weight
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
	//// Don't touch this section -- End


	// Residual forces
	Vector residualMobilityForces(ndof);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
		*state, appliedMobilityForces, appliedBodyForces, knownUdot,
		residualMobilityForces);

	// Extract Joint rotations and velocities
	/* problem below... leave for now, may not be necessary. if add back in, need to update NR above
	osim_double_adouble brach_pos = groundToBrach->getCoordinate();
	osim_double_adouble abrach_pos = brachToAbrach->getCoordinate();

	Vec3 brach_vel = brach->getAngularVelocityInGround(*state);
	Vec3 abrach_vel = abrach->getAngularVelocityInGround(*state);
	*/

	// CasADi may not always request all outputs
	// if res[i] is a null pointer, this means that output i is not required
	if (res[0]) {
		for (int i = 0; i < NU; ++i) {
			res[0][i] = value<T>(residualMobilityForces[i]); // residual torques
		}
		// continue adding as above to get other return variables
		/* problem below... leave for now, may not be necessary. if add back in, need to update NR above
		res[0][ndof + 1] = value<T>(brach_pos);
		res[0][ndof + 2] = value<T>(abrach_pos);
		for (int i = 0; i < 2; ++i) {
			res[0][i + ndof + 2] = value<T>(brach_vel[i]); // brachium velocity (x and y)
		}
		for (int i = 0; i < 2; ++i) {
			res[0][i + ndof + 4] = value<T>(abrach_vel[i]); // antebrachium velocity (x and y)
		}
		*/
	}
	return 0;
}

int main() {

	Recorder x[NX];
	Recorder u[NU];
	Recorder tau[NR];

	for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;

	const Recorder* Recorder_arg[n_in] = { x,u };
	Recorder* Recorder_res[n_out] = { tau };

	F_generic<Recorder>(Recorder_arg, Recorder_res);

	double res[NR];
	for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

	Recorder::stop_recording();

	return 0;
}
