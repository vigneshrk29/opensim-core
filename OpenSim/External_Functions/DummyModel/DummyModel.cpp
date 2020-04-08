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
#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
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
constexpr int ndof = 8;    // degrees of freedom
constexpr int NX = 2 * ndof;  // states
constexpr int NU = ndof;    // controls
constexpr int NR = 23;      // residual forces + GRFs + GRMs (= ndof + 6*nFeet). All in one output
							/// Value
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
	OpenSim::Body* Base;
	OpenSim::Body* Segment1;
	OpenSim::Body* Segment2;
	/// joints
	OpenSim::CustomJoint* ground_base;
	OpenSim::CustomJoint* Joint1;
	OpenSim::CustomJoint* Joint2;
	/// contact elements
	OpenSim::HuntCrossleyForce_smooth* Endpoint_contact;
	/// states
	SimTK::State* state;

	// OpenSim model
	/// Model
	model = new OpenSim::Model();
	/// Bodies - Definition
	Base = new OpenSim::Body("Base", 3.4465, Vec3(0, -0.062356, 0), Inertia(0.0070757, 0.091049, 0.096688, 0, 0, 0));
	Segment1 = new OpenSim::Body("Segment1", 3.6424, Vec3(0, -0.325, 0), Inertia(0.16228, 0.0037128, 0.16448, 0, 0, 0));
	Segment2 = new OpenSim::Body("Segment2", 3.6424, Vec3(0, -0.325, 0), Inertia(0.16228, 0.0037128, 0.16448, 0, 0, 0));
	/// Joints - Transforms
	// Ground-pelvis
	SpatialTransform st_ground_base;
	st_ground_base[0].setCoordinateNames(OpenSim::Array<std::string>("Base_rx", 1, 1));
	st_ground_base[0].setFunction(new LinearFunction());
	st_ground_base[0].setAxis(Vec3(1, 0, 0));
	st_ground_base[1].setCoordinateNames(OpenSim::Array<std::string>("Base_ry", 1, 1));
	st_ground_base[1].setFunction(new LinearFunction());
	st_ground_base[1].setAxis(Vec3(0, 1, 0));
	st_ground_base[2].setCoordinateNames(OpenSim::Array<std::string>("Base_rz", 1, 1));
	st_ground_base[2].setFunction(new LinearFunction());
	st_ground_base[2].setAxis(Vec3(0, 0, 1));
	st_ground_base[3].setCoordinateNames(OpenSim::Array<std::string>("Base_tx", 1, 1));
	st_ground_base[3].setFunction(new LinearFunction());
	st_ground_base[3].setAxis(Vec3(1, 0, 0));
	st_ground_base[4].setCoordinateNames(OpenSim::Array<std::string>("Base_ty", 1, 1));
	st_ground_base[4].setFunction(new LinearFunction());
	st_ground_base[4].setAxis(Vec3(0, 1, 0));
	st_ground_base[5].setCoordinateNames(OpenSim::Array<std::string>("Base_tz", 1, 1));
	st_ground_base[5].setFunction(new LinearFunction());
	st_ground_base[5].setAxis(Vec3(0, 0, 1));
	// Joint1
	SpatialTransform st_joint1;
	st_joint1[2].setCoordinateNames(OpenSim::Array<std::string>("joint_1", 1, 1));
	st_joint1[2].setFunction(new LinearFunction());
	st_joint1[2].setAxis(Vec3(0, 0, 1));
	// Joint2
	SpatialTransform st_joint2;
	st_joint2[2].setCoordinateNames(OpenSim::Array<std::string>("joint_2", 1, 1));
	st_joint2[2].setFunction(new LinearFunction());
	st_joint2[2].setAxis(Vec3(0, 0, 1));
	/// Joints - Definition
	ground_base = new CustomJoint("ground_base", model->getGround(), Vec3(0), Vec3(0), *Base, Vec3(-0.3, 0, 0), Vec3(0), st_ground_base);
	Joint1 = new CustomJoint("Joint1", *Base, Vec3(0, -0.15, 0), Vec3(0), *Segment1, Vec3(0), Vec3(0), st_joint1);
	Joint2 = new CustomJoint("Joint2", *Segment1, Vec3(0, -0.65, 0), Vec3(0), *Segment2, Vec3(0), Vec3(0), st_joint2);
	/// bodies and joints
	model->addBody(Base);		    model->addJoint(ground_base);
	model->addBody(Segment1);		model->addJoint(Joint1);
	model->addBody(Segment2);		model->addJoint(Joint2);
	/// Contact elements
	/// Parameters
	osim_double_adouble radiusSphere = 0.1;
	osim_double_adouble stiffness = 250000;
	osim_double_adouble dissipation = 2;
	osim_double_adouble staticFriction = 0.8;
	osim_double_adouble dynamicFriction = 0.8;
	osim_double_adouble viscousFriction = 0.5;
	osim_double_adouble transitionVelocity = 0.2;
	Vec3 locSphere_contact = Vec3(0, -0.65, 0);

	Vec3 normal = Vec3(0, 1, 0);
	osim_double_adouble offset = 0;
	/// Define left contact bodies
	Endpoint_contact = new HuntCrossleyForce_smooth("Endpoint", "Segment2", locSphere_contact, radiusSphere,
		stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);

	///Add contact components to model
	model->addComponent(Endpoint_contact);
	Endpoint_contact->connectSocket_body_sphere(*Segment2);

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
	/// Contact forces and moments
	/// left
	Array<osim_double_adouble> Force_endpoint = Endpoint_contact->getRecordValues(*state);
	SpatialVec GRF_endpoint;
	GRF_endpoint[0] = Vec3(Force_endpoint[9], Force_endpoint[10], Force_endpoint[11]);	 //GRMs (x,y,z) applied on the contact body by the plane (ground)
	GRF_endpoint[1] = Vec3(Force_endpoint[6], Force_endpoint[7], Force_endpoint[8]);    //GRFs (x,y,z) applied on the contact body by the plane (ground)
	int n_endpoint = model->getBodySet().get("Segment2").getMobilizedBodyIndex();
	appliedBodyForces[n_endpoint] = appliedBodyForces[n_endpoint] + GRF_endpoint;

	/// knownUdot
	Vector knownUdot(ndof);
	knownUdot.setToZero();
	for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];

	// Residual forces
	Vector residualMobilityForces(ndof);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
		*state, appliedMobilityForces, appliedBodyForces, knownUdot,
		residualMobilityForces);

    // Get point of application in the ground

    Vec3 contact_inGround = Segment2->findStationLocationInGround(*state, Vec3(0, -0.65, 0));

	// CasADi may not always request all outputs
	// if res[i] is a null pointer, this means that output i is not required
	if (res[0]) {
		for (int i = 0; i < NU; ++i) {
			res[0][i] = value<T>(residualMobilityForces[i]); // residual torques
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof] = value<T>(GRF_endpoint[1][i]); // GRF (x, y and z)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 3] = value<T>(GRF_endpoint[0][i]); // GRM (x, y and z)
		}
        for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 6] = value<T>(appliedBodyForces[n_endpoint][1][i]); // GRF (x, y and z)
		}
        for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 9] = value<T>(appliedBodyForces[n_endpoint][0][i]); // GRM (x, y and z)
		}
        for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 12] = value<T>(contact_inGround[i]); // contact location
		}
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
