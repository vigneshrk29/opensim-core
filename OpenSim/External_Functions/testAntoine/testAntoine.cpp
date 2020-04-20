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
constexpr int ndof = 18;    // degrees of freedom
constexpr int NX = 2 * ndof;  // states
constexpr int NU = ndof;    // controls
constexpr int NR = 42;      // residual forces + GRFs + GRMs + 3 + 3 + 6 (= ndof + 6*nFeet + 3 for COM position
                            // + 3 for COM velocity + 2*3 points for OCP penetration constraint). All in one output
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
	OpenSim::Body* Body;
	OpenSim::Body* Left_thigh;
	OpenSim::Body* Left_crus;
	OpenSim::Body* Left_pes;
	OpenSim::Body* Left_digits;
	OpenSim::Body* Right_thigh;
	OpenSim::Body* Right_crus;
	OpenSim::Body* Right_pes;
	OpenSim::Body* Right_digits;
	/// joints
	OpenSim::CustomJoint* ground_pelvis;
	OpenSim::CustomJoint* L_hip;
	OpenSim::CustomJoint* L_knee;
	OpenSim::CustomJoint* L_ankle;
	OpenSim::CustomJoint* L_mtp;
	OpenSim::CustomJoint* R_hip;
	OpenSim::CustomJoint* R_knee;
	OpenSim::CustomJoint* R_ankle;
	OpenSim::CustomJoint* R_mtp;
	/// contact elements
	OpenSim::HuntCrossleyForce_smooth* L_pes_contact;
	OpenSim::HuntCrossleyForce_smooth* R_pes_contact;
	/// states
	SimTK::State* state;

	// OpenSim model
	/// Model
	model = new OpenSim::Model();
	/// Bodies - Definition
	Body = new OpenSim::Body("Body", 0.42506325, Vec3(-0.042913, 0, -0.018888), Inertia(0.0002253409375, 0.00084157625, 0.0008987103125, 2.52399375e-005, 2.659759375e-005, 6.413775e-006));
	Left_thigh = new OpenSim::Body("Left_thigh", 0.034807125, Vec3(-0.021695, -0.0054673, 0.0030261), Inertia(7.17028125e-006, 8.82265625e-006, 1.47005e-005, 1.8805e-006, -1.225365625e-006, 2.1577625e-007));
	Left_crus = new OpenSim::Body("Left_crus", 0.019371625, Vec3(-0.025672, -0.00047876, -0.00067702), Inertia(9.51618125e-007, 6.548009375e-006, 6.82638125e-006, 3.963553125e-007, -1.807321875e-007, 9.42809375e-008));
	Left_pes = new OpenSim::Body("Left_pes", 0.00388005, Vec3(-0.023325, -0.0017032, 0.00124), Inertia(4.3570625e-008, 9.937528125e-007, 1.0233653125e-006, 1.0086875e-008, -1.9331875e-008, 2.040471875e-009));
	Left_digits = new OpenSim::Body("Left_digits", 0.002045474375, Vec3(0.015011, 0.0025791, -0.00024952), Inertia(4.92443125e-008, 1.49179375e-007, 1.329793125e-007, 6.8456875e-009, 3.114965625e-009, -1.312746875e-008));
	Right_thigh = new OpenSim::Body("Right_thigh", 0.034807125, Vec3(-0.021695, -0.0054673, -0.0030261), Inertia(7.17028125e-006, 8.82265625e-006, 1.47005e-005, 1.8805e-006, -1.225365625e-006, 2.1577625e-007));
	Right_crus = new OpenSim::Body("Right_crus", 0.019371625, Vec3(-0.025672, -0.00047876, 0.00067702), Inertia(9.51618125e-007, 6.548009375e-006, 6.82638125e-006, 3.963553125e-007, -1.807321875e-007, 9.42809375e-008));
	Right_pes = new OpenSim::Body("Right_pes", 0.00388005, Vec3(-0.023325, -0.0017032, -0.00124), Inertia(4.3570625e-008, 9.937528125e-007, 1.0233653125e-006, 1.0086875e-008, -1.9331875e-008, 2.040471875e-009));
	Right_digits = new OpenSim::Body("Right_digits", 0.002045474375, Vec3(0.015011, 0.0025791, 0.00024952), Inertia(4.92443125e-008, 1.49179375e-007, 1.329793125e-007, 6.8456875e-009, 3.114965625e-009, -1.312746875e-008));
	/// Joints - Transforms
	// Ground-pelvis
	SpatialTransform st_ground_pelvis;
	st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("Pelvis_yaw", 1, 1));
	st_ground_pelvis[0].setFunction(new LinearFunction());
	st_ground_pelvis[0].setAxis(Vec3(0, 0, 1));
	st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("Pelvis_pitch", 1, 1));
	st_ground_pelvis[1].setFunction(new LinearFunction());
	st_ground_pelvis[1].setAxis(Vec3(0, 1, 0));
	st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("Pelvis_roll", 1, 1));
	st_ground_pelvis[2].setFunction(new LinearFunction());
	st_ground_pelvis[2].setAxis(Vec3(1, 0, 0));
	st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("Pelvis_tx", 1, 1));
	st_ground_pelvis[3].setFunction(new LinearFunction());
	st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
	st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("Pelvis_ty", 1, 1));
	st_ground_pelvis[4].setFunction(new LinearFunction());
	st_ground_pelvis[4].setAxis(Vec3(0, 0, 1));
	st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("Pelvis_tz", 1, 1));
	st_ground_pelvis[5].setFunction(new LinearFunction());
	st_ground_pelvis[5].setAxis(Vec3(0, -1, 0));
	// Hip_l
	SpatialTransform st_hip_l;
	st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("L_hip_extension", 1, 1));
	st_hip_l[0].setFunction(new LinearFunction());
	st_hip_l[0].setAxis(Vec3(0, 0, 1));
	st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("L_hip_abduction", 1, 1));
	st_hip_l[1].setFunction(new LinearFunction());
	st_hip_l[1].setAxis(Vec3(0, 1, 0));
	st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("L_hip_rotation", 1, 1));
	st_hip_l[2].setFunction(new LinearFunction());
	st_hip_l[2].setAxis(Vec3(1, 0, 0));
	// Hip_r
	SpatialTransform st_hip_r;
	st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("R_hip_extension", 1, 1));
	st_hip_r[0].setFunction(new LinearFunction());
	st_hip_r[0].setAxis(Vec3(0, 0, 1));
	st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("R_hip_abduction", 1, 1));
	st_hip_r[1].setFunction(new LinearFunction());
	st_hip_r[1].setAxis(Vec3(0, -1, 0));
	st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("R_hip_rotation", 1, 1));
	st_hip_r[2].setFunction(new LinearFunction());
	st_hip_r[2].setAxis(Vec3(-1, 0, 0));
	// Knee_l
	SpatialTransform st_knee_l;
	st_knee_l[2].setCoordinateNames(OpenSim::Array<std::string>("L_knee_extension", 1, 1));
	st_knee_l[2].setFunction(new LinearFunction());
	st_knee_l[2].setAxis(Vec3(0, 0, 1));
	// Knee_r
	SpatialTransform st_knee_r;
	st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("R_knee_extension", 1, 1));
	st_knee_r[2].setFunction(new LinearFunction());
	st_knee_r[2].setAxis(Vec3(0, 0, 1));
	// Ankle_l
	SpatialTransform st_ankle_l;
	st_ankle_l[2].setCoordinateNames(OpenSim::Array<std::string>("L_ankle_extension", 1, 1));
	st_ankle_l[2].setFunction(new LinearFunction());
	st_ankle_l[2].setAxis(Vec3(0, 0, 1));
	// Ankle_r
	SpatialTransform st_ankle_r;
	st_ankle_r[2].setCoordinateNames(OpenSim::Array<std::string>("R_ankle_extension", 1, 1));
	st_ankle_r[2].setFunction(new LinearFunction());
	st_ankle_r[2].setAxis(Vec3(0, 0, 1));
	// MTP_l
	SpatialTransform st_mtp_l;
	st_mtp_l[2].setCoordinateNames(OpenSim::Array<std::string>("L_mtp_angle", 1, 1));
	st_mtp_l[2].setFunction(new LinearFunction());
	st_mtp_l[2].setAxis(Vec3(0, 0, 1));
	// MTP_r
	SpatialTransform st_mtp_r;
	st_mtp_r[2].setCoordinateNames(OpenSim::Array<std::string>("R_mtp_angle", 1, 1));
	st_mtp_r[2].setFunction(new LinearFunction());
	st_mtp_r[2].setAxis(Vec3(0, 0, 1));
	/// Joints - Definition
	ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(-1.57079633, 0, 0), *Body, Vec3(0), Vec3(0), st_ground_pelvis);
	L_hip = new CustomJoint("hip_l", *Body, Vec3(3.2e-007, -0.018927, 4.56e-006), Vec3(1.5712, 2.46e-006, 1.5708), *Left_thigh, Vec3(0), Vec3(0), st_hip_l);
	R_hip = new CustomJoint("hip_r", *Body, Vec3(3.2e-007, 0.018927, 4.56e-006), Vec3(1.5712, 2.46e-006, 1.5708), *Right_thigh, Vec3(0), Vec3(0), st_hip_r);
	L_knee = new CustomJoint("knee_l", *Left_thigh, Vec3(-0.047345, 6.622e-007, 3.6597e-006), Vec3(3.1411, 0.28735, 4.59e-005), *Left_crus, Vec3(0.005, 0, 0), Vec3(0), st_knee_l);
	R_knee = new CustomJoint("knee_r", *Right_thigh, Vec3(-0.047345, 6.622e-007, 3.6597e-006), Vec3(3.1411, -0.28735, 4.59e-005), *Right_crus, Vec3(0.005, 0, 0), Vec3(0), st_knee_r);
	L_ankle = new CustomJoint("ankle_l", *Left_crus, Vec3(-0.069861, -1.6171e-006, 5.4602e-006), Vec3(-3.1415, 0.034637, 5.411e-005), *Left_pes, Vec3(0.004, 0, 0), Vec3(0), st_ankle_l);
	R_ankle = new CustomJoint("ankle_r", *Right_crus, Vec3(-0.069861, -1.6171e-006, 5.4602e-006), Vec3(-3.1415, -0.034637, 5.411e-005), *Right_pes, Vec3(0.004, 0, 0), Vec3(0), st_ankle_r);
	L_mtp = new CustomJoint("mtp_l", *Left_pes, Vec3(-0.043722, 2.0753e-006, 3.1971e-007), Vec3(0.00021799, -0.071172, 1.5708), *Left_digits, Vec3(0.0035, 0, 0), Vec3(0), st_mtp_l);
	R_mtp = new CustomJoint("mtp_r", *Right_pes, Vec3(-0.043722, 2.0753e-006, 3.1971e-007), Vec3(0.00021799, 0.071172, 1.5708), *Right_digits, Vec3(0.0035, 0, 0), Vec3(0), st_mtp_r);
	/// bodies and joints
	model->addBody(Body);		    model->addJoint(ground_pelvis);
	model->addBody(Left_thigh);		model->addJoint(L_hip);
	model->addBody(Right_thigh);	model->addJoint(R_hip);
	model->addBody(Left_crus);		model->addJoint(L_knee);
	model->addBody(Right_crus);		model->addJoint(R_knee);
	model->addBody(Left_pes);		model->addJoint(L_ankle);
	model->addBody(Right_pes);		model->addJoint(R_ankle);
	model->addBody(Left_digits);	model->addJoint(L_mtp);
	model->addBody(Right_digits);	model->addJoint(R_mtp);
	/// Contact elements
	/// Parameters
	osim_double_adouble radiusSphere_pes = 0.012;
	osim_double_adouble stiffness_pes = 250000;
	osim_double_adouble dissipation = 2;
	osim_double_adouble staticFriction = 0.8;
	osim_double_adouble dynamicFriction = 0.8;
	osim_double_adouble viscousFriction = 0.5;
	osim_double_adouble transitionVelocity = 0.2;
	Vec3 locSphere_L_pes_contact = Vec3(0.014, -0.006, -0.0008);
	Vec3 locSphere_R_pes_contact = Vec3(0.014, -0.006, 0.0008);

	Vec3 normal = Vec3(0, 1, 0);
	osim_double_adouble offset = 0;
	/// Define left contact bodies
	L_pes_contact = new HuntCrossleyForce_smooth("L_pes_1", "Left_digits", locSphere_L_pes_contact, radiusSphere_pes,
		stiffness_pes, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
	/// Define right contact bodies
	R_pes_contact = new HuntCrossleyForce_smooth("R_pes_1", "Right_digits", locSphere_R_pes_contact, radiusSphere_pes,
		stiffness_pes, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);

	///Add contact components to model
	model->addComponent(L_pes_contact);
	L_pes_contact->connectSocket_body_sphere(*Left_digits);
	model->addComponent(R_pes_contact);
	R_pes_contact->connectSocket_body_sphere(*Right_digits);

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
	Array<osim_double_adouble> Force_values_pes_1_l = L_pes_contact->getRecordValues(*state);
	SpatialVec GRF_pes_1_l;
	GRF_pes_1_l[0] = Vec3(Force_values_pes_1_l[9], Force_values_pes_1_l[10], Force_values_pes_1_l[11]);	 //GRMs (x,y,z) applied on the contact body by the plane (ground)
	GRF_pes_1_l[1] = Vec3(Force_values_pes_1_l[6], Force_values_pes_1_l[7], Force_values_pes_1_l[8]);    //GRFs (x,y,z) applied on the contact body by the plane (ground)
	int npes_l = model->getBodySet().get("Left_digits").getMobilizedBodyIndex();
	appliedBodyForces[npes_l] = appliedBodyForces[npes_l] + GRF_pes_1_l;

	/// right
	Array<osim_double_adouble> Force_values_pes_1_r = R_pes_contact->getRecordValues(*state);
	SpatialVec GRF_pes_1_r;
	GRF_pes_1_r[0] = Vec3(Force_values_pes_1_r[9], Force_values_pes_1_r[10], Force_values_pes_1_r[11]);
	GRF_pes_1_r[1] = Vec3(Force_values_pes_1_r[6], Force_values_pes_1_r[7], Force_values_pes_1_r[8]);
	int npes_r = model->getBodySet().get("Right_digits").getMobilizedBodyIndex();
	appliedBodyForces[npes_r] = appliedBodyForces[npes_r] + GRF_pes_1_r;

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

	// Ground reaction forces and moments
	SpatialVec GRF_l = GRF_pes_1_l;
	SpatialVec GRF_r = GRF_pes_1_r;

	// Extract position and velocity of whole-body COM
	Vec3 COM_pos = model->calcMassCenterPosition(*state);
	Vec3 COM_vel = model->calcMassCenterVelocity(*state);

	// Extract location of contact spheres (w.r.t. ground) to impose constraints in OCP
	Vec3 L_contact_inGround = Left_digits->findStationLocationInGround(*state, Vec3(0.014, -0.006, -0.0008));
	Vec3 R_contact_inGround = Right_digits->findStationLocationInGround(*state, Vec3(0.014, -0.006, 0.0008));

	// CasADi may not always request all outputs
	// if res[i] is a null pointer, this means that output i is not required
	if (res[0]) {
		for (int i = 0; i < NU; ++i) {
			res[0][i] = value<T>(residualMobilityForces[i]); // residual torques
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof] = value<T>(GRF_l[1][i]); // GRF_l (x, y and z)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 3] = value<T>(GRF_r[1][i]); // GRF_r (x, y and z)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 6] = value<T>(GRF_l[0][i]); // GRM_l (x, y and z)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 9] = value<T>(GRF_r[0][i]); // GRM_r (x, y and z)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 12] = value<T>(COM_pos[i]); // COM_pos (x, y and z)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 15] = value<T>(COM_vel[i]); // COM_vel (x, y and z)
		}
		res[0][ndof + 18] = value<T>(L_contact_inGround[0]);        // L_contact_location_x
		res[0][ndof + 18 + 1] = value<T>(L_contact_inGround[1]);    // L_contact_location_y
		res[0][ndof + 18 + 2] = value<T>(L_contact_inGround[2]);    // L_contact_location_z
		res[0][ndof + 18 + 3] = value<T>(R_contact_inGround[0]);    // R_contact_location_x
		res[0][ndof + 18 + 4] = value<T>(R_contact_inGround[1]);    // R_contact_location_y
		res[0][ndof + 18 + 5] = value<T>(R_contact_inGround[2]);    // R_contact_location_z
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
