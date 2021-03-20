/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
//#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Common/SimmSpline.h>
//#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
//#include <OpenSim/Simulation/Model/MovingPathPoint.h>
//#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
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
    joint accelerations (controls u), contact forces as well as several
    parameters related to the contact models (parameters p), and returns the
    joint torques as well as several variables for use in the
    optimal control problems. F is templatized using type T. F(x,u,p)->(r).
*/

// Inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 3;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 31;        // # degrees of freedom (excluding locked)
constexpr int ndofr = ndof+2;   // # degrees of freedom (including locked)
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NP = 54;          // # parameters
constexpr int NR = ndof+6+6;    // # residual torques + # GRFs + # GRMs

template<typename T>
T value(const Recorder& e) { return e; };
template<>
double value(const Recorder& e) { return e.getValue(); };

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
				idxOSInSimbody[iy] = isv / 2;
				break;
			}
		}
	}
	return idxOSInSimbody;
}

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

template<typename T>
int F_generic(const T** arg, T** res) {
	OpenSim::Model* model;
	model = new OpenSim::Model();;

	// Definition of bodies
	OpenSim::Body* pelvis;
	pelvis = new OpenSim::Body("pelvis", 9.71433360917240484866, Vec3(-0.06827780017111793887, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.08479523605527072849, 0.07184499086005914636, 0.04775918450972933132, 0., 0., 0.));
	model->addBody(pelvis);

	OpenSim::Body* femur_r;
	femur_r = new OpenSim::Body("femur_r", 7.67231915023827859557, Vec3(0.00000000000000000000, -0.17046665777369882089, 0.00000000000000000000), Inertia(0.11105547289013888157, 0.02911162881586165305, 0.11711002817093063566, 0., 0., 0.));
	model->addBody(femur_r);

	OpenSim::Body* tibia_r;
	tibia_r = new OpenSim::Body("tibia_r", 3.05815503574821212496, Vec3(0.00000000000000000000, -0.18048889444253254921, 0.00000000000000000000), Inertia(0.03885270037340381871, 0.00393152325207062493, 0.03939232121192332015, 0., 0., 0.));
	model->addBody(tibia_r);

	OpenSim::Body* talus_r;
	talus_r = new OpenSim::Body("talus_r", 0.08248563818606102771, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00068896770091018207, 0.00068896770091018207, 0.00068896770091018207, 0., 0., 0.));
	model->addBody(talus_r);

	OpenSim::Body* calcn_r;
	calcn_r = new OpenSim::Body("calcn_r", 1.03107047732576262433, Vec3(0.09139243779317623995, 0.02741773133795287129, 0.00000000000000000000), Inertia(0.00096455478127425477, 0.00268697403354970941, 0.00282476757373174674, 0., 0., 0.));
	model->addBody(calcn_r);

	OpenSim::Body* toes_r;
	toes_r = new OpenSim::Body("toes_r", 0.17866389231100815449, Vec3(0.03162178347643897908, 0.00548354626759057443, -0.01599367661380584477), Inertia(0.00006889677009101820, 0.00013779354018203640, 0.00068896770091018207, 0., 0., 0.));
	model->addBody(toes_r);

	OpenSim::Body* femur_l;
	femur_l = new OpenSim::Body("femur_l", 7.67231915023827859557, Vec3(0.00000000000000000000, -0.17046665777369882089, 0.00000000000000000000), Inertia(0.11105547289013888157, 0.02911162881586165305, 0.11711002817093063566, 0., 0., 0.));
	model->addBody(femur_l);

	OpenSim::Body* tibia_l;
	tibia_l = new OpenSim::Body("tibia_l", 3.05815503574821212496, Vec3(0.00000000000000000000, -0.18048889444253254921, 0.00000000000000000000), Inertia(0.03885270037340381871, 0.00393152325207062493, 0.03939232121192332015, 0., 0., 0.));
	model->addBody(tibia_l);

	OpenSim::Body* talus_l;
	talus_l = new OpenSim::Body("talus_l", 0.08248563818606102771, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00068896770091018207, 0.00068896770091018207, 0.00068896770091018207, 0., 0., 0.));
	model->addBody(talus_l);

	OpenSim::Body* calcn_l;
	calcn_l = new OpenSim::Body("calcn_l", 1.03107047732576262433, Vec3(0.09139243779317623995, 0.02741773133795287129, 0.00000000000000000000), Inertia(0.00096455478127425477, 0.00268697403354970941, 0.00282476757373174674, 0., 0., 0.));
	model->addBody(calcn_l);

	OpenSim::Body* toes_l;
	toes_l = new OpenSim::Body("toes_l", 0.17866389231100815449, Vec3(0.03162178347643897908, 0.00548354626759057443, 0.01599367661380584477), Inertia(0.00006889677009101820, 0.00013779354018203640, 0.00068896770091018207, 0., 0., 0.));
	model->addBody(toes_l);

	OpenSim::Body* torso;
	torso = new OpenSim::Body("torso", 22.12809221362184430859, Vec3(-0.02676026676359908804, 0.30650513962530268053, 0.00000000000000000000), Inertia(1.21625073505346970038, 0.62317899649569097331, 1.18069942499527735791, 0., 0., 0.));
	model->addBody(torso);

	OpenSim::Body* humerus_r;
	humerus_r = new OpenSim::Body("humerus_r", 1.67652059613169024388, Vec3(0.00000000000000000000, -0.16903343413648269644, 0.00000000000000000000), Inertia(0.01040408074495897950, 0.00358908561442959542, 0.01167824532974677358, 0., 0., 0.));
	model->addBody(humerus_r);

	OpenSim::Body* ulna_r;
	ulna_r = new OpenSim::Body("ulna_r", 0.50110025198032071003, Vec3(0.00000000000000000000, -0.11842869140885815826, 0.00000000000000000000), Inertia(0.00235897301981347410, 0.00049218275700362165, 0.00255887248908193620, 0., 0., 0.));
	model->addBody(ulna_r);

	OpenSim::Body* radius_r;
	radius_r = new OpenSim::Body("radius_r", 0.50110025198032071003, Vec3(0.00000000000000000000, -0.11842869140885815826, 0.00000000000000000000), Inertia(0.00235897301981347410, 0.00049218275700362165, 0.00255887248908193620, 0., 0., 0.));
	model->addBody(radius_r);

	OpenSim::Body* hand_r;
	hand_r = new OpenSim::Body("hand_r", 0.37737179470122916847, Vec3(0.00000000000000000000, -0.06691061390986266511, 0.00000000000000000000), Inertia(0.00071039970751979049, 0.00043563748880417645, 0.00106719238573600803, 0., 0., 0.));
	model->addBody(hand_r);

	OpenSim::Body* humerus_l;
	humerus_l = new OpenSim::Body("humerus_l", 1.67652059613169024388, Vec3(0.00000000000000000000, -0.16903343413648269644, 0.00000000000000000000), Inertia(0.01040408074495897950, 0.00358908561442959542, 0.01167824532974677358, 0., 0., 0.));
	model->addBody(humerus_l);

	OpenSim::Body* ulna_l;
	ulna_l = new OpenSim::Body("ulna_l", 0.50110025198032071003, Vec3(0.00000000000000000000, -0.11842869140885815826, 0.00000000000000000000), Inertia(0.00235897301981347410, 0.00049218275700362165, 0.00255887248908193620, 0., 0., 0.));
	model->addBody(ulna_l);

	OpenSim::Body* radius_l;
	radius_l = new OpenSim::Body("radius_l", 0.50110025198032071003, Vec3(0.00000000000000000000, -0.11842869140885815826, 0.00000000000000000000), Inertia(0.00235897301981347410, 0.00049218275700362165, 0.00255887248908193620, 0., 0., 0.));
	model->addBody(radius_l);

	OpenSim::Body* hand_l;
	hand_l = new OpenSim::Body("hand_l", 0.37737179470122916847, Vec3(0.00000000000000000000, -0.06691061390986266511, 0.00000000000000000000), Inertia(0.00071039970751979049, 0.00043563748880417645, 0.00106719238573600803, 0., 0., 0.));
	model->addBody(hand_l);

	// Definition of joints
	SpatialTransform st_ground_pelvis;
	st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
	st_ground_pelvis[0].setFunction(new LinearFunction());
	st_ground_pelvis[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_list", 1, 1));
	st_ground_pelvis[1].setFunction(new LinearFunction());
	st_ground_pelvis[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_rotation", 1, 1));
	st_ground_pelvis[2].setFunction(new LinearFunction());
	st_ground_pelvis[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
	st_ground_pelvis[3].setFunction(new LinearFunction());
	st_ground_pelvis[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
	st_ground_pelvis[4].setFunction(new LinearFunction());
	st_ground_pelvis[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tz", 1, 1));
	st_ground_pelvis[5].setFunction(new LinearFunction());
	st_ground_pelvis[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* ground_pelvis;
	ground_pelvis = new OpenSim::CustomJoint("ground_pelvis", model->getGround(), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *pelvis, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ground_pelvis);

	SpatialTransform st_hip_r;
	st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
	st_hip_r[0].setFunction(new LinearFunction());
	st_hip_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
	st_hip_r[1].setFunction(new LinearFunction());
	st_hip_r[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
	st_hip_r[2].setFunction(new LinearFunction());
	st_hip_r[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_hip_r[3].setFunction(new Constant(0));
	st_hip_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_r[4].setFunction(new Constant(0));
	st_hip_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_hip_r[5].setFunction(new Constant(0));
	st_hip_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* hip_r;
	hip_r = new OpenSim::CustomJoint("hip_r", *pelvis, Vec3(-0.06827780017111793887, -0.06383539733113006986, 0.08233069400586875974), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_r);

	SpatialTransform st_knee_r;
	st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[0].setFunction(new LinearFunction());
	st_knee_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_r[1].setFunction(new Constant(0));
	st_knee_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_r[2].setFunction(new Constant(0));
	st_knee_r[2].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_r[3].setFunction(new Constant(0));
	st_knee_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_r[4].setFunction(new Constant(0));
	st_knee_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_r[5].setFunction(new Constant(0));
	st_knee_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* knee_r;
	knee_r = new OpenSim::CustomJoint("knee_r", *femur_r, Vec3(-0.00451221000000000001, -0.39690723999999999450, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_r);

	SpatialTransform st_ankle_r;
	st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
	st_ankle_r[0].setFunction(new LinearFunction());
	st_ankle_r[0].setAxis(Vec3(-0.10501354999999999718, -0.17402244999999999520, 0.97912631999999999444));
	st_ankle_r[1].setFunction(new Constant(0));
	st_ankle_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_r[2].setFunction(new Constant(0));
	st_ankle_r[2].setAxis(Vec3(0.97912631999999999444, -0.00000000000000000000, 0.10501354999999999718));
	st_ankle_r[3].setFunction(new Constant(0));
	st_ankle_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_ankle_r[4].setFunction(new Constant(0));
	st_ankle_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_r[5].setFunction(new Constant(0));
	st_ankle_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* ankle_r;
	ankle_r = new OpenSim::CustomJoint("ankle_r", *tibia_r, Vec3(0.00000000000000000000, -0.41569482919276373734, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_r);

	SpatialTransform st_subtalar_r;
	st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
	st_subtalar_r[0].setFunction(new LinearFunction());
	st_subtalar_r[0].setAxis(Vec3(0.78717961000000002958, 0.60474746000000001445, -0.12094949000000000672));
	st_subtalar_r[1].setFunction(new Constant(0));
	st_subtalar_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_r[2].setFunction(new Constant(0));
	st_subtalar_r[2].setAxis(Vec3(-0.12094949000000000672, 0.00000000000000000000, -0.78717961000000002958));
	st_subtalar_r[3].setFunction(new Constant(0));
	st_subtalar_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_subtalar_r[4].setFunction(new Constant(0));
	st_subtalar_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_r[5].setFunction(new Constant(0));
	st_subtalar_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* subtalar_r;
	subtalar_r = new OpenSim::CustomJoint("subtalar_r", *talus_r, Vec3(-0.04457209191173205215, -0.03833912765423743568, 0.00723828107321955808), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_r);

	OpenSim::PinJoint* mtp_r;
	mtp_r = new OpenSim::PinJoint("mtp_r", *calcn_r, Vec3(0.16340967877419909637, -0.00182784875586352474, 0.00098703832816630335), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

	SpatialTransform st_hip_l;
	st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
	st_hip_l[0].setFunction(new LinearFunction());
	st_hip_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
	st_hip_l[1].setFunction(new LinearFunction());
	st_hip_l[1].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
	st_hip_l[2].setFunction(new LinearFunction());
	st_hip_l[2].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_hip_l[3].setFunction(new Constant(0));
	st_hip_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_hip_l[4].setFunction(new Constant(0));
	st_hip_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_hip_l[5].setFunction(new Constant(0));
	st_hip_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* hip_l;
	hip_l = new OpenSim::CustomJoint("hip_l", *pelvis, Vec3(-0.06827780017111793887, -0.06383539733113006986, -0.08233069400586875974), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_l);

	SpatialTransform st_knee_l;
	st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
	st_knee_l[0].setFunction(new LinearFunction());
	st_knee_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_l[1].setFunction(new Constant(0));
	st_knee_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_l[2].setFunction(new Constant(0));
	st_knee_l[2].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_l[3].setFunction(new Constant(0));
	st_knee_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_l[4].setFunction(new Constant(0));
	st_knee_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_l[5].setFunction(new Constant(0));
	st_knee_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* knee_l;
	knee_l = new OpenSim::CustomJoint("knee_l", *femur_l, Vec3(-0.00451221000000000001, -0.39690723999999999450, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_l);

	SpatialTransform st_ankle_l;
	st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
	st_ankle_l[0].setFunction(new LinearFunction());
	st_ankle_l[0].setAxis(Vec3(0.10501354999999999718, 0.17402244999999999520, 0.97912631999999999444));
	st_ankle_l[1].setFunction(new Constant(0));
	st_ankle_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_l[2].setFunction(new Constant(0));
	st_ankle_l[2].setAxis(Vec3(0.97912631999999999444, 0.00000000000000000000, -0.10501354999999999718));
	st_ankle_l[3].setFunction(new Constant(0));
	st_ankle_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_ankle_l[4].setFunction(new Constant(0));
	st_ankle_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_ankle_l[5].setFunction(new Constant(0));
	st_ankle_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* ankle_l;
	ankle_l = new OpenSim::CustomJoint("ankle_l", *tibia_l, Vec3(0.00000000000000000000, -0.41569482919276373734, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_l);

	SpatialTransform st_subtalar_l;
	st_subtalar_l[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_l", 1, 1));
	st_subtalar_l[0].setFunction(new LinearFunction());
	st_subtalar_l[0].setAxis(Vec3(-0.78717961000000002958, -0.60474746000000001445, -0.12094949000000000672));
	st_subtalar_l[1].setFunction(new Constant(0));
	st_subtalar_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_l[2].setFunction(new Constant(0));
	st_subtalar_l[2].setAxis(Vec3(-0.12094949000000000672, 0.00000000000000000000, 0.78717961000000002958));
	st_subtalar_l[3].setFunction(new Constant(0));
	st_subtalar_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_subtalar_l[4].setFunction(new Constant(0));
	st_subtalar_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_subtalar_l[5].setFunction(new Constant(0));
	st_subtalar_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* subtalar_l;
	subtalar_l = new OpenSim::CustomJoint("subtalar_l", *talus_l, Vec3(-0.04457209191173205215, -0.03833912765423743568, -0.00723828107321955808), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_l);

	OpenSim::PinJoint* mtp_l;
	mtp_l = new OpenSim::PinJoint("mtp_l", *calcn_l, Vec3(0.16340967877419909637, -0.00182784875586352474, -0.00098703832816630335), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

	SpatialTransform st_back;
	st_back[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_extension", 1, 1));
	st_back[0].setFunction(new LinearFunction());
	st_back[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_back[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_bending", 1, 1));
	st_back[1].setFunction(new LinearFunction());
	st_back[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_back[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_rotation", 1, 1));
	st_back[2].setFunction(new LinearFunction());
	st_back[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_back[3].setFunction(new Constant(0));
	st_back[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_back[4].setFunction(new Constant(0));
	st_back[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_back[5].setFunction(new Constant(0));
	st_back[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* back;
	back = new OpenSim::CustomJoint("back", *pelvis, Vec3(-0.09724999260582144200, 0.07870778944761119833, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *torso, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_back);

	SpatialTransform st_acromial_r;
	st_acromial_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
	st_acromial_r[0].setFunction(new LinearFunction());
	st_acromial_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_acromial_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
	st_acromial_r[1].setFunction(new LinearFunction());
	st_acromial_r[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
	st_acromial_r[2].setFunction(new LinearFunction());
	st_acromial_r[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[3].setFunction(new Constant(0));
	st_acromial_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[4].setFunction(new Constant(0));
	st_acromial_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_acromial_r[5].setFunction(new Constant(0));
	st_acromial_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* acromial_r;
	acromial_r = new OpenSim::CustomJoint("acromial_r", *torso, Vec3(0.00281428805463850408, 0.35583331053374983588, 0.15164151166039482876), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_r);

	SpatialTransform st_elbow_r;
	st_elbow_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
	st_elbow_r[0].setFunction(new LinearFunction());
	st_elbow_r[0].setAxis(Vec3(0.22604695999999999123, 0.02226900000000000060, 0.97386183000000003940));
	st_elbow_r[1].setFunction(new Constant(0));
	st_elbow_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_r[2].setFunction(new Constant(0));
	st_elbow_r[2].setAxis(Vec3(0.97386183000000003940, 0.00000000000000000000, -0.22604695999999999123));
	st_elbow_r[3].setFunction(new Constant(0));
	st_elbow_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_elbow_r[4].setFunction(new Constant(0));
	st_elbow_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_r[5].setFunction(new Constant(0));
	st_elbow_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* elbow_r;
	elbow_r = new OpenSim::CustomJoint("elbow_r", *humerus_r, Vec3(0.01350606958146362106, -0.29415878403030548682, -0.00985930748890318301), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_r);

	SpatialTransform st_radioulnar_r;
	st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
	st_radioulnar_r[0].setFunction(new LinearFunction());
	st_radioulnar_r[0].setAxis(Vec3(0.05639803000000000177, 0.99840645999999999560, 0.00195199999999999996));
	st_radioulnar_r[1].setFunction(new Constant(0));
	st_radioulnar_r[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_r[2].setFunction(new Constant(0));
	st_radioulnar_r[2].setAxis(Vec3(0.00195199999999999996, 0.00000000000000000000, -0.05639803000000000177));
	st_radioulnar_r[3].setFunction(new Constant(0));
	st_radioulnar_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_r[4].setFunction(new Constant(0));
	st_radioulnar_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_r[5].setFunction(new Constant(0));
	st_radioulnar_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* radioulnar_r;
	radioulnar_r = new OpenSim::CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.00660999632530503197, -0.01278076738564628244, 0.02562933464440777728), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_r);

	OpenSim::WeldJoint* radius_hand_r;
	radius_hand_r = new OpenSim::WeldJoint("radius_hand_r", *radius_r, Vec3(-0.00864399251876146225, -0.23173898370094603294, 0.01337327932026185183), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

	SpatialTransform st_acromial_l;
	st_acromial_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
	st_acromial_l[0].setFunction(new LinearFunction());
	st_acromial_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_acromial_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
	st_acromial_l[1].setFunction(new LinearFunction());
	st_acromial_l[1].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
	st_acromial_l[2].setFunction(new LinearFunction());
	st_acromial_l[2].setAxis(Vec3(0.00000000000000000000, -1.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[3].setFunction(new Constant(0));
	st_acromial_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[4].setFunction(new Constant(0));
	st_acromial_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_acromial_l[5].setFunction(new Constant(0));
	st_acromial_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* acromial_l;
	acromial_l = new OpenSim::CustomJoint("acromial_l", *torso, Vec3(0.00281428805463850408, 0.35583331053374983588, -0.15164151166039482876), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_l);

	SpatialTransform st_elbow_l;
	st_elbow_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
	st_elbow_l[0].setFunction(new LinearFunction());
	st_elbow_l[0].setAxis(Vec3(-0.22604695999999999123, -0.02226900000000000060, 0.97386183000000003940));
	st_elbow_l[1].setFunction(new Constant(0));
	st_elbow_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_l[2].setFunction(new Constant(0));
	st_elbow_l[2].setAxis(Vec3(0.97386183000000003940, -0.00000000000000000000, 0.22604695999999999123));
	st_elbow_l[3].setFunction(new Constant(0));
	st_elbow_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_elbow_l[4].setFunction(new Constant(0));
	st_elbow_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_elbow_l[5].setFunction(new Constant(0));
	st_elbow_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* elbow_l;
	elbow_l = new OpenSim::CustomJoint("elbow_l", *humerus_l, Vec3(0.01350606958146362106, -0.29415878403030548682, 0.00985930748890318301), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_l);

	SpatialTransform st_radioulnar_l;
	st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
	st_radioulnar_l[0].setFunction(new LinearFunction());
	st_radioulnar_l[0].setAxis(Vec3(-0.05639803000000000177, -0.99840645999999999560, 0.00195199999999999996));
	st_radioulnar_l[1].setFunction(new Constant(0));
	st_radioulnar_l[1].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_l[2].setFunction(new Constant(0));
	st_radioulnar_l[2].setAxis(Vec3(0.00195199999999999996, -0.00000000000000000000, 0.05639803000000000177));
	st_radioulnar_l[3].setFunction(new Constant(0));
	st_radioulnar_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_l[4].setFunction(new Constant(0));
	st_radioulnar_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_radioulnar_l[5].setFunction(new Constant(0));
	st_radioulnar_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* radioulnar_l;
	radioulnar_l = new OpenSim::CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.00660999632530503197, -0.01278076738564628244, -0.02562933464440777728), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_l);

	OpenSim::WeldJoint* radius_hand_l;
	radius_hand_l = new OpenSim::WeldJoint("radius_hand_l", *radius_l, Vec3(-0.00864399251876146225, -0.23173898370094603294, -0.01337327932026185183), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

	model->addJoint(ground_pelvis);
	model->addJoint(hip_l);
	model->addJoint(hip_r);
	model->addJoint(knee_l);
	model->addJoint(knee_r);
	model->addJoint(ankle_l);
	model->addJoint(ankle_r);
	model->addJoint(subtalar_l);
	model->addJoint(subtalar_r);
	model->addJoint(mtp_l);
	model->addJoint(mtp_r);
	model->addJoint(back);
	model->addJoint(acromial_l);
	model->addJoint(acromial_r);
	model->addJoint(elbow_l);
	model->addJoint(elbow_r);
	model->addJoint(radioulnar_l);
	model->addJoint(radioulnar_r);
	model->addJoint(radius_hand_l);
	model->addJoint(radius_hand_r);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);
    std::vector<T> p(arg[2], arg[2] + NP);

    // States and controls
    T ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    T up[NP]; /// contact model parameters - parameters
    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// pro_sup dofs are locked so Qs and Qdots are hard coded
    QsUs[NX] = SimTK::Pi / 2;
    QsUs[NX+1] = 0;
    QsUs[NX+2] = SimTK::Pi / 2;
    QsUs[NX+3] = 0;
    /// Controls
    T ut[NU+2];
    for (int i = 0; i < NU; ++i) ut[i] = u[i];
    /// pro_sup dofs are locked so Qdotdots are hard coded (0)
    /// Need to have a temporary vector to add 0s to the vector before
    /// adjusting for the index difference between OpenSim and Simbody.
    ut[NU] = 0;
    ut[NU+1] = 0;
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
    for (int i = 0; i < ndofr; ++i) ua[i] = ut[indicesOSInSimbody[i]];
    /// Parameters
    for (int i = 0; i < NP; ++i) up[i] = p[i];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    // Compute residual forces
    /// appliedMobilityForces (# mobilities)
    Vector appliedMobilityForces(ndofr);
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
    /// Extract contact forces
    Vec3 AppliedPointForce_s1_l, AppliedPointForce_s2_l;
    Vec3 AppliedPointForce_s3_l, AppliedPointForce_s4_l;
    Vec3 AppliedPointForce_s5_l, AppliedPointForce_s6_l;
    Vec3 AppliedPointForce_s1_r, AppliedPointForce_s2_r;
    Vec3 AppliedPointForce_s3_r, AppliedPointForce_s4_r;
    Vec3 AppliedPointForce_s5_r, AppliedPointForce_s6_r;
    int nc = 3;
    for (int i = 0; i < nc; ++i) {
        AppliedPointForce_s1_l[i]   = up[i];
        AppliedPointForce_s2_l[i]   = up[i + nc];
        AppliedPointForce_s3_l[i]   = up[i + nc + nc];
        AppliedPointForce_s4_l[i]   = up[i + nc + nc + nc];
        AppliedPointForce_s5_l[i]   = up[i + nc + nc + nc + nc];
        AppliedPointForce_s6_l[i]   = up[i + nc + nc + nc + nc + nc];
        AppliedPointForce_s1_r[i]   = up[i + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s2_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s3_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s4_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s5_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
        AppliedPointForce_s6_r[i]   = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
    }
    /// Extract contact sphere locations
    Vec3 locSphere_s1_r, locSphere_s2_r;
    Vec3 locSphere_s3_r, locSphere_s4_r;
    Vec3 locSphere_s5_r, locSphere_s6_r;
    /// Vertical positions are fixed
    locSphere_s1_r[1] = -0.01;
    locSphere_s2_r[1] = -0.01;
    locSphere_s3_r[1] = -0.01;
    locSphere_s4_r[1] = -0.01;
    locSphere_s5_r[1] = -0.01;
    locSphere_s6_r[1] = -0.01;
    int count = 0;
    for (int i = 0; i < nc; i+=2) {
        locSphere_s1_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
        locSphere_s2_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1];
        locSphere_s3_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1];
        locSphere_s4_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1];
        locSphere_s5_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1];
        locSphere_s6_r[i]   = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1];
        ++count;
    }
    Vec3 locSphere_s1_l(locSphere_s1_r[0],locSphere_s1_r[1],-locSphere_s1_r[2]);
    Vec3 locSphere_s2_l(locSphere_s2_r[0],locSphere_s2_r[1],-locSphere_s2_r[2]);
    Vec3 locSphere_s3_l(locSphere_s3_r[0],locSphere_s3_r[1],-locSphere_s3_r[2]);
    Vec3 locSphere_s4_l(locSphere_s4_r[0],locSphere_s4_r[1],-locSphere_s4_r[2]);
    Vec3 locSphere_s5_l(locSphere_s5_r[0],locSphere_s5_r[1],-locSphere_s5_r[2]);
    Vec3 locSphere_s6_l(locSphere_s6_r[0],locSphere_s6_r[1],-locSphere_s6_r[2]);
    /// Extract radii
    osim_double_adouble radius_s1, radius_s2, radius_s3, radius_s4, radius_s5, radius_s6;
    radius_s1 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1];
    radius_s2 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 1];
    radius_s3 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 2];
    radius_s4 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 3];
    radius_s5 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 4];
    radius_s6 =  up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + nc-1 + 5];
    /// Compute contact point positions in body frames
    Vec3 normal = Vec3(0, 1, 0);
    /// sphere 1 left
    Vec3 pos_InGround_HC_s1_l = calcn_l->findStationLocationInGround(*state, locSphere_s1_l);
    Vec3 contactPointpos_InGround_HC_s1_l = pos_InGround_HC_s1_l - radius_s1*normal;
    Vec3 contactPointpos_InGround_HC_s1_l_adj = contactPointpos_InGround_HC_s1_l - 0.5*contactPointpos_InGround_HC_s1_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s1_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_l_adj, *calcn_l);
    /// sphere 2 left
    Vec3 pos_InGround_HC_s2_l = calcn_l->findStationLocationInGround(*state, locSphere_s2_l);
    Vec3 contactPointpos_InGround_HC_s2_l = pos_InGround_HC_s2_l - radius_s2*normal;
    Vec3 contactPointpos_InGround_HC_s2_l_adj = contactPointpos_InGround_HC_s2_l - 0.5*contactPointpos_InGround_HC_s2_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s2_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_l_adj, *calcn_l);
    /// sphere 3 left
    Vec3 pos_InGround_HC_s3_l = calcn_l->findStationLocationInGround(*state, locSphere_s3_l);
    Vec3 contactPointpos_InGround_HC_s3_l = pos_InGround_HC_s3_l - radius_s3*normal;
    Vec3 contactPointpos_InGround_HC_s3_l_adj = contactPointpos_InGround_HC_s3_l - 0.5*contactPointpos_InGround_HC_s3_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s3_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_l_adj, *calcn_l);
    /// sphere 4 left
    Vec3 pos_InGround_HC_s4_l = calcn_l->findStationLocationInGround(*state, locSphere_s4_l);
    Vec3 contactPointpos_InGround_HC_s4_l = pos_InGround_HC_s4_l - radius_s4*normal;
    Vec3 contactPointpos_InGround_HC_s4_l_adj = contactPointpos_InGround_HC_s4_l - 0.5*contactPointpos_InGround_HC_s4_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s4_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_l_adj, *calcn_l);
    /// sphere 5 left
    Vec3 pos_InGround_HC_s5_l = toes_l->findStationLocationInGround(*state, locSphere_s5_l);
    Vec3 contactPointpos_InGround_HC_s5_l = pos_InGround_HC_s5_l - radius_s5*normal;
    Vec3 contactPointpos_InGround_HC_s5_l_adj = contactPointpos_InGround_HC_s5_l - 0.5*contactPointpos_InGround_HC_s5_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s5_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_l_adj, *toes_l);
    /// sphere 6 left
    Vec3 pos_InGround_HC_s6_l = toes_l->findStationLocationInGround(*state, locSphere_s6_l);
    Vec3 contactPointpos_InGround_HC_s6_l = pos_InGround_HC_s6_l - radius_s6*normal;
    Vec3 contactPointpos_InGround_HC_s6_l_adj = contactPointpos_InGround_HC_s6_l - 0.5*contactPointpos_InGround_HC_s6_l[1]*normal;
    Vec3 contactPointPos_InBody_HC_s6_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_l_adj, *toes_l);
    /// sphere 1 right
    Vec3 pos_InGround_HC_s1_r = calcn_r->findStationLocationInGround(*state, locSphere_s1_r);
    Vec3 contactPointpos_InGround_HC_s1_r = pos_InGround_HC_s1_r - radius_s1*normal;
    Vec3 contactPointpos_InGround_HC_s1_r_adj = contactPointpos_InGround_HC_s1_r - 0.5*contactPointpos_InGround_HC_s1_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s1_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_r_adj, *calcn_r);
    /// sphere 2 right
    Vec3 pos_InGround_HC_s2_r = calcn_r->findStationLocationInGround(*state, locSphere_s2_r);
    Vec3 contactPointpos_InGround_HC_s2_r = pos_InGround_HC_s2_r - radius_s2*normal;
    Vec3 contactPointpos_InGround_HC_s2_r_adj = contactPointpos_InGround_HC_s2_r - 0.5*contactPointpos_InGround_HC_s2_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s2_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_r_adj, *calcn_r);
    /// sphere 3 right
    Vec3 pos_InGround_HC_s3_r = calcn_r->findStationLocationInGround(*state, locSphere_s3_r);
    Vec3 contactPointpos_InGround_HC_s3_r = pos_InGround_HC_s3_r - radius_s3*normal;
    Vec3 contactPointpos_InGround_HC_s3_r_adj = contactPointpos_InGround_HC_s3_r - 0.5*contactPointpos_InGround_HC_s3_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s3_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_r_adj, *calcn_r);
    /// sphere 4 right
    Vec3 pos_InGround_HC_s4_r = calcn_r->findStationLocationInGround(*state, locSphere_s4_r);
    Vec3 contactPointpos_InGround_HC_s4_r = pos_InGround_HC_s4_r - radius_s4*normal;
    Vec3 contactPointpos_InGround_HC_s4_r_adj = contactPointpos_InGround_HC_s4_r - 0.5*contactPointpos_InGround_HC_s4_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s4_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_r_adj, *calcn_r);
    /// sphere 5 right
    Vec3 pos_InGround_HC_s5_r = toes_r->findStationLocationInGround(*state, locSphere_s5_r);
    Vec3 contactPointpos_InGround_HC_s5_r = pos_InGround_HC_s5_r - radius_s5*normal;
    Vec3 contactPointpos_InGround_HC_s5_r_adj = contactPointpos_InGround_HC_s5_r - 0.5*contactPointpos_InGround_HC_s5_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s5_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_r_adj, *toes_r);
    /// sphere 6 right
    Vec3 pos_InGround_HC_s6_r = toes_r->findStationLocationInGround(*state, locSphere_s6_r);
    Vec3 contactPointpos_InGround_HC_s6_r = pos_InGround_HC_s6_r - radius_s6*normal;
    Vec3 contactPointpos_InGround_HC_s6_r_adj = contactPointpos_InGround_HC_s6_r - 0.5*contactPointpos_InGround_HC_s6_r[1]*normal;
    Vec3 contactPointPos_InBody_HC_s6_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_r_adj, *toes_r);
    /// Add contact forces to appliedBodyForces
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s1_l, AppliedPointForce_s1_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s2_l, AppliedPointForce_s2_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s3_l, AppliedPointForce_s3_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s4_l, AppliedPointForce_s4_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s5_l, AppliedPointForce_s5_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_l->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s6_l, AppliedPointForce_s6_l, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s1_r, AppliedPointForce_s1_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s2_r, AppliedPointForce_s2_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s3_r, AppliedPointForce_s3_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, calcn_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s4_r, AppliedPointForce_s4_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s5_r, AppliedPointForce_s5_r, appliedBodyForces);
    model->getMatterSubsystem().addInStationForce(*state, toes_r->getMobilizedBodyIndex(), contactPointPos_InBody_HC_s6_r, AppliedPointForce_s6_r, appliedBodyForces);
    /// knownUdot
    Vector knownUdot(ndofr);
    knownUdot.setToZero();
    for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];
    ///  Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Compute contact torques about the ground frame origin
    /// Get transforms
    SimTK::Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_toes_l = toes_l->getMobilizedBody().getBodyTransform(*state);
    SimTK::Transform TR_GB_toes_r = toes_r->getMobilizedBody().getBodyTransform(*state);
    /// Calculate torques
    Vec3 AppliedPointTorque_s1_l, AppliedPointTorque_s2_l, AppliedPointTorque_s3_l, AppliedPointTorque_s4_l, AppliedPointTorque_s5_l, AppliedPointTorque_s6_l;
    Vec3 AppliedPointTorque_s1_r, AppliedPointTorque_s2_r, AppliedPointTorque_s3_r, AppliedPointTorque_s4_r, AppliedPointTorque_s5_r, AppliedPointTorque_s6_r;
    AppliedPointTorque_s1_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s1_l) % AppliedPointForce_s1_l;
    AppliedPointTorque_s2_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s2_l) % AppliedPointForce_s2_l;
    AppliedPointTorque_s3_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s3_l) % AppliedPointForce_s3_l;
    AppliedPointTorque_s4_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s4_l) % AppliedPointForce_s4_l;
    AppliedPointTorque_s5_l = (TR_GB_toes_l*contactPointPos_InBody_HC_s5_l) % AppliedPointForce_s5_l;
    AppliedPointTorque_s6_l = (TR_GB_toes_l*contactPointPos_InBody_HC_s6_l) % AppliedPointForce_s6_l;
    AppliedPointTorque_s1_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s1_r) % AppliedPointForce_s1_r;
    AppliedPointTorque_s2_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s2_r) % AppliedPointForce_s2_r;
    AppliedPointTorque_s3_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s3_r) % AppliedPointForce_s3_r;
    AppliedPointTorque_s4_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s4_r) % AppliedPointForce_s4_r;
    AppliedPointTorque_s5_r = (TR_GB_toes_r*contactPointPos_InBody_HC_s5_r) % AppliedPointForce_s5_r;
    AppliedPointTorque_s6_r = (TR_GB_toes_r*contactPointPos_InBody_HC_s6_r) % AppliedPointForce_s6_r;
    /// Contact torques
    Vec3 MOM_l, MOM_r;
    MOM_l = AppliedPointTorque_s1_l + AppliedPointTorque_s2_l + AppliedPointTorque_s3_l + AppliedPointTorque_s4_l + AppliedPointTorque_s5_l + AppliedPointTorque_s6_l;
    MOM_r = AppliedPointTorque_s1_r + AppliedPointTorque_s2_r + AppliedPointTorque_s3_r + AppliedPointTorque_s4_r + AppliedPointTorque_s5_r + AppliedPointTorque_s6_r;
    /// Contact forces
    Vec3 GRF_r = AppliedPointForce_s1_r + AppliedPointForce_s2_r + AppliedPointForce_s3_r + AppliedPointForce_s4_r + AppliedPointForce_s5_r + AppliedPointForce_s6_r;
    Vec3 GRF_l = AppliedPointForce_s1_l + AppliedPointForce_s2_l + AppliedPointForce_s3_l + AppliedPointForce_s4_l + AppliedPointForce_s5_l + AppliedPointForce_s6_l;

    // Residual forces in OpenSim order
    T res_os[ndofr];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < ndofr; ++i) res_os[i] =
            value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
    // Extract results
    /// Residual forces
    /// We do want to extract the pro_sup torques (last two -> till NU)
    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];
    /// Contact forces
    for (int i = 0; i < nc; ++i) {
        res[0][i + NU] = value<T>(GRF_r[i]);      /// GRF_r
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i + NU + nc] = value<T>(GRF_l[i]); /// GRF_l
    }
    /// Contact torques
    for (int i = 0; i < nc; ++i) {
        res[0][i + NU + nc + nc] = value<T>(MOM_r[i]);        /// GRM_r
    }
    for (int i = 0; i < nc; ++i) {
        res[0][i + NU + nc + nc + nc] = value<T>(MOM_l[i]);   /// GRM_l
    }
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
    Recorder p[NP];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    for (int i = 0; i < NU; ++i) u[i] <<= 0;
    for (int i = 0; i < NP; ++i) p[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x,u,p };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}
