#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
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

constexpr int n_in = 3;
constexpr int n_out = 1;
constexpr int nCoordinates = 33;
constexpr int NX = nCoordinates * 2;
constexpr int NU = nCoordinates;
constexpr int NP = 54;
constexpr int NR = nCoordinates + 29 * 3 + 3 * 4 + 3 * 2 * 12; ;

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
	pelvis = new OpenSim::Body("pelvis", 12.53462401183536201188, Vec3(-0.07781555017967077548, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.10941320781325256095, 0.09270321401297954012, 0.06162475420610236121, 0., 0., 0.));
	model->addBody(pelvis);

	OpenSim::Body* femur_r;
	femur_r = new OpenSim::Body("femur_r", 9.89976664546874829398, Vec3(0.00000000000000000000, -0.20242095410419783108, 0.00000000000000000000), Inertia(0.20205532747833265805, 0.05296595963024254561, 0.21307103988006403927, 0., 0., 0.));
	model->addBody(femur_r);

	OpenSim::Body* tibia_r;
	tibia_r = new OpenSim::Body("tibia_r", 3.94600649773962874889, Vec3(0.00000000000000000000, -0.22009009022360648267, 0.00000000000000000000), Inertia(0.07454517074013587707, 0.00754326132489470218, 0.07558052033374887402, 0., 0., 0.));
	model->addBody(tibia_r);

	OpenSim::Body* talus_r;
	talus_r = new OpenSim::Body("talus_r", 0.10643308153040131891, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00102354949748176402, 0.00102354949748176402, 0.00102354949748176402, 0., 0., 0.));
	model->addBody(talus_r);

	OpenSim::Body* calcn_r;
	calcn_r = new OpenSim::Body("calcn_r", 1.33041351913001637541, Vec3(0.09806546828421421058, 0.02941964048526425970, 0.00000000000000000000), Inertia(0.00143296929647446984, 0.00399184304017887945, 0.00419655293967523269, 0., 0., 0.));
	model->addBody(calcn_r);

	OpenSim::Body* toes_r;
	toes_r = new OpenSim::Body("toes_r", 0.23053405459484924145, Vec3(0.03393065202633811783, 0.00588392809705285281, -0.01716145694973748859), Inertia(0.00010235494974817643, 0.00020470989949635286, 0.00102354949748176402, 0., 0., 0.));
	model->addBody(toes_r);

	OpenSim::Body* femur_l;
	femur_l = new OpenSim::Body("femur_l", 9.89976664546874829398, Vec3(0.00000000000000000000, -0.20242095410419783108, 0.00000000000000000000), Inertia(0.20205532747833265805, 0.05296595963024254561, 0.21307103988006403927, 0., 0., 0.));
	model->addBody(femur_l);

	OpenSim::Body* tibia_l;
	tibia_l = new OpenSim::Body("tibia_l", 3.94600649773962874889, Vec3(0.00000000000000000000, -0.22009009022360648267, 0.00000000000000000000), Inertia(0.07454517074013587707, 0.00754326132489470218, 0.07558052033374887402, 0., 0., 0.));
	model->addBody(tibia_l);

	OpenSim::Body* talus_l;
	talus_l = new OpenSim::Body("talus_l", 0.10643308153040131891, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00102354949748176402, 0.00102354949748176402, 0.00102354949748176402, 0., 0., 0.));
	model->addBody(talus_l);

	OpenSim::Body* calcn_l;
	calcn_l = new OpenSim::Body("calcn_l", 1.33041351913001637541, Vec3(0.09806546828421421058, 0.02941964048526425970, 0.00000000000000000000), Inertia(0.00143296929647446984, 0.00399184304017887945, 0.00419655293967523269, 0., 0., 0.));
	model->addBody(calcn_l);

	OpenSim::Body* toes_l;
	toes_l = new OpenSim::Body("toes_l", 0.23053405459484924145, Vec3(0.03393065202633811783, 0.00588392809705285281, 0.01716145694973748859), Inertia(0.00010235494974817643, 0.00020470989949635286, 0.00102354949748176402, 0., 0., 0.));
	model->addBody(toes_l);

	OpenSim::Body* torso;
	torso = new OpenSim::Body("torso", 28.55237704983463942199, Vec3(-0.03029129157016192953, 0.32310711008172726677, 0.00000000000000000000), Inertia(1.56935578716576729086, 0.80410193096218196551, 1.52348312902616456199, 0., 0., 0.));
	model->addBody(torso);

	OpenSim::Body* humerus_r;
	humerus_r = new OpenSim::Body("humerus_r", 2.16325238210540682360, Vec3(0.00000000000000000000, -0.17217736840018482702, 0.00000000000000000000), Inertia(0.01392864639848861129, 0.00480495159954558507, 0.01563445668485968429, 0., 0., 0.));
	model->addBody(humerus_r);

	OpenSim::Body* ulna_r;
	ulna_r = new OpenSim::Body("ulna_r", 0.64658097029718808457, Vec3(0.00000000000000000000, -0.13420723621273464299, 0.00000000000000000000), Inertia(0.00390894244458061483, 0.00081557273151614451, 0.00424018638569801428, 0., 0., 0.));
	model->addBody(ulna_r);

	OpenSim::Body* radius_r;
	radius_r = new OpenSim::Body("radius_r", 0.64658097029718808457, Vec3(0.00000000000000000000, -0.13420723621273464299, 0.00000000000000000000), Inertia(0.00390894244458061483, 0.00081557273151614451, 0.00424018638569801428, 0., 0., 0.));
	model->addBody(radius_r);

	OpenSim::Body* hand_r;
	hand_r = new OpenSim::Body("hand_r", 0.48693134800158605069, Vec3(0.00000000000000000000, -0.07582527898698333824, 0.00000000000000000000), Inertia(0.00117716970309449998, 0.00072187424618014729, 0.00176839394859487674, 0., 0., 0.));
	model->addBody(hand_r);

	OpenSim::Body* humerus_l;
	humerus_l = new OpenSim::Body("humerus_l", 2.16325238210540682360, Vec3(0.00000000000000000000, -0.17217736840018482702, 0.00000000000000000000), Inertia(0.01392864639848861129, 0.00480495159954558507, 0.01563445668485968429, 0., 0., 0.));
	model->addBody(humerus_l);

	OpenSim::Body* ulna_l;
	ulna_l = new OpenSim::Body("ulna_l", 0.64658097029718808457, Vec3(0.00000000000000000000, -0.13420723621273464299, 0.00000000000000000000), Inertia(0.00390894244458061483, 0.00081557273151614451, 0.00424018638569801428, 0., 0., 0.));
	model->addBody(ulna_l);

	OpenSim::Body* radius_l;
	radius_l = new OpenSim::Body("radius_l", 0.64658097029718808457, Vec3(0.00000000000000000000, -0.13420723621273464299, 0.00000000000000000000), Inertia(0.00390894244458061483, 0.00081557273151614451, 0.00424018638569801428, 0., 0., 0.));
	model->addBody(radius_l);

	OpenSim::Body* hand_l;
	hand_l = new OpenSim::Body("hand_l", 0.48693134800158605069, Vec3(0.00000000000000000000, -0.07582527898698333824, 0.00000000000000000000), Inertia(0.00117716970309449998, 0.00072187424618014729, 0.00176839394859487674, 0., 0., 0.));
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
	hip_r = new OpenSim::CustomJoint("hip_r", *pelvis, Vec3(-0.07781555017967077548, -0.07275258651875868288, 0.10221907923123482731), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_r);

	SpatialTransform st_knee_r;
	st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[0].setFunction(new LinearFunction());
	st_knee_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_r[1].setFunction(new Constant(0));
	st_knee_r[1].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_r[2].setFunction(new Constant(0));
	st_knee_r[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_r[3].setFunction(new Constant(0));
	st_knee_r[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_r[4].setFunction(new Constant(0));
	st_knee_r[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_r[5].setFunction(new Constant(0));
	st_knee_r[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* knee_r;
	knee_r = new OpenSim::CustomJoint("knee_r", *femur_r, Vec3(-0.00535802999999999974, -0.47130824999999998415, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_r);

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
	ankle_r = new OpenSim::CustomJoint("ankle_r", *tibia_r, Vec3(0.00000000000000000000, -0.50690272520702073233, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_r);

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
	subtalar_r = new OpenSim::CustomJoint("subtalar_r", *talus_r, Vec3(-0.04782652888221126941, -0.04113846394522786137, 0.00776678508810976515), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_r);

	OpenSim::PinJoint* mtp_r;
	mtp_r = new OpenSim::PinJoint("mtp_r", *calcn_r, Vec3(0.17534105729217500103, -0.00196130936568428412, 0.00105910705746951341), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	hip_l = new OpenSim::CustomJoint("hip_l", *pelvis, Vec3(-0.07781555017967077548, -0.07275258651875868288, -0.10221907923123482731), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_l);

	SpatialTransform st_knee_l;
	st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
	st_knee_l[0].setFunction(new LinearFunction());
	st_knee_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_l[1].setFunction(new Constant(0));
	st_knee_l[1].setAxis(Vec3(-1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_l[2].setFunction(new Constant(0));
	st_knee_l[2].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_l[3].setFunction(new Constant(0));
	st_knee_l[3].setAxis(Vec3(1.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	st_knee_l[4].setFunction(new Constant(0));
	st_knee_l[4].setAxis(Vec3(0.00000000000000000000, 1.00000000000000000000, 0.00000000000000000000));
	st_knee_l[5].setFunction(new Constant(0));
	st_knee_l[5].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	OpenSim::CustomJoint* knee_l;
	knee_l = new OpenSim::CustomJoint("knee_l", *femur_l, Vec3(-0.00535802999999999974, -0.47130824999999998415, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_l);

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
	ankle_l = new OpenSim::CustomJoint("ankle_l", *tibia_l, Vec3(0.00000000000000000000, -0.50690272520702073233, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_l);

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
	subtalar_l = new OpenSim::CustomJoint("subtalar_l", *talus_l, Vec3(-0.04782652888221126941, -0.04113846394522786137, -0.00776678508810976515), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_l);

	OpenSim::PinJoint* mtp_l;
	mtp_l = new OpenSim::PinJoint("mtp_l", *calcn_l, Vec3(0.17534105729217500103, -0.00196130936568428412, -0.00105910705746951341), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	back = new OpenSim::CustomJoint("back", *pelvis, Vec3(-0.11083487840301058103, 0.08970250834007310881, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *torso, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_back);

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
	acromial_r = new OpenSim::CustomJoint("acromial_r", *torso, Vec3(0.00318563416346202919, 0.37510716061050519698, 0.19965066549714699518), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_r);

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
	elbow_r = new OpenSim::CustomJoint("elbow_r", *humerus_r, Vec3(0.01375727547538649342, -0.29962998494866993626, -0.01004268549804727614), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_r);

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
	radioulnar_r = new OpenSim::CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.00749066233564045763, -0.01448358034780368973, 0.02904399371198305643), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_r);

	OpenSim::WeldJoint* radius_hand_r;
	radius_hand_r = new OpenSim::WeldJoint("radius_hand_r", *radius_r, Vec3(-0.00979565282691082125, -0.26261413645009379358, 0.01515503409960853622), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	acromial_l = new OpenSim::CustomJoint("acromial_l", *torso, Vec3(0.00318563416346202919, 0.37510716061050519698, -0.19965066549714699518), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_l);

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
	elbow_l = new OpenSim::CustomJoint("elbow_l", *humerus_l, Vec3(0.01375727547538649342, -0.29962998494866993626, 0.01004268549804727614), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_l);

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
	radioulnar_l = new OpenSim::CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.00749066233564045763, -0.01448358034780368973, -0.02904399371198305643), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_l);

	OpenSim::WeldJoint* radius_hand_l;
	radius_hand_l = new OpenSim::WeldJoint("radius_hand_l", *radius_l, Vec3(-0.00979565282691082125, -0.26261413645009379358, -0.01515503409960853622), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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

	OpenSim::Station* C7_study;
	C7_study = new Station(*torso, Vec3(-0.08691489959863486769, 0.41085535577464793455, 0.00215944843199364356));
	model->addComponent(C7_study);
	OpenSim::Station* r_shoulder_study;
	r_shoulder_study = new Station(*torso, Vec3(0.00332986663766431101, 0.40784739177318485304, 0.16743157284992438161));
	model->addComponent(r_shoulder_study);
	OpenSim::Station* L_shoulder_study;
	L_shoulder_study = new Station(*torso, Vec3(0.00503539892362113439, 0.40637347082973906964, -0.16609247948493141567));
	model->addComponent(L_shoulder_study);
	OpenSim::Station* r_ASIS_study;
	r_ASIS_study = new Station(*pelvis, Vec3(0.00955096695431409404, 0.02530129165824734905, 0.16085986290878229177));
	model->addComponent(r_ASIS_study);
	OpenSim::Station* L_ASIS_study;
	L_ASIS_study = new Station(*pelvis, Vec3(0.02564119885680868349, 0.02401945067417710078, -0.15211044304060905574));
	model->addComponent(L_ASIS_study);
	OpenSim::Station* r_PSIS_study;
	r_PSIS_study = new Station(*pelvis, Vec3(-0.18082197035786604089, 0.00994948806250772755, 0.05486052708103439368));
	model->addComponent(r_PSIS_study);
	OpenSim::Station* L_PSIS_study;
	L_PSIS_study = new Station(*pelvis, Vec3(-0.17572816583532080426, 0.00761410674148277344, -0.06841303205075810467));
	model->addComponent(L_PSIS_study);
	OpenSim::Station* r_knee_study;
	r_knee_study = new Station(*femur_r, Vec3(0.03345480963930960727, -0.46931425407606564004, 0.06979186145956173259));
	model->addComponent(r_knee_study);
	OpenSim::Station* L_knee_study;
	L_knee_study = new Station(*femur_l, Vec3(0.02192737542091033331, -0.47380601276466693950, -0.06476139111937058435));
	model->addComponent(L_knee_study);
	OpenSim::Station* r_mknee_study;
	r_mknee_study = new Station(*femur_r, Vec3(0.04228724959699819985, -0.46735936272769362798, -0.06045050898334847211));
	model->addComponent(r_mknee_study);
	OpenSim::Station* L_mknee_study;
	L_mknee_study = new Station(*femur_l, Vec3(0.03043107103890826948, -0.47606823408657350516, 0.06052987875803789164));
	model->addComponent(L_mknee_study);
	OpenSim::Station* r_ankle_study;
	r_ankle_study = new Station(*tibia_r, Vec3(-0.02410389745808870499, -0.48499494856937380716, 0.05542652958730454049));
	model->addComponent(r_ankle_study);
	OpenSim::Station* L_ankle_study;
	L_ankle_study = new Station(*tibia_l, Vec3(-0.02819856035310782522, -0.48129688675141424348, -0.05752294834898633719));
	model->addComponent(L_ankle_study);
	OpenSim::Station* r_mankle_study;
	r_mankle_study = new Station(*tibia_r, Vec3(0.01120425688968612610, -0.48112593209006893069, -0.04203190925811695067));
	model->addComponent(r_mankle_study);
	OpenSim::Station* L_mankle_study;
	L_mankle_study = new Station(*tibia_l, Vec3(-0.00064421143949230236, -0.47280968984969473645, 0.04356052822436426442));
	model->addComponent(L_mankle_study);
	OpenSim::Station* r_calc_study;
	r_calc_study = new Station(*calcn_r, Vec3(-0.03390875624527270316, 0.05503903171326335331, -0.01235520996880051792));
	model->addComponent(r_calc_study);
	OpenSim::Station* L_calc_study;
	L_calc_study = new Station(*calcn_l, Vec3(-0.04542548013250082006, 0.03105539759700204350, 0.00988728924325359138));
	model->addComponent(L_calc_study);
	OpenSim::Station* r_toe_study;
	r_toe_study = new Station(*calcn_r, Vec3(0.16626928501547857953, -0.01236271270929323318, 0.00814005250893468091));
	model->addComponent(r_toe_study);
	OpenSim::Station* L_toe_study;
	L_toe_study = new Station(*calcn_l, Vec3(0.17199081736829568978, 0.02109135523138788715, -0.00750206767972391719));
	model->addComponent(L_toe_study);
	OpenSim::Station* r_5meta_study;
	r_5meta_study = new Station(*calcn_r, Vec3(0.12022638178659206254, -0.00704063595497178563, 0.06411494260185424121));
	model->addComponent(r_5meta_study);
	OpenSim::Station* L_5meta_study;
	L_5meta_study = new Station(*calcn_l, Vec3(0.12519991560678261910, 0.01489495250304473384, -0.06422578065261456970));
	model->addComponent(L_5meta_study);
	OpenSim::Station* r_lelbow_study;
	r_lelbow_study = new Station(*humerus_r, Vec3(0.01493164875981434214, -0.30588528623942234930, 0.04766515697291384690));
	model->addComponent(r_lelbow_study);
	OpenSim::Station* L_lelbow_study;
	L_lelbow_study = new Station(*humerus_l, Vec3(0.00911155840919679327, -0.30982118084203524866, -0.04150489742070326282));
	model->addComponent(L_lelbow_study);
	OpenSim::Station* r_melbow_study;
	r_melbow_study = new Station(*humerus_r, Vec3(0.00175195450064971614, -0.31112760282371176856, -0.05327054416199250575));
	model->addComponent(r_melbow_study);
	OpenSim::Station* L_melbow_study;
	L_melbow_study = new Station(*humerus_l, Vec3(0.00603720920044156784, -0.31084672542060975964, 0.05009729438481541619));
	model->addComponent(L_melbow_study);
	OpenSim::Station* r_lwrist_study;
	r_lwrist_study = new Station(*radius_r, Vec3(0.00242404636437115739, -0.25988529106923285994, 0.06133738847660302751));
	model->addComponent(r_lwrist_study);
	OpenSim::Station* L_lwrist_study;
	L_lwrist_study = new Station(*radius_l, Vec3(-0.00738392367045467846, -0.25550629365319510455, -0.05927623282973110141));
	model->addComponent(L_lwrist_study);
	OpenSim::Station* r_mwrist_study;
	r_mwrist_study = new Station(*radius_r, Vec3(-0.02542034525700959113, -0.26368048141040811849, -0.03194439714309565970));
	model->addComponent(r_mwrist_study);
	OpenSim::Station* L_mwrist_study;
	L_mwrist_study = new Station(*radius_l, Vec3(-0.01922375751635518187, -0.26585310860241695430, 0.02885616330348989012));
	model->addComponent(L_mwrist_study);

	// Initialize system.
	SimTK::State* state;
	state = new State(model->initSystem());

	// Read inputs.
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);
	std::vector<T> p(arg[2], arg[2] + NP);

	// States and controls.
	T ua[NU];
	Vector QsUs(NX);
	T up[NP]; /// contact model parameters - parameters
	/// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
	/// Controls
	/// OpenSim and Simbody have different state orders.
	auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
	for (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];
	/// Parameters
	for (int i = 0; i < NP; ++i) up[i] = p[i];

	// Set state variables and realize.
	model->setStateVariableValues(*state, QsUs);
	model->realizeVelocity(*state);

	// Compute residual forces.
	/// Set appliedMobilityForces (# mobilities).
	Vector appliedMobilityForces(nCoordinates);
	appliedMobilityForces.setToZero();
	/// Set appliedBodyForces (# bodies + ground).
	Vector_<SpatialVec> appliedBodyForces;
	int nbodies = model->getBodySet().getSize() + 1;
	appliedBodyForces.resize(nbodies);
	appliedBodyForces.setToZero();
	/// Set gravity.
	Vec3 gravity(0);
	gravity[1] = -9.80664999999999942304;
	/// Add weights to appliedBodyForces.
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
		AppliedPointForce_s1_l[i] = up[i];
		AppliedPointForce_s2_l[i] = up[i + nc];
		AppliedPointForce_s3_l[i] = up[i + nc + nc];
		AppliedPointForce_s4_l[i] = up[i + nc + nc + nc];
		AppliedPointForce_s5_l[i] = up[i + nc + nc + nc + nc];
		AppliedPointForce_s6_l[i] = up[i + nc + nc + nc + nc + nc];
		AppliedPointForce_s1_r[i] = up[i + nc + nc + nc + nc + nc + nc];
		AppliedPointForce_s2_r[i] = up[i + nc + nc + nc + nc + nc + nc + nc];
		AppliedPointForce_s3_r[i] = up[i + nc + nc + nc + nc + nc + nc + nc + nc];
		AppliedPointForce_s4_r[i] = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc];
		AppliedPointForce_s5_r[i] = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
		AppliedPointForce_s6_r[i] = up[i + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
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
	for (int i = 0; i < nc; i += 2) {
		locSphere_s1_r[i] = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc];
		locSphere_s2_r[i] = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1];
		locSphere_s3_r[i] = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1];
		locSphere_s4_r[i] = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1];
		locSphere_s5_r[i] = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1];
		locSphere_s6_r[i] = up[count + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1];
		++count;
	}
	Vec3 locSphere_s1_l(locSphere_s1_r[0], locSphere_s1_r[1], -locSphere_s1_r[2]);
	Vec3 locSphere_s2_l(locSphere_s2_r[0], locSphere_s2_r[1], -locSphere_s2_r[2]);
	Vec3 locSphere_s3_l(locSphere_s3_r[0], locSphere_s3_r[1], -locSphere_s3_r[2]);
	Vec3 locSphere_s4_l(locSphere_s4_r[0], locSphere_s4_r[1], -locSphere_s4_r[2]);
	Vec3 locSphere_s5_l(locSphere_s5_r[0], locSphere_s5_r[1], -locSphere_s5_r[2]);
	Vec3 locSphere_s6_l(locSphere_s6_r[0], locSphere_s6_r[1], -locSphere_s6_r[2]);
	/// Extract radii
	osim_double_adouble radius_s1, radius_s2, radius_s3, radius_s4, radius_s5, radius_s6;
	radius_s1 = up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1];
	radius_s2 = up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + 1];
	radius_s3 = up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + 2];
	radius_s4 = up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + 3];
	radius_s5 = up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + 4];
	radius_s6 = up[nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + nc - 1 + 5];
	/// Compute contact point positions in body frames
	Vec3 normal = Vec3(0, 1, 0);
	/// sphere 1 left
	Vec3 pos_InGround_HC_s1_l = calcn_l->findStationLocationInGround(*state, locSphere_s1_l);
	Vec3 contactPointpos_InGround_HC_s1_l = pos_InGround_HC_s1_l - radius_s1 * normal;
	Vec3 contactPointpos_InGround_HC_s1_l_adj = contactPointpos_InGround_HC_s1_l - 0.5*contactPointpos_InGround_HC_s1_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s1_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_l_adj, *calcn_l);
	/// sphere 2 left
	Vec3 pos_InGround_HC_s2_l = calcn_l->findStationLocationInGround(*state, locSphere_s2_l);
	Vec3 contactPointpos_InGround_HC_s2_l = pos_InGround_HC_s2_l - radius_s2 * normal;
	Vec3 contactPointpos_InGround_HC_s2_l_adj = contactPointpos_InGround_HC_s2_l - 0.5*contactPointpos_InGround_HC_s2_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s2_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_l_adj, *calcn_l);
	/// sphere 3 left
	Vec3 pos_InGround_HC_s3_l = calcn_l->findStationLocationInGround(*state, locSphere_s3_l);
	Vec3 contactPointpos_InGround_HC_s3_l = pos_InGround_HC_s3_l - radius_s3 * normal;
	Vec3 contactPointpos_InGround_HC_s3_l_adj = contactPointpos_InGround_HC_s3_l - 0.5*contactPointpos_InGround_HC_s3_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s3_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_l_adj, *calcn_l);
	/// sphere 4 left
	Vec3 pos_InGround_HC_s4_l = calcn_l->findStationLocationInGround(*state, locSphere_s4_l);
	Vec3 contactPointpos_InGround_HC_s4_l = pos_InGround_HC_s4_l - radius_s4 * normal;
	Vec3 contactPointpos_InGround_HC_s4_l_adj = contactPointpos_InGround_HC_s4_l - 0.5*contactPointpos_InGround_HC_s4_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s4_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_l_adj, *calcn_l);
	/// sphere 5 left
	Vec3 pos_InGround_HC_s5_l = toes_l->findStationLocationInGround(*state, locSphere_s5_l);
	Vec3 contactPointpos_InGround_HC_s5_l = pos_InGround_HC_s5_l - radius_s5 * normal;
	Vec3 contactPointpos_InGround_HC_s5_l_adj = contactPointpos_InGround_HC_s5_l - 0.5*contactPointpos_InGround_HC_s5_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s5_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_l_adj, *toes_l);
	/// sphere 6 left
	Vec3 pos_InGround_HC_s6_l = toes_l->findStationLocationInGround(*state, locSphere_s6_l);
	Vec3 contactPointpos_InGround_HC_s6_l = pos_InGround_HC_s6_l - radius_s6 * normal;
	Vec3 contactPointpos_InGround_HC_s6_l_adj = contactPointpos_InGround_HC_s6_l - 0.5*contactPointpos_InGround_HC_s6_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s6_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_l_adj, *toes_l);
	/// sphere 1 right
	Vec3 pos_InGround_HC_s1_r = calcn_r->findStationLocationInGround(*state, locSphere_s1_r);
	Vec3 contactPointpos_InGround_HC_s1_r = pos_InGround_HC_s1_r - radius_s1 * normal;
	Vec3 contactPointpos_InGround_HC_s1_r_adj = contactPointpos_InGround_HC_s1_r - 0.5*contactPointpos_InGround_HC_s1_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s1_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_r_adj, *calcn_r);
	/// sphere 2 right
	Vec3 pos_InGround_HC_s2_r = calcn_r->findStationLocationInGround(*state, locSphere_s2_r);
	Vec3 contactPointpos_InGround_HC_s2_r = pos_InGround_HC_s2_r - radius_s2 * normal;
	Vec3 contactPointpos_InGround_HC_s2_r_adj = contactPointpos_InGround_HC_s2_r - 0.5*contactPointpos_InGround_HC_s2_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s2_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_r_adj, *calcn_r);
	/// sphere 3 right
	Vec3 pos_InGround_HC_s3_r = calcn_r->findStationLocationInGround(*state, locSphere_s3_r);
	Vec3 contactPointpos_InGround_HC_s3_r = pos_InGround_HC_s3_r - radius_s3 * normal;
	Vec3 contactPointpos_InGround_HC_s3_r_adj = contactPointpos_InGround_HC_s3_r - 0.5*contactPointpos_InGround_HC_s3_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s3_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_r_adj, *calcn_r);
	/// sphere 4 right
	Vec3 pos_InGround_HC_s4_r = calcn_r->findStationLocationInGround(*state, locSphere_s4_r);
	Vec3 contactPointpos_InGround_HC_s4_r = pos_InGround_HC_s4_r - radius_s4 * normal;
	Vec3 contactPointpos_InGround_HC_s4_r_adj = contactPointpos_InGround_HC_s4_r - 0.5*contactPointpos_InGround_HC_s4_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s4_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_r_adj, *calcn_r);
	/// sphere 5 right
	Vec3 pos_InGround_HC_s5_r = toes_r->findStationLocationInGround(*state, locSphere_s5_r);
	Vec3 contactPointpos_InGround_HC_s5_r = pos_InGround_HC_s5_r - radius_s5 * normal;
	Vec3 contactPointpos_InGround_HC_s5_r_adj = contactPointpos_InGround_HC_s5_r - 0.5*contactPointpos_InGround_HC_s5_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s5_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_r_adj, *toes_r);
	/// sphere 6 right
	Vec3 pos_InGround_HC_s6_r = toes_r->findStationLocationInGround(*state, locSphere_s6_r);
	Vec3 contactPointpos_InGround_HC_s6_r = pos_InGround_HC_s6_r - radius_s6 * normal;
	Vec3 contactPointpos_InGround_HC_s6_r_adj = contactPointpos_InGround_HC_s6_r - 0.5*contactPointpos_InGround_HC_s6_r[1] * normal;
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

	/// knownUdot.
	Vector knownUdot(nCoordinates);
	knownUdot.setToZero();
	for (int i = 0; i < nCoordinates; ++i) knownUdot[i] = ua[i];
	/// Calculate residual forces.
	Vector residualMobilityForces(nCoordinates);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
		appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);

	/// Marker positions.
	Vec3 C7_study_location = C7_study->getLocationInGround(*state);
	Vec3 r_shoulder_study_location = r_shoulder_study->getLocationInGround(*state);
	Vec3 L_shoulder_study_location = L_shoulder_study->getLocationInGround(*state);
	Vec3 r_ASIS_study_location = r_ASIS_study->getLocationInGround(*state);
	Vec3 L_ASIS_study_location = L_ASIS_study->getLocationInGround(*state);
	Vec3 r_PSIS_study_location = r_PSIS_study->getLocationInGround(*state);
	Vec3 L_PSIS_study_location = L_PSIS_study->getLocationInGround(*state);
	Vec3 r_knee_study_location = r_knee_study->getLocationInGround(*state);
	Vec3 L_knee_study_location = L_knee_study->getLocationInGround(*state);
	Vec3 r_mknee_study_location = r_mknee_study->getLocationInGround(*state);
	Vec3 L_mknee_study_location = L_mknee_study->getLocationInGround(*state);
	Vec3 r_ankle_study_location = r_ankle_study->getLocationInGround(*state);
	Vec3 L_ankle_study_location = L_ankle_study->getLocationInGround(*state);
	Vec3 r_mankle_study_location = r_mankle_study->getLocationInGround(*state);
	Vec3 L_mankle_study_location = L_mankle_study->getLocationInGround(*state);
	Vec3 r_calc_study_location = r_calc_study->getLocationInGround(*state);
	Vec3 L_calc_study_location = L_calc_study->getLocationInGround(*state);
	Vec3 r_toe_study_location = r_toe_study->getLocationInGround(*state);
	Vec3 L_toe_study_location = L_toe_study->getLocationInGround(*state);
	Vec3 r_5meta_study_location = r_5meta_study->getLocationInGround(*state);
	Vec3 L_5meta_study_location = L_5meta_study->getLocationInGround(*state);
	Vec3 r_lelbow_study_location = r_lelbow_study->getLocationInGround(*state);
	Vec3 L_lelbow_study_location = L_lelbow_study->getLocationInGround(*state);
	Vec3 r_melbow_study_location = r_melbow_study->getLocationInGround(*state);
	Vec3 L_melbow_study_location = L_melbow_study->getLocationInGround(*state);
	Vec3 r_lwrist_study_location = r_lwrist_study->getLocationInGround(*state);
	Vec3 L_lwrist_study_location = L_lwrist_study->getLocationInGround(*state);
	Vec3 r_mwrist_study_location = r_mwrist_study->getLocationInGround(*state);
	Vec3 L_mwrist_study_location = L_mwrist_study->getLocationInGround(*state);

	/// Ground reaction forces and moments.
	Vec3 GRF_r, GRF_l, GRM_r, GRM_l;

	SimTK::Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);
	SimTK::Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
	SimTK::Transform TR_GB_toes_r = toes_r->getMobilizedBody().getBodyTransform(*state);
	SimTK::Transform TR_GB_toes_l = toes_l->getMobilizedBody().getBodyTransform(*state);

	Vec3 GRM_0 = (TR_GB_calcn_r*contactPointPos_InBody_HC_s1_r) % AppliedPointForce_s1_r;
	GRF_r += AppliedPointForce_s1_r;
	GRM_r += GRM_0;

	Vec3 GRM_1 = (TR_GB_calcn_r*contactPointPos_InBody_HC_s2_r) % AppliedPointForce_s2_r;
	GRF_r += AppliedPointForce_s2_r;
	GRM_r += GRM_1;

	Vec3 GRM_2 = (TR_GB_calcn_r*contactPointPos_InBody_HC_s3_r) % AppliedPointForce_s3_r;
	GRF_r += AppliedPointForce_s3_r;
	GRM_r += GRM_2;

	Vec3 GRM_3 = (TR_GB_calcn_r*contactPointPos_InBody_HC_s4_r) % AppliedPointForce_s4_r;
	GRF_r += AppliedPointForce_s4_r;
	GRM_r += GRM_3;

	Vec3 GRM_4 = (TR_GB_toes_r*contactPointPos_InBody_HC_s5_r) % AppliedPointForce_s5_r;
	GRF_r += AppliedPointForce_s5_r;
	GRM_r += GRM_4;

	Vec3 GRM_5 = (TR_GB_toes_r*contactPointPos_InBody_HC_s6_r) % AppliedPointForce_s6_r;
	GRF_r += AppliedPointForce_s6_r;
	GRM_r += GRM_5;

	Vec3 GRM_6 = (TR_GB_calcn_l*contactPointPos_InBody_HC_s1_l) % AppliedPointForce_s1_l;
	GRF_l += AppliedPointForce_s1_l;
	GRM_l += GRM_6;

	Vec3 GRM_7 = (TR_GB_calcn_l*contactPointPos_InBody_HC_s2_l) % AppliedPointForce_s2_l;
	GRF_l += AppliedPointForce_s2_l;
	GRM_l += GRM_7;

	Vec3 GRM_8 = (TR_GB_calcn_l*contactPointPos_InBody_HC_s3_l) % AppliedPointForce_s3_l;
	GRF_l += AppliedPointForce_s3_l;
	GRM_l += GRM_8;

	Vec3 GRM_9 = (TR_GB_calcn_l*contactPointPos_InBody_HC_s4_l) % AppliedPointForce_s4_l;
	GRF_l += AppliedPointForce_s4_l;
	GRM_l += GRM_9;

	Vec3 GRM_10 = (TR_GB_toes_l*contactPointPos_InBody_HC_s5_l) % AppliedPointForce_s5_l;
	GRF_l += AppliedPointForce_s5_l;
	GRM_l += GRM_10;

	Vec3 GRM_11 = (TR_GB_toes_l*contactPointPos_InBody_HC_s6_l) % AppliedPointForce_s6_l;
	GRF_l += AppliedPointForce_s6_l;
	GRM_l += GRM_11;

	/// Residual forces.
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
	for (int i = 0; i < NU; ++i) res[0][i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
	/// Marker positions.
	for (int i = 0; i < nc; ++i) res[0][i + NU + 0 * nc] = value<T>(C7_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 1 * nc] = value<T>(r_shoulder_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 2 * nc] = value<T>(L_shoulder_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 3 * nc] = value<T>(r_ASIS_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 4 * nc] = value<T>(L_ASIS_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 5 * nc] = value<T>(r_PSIS_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 6 * nc] = value<T>(L_PSIS_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 7 * nc] = value<T>(r_knee_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 8 * nc] = value<T>(L_knee_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 9 * nc] = value<T>(r_mknee_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 10 * nc] = value<T>(L_mknee_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 11 * nc] = value<T>(r_ankle_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 12 * nc] = value<T>(L_ankle_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 13 * nc] = value<T>(r_mankle_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 14 * nc] = value<T>(L_mankle_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 15 * nc] = value<T>(r_calc_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 16 * nc] = value<T>(L_calc_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 17 * nc] = value<T>(r_toe_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 18 * nc] = value<T>(L_toe_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 19 * nc] = value<T>(r_5meta_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 20 * nc] = value<T>(L_5meta_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 21 * nc] = value<T>(r_lelbow_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 22 * nc] = value<T>(L_lelbow_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 23 * nc] = value<T>(r_melbow_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 24 * nc] = value<T>(L_melbow_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 25 * nc] = value<T>(r_lwrist_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 26 * nc] = value<T>(L_lwrist_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 27 * nc] = value<T>(r_mwrist_study_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 28 * nc] = value<T>(L_mwrist_study_location[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 29 * nc] = value<T>(GRF_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 30 * nc] = value<T>(GRF_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 31 * nc] = value<T>(GRM_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 32 * nc] = value<T>(GRM_l[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 33 * nc] = value<T>(AppliedPointForce_s1_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 34 * nc] = value<T>(contactPointpos_InGround_HC_s1_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 35 * nc] = value<T>(AppliedPointForce_s2_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 36 * nc] = value<T>(contactPointpos_InGround_HC_s2_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 37 * nc] = value<T>(AppliedPointForce_s3_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 38 * nc] = value<T>(contactPointpos_InGround_HC_s3_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 39 * nc] = value<T>(AppliedPointForce_s4_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 40 * nc] = value<T>(contactPointpos_InGround_HC_s4_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 41 * nc] = value<T>(AppliedPointForce_s5_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 42 * nc] = value<T>(contactPointpos_InGround_HC_s5_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 43 * nc] = value<T>(AppliedPointForce_s6_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 44 * nc] = value<T>(contactPointpos_InGround_HC_s6_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 45 * nc] = value<T>(AppliedPointForce_s1_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 46 * nc] = value<T>(contactPointpos_InGround_HC_s1_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 47 * nc] = value<T>(AppliedPointForce_s2_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 48 * nc] = value<T>(contactPointpos_InGround_HC_s2_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 49 * nc] = value<T>(AppliedPointForce_s3_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 50 * nc] = value<T>(contactPointpos_InGround_HC_s3_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 51 * nc] = value<T>(AppliedPointForce_s4_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 52 * nc] = value<T>(contactPointpos_InGround_HC_s4_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 53 * nc] = value<T>(AppliedPointForce_s5_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 54 * nc] = value<T>(contactPointpos_InGround_HC_s5_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 55 * nc] = value<T>(AppliedPointForce_s6_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 56 * nc] = value<T>(contactPointpos_InGround_HC_s6_l_adj[i]);

	return 0;
}

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
