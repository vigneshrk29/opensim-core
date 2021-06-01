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
constexpr int nCoordinates = 35;
constexpr int NX = nCoordinates * 2;
constexpr int NU = nCoordinates;
constexpr int NP = 54;
constexpr int NR = nCoordinates + 3 * 4 + 3 * 2 * 12; ;

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
	pelvis = new OpenSim::Body("pelvis", 11.92356109125838514728, Vec3(-0.07145762715894608819, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.10407931393235647644, 0.08818393232984676755, 0.05862054743855486111, 0., 0., 0.));
	model->addBody(pelvis);

	OpenSim::Body* femur_r;
	femur_r = new OpenSim::Body("femur_r", 9.41715302150214306209, Vec3(0.00000000000000000000, -0.21050348397114138077, 0.00000000000000000000), Inertia(0.13556634373095846624, 0.03553680855083377055, 0.14295718995378145499, 0., 0., 0.));
	model->addBody(femur_r);

	OpenSim::Body* tibia_r;
	tibia_r = new OpenSim::Body("tibia_r", 3.75363868097482100694, Vec3(0.00000000000000000000, -0.24531293849250462857, 0.00000000000000000000), Inertia(0.05102721227812029303, 0.00516346790909550605, 0.05173592355976085189, 0., 0., 0.));
	model->addBody(tibia_r);

	OpenSim::Body* talus_r;
	talus_r = new OpenSim::Body("talus_r", 0.10124446880579424102, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00101244468805794237, 0.00101244468805794237, 0.00101244468805794237, 0., 0., 0.));
	model->addBody(talus_r);

	OpenSim::Body* calcn_r;
	calcn_r = new OpenSim::Body("calcn_r", 1.26555586007242792945, Vec3(0.10265105136586288870, 0.02999999999999999889, 0.00000000000000000000), Inertia(0.00141742256328111923, 0.00394853428342597460, 0.00415102322103756347, 0., 0., 0.));
	model->addBody(calcn_r);

	OpenSim::Body* toes_r;
	toes_r = new OpenSim::Body("toes_r", 0.21929551943335029440, Vec3(0.03551726377258855594, 0.00600000000000000012, -0.02057313799156976555), Inertia(0.00010124446880579423, 0.00020248893761158846, 0.00101244468805794237, 0., 0., 0.));
	model->addBody(toes_r);

	OpenSim::Body* femur_l;
	femur_l = new OpenSim::Body("femur_l", 9.41715302150214306209, Vec3(0.00000000000000000000, -0.21070155502705145700, 0.00000000000000000000), Inertia(0.13556634373095846624, 0.03553680855083377055, 0.14295718995378145499, 0., 0., 0.));
	model->addBody(femur_l);

	OpenSim::Body* tibia_l;
	tibia_l = new OpenSim::Body("tibia_l", 3.75363868097482100694, Vec3(0.00000000000000000000, -0.24398880428036534562, 0.00000000000000000000), Inertia(0.05102721227812029303, 0.00516346790909550605, 0.05173592355976085189, 0., 0., 0.));
	model->addBody(tibia_l);

	OpenSim::Body* talus_l;
	talus_l = new OpenSim::Body("talus_l", 0.10124446880579424102, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00101244468805794237, 0.00101244468805794237, 0.00101244468805794237, 0., 0., 0.));
	model->addBody(talus_l);

	OpenSim::Body* calcn_l;
	calcn_l = new OpenSim::Body("calcn_l", 1.26555586007242792945, Vec3(0.10857438534924904183, 0.02999999999999999889, 0.00000000000000000000), Inertia(0.00141742256328111923, 0.00394853428342597460, 0.00415102322103756347, 0., 0., 0.));
	model->addBody(calcn_l);

	OpenSim::Body* toes_l;
	toes_l = new OpenSim::Body("toes_l", 0.21929551943335029440, Vec3(0.03756673733084016364, 0.00600000000000000012, 0.02032169196515020565), Inertia(0.00010124446880579423, 0.00020248893761158846, 0.00101244468805794237, 0., 0., 0.));
	model->addBody(toes_l);

	OpenSim::Body* torso;
	torso = new OpenSim::Body("torso", 27.16044866865519225030, Vec3(-0.03707472427661409609, 0.32334285437141946362, 0.00000000000000000000), Inertia(1.49284969254143584649, 0.76490196182777536738, 1.44921332648613865324, 0., 0., 0.));
	model->addBody(torso);

	OpenSim::Body* humerus_r;
	humerus_r = new OpenSim::Body("humerus_r", 2.05779382847776792786, Vec3(0.00000000000000000000, -0.19874728826216173094, 0.00000000000000000000), Inertia(0.01765443711354647630, 0.00609023399840323380, 0.01981653668638412233, 0., 0., 0.));
	model->addBody(humerus_r);

	OpenSim::Body* ulna_r;
	ulna_r = new OpenSim::Body("ulna_r", 0.61506014799520003500, Vec3(0.00000000000000000000, -0.14423013458948741317, 0.00000000000000000000), Inertia(0.00429451473863476287, 0.00089601961798659133, 0.00465843209157106464, 0., 0., 0.));
	model->addBody(ulna_r);

	OpenSim::Body* radius_r;
	radius_r = new OpenSim::Body("radius_r", 0.61506014799520003500, Vec3(0.00000000000000000000, -0.14423013458948741317, 0.00000000000000000000), Inertia(0.00429451473863476287, 0.00089601961798659133, 0.00465843209157106464, 0., 0., 0.));
	model->addBody(radius_r);

	OpenSim::Body* hand_r;
	hand_r = new OpenSim::Body("hand_r", 0.46319344478650864572, Vec3(0.00000000000000000000, -0.08148808143431773843, 0.00000000000000000000), Inertia(0.00129328397935928735, 0.00079307885281337451, 0.00194282570890296521, 0., 0., 0.));
	model->addBody(hand_r);

	OpenSim::Body* humerus_l;
	humerus_l = new OpenSim::Body("humerus_l", 2.05779382847776792786, Vec3(0.00000000000000000000, -0.19512614110079826157, 0.00000000000000000000), Inertia(0.01701697509927250010, 0.00587032934740515377, 0.01910100611971747575, 0., 0., 0.));
	model->addBody(humerus_l);

	OpenSim::Body* ulna_l;
	ulna_l = new OpenSim::Body("ulna_l", 0.61506014799520003500, Vec3(0.00000000000000000000, -0.13592145359085333611, 0.00000000000000000000), Inertia(0.00381397727260398027, 0.00079575893128604322, 0.00413717386120073926, 0., 0., 0.));
	model->addBody(ulna_l);

	OpenSim::Body* radius_l;
	radius_l = new OpenSim::Body("radius_l", 0.61506014799520003500, Vec3(0.00000000000000000000, -0.13592145359085333611, 0.00000000000000000000), Inertia(0.00381397727260398027, 0.00079575893128604322, 0.00413717386120073926, 0., 0., 0.));
	model->addBody(radius_l);

	OpenSim::Body* hand_l;
	hand_l = new OpenSim::Body("hand_l", 0.46319344478650864572, Vec3(0.00000000000000000000, -0.07679378869337613600, 0.00000000000000000000), Inertia(0.00114857114353907855, 0.00070433678869492812, 0.00172543198693090281, 0., 0., 0.));
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
	hip_r = new OpenSim::CustomJoint("hip_r", *pelvis, Vec3(-0.07145762715894608819, -0.06745375677175843854, 0.10280162834777305247), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_r);

	SpatialTransform st_knee_r;
	st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
	st_knee_r[0].setFunction(new LinearFunction());
	st_knee_r[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_r[1].setCoordinateNames(OpenSim::Array<std::string>("knee_adduction_r", 1, 1));
	st_knee_r[1].setFunction(new LinearFunction());
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
	knee_r = new OpenSim::CustomJoint("knee_r", *femur_r, Vec3(-0.00519042000000000028, -0.49012727000000000377, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_r);

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
	ankle_r = new OpenSim::CustomJoint("ankle_r", *tibia_r, Vec3(0.00000000000000000000, -0.56499498420876803628, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_r);

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
	subtalar_r = new OpenSim::CustomJoint("subtalar_r", *talus_r, Vec3(-0.05006291775113132919, -0.04195000000000000118, 0.00931081445104185897), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_r);

	OpenSim::PinJoint* mtp_r;
	mtp_r = new OpenSim::PinJoint("mtp_r", *calcn_r, Vec3(0.18354007984216280658, -0.00200000000000000004, 0.00126965651605116261), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	hip_l = new OpenSim::CustomJoint("hip_l", *pelvis, Vec3(-0.07145762715894608819, -0.06745375677175843854, -0.10280162834777305247), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_l);

	SpatialTransform st_knee_l;
	st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
	st_knee_l[0].setFunction(new LinearFunction());
	st_knee_l[0].setAxis(Vec3(0.00000000000000000000, 0.00000000000000000000, 1.00000000000000000000));
	st_knee_l[1].setCoordinateNames(OpenSim::Array<std::string>("knee_adduction_l", 1, 1));
	st_knee_l[1].setFunction(new LinearFunction());
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
	knee_l = new OpenSim::CustomJoint("knee_l", *femur_l, Vec3(-0.00519042000000000028, -0.49012727000000000377, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_l);

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
	ankle_l = new OpenSim::CustomJoint("ankle_l", *tibia_l, Vec3(0.00000000000000000000, -0.56194529105815260728, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_l);

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
	subtalar_l = new OpenSim::CustomJoint("subtalar_l", *talus_l, Vec3(-0.05295172773482875739, -0.04195000000000000118, -0.00919701716365654842), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_l);

	OpenSim::PinJoint* mtp_l;
	mtp_l = new OpenSim::PinJoint("mtp_l", *calcn_l, Vec3(0.19413100100445726537, -0.00200000000000000004, -0.00125413870413498395), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	back = new OpenSim::CustomJoint("back", *pelvis, Vec3(-0.10177910968749465104, 0.08316915547501230188, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *torso, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_back);

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
	acromial_r = new OpenSim::CustomJoint("acromial_r", *torso, Vec3(0.00389902516975724903, 0.37538084499681972916, 0.20940702917534731453), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_r);

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
	elbow_r = new OpenSim::CustomJoint("elbow_r", *humerus_r, Vec3(0.01588025894468063359, -0.34586802867244059723, -0.01159244404855528515), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_r);

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
	radioulnar_r = new OpenSim::CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.00805008185342030204, -0.01556524671732389779, 0.03121306451356647971), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_r);

	OpenSim::WeldJoint* radius_hand_r;
	radius_hand_r = new OpenSim::WeldJoint("radius_hand_r", *radius_r, Vec3(-0.01052721422098088118, -0.28222675106176564608, 0.01628684614613502488), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	acromial_l = new OpenSim::CustomJoint("acromial_l", *torso, Vec3(0.00389902516975724903, 0.37538084499681972916, -0.20940702917534731453), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_l);

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
	elbow_l = new OpenSim::CustomJoint("elbow_l", *humerus_l, Vec3(0.01559092289837747843, -0.33956636266640416455, 0.01138123137628818497), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_l);

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
	radioulnar_l = new OpenSim::CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.00758633991541730369, -0.01466857786232092509, -0.02941497012246610879), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_l);

	OpenSim::WeldJoint* radius_hand_l;
	radius_hand_l = new OpenSim::WeldJoint("radius_hand_l", *radius_l, Vec3(-0.00992077185014508985, -0.26596848401842310139, -0.01534860803461119381), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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

	/// Ground reaction forces and moments.
	Vec3 GRF_r(0), GRF_l(0), GRM_r(0), GRM_l(0);

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

	for (int i = 0; i < nc; ++i) res[0][i + NU + 0 * nc] = value<T>(GRF_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 1 * nc] = value<T>(GRF_l[i]);
	
	for (int i = 0; i < nc; ++i) res[0][i + NU + 2 * nc] = value<T>(GRM_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 3 * nc] = value<T>(GRM_l[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 4 * nc] = value<T>(AppliedPointForce_s1_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 5 * nc] = value<T>(contactPointpos_InGround_HC_s1_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 6 * nc] = value<T>(AppliedPointForce_s2_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 7 * nc] = value<T>(contactPointpos_InGround_HC_s2_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 8 * nc] = value<T>(AppliedPointForce_s3_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 9 * nc] = value<T>(contactPointpos_InGround_HC_s3_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 10 * nc] = value<T>(AppliedPointForce_s4_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 11 * nc] = value<T>(contactPointpos_InGround_HC_s4_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 12 * nc] = value<T>(AppliedPointForce_s5_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 13 * nc] = value<T>(contactPointpos_InGround_HC_s5_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 14 * nc] = value<T>(AppliedPointForce_s6_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 15 * nc] = value<T>(contactPointpos_InGround_HC_s6_r_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 16 * nc] = value<T>(AppliedPointForce_s1_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 17 * nc] = value<T>(contactPointpos_InGround_HC_s1_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 18 * nc] = value<T>(AppliedPointForce_s2_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 19 * nc] = value<T>(contactPointpos_InGround_HC_s2_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 20 * nc] = value<T>(AppliedPointForce_s3_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 21 * nc] = value<T>(contactPointpos_InGround_HC_s3_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 22 * nc] = value<T>(AppliedPointForce_s4_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 23 * nc] = value<T>(contactPointpos_InGround_HC_s4_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 24 * nc] = value<T>(AppliedPointForce_s5_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 25 * nc] = value<T>(contactPointpos_InGround_HC_s5_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 26 * nc] = value<T>(AppliedPointForce_s6_l[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 27 * nc] = value<T>(contactPointpos_InGround_HC_s6_l_adj[i]);

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
