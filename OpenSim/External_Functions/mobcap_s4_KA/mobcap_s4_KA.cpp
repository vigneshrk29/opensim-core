#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
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

constexpr int n_in = 2;
constexpr int n_out = 1;
constexpr int nCoordinates = 33;
constexpr int NX = nCoordinates * 2;
constexpr int NU = nCoordinates;
constexpr int NR = nCoordinates;

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
	pelvis = new OpenSim::Body("pelvis", 13.03600897230877464494, Vec3(-0.08403627878184083411, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.11378973612578265062, 0.09641134257349871783, 0.06408974437434644678, 0., 0., 0.));
	model->addBody(pelvis);

	OpenSim::Body* femur_r;
	femur_r = new OpenSim::Body("femur_r", 10.29575731128749538357, Vec3(0.00000000000000000000, -0.18950135017324720721, 0.00000000000000000000), Inertia(0.18416933804266835684, 0.04827739929273831138, 0.19420993675597292283, 0., 0., 0.));
	model->addBody(femur_r);

	OpenSim::Body* tibia_r;
	tibia_r = new OpenSim::Body("tibia_r", 4.10384675764921347252, Vec3(0.00000000000000000000, -0.20438279111996127124, 0.00000000000000000000), Inertia(0.06685602132009388354, 0.00676519263358092877, 0.06778457717176185338, 0., 0., 0.));
	model->addBody(tibia_r);

	OpenSim::Body* talus_r;
	talus_r = new OpenSim::Body("talus_r", 0.11069040479161736112, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00118617421024884710, 0.00118617421024884710, 0.00118617421024884710, 0., 0., 0.));
	model->addBody(talus_r);

	OpenSim::Body* calcn_r;
	calcn_r = new OpenSim::Body("calcn_r", 1.38363005989521692385, Vec3(0.10351880521960737935, 0.03105564156588221103, 0.00000000000000000000), Inertia(0.00166064389434838599, 0.00462607941997050364, 0.00486331426202027406, 0., 0., 0.));
	model->addBody(calcn_r);

	OpenSim::Body* toes_r;
	toes_r = new OpenSim::Body("toes_r", 0.23975541677864317669, Vec3(0.03581750660598415070, 0.00621112831317644203, -0.01811579091343129208), Inertia(0.00011861742102488473, 0.00023723484204976947, 0.00118617421024884710, 0., 0., 0.));
	model->addBody(toes_r);

	OpenSim::Body* femur_l;
	femur_l = new OpenSim::Body("femur_l", 10.29575731128749538357, Vec3(0.00000000000000000000, -0.18950135017324720721, 0.00000000000000000000), Inertia(0.18416933804266835684, 0.04827739929273831138, 0.19420993675597292283, 0., 0., 0.));
	model->addBody(femur_l);

	OpenSim::Body* tibia_l;
	tibia_l = new OpenSim::Body("tibia_l", 4.10384675764921347252, Vec3(0.00000000000000000000, -0.20438279111996127124, 0.00000000000000000000), Inertia(0.06685602132009388354, 0.00676519263358092877, 0.06778457717176185338, 0., 0., 0.));
	model->addBody(tibia_l);

	OpenSim::Body* talus_l;
	talus_l = new OpenSim::Body("talus_l", 0.11069040479161736112, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Inertia(0.00118617421024884710, 0.00118617421024884710, 0.00118617421024884710, 0., 0., 0.));
	model->addBody(talus_l);

	OpenSim::Body* calcn_l;
	calcn_l = new OpenSim::Body("calcn_l", 1.38363005989521692385, Vec3(0.10351880521960737935, 0.03105564156588221103, 0.00000000000000000000), Inertia(0.00166064389434838599, 0.00462607941997050364, 0.00486331426202027406, 0., 0., 0.));
	model->addBody(calcn_l);

	OpenSim::Body* toes_l;
	toes_l = new OpenSim::Body("toes_l", 0.23975541677864317669, Vec3(0.03581750660598415070, 0.00621112831317644203, 0.01811579091343129208), Inertia(0.00011861742102488473, 0.00023723484204976947, 0.00118617421024884710, 0., 0., 0.));
	model->addBody(toes_l);

	OpenSim::Body* torso;
	torso = new OpenSim::Body("torso", 29.69447213182802158826, Vec3(-0.03020373518498941506, 0.32217317530655381619, 0.00000000000000000000), Inertia(1.63213001865239792920, 0.83626600820066909758, 1.58442245418721094907, 0., 0., 0.));
	model->addBody(torso);

	OpenSim::Body* humerus_r;
	humerus_r = new OpenSim::Body("humerus_r", 2.24978247738962311431, Vec3(0.00000000000000000000, -0.17705033960528174242, 0.00000000000000000000), Inertia(0.01531735047257441577, 0.00528401149317588806, 0.01719323225236484010, 0., 0., 0.));
	model->addBody(humerus_r);

	OpenSim::Body* ulna_r;
	ulna_r = new OpenSim::Body("ulna_r", 0.67244420910907554134, Vec3(0.00000000000000000000, -0.13347384756586416765, 0.00000000000000000000), Inertia(0.00402099107446640767, 0.00083895087239035774, 0.00436173002101977302, 0., 0., 0.));
	model->addBody(ulna_r);

	OpenSim::Body* radius_r;
	radius_r = new OpenSim::Body("radius_r", 0.67244420910907554134, Vec3(0.00000000000000000000, -0.13347384756586416765, 0.00000000000000000000), Inertia(0.00402099107446640767, 0.00083895087239035774, 0.00436173002101977302, 0., 0., 0.));
	model->addBody(radius_r);

	OpenSim::Body* hand_r;
	hand_r = new OpenSim::Body("hand_r", 0.50640860192164938169, Vec3(0.00000000000000000000, -0.07541092428954591920, 0.00000000000000000000), Inertia(0.00121091290966375250, 0.00074256654886331022, 0.00181908441586258818, 0., 0., 0.));
	model->addBody(hand_r);

	OpenSim::Body* humerus_l;
	humerus_l = new OpenSim::Body("humerus_l", 2.24978247738962311431, Vec3(0.00000000000000000000, -0.17705033960528174242, 0.00000000000000000000), Inertia(0.01531735047257441577, 0.00528401149317588806, 0.01719323225236484010, 0., 0., 0.));
	model->addBody(humerus_l);

	OpenSim::Body* ulna_l;
	ulna_l = new OpenSim::Body("ulna_l", 0.67244420910907554134, Vec3(0.00000000000000000000, -0.13347384756586416765, 0.00000000000000000000), Inertia(0.00402099107446640767, 0.00083895087239035774, 0.00436173002101977302, 0., 0., 0.));
	model->addBody(ulna_l);

	OpenSim::Body* radius_l;
	radius_l = new OpenSim::Body("radius_l", 0.67244420910907554134, Vec3(0.00000000000000000000, -0.13347384756586416765, 0.00000000000000000000), Inertia(0.00402099107446640767, 0.00083895087239035774, 0.00436173002101977302, 0., 0., 0.));
	model->addBody(radius_l);

	OpenSim::Body* hand_l;
	hand_l = new OpenSim::Body("hand_l", 0.50640860192164938169, Vec3(0.00000000000000000000, -0.07541092428954591920, 0.00000000000000000000), Inertia(0.00121091290966375250, 0.00074256654886331022, 0.00181908441586258818, 0., 0., 0.));
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
	hip_r = new OpenSim::CustomJoint("hip_r", *pelvis, Vec3(-0.08403627878184083411, -0.07856857181725147854, 0.10538323133391283570), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_r);

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
	knee_r = new OpenSim::CustomJoint("knee_r", *femur_r, Vec3(-0.00501605000000000020, -0.44122679999999997458, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_r);

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
	ankle_r = new OpenSim::CustomJoint("ankle_r", *tibia_r, Vec3(0.00000000000000000000, -0.47072629984779507106, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_r);

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
	subtalar_r = new OpenSim::CustomJoint("subtalar_r", *talus_r, Vec3(-0.05048612130560251737, -0.04342613878962529178, 0.00819868937339290373), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_r);

	OpenSim::WeldJoint* mtp_r;
	mtp_r = new OpenSim::WeldJoint("mtp_r", *calcn_r, Vec3(0.18509162373265797275, -0.00207037610439214734, 0.00111800309637175954), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	hip_l = new OpenSim::CustomJoint("hip_l", *pelvis, Vec3(-0.08403627878184083411, -0.07856857181725147854, -0.10538323133391283570), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *femur_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_hip_l);

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
	knee_l = new OpenSim::CustomJoint("knee_l", *femur_l, Vec3(-0.00501605000000000020, -0.44122679999999997458, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *tibia_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_knee_l);

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
	ankle_l = new OpenSim::CustomJoint("ankle_l", *tibia_l, Vec3(0.00000000000000000000, -0.47072629984779507106, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *talus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_ankle_l);

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
	subtalar_l = new OpenSim::CustomJoint("subtalar_l", *talus_l, Vec3(-0.05048612130560251737, -0.04342613878962529178, -0.00819868937339290373), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *calcn_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_subtalar_l);

	OpenSim::WeldJoint* mtp_l;
	mtp_l = new OpenSim::WeldJoint("mtp_l", *calcn_l, Vec3(0.18509162373265797275, -0.00207037610439214734, -0.00111800309637175954), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *toes_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	back = new OpenSim::CustomJoint("back", *pelvis, Vec3(-0.11969523724655406482, 0.09687350382913760549, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *torso, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_back);

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
	acromial_r = new OpenSim::CustomJoint("acromial_r", *torso, Vec3(0.00317642615028805369, 0.37402292070745229990, 0.18551486296816507848), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_r);

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
	elbow_r = new OpenSim::CustomJoint("elbow_r", *humerus_r, Vec3(0.01414663447114213046, -0.30811012552931160702, -0.01032691401024108009), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_r);

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
	radioulnar_r = new OpenSim::CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.00744972887430465393, -0.01440443339796054940, 0.02888527995071922885), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_r);

	OpenSim::WeldJoint* radius_hand_r;
	radius_hand_r = new OpenSim::WeldJoint("radius_hand_r", *radius_r, Vec3(-0.00974212351824855528, -0.26117905566298255104, 0.01507221792467464475), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_r, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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
	acromial_l = new OpenSim::CustomJoint("acromial_l", *torso, Vec3(0.00317642615028805369, 0.37402292070745229990, -0.18551486296816507848), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *humerus_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_acromial_l);

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
	elbow_l = new OpenSim::CustomJoint("elbow_l", *humerus_l, Vec3(0.01414663447114213046, -0.30811012552931160702, 0.01032691401024108009), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *ulna_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_elbow_l);

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
	radioulnar_l = new OpenSim::CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.00744972887430465393, -0.01440443339796054940, -0.02888527995071922885), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *radius_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), st_radioulnar_l);

	OpenSim::WeldJoint* radius_hand_l;
	radius_hand_l = new OpenSim::WeldJoint("radius_hand_l", *radius_l, Vec3(-0.00974212351824855528, -0.26117905566298255104, -0.01507221792467464475), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), *hand_l, Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000), Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));

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

	// Definition of contacts
	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s1_r;
	SmoothSphereHalfSpaceForce_s1_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s1_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s1_r_location(-0.00047744940253422494, -0.01000000000000000021, -0.00566025373492130661);
	SmoothSphereHalfSpaceForce_s1_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s1_r_location);
	double SmoothSphereHalfSpaceForce_s1_r_radius = (0.03232000000000000151);
	SmoothSphereHalfSpaceForce_s1_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s1_r_radius);
	SmoothSphereHalfSpaceForce_s1_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s1_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s1_r->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s1_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s1_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s1_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s1_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s1_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s1_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_s1_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s1_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s2_r;
	SmoothSphereHalfSpaceForce_s2_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s2_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s2_r_location(0.06796110303675625963, -0.01000000000000000021, 0.02265483369730270091);
	SmoothSphereHalfSpaceForce_s2_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s2_r_location);
	double SmoothSphereHalfSpaceForce_s2_r_radius = (0.03232000000000000151);
	SmoothSphereHalfSpaceForce_s2_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s2_r_radius);
	SmoothSphereHalfSpaceForce_s2_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s2_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s2_r->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s2_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s2_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s2_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s2_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s2_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s2_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_s2_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s2_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s3_r;
	SmoothSphereHalfSpaceForce_s3_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s3_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s3_r_location(0.18689303335107970705, -0.01000000000000000021, 0.02399366742712679898);
	SmoothSphereHalfSpaceForce_s3_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s3_r_location);
	double SmoothSphereHalfSpaceForce_s3_r_radius = (0.02337399999999999894);
	SmoothSphereHalfSpaceForce_s3_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s3_r_radius);
	SmoothSphereHalfSpaceForce_s3_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s3_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s3_r->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s3_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s3_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s3_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s3_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s3_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s3_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_s3_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s3_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s4_r;
	SmoothSphereHalfSpaceForce_s4_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s4_r", *calcn_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s4_r_location(0.18689303335107970705, -0.01000000000000000021, -0.01132685050612604327);
	SmoothSphereHalfSpaceForce_s4_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s4_r_location);
	double SmoothSphereHalfSpaceForce_s4_r_radius = (0.02050799999999999845);
	SmoothSphereHalfSpaceForce_s4_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s4_r_radius);
	SmoothSphereHalfSpaceForce_s4_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s4_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s4_r->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s4_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s4_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s4_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s4_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s4_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s4_r->connectSocket_sphere_frame(*calcn_r);
	SmoothSphereHalfSpaceForce_s4_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s4_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s5_r;
	SmoothSphereHalfSpaceForce_s5_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s5_r", *toes_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s5_r_location(0.06020674118026236865, -0.01000000000000000021, -0.00387072462345845270);
	SmoothSphereHalfSpaceForce_s5_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s5_r_location);
	double SmoothSphereHalfSpaceForce_s5_r_radius = (0.01624400000000000149);
	SmoothSphereHalfSpaceForce_s5_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s5_r_radius);
	SmoothSphereHalfSpaceForce_s5_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s5_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s5_r->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s5_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s5_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s5_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s5_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s5_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s5_r->connectSocket_sphere_frame(*toes_r);
	SmoothSphereHalfSpaceForce_s5_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s5_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s6_r;
	SmoothSphereHalfSpaceForce_s6_r = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s6_r", *toes_r, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s6_r_location(0.00000196871988646977, -0.01000000000000000021, 0.02525208051835740181);
	SmoothSphereHalfSpaceForce_s6_r->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s6_r_location);
	double SmoothSphereHalfSpaceForce_s6_r_radius = (0.01841399999999999981);
	SmoothSphereHalfSpaceForce_s6_r->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s6_r_radius);
	SmoothSphereHalfSpaceForce_s6_r->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s6_r->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s6_r->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s6_r->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s6_r->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s6_r->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s6_r->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s6_r->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s6_r->connectSocket_sphere_frame(*toes_r);
	SmoothSphereHalfSpaceForce_s6_r->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s6_r);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s1_l;
	SmoothSphereHalfSpaceForce_s1_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s1_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s1_l_location(-0.00047744940253422494, -0.01000000000000000021, 0.00566025373492130661);
	SmoothSphereHalfSpaceForce_s1_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s1_l_location);
	double SmoothSphereHalfSpaceForce_s1_l_radius = (0.03232000000000000151);
	SmoothSphereHalfSpaceForce_s1_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s1_l_radius);
	SmoothSphereHalfSpaceForce_s1_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s1_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s1_l->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s1_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s1_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s1_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s1_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s1_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s1_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_s1_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s1_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s2_l;
	SmoothSphereHalfSpaceForce_s2_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s2_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s2_l_location(0.06796110303675625963, -0.01000000000000000021, -0.02265483369730270091);
	SmoothSphereHalfSpaceForce_s2_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s2_l_location);
	double SmoothSphereHalfSpaceForce_s2_l_radius = (0.03232000000000000151);
	SmoothSphereHalfSpaceForce_s2_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s2_l_radius);
	SmoothSphereHalfSpaceForce_s2_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s2_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s2_l->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s2_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s2_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s2_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s2_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s2_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s2_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_s2_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s2_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s3_l;
	SmoothSphereHalfSpaceForce_s3_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s3_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s3_l_location(0.18689303335107970705, -0.01000000000000000021, -0.02399366742712679898);
	SmoothSphereHalfSpaceForce_s3_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s3_l_location);
	double SmoothSphereHalfSpaceForce_s3_l_radius = (0.02337399999999999894);
	SmoothSphereHalfSpaceForce_s3_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s3_l_radius);
	SmoothSphereHalfSpaceForce_s3_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s3_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s3_l->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s3_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s3_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s3_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s3_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s3_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s3_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_s3_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s3_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s4_l;
	SmoothSphereHalfSpaceForce_s4_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s4_l", *calcn_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s4_l_location(0.18689303335107970705, -0.01000000000000000021, 0.01132685050612604327);
	SmoothSphereHalfSpaceForce_s4_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s4_l_location);
	double SmoothSphereHalfSpaceForce_s4_l_radius = (0.02050799999999999845);
	SmoothSphereHalfSpaceForce_s4_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s4_l_radius);
	SmoothSphereHalfSpaceForce_s4_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s4_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s4_l->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s4_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s4_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s4_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s4_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s4_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s4_l->connectSocket_sphere_frame(*calcn_l);
	SmoothSphereHalfSpaceForce_s4_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s4_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s5_l;
	SmoothSphereHalfSpaceForce_s5_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s5_l", *toes_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s5_l_location(0.06020674118026236865, -0.01000000000000000021, 0.00387072462345845270);
	SmoothSphereHalfSpaceForce_s5_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s5_l_location);
	double SmoothSphereHalfSpaceForce_s5_l_radius = (0.01624400000000000149);
	SmoothSphereHalfSpaceForce_s5_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s5_l_radius);
	SmoothSphereHalfSpaceForce_s5_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s5_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s5_l->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s5_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s5_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s5_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s5_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s5_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s5_l->connectSocket_sphere_frame(*toes_l);
	SmoothSphereHalfSpaceForce_s5_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s5_l);

	OpenSim::SmoothSphereHalfSpaceForce* SmoothSphereHalfSpaceForce_s6_l;
	SmoothSphereHalfSpaceForce_s6_l = new SmoothSphereHalfSpaceForce("SmoothSphereHalfSpaceForce_s6_l", *toes_l, model->getGround());
	Vec3 SmoothSphereHalfSpaceForce_s6_l_location(0.00000196871988646977, -0.01000000000000000021, -0.02525208051835740181);
	SmoothSphereHalfSpaceForce_s6_l->set_contact_sphere_location(SmoothSphereHalfSpaceForce_s6_l_location);
	double SmoothSphereHalfSpaceForce_s6_l_radius = (0.01841399999999999981);
	SmoothSphereHalfSpaceForce_s6_l->set_contact_sphere_radius(SmoothSphereHalfSpaceForce_s6_l_radius);
	SmoothSphereHalfSpaceForce_s6_l->set_contact_half_space_location(Vec3(0.00000000000000000000, 0.00000000000000000000, 0.00000000000000000000));
	SmoothSphereHalfSpaceForce_s6_l->set_contact_half_space_orientation(Vec3(0.00000000000000000000, 0.00000000000000000000, -1.57079632679489655800));
	SmoothSphereHalfSpaceForce_s6_l->set_stiffness(1000000.00000000000000000000);
	SmoothSphereHalfSpaceForce_s6_l->set_dissipation(2.00000000000000000000);
	SmoothSphereHalfSpaceForce_s6_l->set_static_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s6_l->set_dynamic_friction(0.80000000000000004441);
	SmoothSphereHalfSpaceForce_s6_l->set_viscous_friction(0.50000000000000000000);
	SmoothSphereHalfSpaceForce_s6_l->set_transition_velocity(0.20000000000000001110);
	SmoothSphereHalfSpaceForce_s6_l->connectSocket_sphere_frame(*toes_l);
	SmoothSphereHalfSpaceForce_s6_l->connectSocket_half_space_frame(model->getGround());
	model->addComponent(SmoothSphereHalfSpaceForce_s6_l);

	// Initialize system.
	SimTK::State* state;
	state = new State(model->initSystem());

	// Read inputs.
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);

	// States and controls.
	T ua[NU];
	Vector QsUs(NX);
	/// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
	/// Controls
	/// OpenSim and Simbody have different state orders.
	auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
	for (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];

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
	/// Add contact forces to appliedBodyForces.
	Array<osim_double_adouble> Force_0 = SmoothSphereHalfSpaceForce_s1_r->getRecordValues(*state);
	SpatialVec GRF_0;
	GRF_0[0] = Vec3(Force_0[3], Force_0[4], Force_0[5]);
	GRF_0[1] = Vec3(Force_0[0], Force_0[1], Force_0[2]);
	int c_idx_0 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_0] += GRF_0;

	Array<osim_double_adouble> Force_1 = SmoothSphereHalfSpaceForce_s2_r->getRecordValues(*state);
	SpatialVec GRF_1;
	GRF_1[0] = Vec3(Force_1[3], Force_1[4], Force_1[5]);
	GRF_1[1] = Vec3(Force_1[0], Force_1[1], Force_1[2]);
	int c_idx_1 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_1] += GRF_1;

	Array<osim_double_adouble> Force_2 = SmoothSphereHalfSpaceForce_s3_r->getRecordValues(*state);
	SpatialVec GRF_2;
	GRF_2[0] = Vec3(Force_2[3], Force_2[4], Force_2[5]);
	GRF_2[1] = Vec3(Force_2[0], Force_2[1], Force_2[2]);
	int c_idx_2 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_2] += GRF_2;

	Array<osim_double_adouble> Force_3 = SmoothSphereHalfSpaceForce_s4_r->getRecordValues(*state);
	SpatialVec GRF_3;
	GRF_3[0] = Vec3(Force_3[3], Force_3[4], Force_3[5]);
	GRF_3[1] = Vec3(Force_3[0], Force_3[1], Force_3[2]);
	int c_idx_3 = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_3] += GRF_3;

	Array<osim_double_adouble> Force_4 = SmoothSphereHalfSpaceForce_s5_r->getRecordValues(*state);
	SpatialVec GRF_4;
	GRF_4[0] = Vec3(Force_4[3], Force_4[4], Force_4[5]);
	GRF_4[1] = Vec3(Force_4[0], Force_4[1], Force_4[2]);
	int c_idx_4 = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_4] += GRF_4;

	Array<osim_double_adouble> Force_5 = SmoothSphereHalfSpaceForce_s6_r->getRecordValues(*state);
	SpatialVec GRF_5;
	GRF_5[0] = Vec3(Force_5[3], Force_5[4], Force_5[5]);
	GRF_5[1] = Vec3(Force_5[0], Force_5[1], Force_5[2]);
	int c_idx_5 = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_5] += GRF_5;

	Array<osim_double_adouble> Force_6 = SmoothSphereHalfSpaceForce_s1_l->getRecordValues(*state);
	SpatialVec GRF_6;
	GRF_6[0] = Vec3(Force_6[3], Force_6[4], Force_6[5]);
	GRF_6[1] = Vec3(Force_6[0], Force_6[1], Force_6[2]);
	int c_idx_6 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_6] += GRF_6;

	Array<osim_double_adouble> Force_7 = SmoothSphereHalfSpaceForce_s2_l->getRecordValues(*state);
	SpatialVec GRF_7;
	GRF_7[0] = Vec3(Force_7[3], Force_7[4], Force_7[5]);
	GRF_7[1] = Vec3(Force_7[0], Force_7[1], Force_7[2]);
	int c_idx_7 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_7] += GRF_7;

	Array<osim_double_adouble> Force_8 = SmoothSphereHalfSpaceForce_s3_l->getRecordValues(*state);
	SpatialVec GRF_8;
	GRF_8[0] = Vec3(Force_8[3], Force_8[4], Force_8[5]);
	GRF_8[1] = Vec3(Force_8[0], Force_8[1], Force_8[2]);
	int c_idx_8 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_8] += GRF_8;

	Array<osim_double_adouble> Force_9 = SmoothSphereHalfSpaceForce_s4_l->getRecordValues(*state);
	SpatialVec GRF_9;
	GRF_9[0] = Vec3(Force_9[3], Force_9[4], Force_9[5]);
	GRF_9[1] = Vec3(Force_9[0], Force_9[1], Force_9[2]);
	int c_idx_9 = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_9] += GRF_9;

	Array<osim_double_adouble> Force_10 = SmoothSphereHalfSpaceForce_s5_l->getRecordValues(*state);
	SpatialVec GRF_10;
	GRF_10[0] = Vec3(Force_10[3], Force_10[4], Force_10[5]);
	GRF_10[1] = Vec3(Force_10[0], Force_10[1], Force_10[2]);
	int c_idx_10 = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_10] += GRF_10;

	Array<osim_double_adouble> Force_11 = SmoothSphereHalfSpaceForce_s6_l->getRecordValues(*state);
	SpatialVec GRF_11;
	GRF_11[0] = Vec3(Force_11[3], Force_11[4], Force_11[5]);
	GRF_11[1] = Vec3(Force_11[0], Force_11[1], Force_11[2]);
	int c_idx_11 = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
	appliedBodyForces[c_idx_11] += GRF_11;

	/// knownUdot.
	Vector knownUdot(nCoordinates);
	knownUdot.setToZero();
	for (int i = 0; i < nCoordinates; ++i) knownUdot[i] = ua[i];
	/// Calculate residual forces.
	Vector residualMobilityForces(nCoordinates);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
		appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);

	/// Residual forces.
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
	for (int i = 0; i < NU; ++i) res[0][i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);

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

