/*  This code describes the OpenSim model and the skeleton dynamics
	Author: Antoine Falisse
	Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SliderJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Common/SimmSpline.h>
//#include <OpenSim/Common/MultiplierFunction.h>
#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
#include "SimTKcommon/internal/recorder.h"
//#include <OpenSim/Simulation/Model/PhysicalOffsetFrame.h>

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
constexpr int ndof = 25;        // # degrees of freedom
constexpr int ndof_treadmill = ndof + 1;  // # degrees of freedom including treadmill
constexpr int NX = ndof * 2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NR = ndof + 44 * 3;        // # residual torques

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
    OpenSim::Body* femur_l;
    OpenSim::Body* tibia_r;
    OpenSim::Body* tibia_l;
    OpenSim::Body* talus_r;
    OpenSim::Body* talus_l;
    OpenSim::Body* calcn_r;
    OpenSim::Body* calcn_l;
    OpenSim::Body* toes_r;
    OpenSim::Body* toes_l;
    OpenSim::Body* torso;
    //OpenSim::Body* humerus_r;
    //OpenSim::Body* humerus_l;
    //OpenSim::Body* ulna_r;
    //OpenSim::Body* ulna_l;
    //OpenSim::Body* radius_r;
    //OpenSim::Body* radius_l;
    //OpenSim::Body* hand_r;
    //OpenSim::Body* hand_l;
    OpenSim::Body* treadmill;

    /// Joints
    OpenSim::CustomJoint* ground_pelvis;
    OpenSim::CustomJoint* hip_r;
    OpenSim::CustomJoint* hip_l;
    OpenSim::CustomJoint* knee_r;
    OpenSim::CustomJoint* knee_l;
    OpenSim::CustomJoint* ankle_r;
    OpenSim::CustomJoint* ankle_l;
    OpenSim::CustomJoint* subtalar_r;
    OpenSim::CustomJoint* subtalar_l;
    OpenSim::PinJoint* mtp_r;
    OpenSim::PinJoint* mtp_l;
    OpenSim::CustomJoint* back;
    //OpenSim::CustomJoint* shoulder_r;
    //OpenSim::CustomJoint* shoulder_l;
    //OpenSim::CustomJoint* elbow_r;
    //OpenSim::CustomJoint* elbow_l;
    //OpenSim::CustomJoint* radioulnar_r;
    //OpenSim::CustomJoint* radioulnar_l;
    //OpenSim::WeldJoint* radius_hand_r;
    //OpenSim::WeldJoint* radius_hand_l;
    OpenSim::SliderJoint* ground_treadmill;

    /// Contact elements
    OpenSim::SmoothSphereHalfSpaceForce* HC_1_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_2_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_3_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_4_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_5_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_6_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_1_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_2_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_3_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_4_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_5_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_6_l;

    /// OpenPose markers
    OpenSim::Station* Neck;
    OpenSim::Station* RShoulder;
    OpenSim::Station* LShoulder;
    OpenSim::Station* MidHip;
    OpenSim::Station* RHip;
    OpenSim::Station* LHip;
    OpenSim::Station* RKnee;
    OpenSim::Station* LKnee;
    OpenSim::Station* RAnkle;
    OpenSim::Station* LAnkle;
    OpenSim::Station* RHeel;
    OpenSim::Station* LHeel;
    OpenSim::Station* RBigToe;
    OpenSim::Station* LBigToe;
    OpenSim::Station* RSmallToe;
    OpenSim::Station* LSmallToe;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
    pelvis = new OpenSim::Body("pelvis", 8.883914768388312, Vec3(-0.06793666928317732, 0, 0), Inertia(0.077546611037642738, 0.065703402931699253, 0.043676544543575044, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 7.016459609975974, Vec3(0, -0.15838211714949457, 0), Inertia(0.087672770480383189, 0.022982182553110159, 0.092452540641001543, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 7.016459609975974, Vec3(0, -0.15838211714949457, 0), Inertia(0.087672770480383189, 0.022982182553110159, 0.092452540641001543, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 2.7967321052729619, Vec3(0, -0.1878651555359227, 0), Inertia(0.038494978516929027, 0.0038953252070701998, 0.039029630996330819, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 2.7967321052729619, Vec3(0, -0.1878651555359227, 0), Inertia(0.038494978516929027, 0.0038953252070701998, 0.039029630996330819, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.075434446534671934, Vec3(0, 0, 0), Inertia(0.00075434446534671927, 0.00075434446534671927, 0.00075434446534671927, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.075434446534671934, Vec3(0, 0, 0), Inertia(0.00075434446534671927, 0.00075434446534671927, 0.00075434446534671927, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 0.94293058168339905, Vec3(0.090067006194916144, 0.027020101858474841, 0), Inertia(0.001056082251485407, 0.0029419434148522049, 0.0030928123079215493, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 0.94293058168339905, Vec3(0.090067006194916144, 0.027020101858474841, 0), Inertia(0.001056082251485407, 0.0029419434148522049, 0.0030928123079215493, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.16339101119409938, Vec3(0.031163184143440985, 0.0054040203716949689, 0.014995605804549063), Inertia(7.5434446534671933e-05, 0.00015086889306934387, 0.00075434446534671927, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.16339101119409938, Vec3(0.031163184143440985, 0.0054040203716949689, -0.014995605804549063), Inertia(7.5434446534671933e-05, 0.00015086889306934387, 0.00075434446534671927, 0, 0, 0));
    torso = new OpenSim::Body("torso", 25.826189722289492, Vec3(-0.027404477114772164, 0.30851856360522822, 0), Inertia(1.1122809141537375, 0.56990724356944633, 1.079768667697294, 0, 0, 0));
    //humerus_l = new OpenSim::Body("humerus_l", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
    //humerus_r = new OpenSim::Body("humerus_r", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
    //ulna_l = new OpenSim::Body("ulna_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //ulna_r = new OpenSim::Body("ulna_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //radius_l = new OpenSim::Body("radius_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //radius_r = new OpenSim::Body("radius_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //hand_l = new OpenSim::Body("hand_l", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
    //hand_r = new OpenSim::Body("hand_r", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0)); 
    treadmill = new OpenSim::Body("treadmill", 1., Vec3(0), Inertia(1, 1, 1, 0, 0, 0));

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
    /// Hip_l transform
    SpatialTransform st_hip_l;
    st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
    st_hip_l[0].setFunction(new LinearFunction());
    st_hip_l[0].setAxis(Vec3(0, 0, 1));
    st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
    st_hip_l[1].setFunction(new LinearFunction());
    st_hip_l[1].setAxis(Vec3(-1, 0, 0));
    st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
    st_hip_l[2].setFunction(new LinearFunction());
    st_hip_l[2].setAxis(Vec3(0, -1, 0));
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
    /// Knee_l transform
    SpatialTransform st_knee_l;
    st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    st_knee_l[0].setFunction(new LinearFunction());
    st_knee_l[0].setAxis(Vec3(0, 0, 1));
    st_knee_l[1].setCoordinateNames(OpenSim::Array<std::string>("knee_adduction_l", 1, 1));
    st_knee_l[1].setFunction(new LinearFunction());
    st_knee_l[1].setAxis(Vec3(-1, 0, 0));
    st_knee_l[2].setFunction(new Constant(0));
    st_knee_l[2].setAxis(Vec3(0, -1, 0));
    /// Knee_r transform
    SpatialTransform st_knee_r;
    st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[0].setFunction(new LinearFunction());
    st_knee_r[0].setAxis(Vec3(0, 0, 1));
    st_knee_r[1].setCoordinateNames(OpenSim::Array<std::string>("knee_adduction_r", 1, 1));
    st_knee_r[1].setFunction(new LinearFunction());
    st_knee_r[1].setAxis(Vec3(1, 0, 0));
    st_knee_r[2].setFunction(new Constant(0));
    st_knee_r[2].setAxis(Vec3(0, 1, 0));
    /// Ankle_l transform
    SpatialTransform st_ankle_l;
    st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
    st_ankle_l[0].setFunction(new LinearFunction());
    st_ankle_l[0].setAxis(Vec3(0.10501355, 0.17402245, 0.97912632));
    /// Ankle_r transform
    SpatialTransform st_ankle_r;
    st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
    st_ankle_r[0].setFunction(new LinearFunction());
    st_ankle_r[0].setAxis(Vec3(-0.10501355, -0.17402245, 0.97912632));
    /// Subtalar_l transform
    SpatialTransform st_subtalar_l;
    st_subtalar_l[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_l", 1, 1));
    st_subtalar_l[0].setFunction(new LinearFunction());
    st_subtalar_l[0].setAxis(Vec3(-0.78717961, -0.60474746, -0.12094949));
    /// Subtalar_r transform
    SpatialTransform st_subtalar_r;
    st_subtalar_r[0].setCoordinateNames(OpenSim::Array<std::string>("subtalar_angle_r", 1, 1));
    st_subtalar_r[0].setFunction(new LinearFunction());
    st_subtalar_r[0].setAxis(Vec3(0.78717961, 0.60474746, -0.12094949));
    /// Back transform
    SpatialTransform st_back;
    st_back[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_extension", 1, 1));
    st_back[0].setFunction(new LinearFunction());
    st_back[0].setAxis(Vec3(0, 0, 1));
    st_back[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_bending", 1, 1));
    st_back[1].setFunction(new LinearFunction());
    st_back[1].setAxis(Vec3(1, 0, 0));
    st_back[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_rotation", 1, 1));
    st_back[2].setFunction(new LinearFunction());
    st_back[2].setAxis(Vec3(0, 1, 0));
    ///// Shoulder_l transform
    //SpatialTransform st_sho_l;
    //st_sho_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
    //st_sho_l[0].setFunction(new LinearFunction());
    //st_sho_l[0].setAxis(Vec3(0, 0, 1));
    //st_sho_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
    //st_sho_l[1].setFunction(new LinearFunction());
    //st_sho_l[1].setAxis(Vec3(-1, 0, 0));
    //st_sho_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
    //st_sho_l[2].setFunction(new LinearFunction());
    //st_sho_l[2].setAxis(Vec3(0, -1, 0));
    ///// Shoulder_r transform
    //SpatialTransform st_sho_r;
    //st_sho_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
    //st_sho_r[0].setFunction(new LinearFunction());
    //st_sho_r[0].setAxis(Vec3(0, 0, 1));
    //st_sho_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
    //st_sho_r[1].setFunction(new LinearFunction());
    //st_sho_r[1].setAxis(Vec3(1, 0, 0));
    //st_sho_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
    //st_sho_r[2].setFunction(new LinearFunction());
    //st_sho_r[2].setAxis(Vec3(0, 1, 0));
    ///// Elbow_l transform
    //SpatialTransform st_elb_l;
    //st_elb_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
    //st_elb_l[0].setFunction(new LinearFunction());
    //st_elb_l[0].setAxis(Vec3(-0.22604696, -0.022269, 0.97386183));
    ///// Elbow_r transform
    //SpatialTransform st_elb_r;
    //st_elb_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
    //st_elb_r[0].setFunction(new LinearFunction());
    //st_elb_r[0].setAxis(Vec3(0.22604696, 0.022269, 0.97386183));
    ///// Radioulnar_l transform
    //SpatialTransform st_radioulnar_l;
    //st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
    //st_radioulnar_l[0].setFunction(new LinearFunction());
    //st_radioulnar_l[0].setAxis(Vec3(-0.05639803, -0.99840646, 0.001952));
    ///// Radioulnar_r transform
    //SpatialTransform st_radioulnar_r;
    //st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
    //st_radioulnar_r[0].setFunction(new LinearFunction());
    //st_radioulnar_r[0].setAxis(Vec3(0.05639803, 0.99840646, 0.001952));
    /// Joint specifications
    ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.06793666928317732, -0.063728365794704964, -0.080236377441942106), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.06793666928317732, -0.063728365794704964, 0.080236377441942106), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *femur_l, Vec3(-0.00419233, -0.36877012, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(-0.00419233, -0.36877012, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.43268353979885782, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.43268353979885782, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
    subtalar_l = new CustomJoint("subtalar_l", *talus_l, Vec3(-0.043925678921260605, -0.037783109098767323, -0.0067865827412587751), Vec3(0), *calcn_l, Vec3(0), Vec3(0), st_subtalar_l);
    subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.043925678921260605, -0.037783109098767323, 0.0067865827412587751), Vec3(0), *calcn_r, Vec3(0), Vec3(0), st_subtalar_r);
    mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.16103980707651006, -0.0018013401238983228, -0.00092544310108074216), Vec3(0), *toes_l, Vec3(0), Vec3(0));
    mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.16103980707651006, -0.0018013401238983228, 0.00092544310108074216), Vec3(0), *toes_r, Vec3(0), Vec3(0));
    back = new CustomJoint("back", *pelvis, Vec3(-0.09676411028028227, 0.078575821668206564, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
    //shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, -0.16500992260263453), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
    //shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, 0.16500992260263453), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
    //elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.014421812854093517, -0.31410344120358441, 0.010527791717515771), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
    //elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.014421812854093517, -0.31410344120358441, -0.010527791717515771), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
    //radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.0078589566368657427, -0.015195696294888168, -0.030472003264362887), Vec3(0), *radius_l, Vec3(0), Vec3(0),st_radioulnar_l);
    //radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.0078589566368657427, -0.015195696294888168, 0.030472003264362887), Vec3(0), *radius_r, Vec3(0), Vec3(0),st_radioulnar_r);
    //radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.010277276874462305, -0.27552611746618899, -0.015900163494535866), Vec3(0), *hand_l, Vec3(0), Vec3(0));
    //radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.010277276874462305, -0.27552611746618899, 0.015900163494535866), Vec3(0), *hand_r, Vec3(0), Vec3(0));
    ground_treadmill = new SliderJoint("ground_treadmill", model->getGround(), Vec3(0), Vec3(0), *treadmill, Vec3(0), Vec3(0));

    /// Add bodies and joints to model
    model->addBody(pelvis);		    model->addJoint(ground_pelvis);
    model->addBody(femur_l);		model->addJoint(hip_l);
    model->addBody(femur_r);		model->addJoint(hip_r);
    model->addBody(tibia_l);		model->addJoint(knee_l);
    model->addBody(tibia_r);		model->addJoint(knee_r);
    model->addBody(talus_l);		model->addJoint(ankle_l);
    model->addBody(talus_r);		model->addJoint(ankle_r);
    model->addBody(calcn_l);		model->addJoint(subtalar_l);
    model->addBody(calcn_r);		model->addJoint(subtalar_r);
    model->addBody(toes_l);		    model->addJoint(mtp_l);
    model->addBody(toes_r);		    model->addJoint(mtp_r);
    model->addBody(torso);          model->addJoint(back);
    //model->addBody(humerus_l);      model->addJoint(shoulder_l);
    //model->addBody(humerus_r);      model->addJoint(shoulder_r);
    //model->addBody(ulna_l);         model->addJoint(elbow_l);
    //model->addBody(ulna_r);         model->addJoint(elbow_r);
    //model->addBody(radius_l);       model->addJoint(radioulnar_l);
    //model->addBody(radius_r);       model->addJoint(radioulnar_r);
    //model->addBody(hand_l);         model->addJoint(radius_hand_l);
    //model->addBody(hand_r);         model->addJoint(radius_hand_r);
    model->addBody(treadmill);		model->addJoint(ground_treadmill);
    /// Contact elements
    /// Parameters
    osim_double_adouble radiusSphere_1 = 0.03232;
    osim_double_adouble radiusSphere_2 = 0.03232;
    osim_double_adouble radiusSphere_3 = 0.023374;
    osim_double_adouble radiusSphere_4 = 0.020508;
    osim_double_adouble radiusSphere_5 = 0.016244;
    osim_double_adouble radiusSphere_6 = 0.018414;
    osim_double_adouble stiffness = 1000000;
    osim_double_adouble dissipation = 2.0;
    osim_double_adouble staticFriction = 0.8;
    osim_double_adouble dynamicFriction = 0.8;
    osim_double_adouble viscousFriction = 0.5;
    osim_double_adouble transitionVelocity = 0.2;
    Vec3 halfSpaceLocation(0);
    Vec3 halfSpaceOrientation(0, 0, -0.5*SimTK::Pi);
    Vec3 locSphere_1_r(-0.00041541, -0.00840834, -0.00468536);
    Vec3 locSphere_2_r(0.05912989, -0.00840834, 0.01875287);
    Vec3 locSphere_3_r(0.1626072, -0.00840834, 0.01986111);
    Vec3 locSphere_4_r(0.1626072, -0.00840834, -0.00937597);
    Vec3 locSphere_5_r(0.05238317, -0.00840834, -0.00320405);
    Vec3 locSphere_6_r(1.71289440e-06, -8.40833608e-03, 2.09027808e-02);
    Vec3 locSphere_1_l(locSphere_1_r[0], locSphere_1_r[1], -locSphere_1_r[2]);
    Vec3 locSphere_2_l(locSphere_2_r[0], locSphere_2_r[1], -locSphere_2_r[2]);
    Vec3 locSphere_3_l(locSphere_3_r[0], locSphere_3_r[1], -locSphere_3_r[2]);
    Vec3 locSphere_4_l(locSphere_4_r[0], locSphere_4_r[1], -locSphere_4_r[2]);
    Vec3 locSphere_5_l(locSphere_5_r[0], locSphere_5_r[1], -locSphere_5_r[2]);
    Vec3 locSphere_6_l(locSphere_6_r[0], locSphere_6_r[1], -locSphere_6_r[2]);
    /// Left foot contact shere specifications
    HC_1_l = new SmoothSphereHalfSpaceForce("sphere_1_l", *calcn_l, *treadmill);
    HC_1_l->set_contact_sphere_location(locSphere_1_l);
    HC_1_l->set_contact_sphere_radius(radiusSphere_1);
    HC_1_l->set_contact_half_space_location(halfSpaceLocation);
    HC_1_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_1_l->set_stiffness(stiffness);
    HC_1_l->set_dissipation(dissipation);
    HC_1_l->set_static_friction(staticFriction);
    HC_1_l->set_dynamic_friction(dynamicFriction);
    HC_1_l->set_viscous_friction(viscousFriction);
    HC_1_l->set_transition_velocity(transitionVelocity);
    HC_1_l->connectSocket_sphere_frame(*calcn_l);
    HC_1_l->connectSocket_half_space_frame(*treadmill);

    HC_2_l = new SmoothSphereHalfSpaceForce("sphere_2_l", *calcn_l, *treadmill);
    HC_2_l->set_contact_sphere_location(locSphere_2_l);
    HC_2_l->set_contact_sphere_radius(radiusSphere_2);
    HC_2_l->set_contact_half_space_location(halfSpaceLocation);
    HC_2_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_2_l->set_stiffness(stiffness);
    HC_2_l->set_dissipation(dissipation);
    HC_2_l->set_static_friction(staticFriction);
    HC_2_l->set_dynamic_friction(dynamicFriction);
    HC_2_l->set_viscous_friction(viscousFriction);
    HC_2_l->set_transition_velocity(transitionVelocity);
    HC_2_l->connectSocket_sphere_frame(*calcn_l);
    HC_2_l->connectSocket_half_space_frame(*treadmill);

    HC_3_l = new SmoothSphereHalfSpaceForce("sphere_3_l", *calcn_l, *treadmill);
    HC_3_l->set_contact_sphere_location(locSphere_3_l);
    HC_3_l->set_contact_sphere_radius(radiusSphere_3);
    HC_3_l->set_contact_half_space_location(halfSpaceLocation);
    HC_3_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_3_l->set_stiffness(stiffness);
    HC_3_l->set_dissipation(dissipation);
    HC_3_l->set_static_friction(staticFriction);
    HC_3_l->set_dynamic_friction(dynamicFriction);
    HC_3_l->set_viscous_friction(viscousFriction);
    HC_3_l->set_transition_velocity(transitionVelocity);
    HC_3_l->connectSocket_sphere_frame(*calcn_l);
    HC_3_l->connectSocket_half_space_frame(*treadmill);

    HC_4_l = new SmoothSphereHalfSpaceForce("sphere_4_l", *calcn_l, *treadmill);
    HC_4_l->set_contact_sphere_location(locSphere_4_l);
    HC_4_l->set_contact_sphere_radius(radiusSphere_4);
    HC_4_l->set_contact_half_space_location(halfSpaceLocation);
    HC_4_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_4_l->set_stiffness(stiffness);
    HC_4_l->set_dissipation(dissipation);
    HC_4_l->set_static_friction(staticFriction);
    HC_4_l->set_dynamic_friction(dynamicFriction);
    HC_4_l->set_viscous_friction(viscousFriction);
    HC_4_l->set_transition_velocity(transitionVelocity);
    HC_4_l->connectSocket_sphere_frame(*calcn_l);
    HC_4_l->connectSocket_half_space_frame(*treadmill);

    HC_5_l = new SmoothSphereHalfSpaceForce("sphere_5_l", *toes_l, *treadmill);
    HC_5_l->set_contact_sphere_location(locSphere_5_l);
    HC_5_l->set_contact_sphere_radius(radiusSphere_5);
    HC_5_l->set_contact_half_space_location(halfSpaceLocation);
    HC_5_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_5_l->set_stiffness(stiffness);
    HC_5_l->set_dissipation(dissipation);
    HC_5_l->set_static_friction(staticFriction);
    HC_5_l->set_dynamic_friction(dynamicFriction);
    HC_5_l->set_viscous_friction(viscousFriction);
    HC_5_l->set_transition_velocity(transitionVelocity);
    HC_5_l->connectSocket_sphere_frame(*toes_l);
    HC_5_l->connectSocket_half_space_frame(*treadmill);

    HC_6_l = new SmoothSphereHalfSpaceForce("sphere_6_l", *toes_l, *treadmill);
    HC_6_l->set_contact_sphere_location(locSphere_6_l);
    HC_6_l->set_contact_sphere_radius(radiusSphere_6);
    HC_6_l->set_contact_half_space_location(halfSpaceLocation);
    HC_6_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_6_l->set_stiffness(stiffness);
    HC_6_l->set_dissipation(dissipation);
    HC_6_l->set_static_friction(staticFriction);
    HC_6_l->set_dynamic_friction(dynamicFriction);
    HC_6_l->set_viscous_friction(viscousFriction);
    HC_6_l->set_transition_velocity(transitionVelocity);
    HC_6_l->connectSocket_sphere_frame(*toes_l);
    HC_6_l->connectSocket_half_space_frame(*treadmill);

    /// Add left foot contact spheres to model
    model->addComponent(HC_1_l);
    model->addComponent(HC_2_l);
    model->addComponent(HC_3_l);
    model->addComponent(HC_4_l);
    model->addComponent(HC_5_l);
    model->addComponent(HC_6_l);
    /// Right foot contact shere specifications
    HC_1_r = new SmoothSphereHalfSpaceForce("sphere_1_r", *calcn_r, *treadmill);
    HC_1_r->set_contact_sphere_location(locSphere_1_r);
    HC_1_r->set_contact_sphere_radius(radiusSphere_1);
    HC_1_r->set_contact_half_space_location(halfSpaceLocation);
    HC_1_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_1_r->set_stiffness(stiffness);
    HC_1_r->set_dissipation(dissipation);
    HC_1_r->set_static_friction(staticFriction);
    HC_1_r->set_dynamic_friction(dynamicFriction);
    HC_1_r->set_viscous_friction(viscousFriction);
    HC_1_r->set_transition_velocity(transitionVelocity);
    HC_1_r->connectSocket_sphere_frame(*calcn_r);
    HC_1_r->connectSocket_half_space_frame(*treadmill);

    HC_2_r = new SmoothSphereHalfSpaceForce("sphere_2_r", *calcn_r, *treadmill);
    HC_2_r->set_contact_sphere_location(locSphere_2_r);
    HC_2_r->set_contact_sphere_radius(radiusSphere_2);
    HC_2_r->set_contact_half_space_location(halfSpaceLocation);
    HC_2_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_2_r->set_stiffness(stiffness);
    HC_2_r->set_dissipation(dissipation);
    HC_2_r->set_static_friction(staticFriction);
    HC_2_r->set_dynamic_friction(dynamicFriction);
    HC_2_r->set_viscous_friction(viscousFriction);
    HC_2_r->set_transition_velocity(transitionVelocity);
    HC_2_r->connectSocket_sphere_frame(*calcn_r);
    HC_2_r->connectSocket_half_space_frame(*treadmill);

    HC_3_r = new SmoothSphereHalfSpaceForce("sphere_3_r", *calcn_r, *treadmill);
    HC_3_r->set_contact_sphere_location(locSphere_3_r);
    HC_3_r->set_contact_sphere_radius(radiusSphere_3);
    HC_3_r->set_contact_half_space_location(halfSpaceLocation);
    HC_3_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_3_r->set_stiffness(stiffness);
    HC_3_r->set_dissipation(dissipation);
    HC_3_r->set_static_friction(staticFriction);
    HC_3_r->set_dynamic_friction(dynamicFriction);
    HC_3_r->set_viscous_friction(viscousFriction);
    HC_3_r->set_transition_velocity(transitionVelocity);
    HC_3_r->connectSocket_sphere_frame(*calcn_r);
    HC_3_r->connectSocket_half_space_frame(*treadmill);

    HC_4_r = new SmoothSphereHalfSpaceForce("sphere_4_r", *calcn_r, *treadmill);
    HC_4_r->set_contact_sphere_location(locSphere_4_r);
    HC_4_r->set_contact_sphere_radius(radiusSphere_4);
    HC_4_r->set_contact_half_space_location(halfSpaceLocation);
    HC_4_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_4_r->set_stiffness(stiffness);
    HC_4_r->set_dissipation(dissipation);
    HC_4_r->set_static_friction(staticFriction);
    HC_4_r->set_dynamic_friction(dynamicFriction);
    HC_4_r->set_viscous_friction(viscousFriction);
    HC_4_r->set_transition_velocity(transitionVelocity);
    HC_4_r->connectSocket_sphere_frame(*calcn_r);
    HC_4_r->connectSocket_half_space_frame(*treadmill);

    HC_5_r = new SmoothSphereHalfSpaceForce("sphere_5_r", *toes_r, *treadmill);
    HC_5_r->set_contact_sphere_location(locSphere_5_r);
    HC_5_r->set_contact_sphere_radius(radiusSphere_5);
    HC_5_r->set_contact_half_space_location(halfSpaceLocation);
    HC_5_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_5_r->set_stiffness(stiffness);
    HC_5_r->set_dissipation(dissipation);
    HC_5_r->set_static_friction(staticFriction);
    HC_5_r->set_dynamic_friction(dynamicFriction);
    HC_5_r->set_viscous_friction(viscousFriction);
    HC_5_r->set_transition_velocity(transitionVelocity);
    HC_5_r->connectSocket_sphere_frame(*toes_r);
    HC_5_r->connectSocket_half_space_frame(*treadmill);

    HC_6_r = new SmoothSphereHalfSpaceForce("sphere_6_r", *toes_r, *treadmill);
    HC_6_r->set_contact_sphere_location(locSphere_6_r);
    HC_6_r->set_contact_sphere_radius(radiusSphere_6);
    HC_6_r->set_contact_half_space_location(halfSpaceLocation);
    HC_6_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_6_r->set_stiffness(stiffness);
    HC_6_r->set_dissipation(dissipation);
    HC_6_r->set_static_friction(staticFriction);
    HC_6_r->set_dynamic_friction(dynamicFriction);
    HC_6_r->set_viscous_friction(viscousFriction);
    HC_6_r->set_transition_velocity(transitionVelocity);
    HC_6_r->connectSocket_sphere_frame(*toes_r);
    HC_6_r->connectSocket_half_space_frame(*treadmill);

    /// Add right foot contact spheres to model
    model->addComponent(HC_1_r);
    model->addComponent(HC_2_r);
    model->addComponent(HC_3_r);
    model->addComponent(HC_4_r);
    model->addComponent(HC_5_r);
    model->addComponent(HC_6_r);

    /// OpenPose markers
    Neck = new Station(*torso, Vec3(0.0048747666408733847, 0.34228202149147169, -0.0019243593351921939));
    RShoulder = new Station(*torso, Vec3(0.0088341996212372154, 0.33298469487384974, 0.13313164879858663));
    LShoulder = new Station(*torso, Vec3(0.00091519173114446017, 0.35157948104656622, -0.13698038315231087));
    MidHip = new Station(*pelvis, Vec3(-0.069751070446270047, -0.066830276621380191, -0.0016527458804501682));
    RHip = new Station(*femur_r, Vec3(-0.0070326417095932792, -0.0010089117897051914, 0.012413416619825535));
    LHip = new Station(*femur_l, Vec3(0.0035089505683221223, -0.0054289827493331222, -0.015568962706259515));
    RKnee = new Station(*femur_r, Vec3(-0.0058124680124881367, -0.37856583725265325, 0.0085423690084551751));
    LKnee = new Station(*femur_l, Vec3(-0.0041224905170602621, -0.37537948153785006, -0.0032743144905940103));
    RAnkle = new Station(*tibia_r, Vec3(0.012131639524641025, -0.39579379282818894, -0.0040043196701304851));
    LAnkle = new Station(*tibia_l, Vec3(0.0085935346092369524, -0.39636413839556778, 0.0039811015371158387));
    RHeel = new Station(*calcn_r, Vec3(0.011226112052901249, 0.019587852073327437, -0.0034214978378686034));
    LHeel = new Station(*calcn_l, Vec3(0.0087513280974961249, 0.02122660957545356, 0.0017464646427649555));
    RBigToe = new Station(*toes_r, Vec3(0.04262457758589433, 0.0018780785029475638, -0.016949830751231909));
    LBigToe = new Station(*toes_l, Vec3(0.049023592065012833, -0.003273494125738248, 0.014207664476521664));
    RSmallToe = new Station(*toes_r, Vec3(-0.003469715746162394, 0.004107514024945462, 0.041646317573895697));
    LSmallToe = new Station(*toes_l, Vec3(-0.0013574633323293339, 0.0083213796177722121, -0.042873568766542691));

    model->addComponent(Neck);
    model->addComponent(RShoulder);
    model->addComponent(LShoulder);
    model->addComponent(MidHip);
    model->addComponent(RHip);
    model->addComponent(LHip);
    model->addComponent(RKnee);
    model->addComponent(LKnee);
    model->addComponent(RAnkle);
    model->addComponent(LAnkle);
    model->addComponent(RHeel);
    model->addComponent(LHeel);
    model->addComponent(RBigToe);
    model->addComponent(LBigToe);
    model->addComponent(RSmallToe);
    model->addComponent(LSmallToe);

	// Initialize system and state
	SimTK::State* state;
	state = new State(model->initSystem());

	// Read inputs
	std::vector<T> x(arg[0], arg[0] + NX);
	std::vector<T> u(arg[1], arg[1] + NU);

	// States and controls
	T ua[NU + 1]; /// joint accelerations (Qdotdots) - controls
	Vector QsUs(NX + 2); /// joint positions (Qs) and velocities (Us) - states

	// Assign inputs to model variables
	/// States
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
	// Warning: we impose the treadmill coordinate value and speed.
	// We know that treadmill is last in the state vector.
	QsUs[NX] = 0;
	QsUs[NX + 1] = -1.22; // Treadmill speed.
	/// Controls
	T ut[NU + 1];
	for (int i = 0; i < NU; ++i) ut[i] = u[i];
	// Warning: we impose the treadmill coordinate speed derivative.
	// We know that treadmill is last in the state vector.
	ut[NU] = 0;
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
	for (int i = 0; i < ndof_treadmill; ++i) ua[i] = ut[indicesOSInSimbody[i]];

	// Set state variables and realize
	model->setStateVariableValues(*state, QsUs);
	model->realizeVelocity(*state);

	// Compute residual forces
	/// appliedMobilityForces (# mobilities)
	Vector appliedMobilityForces(ndof_treadmill);
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
	/// Add contact forces to appliedBodyForces
	/// Right foot
	Array<osim_double_adouble> Force_values_1_r = HC_1_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_2_r = HC_2_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_3_r = HC_3_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_4_r = HC_4_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_5_r = HC_5_r->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_6_r = HC_6_r->getRecordValues(*state);
	SpatialVec GRF_1_r;
	GRF_1_r[0] = Vec3(Force_values_1_r[3], Force_values_1_r[4], Force_values_1_r[5]);
	GRF_1_r[1] = Vec3(Force_values_1_r[0], Force_values_1_r[1], Force_values_1_r[2]);
	SpatialVec GRF_2_r;
	GRF_2_r[0] = Vec3(Force_values_2_r[3], Force_values_2_r[4], Force_values_2_r[5]);
	GRF_2_r[1] = Vec3(Force_values_2_r[0], Force_values_2_r[1], Force_values_2_r[2]);
	SpatialVec GRF_3_r;
	GRF_3_r[0] = Vec3(Force_values_3_r[3], Force_values_3_r[4], Force_values_3_r[5]);
	GRF_3_r[1] = Vec3(Force_values_3_r[0], Force_values_3_r[1], Force_values_3_r[2]);
	SpatialVec GRF_4_r;
	GRF_4_r[0] = Vec3(Force_values_4_r[3], Force_values_4_r[4], Force_values_4_r[5]);
	GRF_4_r[1] = Vec3(Force_values_4_r[0], Force_values_4_r[1], Force_values_4_r[2]);
	SpatialVec GRF_5_r;
	GRF_5_r[0] = Vec3(Force_values_5_r[3], Force_values_5_r[4], Force_values_5_r[5]);
	GRF_5_r[1] = Vec3(Force_values_5_r[0], Force_values_5_r[1], Force_values_5_r[2]);
	SpatialVec GRF_6_r;
	GRF_6_r[0] = Vec3(Force_values_6_r[3], Force_values_6_r[4], Force_values_6_r[5]);
	GRF_6_r[1] = Vec3(Force_values_6_r[0], Force_values_6_r[1], Force_values_6_r[2]);
	int ncalcn_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
	int ntoes_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
	appliedBodyForces[ncalcn_r] = appliedBodyForces[ncalcn_r] + GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r;
	appliedBodyForces[ntoes_r] = appliedBodyForces[ntoes_r] + GRF_5_r + GRF_6_r;
	/// Left foot
	Array<osim_double_adouble> Force_values_1_l = HC_1_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_2_l = HC_2_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_3_l = HC_3_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_4_l = HC_4_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_5_l = HC_5_l->getRecordValues(*state);
	Array<osim_double_adouble> Force_values_6_l = HC_6_l->getRecordValues(*state);
	SpatialVec GRF_1_l;
	GRF_1_l[0] = Vec3(Force_values_1_l[3], Force_values_1_l[4], Force_values_1_l[5]);
	GRF_1_l[1] = Vec3(Force_values_1_l[0], Force_values_1_l[1], Force_values_1_l[2]);
	SpatialVec GRF_2_l;
	GRF_2_l[0] = Vec3(Force_values_2_l[3], Force_values_2_l[4], Force_values_2_l[5]);
	GRF_2_l[1] = Vec3(Force_values_2_l[0], Force_values_2_l[1], Force_values_2_l[2]);
	SpatialVec GRF_3_l;
	GRF_3_l[0] = Vec3(Force_values_3_l[3], Force_values_3_l[4], Force_values_3_l[5]);
	GRF_3_l[1] = Vec3(Force_values_3_l[0], Force_values_3_l[1], Force_values_3_l[2]);
	SpatialVec GRF_4_l;
	GRF_4_l[0] = Vec3(Force_values_4_l[3], Force_values_4_l[4], Force_values_4_l[5]);
	GRF_4_l[1] = Vec3(Force_values_4_l[0], Force_values_4_l[1], Force_values_4_l[2]);
	SpatialVec GRF_5_l;
	GRF_5_l[0] = Vec3(Force_values_5_l[3], Force_values_5_l[4], Force_values_5_l[5]);
	GRF_5_l[1] = Vec3(Force_values_5_l[0], Force_values_5_l[1], Force_values_5_l[2]);
	SpatialVec GRF_6_l;
	GRF_6_l[0] = Vec3(Force_values_6_l[3], Force_values_6_l[4], Force_values_6_l[5]);
	GRF_6_l[1] = Vec3(Force_values_6_l[0], Force_values_6_l[1], Force_values_6_l[2]);
	int ncalcn_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
	int ntoes_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
	appliedBodyForces[ncalcn_l] = appliedBodyForces[ncalcn_l] + GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l;
	appliedBodyForces[ntoes_l] = appliedBodyForces[ntoes_l] + GRF_5_l + GRF_6_l;
	/// knownUdot
	Vector knownUdot(ndof_treadmill);
	knownUdot.setToZero();
	for (int i = 0; i < ndof_treadmill; ++i) knownUdot[i] = ua[i];
	/// Calculate residual forces
	Vector residualMobilityForces(ndof_treadmill);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
		appliedMobilityForces, appliedBodyForces, knownUdot,
		residualMobilityForces);

	/// Get location OpenPose markers
	Vec3 Neck_location = Neck->getLocationInGround(*state);
	Vec3 RShoulder_location = RShoulder->getLocationInGround(*state);
	Vec3 LShoulder_location = LShoulder->getLocationInGround(*state);
	Vec3 MidHip_location = MidHip->getLocationInGround(*state);
	Vec3 RHip_location = RHip->getLocationInGround(*state);
	Vec3 LHip_location = LHip->getLocationInGround(*state);
	Vec3 RKnee_location = RKnee->getLocationInGround(*state);
	Vec3 LKnee_location = LKnee->getLocationInGround(*state);
	Vec3 RAnkle_location = RAnkle->getLocationInGround(*state);
	Vec3 LAnkle_location = LAnkle->getLocationInGround(*state);
	Vec3 RHeel_location = RHeel->getLocationInGround(*state);
	Vec3 LHeel_location = LHeel->getLocationInGround(*state);
	Vec3 RBigToe_location = RBigToe->getLocationInGround(*state);
	Vec3 LBigToe_location = LBigToe->getLocationInGround(*state);
	Vec3 RSmallToe_location = RSmallToe->getLocationInGround(*state);
	Vec3 LSmallToe_location = LSmallToe->getLocationInGround(*state);

	SpatialVec GRF_r = GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r + GRF_5_r + GRF_6_r;
	SpatialVec GRF_l = GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l + GRF_5_l + GRF_6_l;

	/// Extract ground reaction forces and moments: method 1.
	/*Vec3 calcn_or_l = calcn_l->getPositionInGround(*state);
	Vec3 calcn_or_r = calcn_r->getPositionInGround(*state);
	Vec3 toes_or_l = toes_l->getPositionInGround(*state);
	Vec3 toes_or_r = toes_r->getPositionInGround(*state);

	SpatialVec GRF_r = GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r + GRF_5_r + GRF_6_r;
	SpatialVec GRF_l = GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l + GRF_5_l + GRF_6_l;

	SpatialVec GRF_calcn_l = GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l;
	SpatialVec GRF_toes_l = GRF_5_l + GRF_6_l;
	SpatialVec GRF_calcn_r = GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r;
	SpatialVec GRF_toes_r = GRF_5_r + GRF_6_r;

	Vec3 GRM_calcn_l = GRF_calcn_l[0] + cross(calcn_or_l, GRF_calcn_l[1]);
	Vec3 GRM_toes_l = GRF_toes_l[0] + cross(toes_or_l, GRF_toes_l[1]);
	Vec3 GRM_calcn_r = GRF_calcn_r[0] + cross(calcn_or_r, GRF_calcn_r[1]);
	Vec3 GRM_toes_r = GRF_toes_r[0] + cross(toes_or_r, GRF_toes_r[1]);

	Vec3 GRM_l = GRM_calcn_l + GRM_toes_l;
	Vec3 GRM_r = GRM_calcn_r + GRM_toes_r;*/

	/// Extract ground reaction forces and moments: method 2.
	/// Compute contact point positions in body frames.
	Vec3 normal(0, 1, 0); // TODO
	/// sphere 1 left
	Vec3 pos_InGround_HC_s1_l = calcn_l->findStationLocationInGround(*state, locSphere_1_l);
	Vec3 contactPointpos_InGround_HC_s1_l = pos_InGround_HC_s1_l - radiusSphere_1 * normal;
	Vec3 contactPointpos_InGround_HC_s1_l_adj = contactPointpos_InGround_HC_s1_l - 0.5*contactPointpos_InGround_HC_s1_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s1_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_l_adj, *calcn_l);
	/// sphere 2 left
	Vec3 pos_InGround_HC_s2_l = calcn_l->findStationLocationInGround(*state, locSphere_2_l);
	Vec3 contactPointpos_InGround_HC_s2_l = pos_InGround_HC_s2_l - radiusSphere_2 * normal;
	Vec3 contactPointpos_InGround_HC_s2_l_adj = contactPointpos_InGround_HC_s2_l - 0.5*contactPointpos_InGround_HC_s2_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s2_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_l_adj, *calcn_l);
	/// sphere 3 left
	Vec3 pos_InGround_HC_s3_l = calcn_l->findStationLocationInGround(*state, locSphere_3_l);
	Vec3 contactPointpos_InGround_HC_s3_l = pos_InGround_HC_s3_l - radiusSphere_3 * normal;
	Vec3 contactPointpos_InGround_HC_s3_l_adj = contactPointpos_InGround_HC_s3_l - 0.5*contactPointpos_InGround_HC_s3_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s3_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_l_adj, *calcn_l);
	/// sphere 4 left
	Vec3 pos_InGround_HC_s4_l = calcn_l->findStationLocationInGround(*state, locSphere_4_l);
	Vec3 contactPointpos_InGround_HC_s4_l = pos_InGround_HC_s4_l - radiusSphere_4 * normal;
	Vec3 contactPointpos_InGround_HC_s4_l_adj = contactPointpos_InGround_HC_s4_l - 0.5*contactPointpos_InGround_HC_s4_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s4_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_l_adj, *calcn_l);
	/// sphere 5 left
	Vec3 pos_InGround_HC_s5_l = toes_l->findStationLocationInGround(*state, locSphere_5_l);
	Vec3 contactPointpos_InGround_HC_s5_l = pos_InGround_HC_s5_l - radiusSphere_5 * normal;
	Vec3 contactPointpos_InGround_HC_s5_l_adj = contactPointpos_InGround_HC_s5_l - 0.5*contactPointpos_InGround_HC_s5_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s5_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_l_adj, *toes_l);
	/// sphere 6 left
	Vec3 pos_InGround_HC_s6_l = toes_l->findStationLocationInGround(*state, locSphere_6_l);
	Vec3 contactPointpos_InGround_HC_s6_l = pos_InGround_HC_s6_l - radiusSphere_6 * normal;
	Vec3 contactPointpos_InGround_HC_s6_l_adj = contactPointpos_InGround_HC_s6_l - 0.5*contactPointpos_InGround_HC_s6_l[1] * normal;
	Vec3 contactPointPos_InBody_HC_s6_l = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_l_adj, *toes_l);
	/// sphere 1 right
	Vec3 pos_InGround_HC_s1_r = calcn_r->findStationLocationInGround(*state, locSphere_1_r);
	Vec3 contactPointpos_InGround_HC_s1_r = pos_InGround_HC_s1_r - radiusSphere_1 * normal;
	Vec3 contactPointpos_InGround_HC_s1_r_adj = contactPointpos_InGround_HC_s1_r - 0.5*contactPointpos_InGround_HC_s1_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s1_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s1_r_adj, *calcn_r);
	/// sphere 2 right
	Vec3 pos_InGround_HC_s2_r = calcn_r->findStationLocationInGround(*state, locSphere_2_r);
	Vec3 contactPointpos_InGround_HC_s2_r = pos_InGround_HC_s2_r - radiusSphere_2 * normal;
	Vec3 contactPointpos_InGround_HC_s2_r_adj = contactPointpos_InGround_HC_s2_r - 0.5*contactPointpos_InGround_HC_s2_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s2_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s2_r_adj, *calcn_r);
	/// sphere 3 right
	Vec3 pos_InGround_HC_s3_r = calcn_r->findStationLocationInGround(*state, locSphere_3_r);
	Vec3 contactPointpos_InGround_HC_s3_r = pos_InGround_HC_s3_r - radiusSphere_3 * normal;
	Vec3 contactPointpos_InGround_HC_s3_r_adj = contactPointpos_InGround_HC_s3_r - 0.5*contactPointpos_InGround_HC_s3_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s3_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s3_r_adj, *calcn_r);
	/// sphere 4 right
	Vec3 pos_InGround_HC_s4_r = calcn_r->findStationLocationInGround(*state, locSphere_4_r);
	Vec3 contactPointpos_InGround_HC_s4_r = pos_InGround_HC_s4_r - radiusSphere_4 * normal;
	Vec3 contactPointpos_InGround_HC_s4_r_adj = contactPointpos_InGround_HC_s4_r - 0.5*contactPointpos_InGround_HC_s4_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s4_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s4_r_adj, *calcn_r);
	/// sphere 5 right
	Vec3 pos_InGround_HC_s5_r = toes_r->findStationLocationInGround(*state, locSphere_5_r);
	Vec3 contactPointpos_InGround_HC_s5_r = pos_InGround_HC_s5_r - radiusSphere_5 * normal;
	Vec3 contactPointpos_InGround_HC_s5_r_adj = contactPointpos_InGround_HC_s5_r - 0.5*contactPointpos_InGround_HC_s5_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s5_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s5_r_adj, *toes_r);
	/// sphere 6 right
	Vec3 pos_InGround_HC_s6_r = toes_r->findStationLocationInGround(*state, locSphere_6_r);
	Vec3 contactPointpos_InGround_HC_s6_r = pos_InGround_HC_s6_r - radiusSphere_6 * normal;
	Vec3 contactPointpos_InGround_HC_s6_r_adj = contactPointpos_InGround_HC_s6_r - 0.5*contactPointpos_InGround_HC_s6_r[1] * normal;
	Vec3 contactPointPos_InBody_HC_s6_r = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_HC_s6_r_adj, *toes_r);
	/// Get transforms
	SimTK::Transform TR_GB_calcn_l = calcn_l->getMobilizedBody().getBodyTransform(*state);
	SimTK::Transform TR_GB_calcn_r = calcn_r->getMobilizedBody().getBodyTransform(*state);
	SimTK::Transform TR_GB_toes_l = toes_l->getMobilizedBody().getBodyTransform(*state);
	SimTK::Transform TR_GB_toes_r = toes_r->getMobilizedBody().getBodyTransform(*state);
	/// Calculate torques about ground origin
	Vec3 AppliedPointTorque_s1_l, AppliedPointTorque_s2_l, AppliedPointTorque_s3_l, AppliedPointTorque_s4_l, AppliedPointTorque_s5_l, AppliedPointTorque_s6_l;
	Vec3 AppliedPointTorque_s1_r, AppliedPointTorque_s2_r, AppliedPointTorque_s3_r, AppliedPointTorque_s4_r, AppliedPointTorque_s5_r, AppliedPointTorque_s6_r;
	AppliedPointTorque_s1_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s1_l) % GRF_1_l[1];
	AppliedPointTorque_s2_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s2_l) % GRF_2_l[1];
	AppliedPointTorque_s3_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s3_l) % GRF_3_l[1];
	AppliedPointTorque_s4_l = (TR_GB_calcn_l*contactPointPos_InBody_HC_s4_l) % GRF_4_l[1];
	AppliedPointTorque_s5_l = (TR_GB_toes_l*contactPointPos_InBody_HC_s5_l) % GRF_5_l[1];
	AppliedPointTorque_s6_l = (TR_GB_toes_l*contactPointPos_InBody_HC_s6_l) % GRF_6_l[1];
	AppliedPointTorque_s1_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s1_r) % GRF_1_r[1];
	AppliedPointTorque_s2_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s2_r) % GRF_2_r[1];
	AppliedPointTorque_s3_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s3_r) % GRF_3_r[1];
	AppliedPointTorque_s4_r = (TR_GB_calcn_r*contactPointPos_InBody_HC_s4_r) % GRF_4_r[1];
	AppliedPointTorque_s5_r = (TR_GB_toes_r*contactPointPos_InBody_HC_s5_r) % GRF_5_r[1];
	AppliedPointTorque_s6_r = (TR_GB_toes_r*contactPointPos_InBody_HC_s6_r) % GRF_6_r[1];
	/// Sum up torques
	Vec3 GRM2_l, GRM2_r;
	GRM2_l = AppliedPointTorque_s1_l + AppliedPointTorque_s2_l + AppliedPointTorque_s3_l + AppliedPointTorque_s4_l + AppliedPointTorque_s5_l + AppliedPointTorque_s6_l;
	GRM2_r = AppliedPointTorque_s1_r + AppliedPointTorque_s2_r + AppliedPointTorque_s3_r + AppliedPointTorque_s4_r + AppliedPointTorque_s5_r + AppliedPointTorque_s6_r;

	// Residual forces in OpenSim order
	T res_os[ndof];
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
	for (int i = 0; i < ndof; ++i) res_os[i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
	for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];

	int nc = 3;
	for (int i = 0; i < nc; ++i) res[0][i + NU + 0 * nc] = value<T>(Neck_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 1 * nc] = value<T>(RShoulder_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 2 * nc] = value<T>(LShoulder_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 3 * nc] = value<T>(MidHip_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 4 * nc] = value<T>(RHip_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 5 * nc] = value<T>(LHip_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 6 * nc] = value<T>(RKnee_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 7 * nc] = value<T>(LKnee_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 8 * nc] = value<T>(RAnkle_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 9 * nc] = value<T>(LAnkle_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 10 * nc] = value<T>(RHeel_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 11 * nc] = value<T>(LHeel_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 12 * nc] = value<T>(RBigToe_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 13 * nc] = value<T>(LBigToe_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 14 * nc] = value<T>(RSmallToe_location[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 15 * nc] = value<T>(LSmallToe_location[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 16 * nc] = value<T>(GRF_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 17 * nc] = value<T>(GRF_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 18 * nc] = value<T>(GRM2_r[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 19 * nc] = value<T>(GRM2_l[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 20 * nc] = value<T>(GRF_1_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 21 * nc] = value<T>(contactPointpos_InGround_HC_s1_r_adj[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 22 * nc] = value<T>(GRF_1_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 23 * nc] = value<T>(contactPointpos_InGround_HC_s1_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 24 * nc] = value<T>(GRF_2_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 25 * nc] = value<T>(contactPointpos_InGround_HC_s2_r_adj[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 26 * nc] = value<T>(GRF_2_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 27 * nc] = value<T>(contactPointpos_InGround_HC_s2_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 28 * nc] = value<T>(GRF_3_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 29 * nc] = value<T>(contactPointpos_InGround_HC_s3_r_adj[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 30 * nc] = value<T>(GRF_3_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 31 * nc] = value<T>(contactPointpos_InGround_HC_s3_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 32 * nc] = value<T>(GRF_4_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 33 * nc] = value<T>(contactPointpos_InGround_HC_s4_r_adj[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 34 * nc] = value<T>(GRF_4_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 35 * nc] = value<T>(contactPointpos_InGround_HC_s4_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 36 * nc] = value<T>(GRF_5_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 37 * nc] = value<T>(contactPointpos_InGround_HC_s5_r_adj[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 38 * nc] = value<T>(GRF_5_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 39 * nc] = value<T>(contactPointpos_InGround_HC_s5_l_adj[i]);

	for (int i = 0; i < nc; ++i) res[0][i + NU + 40 * nc] = value<T>(GRF_6_r[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 41 * nc] = value<T>(contactPointpos_InGround_HC_s6_r_adj[i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 42 * nc] = value<T>(GRF_6_l[1][i]);
	for (int i = 0; i < nc; ++i) res[0][i + NU + 43 * nc] = value<T>(contactPointpos_InGround_HC_s6_l_adj[i]);

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

	for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;

	const Recorder* Recorder_arg[n_in] = { x,u };
	Recorder* Recorder_res[n_out] = { tau };

	F_generic<Recorder>(Recorder_arg, Recorder_res);

	double res[NR];
	for (int i = 0; i < NR; ++i) {
		Recorder_res[0][i] >>= res[i];
		//std::cout << Recorder_res[0][i] << std::endl;
	}

	Recorder::stop_recording();

	return 0;

}
