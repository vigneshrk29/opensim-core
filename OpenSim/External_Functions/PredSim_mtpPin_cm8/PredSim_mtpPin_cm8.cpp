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
//#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Common/SimmSpline.h>
//#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
//#include <OpenSim/Simulation/Model/MovingPathPoint.h>
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
constexpr int ndof = 31;        // # degrees of freedom (excluding locked)
constexpr int ndofr = ndof+2;   // # degrees of freedom (including locked)
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NR = ndof+5*4;    // # residual torques + # joint origins

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
    OpenSim::Body* humerus_r;
    OpenSim::Body* humerus_l;
    OpenSim::Body* ulna_r;
    OpenSim::Body* ulna_l;
    OpenSim::Body* radius_r;
    OpenSim::Body* radius_l;
    OpenSim::Body* hand_r;
    OpenSim::Body* hand_l;
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
    OpenSim::CustomJoint* shoulder_r;
    OpenSim::CustomJoint* shoulder_l;
    OpenSim::CustomJoint* elbow_r;
    OpenSim::CustomJoint* elbow_l;
    OpenSim::CustomJoint* radioulnar_r;
    OpenSim::CustomJoint* radioulnar_l;
    OpenSim::WeldJoint* radius_hand_r;
    OpenSim::WeldJoint* radius_hand_l;
    /// Contact elements
    OpenSim::SmoothSphereHalfSpaceForce* HC_1_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_2_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_3_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_4_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_5_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_6_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_7_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_8_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_9_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_10_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_11_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_12_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_13_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_14_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_15_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_16_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_17_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_18_r;
    OpenSim::SmoothSphereHalfSpaceForce* HC_1_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_2_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_3_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_4_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_5_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_6_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_7_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_8_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_9_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_10_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_11_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_12_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_13_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_14_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_15_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_16_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_17_l;
    OpenSim::SmoothSphereHalfSpaceForce* HC_18_l;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
    pelvis = new OpenSim::Body("pelvis", 8.8425916618972398, Vec3(-0.0682778, 0, 0), Inertia(0.074179899999999993, 0.074179899999999993, 0.040545600000000001, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 6.9838228822256099, Vec3(0, -0.17046700000000001, 0), Inertia(0.10109, 0.026499200000000001, 0.106601, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 6.9838228822256099, Vec3(0, -0.17046700000000001, 0), Inertia(0.10109, 0.026499200000000001, 0.106601, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 2.7837232390663198, Vec3(0, -0.18048900000000001, 0), Inertia(0.035366099999999998, 0.0035787200000000001, 0.035857300000000002, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 2.7837232390663198, Vec3(0, -0.18048900000000001, 0), Inertia(0.035366099999999998, 0.0035787200000000001, 0.035857300000000002, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.075083566798821805, Vec3(0, 0, 0), Inertia(0.000627141, 0.000627141, 0.000627141, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.075083566798821805, Vec3(0, 0, 0), Inertia(0.000627141, 0.000627141, 0.000627141, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 0.93854458498527304, Vec3(0.091392399999999999, 0.0274177, 0), Inertia(0.00087799799999999995, 0.0024458499999999998, 0.00257128, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 0.93854458498527304, Vec3(0.091392399999999999, 0.0274177, 0), Inertia(0.00087799799999999995, 0.0024458499999999998, 0.00257128, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.162631005686248, Vec3(0.031621799999999999, 0.0054835500000000002, 0.0159937), Inertia(6.2714100000000003e-05, 0.00012542799999999999, 6.2714100000000003e-05, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.162631005686248, Vec3(0.031621799999999999, 0.0054835500000000002, -0.0159937), Inertia(6.2714100000000003e-05, 0.00012542799999999999, 6.2714100000000003e-05, 0, 0, 0));
    torso = new OpenSim::Body("torso", 25.706060430645401, Vec3(-0.026760300000000001, 0.30650500000000003, 0), Inertia(0.98116599999999998, 0.45135399999999998, 0.98116599999999998, 0, 0, 0));
    humerus_l = new OpenSim::Body("humerus_l", 1.5261185453261299, Vec3(0, -0.16903299999999999, 0), Inertia(0.0094704400000000001, 0.0032670099999999999, 0.010630300000000001, 0, 0, 0));
    humerus_r = new OpenSim::Body("humerus_r", 1.5261185453261299, Vec3(0, -0.16903299999999999, 0), Inertia(0.0094704400000000001, 0.0032670099999999999, 0.010630300000000001, 0, 0, 0));
    ulna_l = new OpenSim::Body("ulna_l", 0.45613266830284299, Vec3(0, -0.11827500000000001, 0), Inertia(0.0021417200000000002, 0.000446854, 0.00232321, 0, 0, 0));
    ulna_r = new OpenSim::Body("ulna_r", 0.45613266830284299, Vec3(0, -0.11827500000000001, 0), Inertia(0.0021417200000000002, 0.000446854, 0.00232321, 0, 0, 0));
    radius_l = new OpenSim::Body("radius_l", 0.456132668302843, Vec3(0, -0.11827500000000001, 0), Inertia(0.0021417200000000002, 0.000446854, 0.00232321, 0, 0, 0));
    radius_r = new OpenSim::Body("radius_r", 0.456132668302843, Vec3(0, -0.11827500000000001, 0), Inertia(0.0021417200000000002, 0.000446854, 0.00232321, 0, 0, 0));
    hand_l = new OpenSim::Body("hand_l", 0.34350731810460999, Vec3(0, -0.066823900000000006, 0), Inertia(0.00064497399999999998, 0.00039551700000000001, 0.00096890799999999999, 0, 0, 0));
    hand_r = new OpenSim::Body("hand_r", 0.34350731810460999, Vec3(0, -0.066823900000000006, 0), Inertia(0.00064497399999999998, 0.00039551700000000001, 0.00096890799999999999, 0, 0, 0));
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
    st_knee_l[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    st_knee_l[2].setFunction(new LinearFunction());
    st_knee_l[2].setAxis(Vec3(0, 0, 1));
    /// Knee_r transform
    SpatialTransform st_knee_r;
    st_knee_r[2].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[2].setFunction(new LinearFunction());
    st_knee_r[2].setAxis(Vec3(0, 0, 1));
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
    /// Shoulder_l transform
    SpatialTransform st_sho_l;
    st_sho_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
    st_sho_l[0].setFunction(new LinearFunction());
    st_sho_l[0].setAxis(Vec3(0, 0, 1));
    st_sho_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
    st_sho_l[1].setFunction(new LinearFunction());
    st_sho_l[1].setAxis(Vec3(-1, 0, 0));
    st_sho_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
    st_sho_l[2].setFunction(new LinearFunction());
    st_sho_l[2].setAxis(Vec3(0, -1, 0));
    /// Shoulder_r transform
    SpatialTransform st_sho_r;
    st_sho_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
    st_sho_r[0].setFunction(new LinearFunction());
    st_sho_r[0].setAxis(Vec3(0, 0, 1));
    st_sho_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
    st_sho_r[1].setFunction(new LinearFunction());
    st_sho_r[1].setAxis(Vec3(1, 0, 0));
    st_sho_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
    st_sho_r[2].setFunction(new LinearFunction());
    st_sho_r[2].setAxis(Vec3(0, 1, 0));
    /// Elbow_l transform
    SpatialTransform st_elb_l;
    st_elb_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
    st_elb_l[0].setFunction(new LinearFunction());
    st_elb_l[0].setAxis(Vec3(-0.22604696, -0.022269, 0.97386183));
    /// Elbow_r transform
    SpatialTransform st_elb_r;
    st_elb_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
    st_elb_r[0].setFunction(new LinearFunction());
    st_elb_r[0].setAxis(Vec3(0.22604696, 0.022269, 0.97386183));
    /// Radioulnar_l transform
    SpatialTransform st_radioulnar_l;
    st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
    st_radioulnar_l[0].setFunction(new LinearFunction());
    st_radioulnar_l[0].setAxis(Vec3(-0.05639803, -0.99840646, 0.001952));
    /// Radioulnar_r transform
    SpatialTransform st_radioulnar_r;
    st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
    st_radioulnar_r[0].setFunction(new LinearFunction());
    st_radioulnar_r[0].setAxis(Vec3(0.05639803, 0.99840646, 0.001952));
    /// Joint specifications
    ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.0682778, -0.0638354, -0.082330700000000007), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.0682778, -0.0638354, 0.082330700000000007), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *femur_l, Vec3(-0.00451221, -0.39690700000000001, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(-0.00451221, -0.39690700000000001, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.41569499999999998, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.41569499999999998, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
    subtalar_l = new CustomJoint("subtalar_l", *talus_l, Vec3(-0.044572100000000003, -0.038339100000000001, -0.0072382799999999997), Vec3(0), *calcn_l, Vec3(0), Vec3(0),st_subtalar_l);
    subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.044572100000000003, -0.038339100000000001, 0.0072382799999999997), Vec3(0), *calcn_r, Vec3(0), Vec3(0),st_subtalar_r);
    mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.16341, -0.00182785 ,-0.00098703799999999998), Vec3(0), *toes_l, Vec3(0), Vec3(0));
    mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.16341, -0.00182785, 0.00098703799999999998), Vec3(0), *toes_r, Vec3(0), Vec3(0));
    back = new CustomJoint("back", *pelvis, Vec3(-0.097250000000000003, 0.078707799999999994, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
    shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.0028142900000000001, 0.35583300000000001, -0.151642), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
    shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.0028142900000000001, 0.35583300000000001, 0.151642), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
    elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.0135061, -0.294159, 0.0098593099999999996), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
    elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.0135061, -0.294159, -0.0098593099999999996), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
    radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.0066014300000000001, -0.0127642, -0.0255961), Vec3(0), *radius_l, Vec3(0), Vec3(0),st_radioulnar_l);
    radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.0066014300000000001, -0.0127642, 0.0255961), Vec3(0), *radius_r, Vec3(0), Vec3(0),st_radioulnar_r);
    radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.0086327899999999996, -0.23143900000000001, -0.0133559), Vec3(0), *hand_l, Vec3(0), Vec3(0));
    radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.0086327899999999996, -0.23143900000000001, 0.0133559), Vec3(0), *hand_r, Vec3(0), Vec3(0));
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
    model->addBody(humerus_l);      model->addJoint(shoulder_l);
    model->addBody(humerus_r);      model->addJoint(shoulder_r);
    model->addBody(ulna_l);         model->addJoint(elbow_l);
    model->addBody(ulna_r);         model->addJoint(elbow_r);
    model->addBody(radius_l);       model->addJoint(radioulnar_l);
    model->addBody(radius_r);       model->addJoint(radioulnar_r);
    model->addBody(hand_l);         model->addJoint(radius_hand_l);
    model->addBody(hand_r);         model->addJoint(radius_hand_r);
    /// Contact elements    
    /// Parameters
    osim_double_adouble radiusSphere = 0.005;
    osim_double_adouble stiffness = 500000;
    osim_double_adouble dissipation = 2.0;
    osim_double_adouble staticFriction = 0.8;
    osim_double_adouble dynamicFriction = 0.8;
    osim_double_adouble viscousFriction = 0.5;
    osim_double_adouble transitionVelocity = 0.2;
    Vec3 halfSpaceLocation(0);
    Vec3 halfSpaceOrientation(0, 0, -0.5*SimTK::Pi);
    Vec3 locSphere_1_r(0.005  , -0.015, -0.005);
    Vec3 locSphere_2_r(0.02   , -0.015, 0.005 );
    Vec3 locSphere_3_r(0.035  , -0.015, 0.01  );
    Vec3 locSphere_4_r(0.02   , -0.015, -0.015);
    Vec3 locSphere_5_r(0.04   , -0.015, -0.015);
    Vec3 locSphere_6_r(0.05   , -0.015, 0     );
    Vec3 locSphere_7_r(0.06   , -0.015, 0.02  );
    Vec3 locSphere_8_r(0.09   , -0.015, -0.01 );
    Vec3 locSphere_9_r(0.11   , -0.015, 0.03  );
    Vec3 locSphere_10_r(0.14  , -0.015, -0.02 );
    Vec3 locSphere_11_r(-0.02 , -0.01, 0.04   );
    Vec3 locSphere_12_r(-0.01 , -0.01, 0.025  );
    Vec3 locSphere_13_r(-0.005, -0.01, 0.015  );
    Vec3 locSphere_14_r(0     , -0.01, 0      );
    Vec3 locSphere_15_r(0.015 , -0.01, -0.02  );
    Vec3 locSphere_16_r(0.005 , -0.01, 0.045  );
    Vec3 locSphere_17_r(0.035 , -0.01, 0.025  );
    Vec3 locSphere_18_r(0.05  , -0.01, -0.015 );
    Vec3 locSphere_1_l(locSphere_1_r[0], locSphere_1_r[1], -locSphere_1_r[2]);
    Vec3 locSphere_2_l(locSphere_2_r[0], locSphere_2_r[1], -locSphere_2_r[2]);
    Vec3 locSphere_3_l(locSphere_3_r[0], locSphere_3_r[1], -locSphere_3_r[2]);
    Vec3 locSphere_4_l(locSphere_4_r[0], locSphere_4_r[1], -locSphere_4_r[2]);
    Vec3 locSphere_5_l(locSphere_5_r[0], locSphere_5_r[1], -locSphere_5_r[2]);
    Vec3 locSphere_6_l(locSphere_6_r[0], locSphere_6_r[1], -locSphere_6_r[2]);
    Vec3 locSphere_7_l(locSphere_7_r[0], locSphere_7_r[1], -locSphere_7_r[2]);
    Vec3 locSphere_8_l(locSphere_8_r[0], locSphere_8_r[1], -locSphere_8_r[2]);
    Vec3 locSphere_9_l(locSphere_9_r[0], locSphere_9_r[1], -locSphere_9_r[2]);
    Vec3 locSphere_10_l(locSphere_10_r[0], locSphere_10_r[1], -locSphere_10_r[2]);
    Vec3 locSphere_11_l(locSphere_11_r[0], locSphere_11_r[1], -locSphere_11_r[2]);
    Vec3 locSphere_12_l(locSphere_12_r[0], locSphere_12_r[1], -locSphere_12_r[2]);
    Vec3 locSphere_13_l(locSphere_13_r[0], locSphere_13_r[1], -locSphere_13_r[2]);
    Vec3 locSphere_14_l(locSphere_14_r[0], locSphere_14_r[1], -locSphere_14_r[2]);
    Vec3 locSphere_15_l(locSphere_15_r[0], locSphere_15_r[1], -locSphere_15_r[2]);
    Vec3 locSphere_16_l(locSphere_16_r[0], locSphere_16_r[1], -locSphere_16_r[2]);
    Vec3 locSphere_17_l(locSphere_17_r[0], locSphere_17_r[1], -locSphere_17_r[2]);
    Vec3 locSphere_18_l(locSphere_18_r[0], locSphere_18_r[1], -locSphere_18_r[2]);


    /// Left foot contact shere specifications
    HC_1_l = new SmoothSphereHalfSpaceForce("sphere_1_l", *calcn_l, model->getGround());
    HC_1_l->set_contact_sphere_location(locSphere_1_l);
    HC_1_l->set_contact_sphere_radius(radiusSphere);
    HC_1_l->set_contact_half_space_location(halfSpaceLocation);
    HC_1_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_1_l->set_stiffness(stiffness);
    HC_1_l->set_dissipation(dissipation);
    HC_1_l->set_static_friction(staticFriction);
    HC_1_l->set_dynamic_friction(dynamicFriction);
    HC_1_l->set_viscous_friction(viscousFriction);
    HC_1_l->set_transition_velocity(transitionVelocity);
    HC_1_l->connectSocket_sphere_frame(*calcn_l);
    HC_1_l->connectSocket_half_space_frame(model->getGround());

    HC_2_l = new SmoothSphereHalfSpaceForce("sphere_2_l", *calcn_l, model->getGround());
    HC_2_l->set_contact_sphere_location(locSphere_2_l);
    HC_2_l->set_contact_sphere_radius(radiusSphere);
    HC_2_l->set_contact_half_space_location(halfSpaceLocation);
    HC_2_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_2_l->set_stiffness(stiffness);
    HC_2_l->set_dissipation(dissipation);
    HC_2_l->set_static_friction(staticFriction);
    HC_2_l->set_dynamic_friction(dynamicFriction);
    HC_2_l->set_viscous_friction(viscousFriction);
    HC_2_l->set_transition_velocity(transitionVelocity);
    HC_2_l->connectSocket_sphere_frame(*calcn_l);
    HC_2_l->connectSocket_half_space_frame(model->getGround());

    HC_3_l = new SmoothSphereHalfSpaceForce("sphere_3_l", *calcn_l, model->getGround());
    HC_3_l->set_contact_sphere_location(locSphere_3_l);
    HC_3_l->set_contact_sphere_radius(radiusSphere);
    HC_3_l->set_contact_half_space_location(halfSpaceLocation);
    HC_3_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_3_l->set_stiffness(stiffness);
    HC_3_l->set_dissipation(dissipation);
    HC_3_l->set_static_friction(staticFriction);
    HC_3_l->set_dynamic_friction(dynamicFriction);
    HC_3_l->set_viscous_friction(viscousFriction);
    HC_3_l->set_transition_velocity(transitionVelocity);
    HC_3_l->connectSocket_sphere_frame(*calcn_l);
    HC_3_l->connectSocket_half_space_frame(model->getGround());

    HC_4_l = new SmoothSphereHalfSpaceForce("sphere_4_l", *calcn_l, model->getGround());
    HC_4_l->set_contact_sphere_location(locSphere_4_l);
    HC_4_l->set_contact_sphere_radius(radiusSphere);
    HC_4_l->set_contact_half_space_location(halfSpaceLocation);
    HC_4_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_4_l->set_stiffness(stiffness);
    HC_4_l->set_dissipation(dissipation);
    HC_4_l->set_static_friction(staticFriction);
    HC_4_l->set_dynamic_friction(dynamicFriction);
    HC_4_l->set_viscous_friction(viscousFriction);
    HC_4_l->set_transition_velocity(transitionVelocity);
    HC_4_l->connectSocket_sphere_frame(*calcn_l);
    HC_4_l->connectSocket_half_space_frame(model->getGround());

    HC_5_l = new SmoothSphereHalfSpaceForce("sphere_5_l", *calcn_l, model->getGround());
    HC_5_l->set_contact_sphere_location(locSphere_5_l);
    HC_5_l->set_contact_sphere_radius(radiusSphere);
    HC_5_l->set_contact_half_space_location(halfSpaceLocation);
    HC_5_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_5_l->set_stiffness(stiffness);
    HC_5_l->set_dissipation(dissipation);
    HC_5_l->set_static_friction(staticFriction);
    HC_5_l->set_dynamic_friction(dynamicFriction);
    HC_5_l->set_viscous_friction(viscousFriction);
    HC_5_l->set_transition_velocity(transitionVelocity);
    HC_5_l->connectSocket_sphere_frame(*calcn_l);
    HC_5_l->connectSocket_half_space_frame(model->getGround());

    HC_6_l = new SmoothSphereHalfSpaceForce("sphere_6_l", *calcn_l, model->getGround());
    HC_6_l->set_contact_sphere_location(locSphere_6_l);
    HC_6_l->set_contact_sphere_radius(radiusSphere);
    HC_6_l->set_contact_half_space_location(halfSpaceLocation);
    HC_6_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_6_l->set_stiffness(stiffness);
    HC_6_l->set_dissipation(dissipation);
    HC_6_l->set_static_friction(staticFriction);
    HC_6_l->set_dynamic_friction(dynamicFriction);
    HC_6_l->set_viscous_friction(viscousFriction);
    HC_6_l->set_transition_velocity(transitionVelocity);
    HC_6_l->connectSocket_sphere_frame(*calcn_l);
    HC_6_l->connectSocket_half_space_frame(model->getGround());

    HC_7_l = new SmoothSphereHalfSpaceForce("sphere_7_l", *calcn_l, model->getGround());
    HC_7_l->set_contact_sphere_location(locSphere_7_l);
    HC_7_l->set_contact_sphere_radius(radiusSphere);
    HC_7_l->set_contact_half_space_location(halfSpaceLocation);
    HC_7_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_7_l->set_stiffness(stiffness);
    HC_7_l->set_dissipation(dissipation);
    HC_7_l->set_static_friction(staticFriction);
    HC_7_l->set_dynamic_friction(dynamicFriction);
    HC_7_l->set_viscous_friction(viscousFriction);
    HC_7_l->set_transition_velocity(transitionVelocity);
    HC_7_l->connectSocket_sphere_frame(*calcn_l);
    HC_7_l->connectSocket_half_space_frame(model->getGround());

    HC_8_l = new SmoothSphereHalfSpaceForce("sphere_8_l", *calcn_l, model->getGround());
    HC_8_l->set_contact_sphere_location(locSphere_8_l);
    HC_8_l->set_contact_sphere_radius(radiusSphere);
    HC_8_l->set_contact_half_space_location(halfSpaceLocation);
    HC_8_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_8_l->set_stiffness(stiffness);
    HC_8_l->set_dissipation(dissipation);
    HC_8_l->set_static_friction(staticFriction);
    HC_8_l->set_dynamic_friction(dynamicFriction);
    HC_8_l->set_viscous_friction(viscousFriction);
    HC_8_l->set_transition_velocity(transitionVelocity);
    HC_8_l->connectSocket_sphere_frame(*calcn_l);
    HC_8_l->connectSocket_half_space_frame(model->getGround());

    HC_9_l = new SmoothSphereHalfSpaceForce("sphere_9_l", *calcn_l, model->getGround());
    HC_9_l->set_contact_sphere_location(locSphere_9_l);
    HC_9_l->set_contact_sphere_radius(radiusSphere);
    HC_9_l->set_contact_half_space_location(halfSpaceLocation);
    HC_9_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_9_l->set_stiffness(stiffness);
    HC_9_l->set_dissipation(dissipation);
    HC_9_l->set_static_friction(staticFriction);
    HC_9_l->set_dynamic_friction(dynamicFriction);
    HC_9_l->set_viscous_friction(viscousFriction);
    HC_9_l->set_transition_velocity(transitionVelocity);
    HC_9_l->connectSocket_sphere_frame(*calcn_l);
    HC_9_l->connectSocket_half_space_frame(model->getGround());

    HC_10_l = new SmoothSphereHalfSpaceForce("sphere_10_l", *calcn_l, model->getGround());
    HC_10_l->set_contact_sphere_location(locSphere_10_l);
    HC_10_l->set_contact_sphere_radius(radiusSphere);
    HC_10_l->set_contact_half_space_location(halfSpaceLocation);
    HC_10_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_10_l->set_stiffness(stiffness);
    HC_10_l->set_dissipation(dissipation);
    HC_10_l->set_static_friction(staticFriction);
    HC_10_l->set_dynamic_friction(dynamicFriction);
    HC_10_l->set_viscous_friction(viscousFriction);
    HC_10_l->set_transition_velocity(transitionVelocity);
    HC_10_l->connectSocket_sphere_frame(*calcn_l);
    HC_10_l->connectSocket_half_space_frame(model->getGround());

    HC_11_l = new SmoothSphereHalfSpaceForce("sphere_11_l", *toes_l, model->getGround());
    HC_11_l->set_contact_sphere_location(locSphere_11_l);
    HC_11_l->set_contact_sphere_radius(radiusSphere);
    HC_11_l->set_contact_half_space_location(halfSpaceLocation);
    HC_11_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_11_l->set_stiffness(stiffness);
    HC_11_l->set_dissipation(dissipation);
    HC_11_l->set_static_friction(staticFriction);
    HC_11_l->set_dynamic_friction(dynamicFriction);
    HC_11_l->set_viscous_friction(viscousFriction);
    HC_11_l->set_transition_velocity(transitionVelocity);
    HC_11_l->connectSocket_sphere_frame(*toes_l);
    HC_11_l->connectSocket_half_space_frame(model->getGround());

    HC_12_l = new SmoothSphereHalfSpaceForce("sphere_12_l", *toes_l, model->getGround());
    HC_12_l->set_contact_sphere_location(locSphere_12_l);
    HC_12_l->set_contact_sphere_radius(radiusSphere);
    HC_12_l->set_contact_half_space_location(halfSpaceLocation);
    HC_12_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_12_l->set_stiffness(stiffness);
    HC_12_l->set_dissipation(dissipation);
    HC_12_l->set_static_friction(staticFriction);
    HC_12_l->set_dynamic_friction(dynamicFriction);
    HC_12_l->set_viscous_friction(viscousFriction);
    HC_12_l->set_transition_velocity(transitionVelocity);
    HC_12_l->connectSocket_sphere_frame(*toes_l);
    HC_12_l->connectSocket_half_space_frame(model->getGround());

    HC_13_l = new SmoothSphereHalfSpaceForce("sphere_13_l", *toes_l, model->getGround());
    HC_13_l->set_contact_sphere_location(locSphere_13_l);
    HC_13_l->set_contact_sphere_radius(radiusSphere);
    HC_13_l->set_contact_half_space_location(halfSpaceLocation);
    HC_13_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_13_l->set_stiffness(stiffness);
    HC_13_l->set_dissipation(dissipation);
    HC_13_l->set_static_friction(staticFriction);
    HC_13_l->set_dynamic_friction(dynamicFriction);
    HC_13_l->set_viscous_friction(viscousFriction);
    HC_13_l->set_transition_velocity(transitionVelocity);
    HC_13_l->connectSocket_sphere_frame(*toes_l);
    HC_13_l->connectSocket_half_space_frame(model->getGround());

    HC_14_l = new SmoothSphereHalfSpaceForce("sphere_14_l", *toes_l, model->getGround());
    HC_14_l->set_contact_sphere_location(locSphere_14_l);
    HC_14_l->set_contact_sphere_radius(radiusSphere);
    HC_14_l->set_contact_half_space_location(halfSpaceLocation);
    HC_14_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_14_l->set_stiffness(stiffness);
    HC_14_l->set_dissipation(dissipation);
    HC_14_l->set_static_friction(staticFriction);
    HC_14_l->set_dynamic_friction(dynamicFriction);
    HC_14_l->set_viscous_friction(viscousFriction);
    HC_14_l->set_transition_velocity(transitionVelocity);
    HC_14_l->connectSocket_sphere_frame(*toes_l);
    HC_14_l->connectSocket_half_space_frame(model->getGround());

    HC_15_l = new SmoothSphereHalfSpaceForce("sphere_15_l", *toes_l, model->getGround());
    HC_15_l->set_contact_sphere_location(locSphere_15_l);
    HC_15_l->set_contact_sphere_radius(radiusSphere);
    HC_15_l->set_contact_half_space_location(halfSpaceLocation);
    HC_15_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_15_l->set_stiffness(stiffness);
    HC_15_l->set_dissipation(dissipation);
    HC_15_l->set_static_friction(staticFriction);
    HC_15_l->set_dynamic_friction(dynamicFriction);
    HC_15_l->set_viscous_friction(viscousFriction);
    HC_15_l->set_transition_velocity(transitionVelocity);
    HC_15_l->connectSocket_sphere_frame(*toes_l);
    HC_15_l->connectSocket_half_space_frame(model->getGround());

    HC_16_l = new SmoothSphereHalfSpaceForce("sphere_16_l", *toes_l, model->getGround());
    HC_16_l->set_contact_sphere_location(locSphere_16_l);
    HC_16_l->set_contact_sphere_radius(radiusSphere);
    HC_16_l->set_contact_half_space_location(halfSpaceLocation);
    HC_16_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_16_l->set_stiffness(stiffness);
    HC_16_l->set_dissipation(dissipation);
    HC_16_l->set_static_friction(staticFriction);
    HC_16_l->set_dynamic_friction(dynamicFriction);
    HC_16_l->set_viscous_friction(viscousFriction);
    HC_16_l->set_transition_velocity(transitionVelocity);
    HC_16_l->connectSocket_sphere_frame(*toes_l);
    HC_16_l->connectSocket_half_space_frame(model->getGround());

    HC_17_l = new SmoothSphereHalfSpaceForce("sphere_17_l", *toes_l, model->getGround());
    HC_17_l->set_contact_sphere_location(locSphere_17_l);
    HC_17_l->set_contact_sphere_radius(radiusSphere);
    HC_17_l->set_contact_half_space_location(halfSpaceLocation);
    HC_17_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_17_l->set_stiffness(stiffness);
    HC_17_l->set_dissipation(dissipation);
    HC_17_l->set_static_friction(staticFriction);
    HC_17_l->set_dynamic_friction(dynamicFriction);
    HC_17_l->set_viscous_friction(viscousFriction);
    HC_17_l->set_transition_velocity(transitionVelocity);
    HC_17_l->connectSocket_sphere_frame(*toes_l);
    HC_17_l->connectSocket_half_space_frame(model->getGround());

    HC_18_l = new SmoothSphereHalfSpaceForce("sphere_18_l", *toes_l, model->getGround());
    HC_18_l->set_contact_sphere_location(locSphere_18_l);
    HC_18_l->set_contact_sphere_radius(radiusSphere);
    HC_18_l->set_contact_half_space_location(halfSpaceLocation);
    HC_18_l->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_18_l->set_stiffness(stiffness);
    HC_18_l->set_dissipation(dissipation);
    HC_18_l->set_static_friction(staticFriction);
    HC_18_l->set_dynamic_friction(dynamicFriction);
    HC_18_l->set_viscous_friction(viscousFriction);
    HC_18_l->set_transition_velocity(transitionVelocity);
    HC_18_l->connectSocket_sphere_frame(*toes_l);
    HC_18_l->connectSocket_half_space_frame(model->getGround());

    /// Add left foot contact spheres to model
    model->addComponent(HC_1_l);
    model->addComponent(HC_2_l);
    model->addComponent(HC_3_l);
    model->addComponent(HC_4_l);
    model->addComponent(HC_5_l);
    model->addComponent(HC_6_l);
    model->addComponent(HC_7_l);
    model->addComponent(HC_8_l);
    model->addComponent(HC_9_l);
    model->addComponent(HC_10_l);
    model->addComponent(HC_11_l);
    model->addComponent(HC_12_l);
    model->addComponent(HC_13_l);
    model->addComponent(HC_14_l);
    model->addComponent(HC_15_l);
    model->addComponent(HC_16_l);
    model->addComponent(HC_17_l);
    model->addComponent(HC_18_l);

    /// Right foot contact shere specifications
    HC_1_r = new SmoothSphereHalfSpaceForce("sphere_1_r", *calcn_r, model->getGround());
    HC_1_r->set_contact_sphere_location(locSphere_1_r);
    HC_1_r->set_contact_sphere_radius(radiusSphere);
    HC_1_r->set_contact_half_space_location(halfSpaceLocation);
    HC_1_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_1_r->set_stiffness(stiffness);
    HC_1_r->set_dissipation(dissipation);
    HC_1_r->set_static_friction(staticFriction);
    HC_1_r->set_dynamic_friction(dynamicFriction);
    HC_1_r->set_viscous_friction(viscousFriction);
    HC_1_r->set_transition_velocity(transitionVelocity);
    HC_1_r->connectSocket_sphere_frame(*calcn_r);
    HC_1_r->connectSocket_half_space_frame(model->getGround());

    HC_2_r = new SmoothSphereHalfSpaceForce("sphere_2_r", *calcn_r, model->getGround());
    HC_2_r->set_contact_sphere_location(locSphere_2_r);
    HC_2_r->set_contact_sphere_radius(radiusSphere);
    HC_2_r->set_contact_half_space_location(halfSpaceLocation);
    HC_2_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_2_r->set_stiffness(stiffness);
    HC_2_r->set_dissipation(dissipation);
    HC_2_r->set_static_friction(staticFriction);
    HC_2_r->set_dynamic_friction(dynamicFriction);
    HC_2_r->set_viscous_friction(viscousFriction);
    HC_2_r->set_transition_velocity(transitionVelocity);
    HC_2_r->connectSocket_sphere_frame(*calcn_r);
    HC_2_r->connectSocket_half_space_frame(model->getGround());

    HC_3_r = new SmoothSphereHalfSpaceForce("sphere_3_r", *calcn_r, model->getGround());
    HC_3_r->set_contact_sphere_location(locSphere_3_r);
    HC_3_r->set_contact_sphere_radius(radiusSphere);
    HC_3_r->set_contact_half_space_location(halfSpaceLocation);
    HC_3_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_3_r->set_stiffness(stiffness);
    HC_3_r->set_dissipation(dissipation);
    HC_3_r->set_static_friction(staticFriction);
    HC_3_r->set_dynamic_friction(dynamicFriction);
    HC_3_r->set_viscous_friction(viscousFriction);
    HC_3_r->set_transition_velocity(transitionVelocity);
    HC_3_r->connectSocket_sphere_frame(*calcn_r);
    HC_3_r->connectSocket_half_space_frame(model->getGround());

    HC_4_r = new SmoothSphereHalfSpaceForce("sphere_4_r", *calcn_r, model->getGround());
    HC_4_r->set_contact_sphere_location(locSphere_4_r);
    HC_4_r->set_contact_sphere_radius(radiusSphere);
    HC_4_r->set_contact_half_space_location(halfSpaceLocation);
    HC_4_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_4_r->set_stiffness(stiffness);
    HC_4_r->set_dissipation(dissipation);
    HC_4_r->set_static_friction(staticFriction);
    HC_4_r->set_dynamic_friction(dynamicFriction);
    HC_4_r->set_viscous_friction(viscousFriction);
    HC_4_r->set_transition_velocity(transitionVelocity);
    HC_4_r->connectSocket_sphere_frame(*calcn_r);
    HC_4_r->connectSocket_half_space_frame(model->getGround());

    HC_5_r = new SmoothSphereHalfSpaceForce("sphere_5_r", *calcn_r, model->getGround());
    HC_5_r->set_contact_sphere_location(locSphere_5_r);
    HC_5_r->set_contact_sphere_radius(radiusSphere);
    HC_5_r->set_contact_half_space_location(halfSpaceLocation);
    HC_5_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_5_r->set_stiffness(stiffness);
    HC_5_r->set_dissipation(dissipation);
    HC_5_r->set_static_friction(staticFriction);
    HC_5_r->set_dynamic_friction(dynamicFriction);
    HC_5_r->set_viscous_friction(viscousFriction);
    HC_5_r->set_transition_velocity(transitionVelocity);
    HC_5_r->connectSocket_sphere_frame(*calcn_r);
    HC_5_r->connectSocket_half_space_frame(model->getGround());

    HC_6_r = new SmoothSphereHalfSpaceForce("sphere_6_r", *calcn_r, model->getGround());
    HC_6_r->set_contact_sphere_location(locSphere_6_r);
    HC_6_r->set_contact_sphere_radius(radiusSphere);
    HC_6_r->set_contact_half_space_location(halfSpaceLocation);
    HC_6_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_6_r->set_stiffness(stiffness);
    HC_6_r->set_dissipation(dissipation);
    HC_6_r->set_static_friction(staticFriction);
    HC_6_r->set_dynamic_friction(dynamicFriction);
    HC_6_r->set_viscous_friction(viscousFriction);
    HC_6_r->set_transition_velocity(transitionVelocity);
    HC_6_r->connectSocket_sphere_frame(*calcn_r);
    HC_6_r->connectSocket_half_space_frame(model->getGround());

    HC_7_r = new SmoothSphereHalfSpaceForce("sphere_7_r", *calcn_r, model->getGround());
    HC_7_r->set_contact_sphere_location(locSphere_7_r);
    HC_7_r->set_contact_sphere_radius(radiusSphere);
    HC_7_r->set_contact_half_space_location(halfSpaceLocation);
    HC_7_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_7_r->set_stiffness(stiffness);
    HC_7_r->set_dissipation(dissipation);
    HC_7_r->set_static_friction(staticFriction);
    HC_7_r->set_dynamic_friction(dynamicFriction);
    HC_7_r->set_viscous_friction(viscousFriction);
    HC_7_r->set_transition_velocity(transitionVelocity);
    HC_7_r->connectSocket_sphere_frame(*calcn_r);
    HC_7_r->connectSocket_half_space_frame(model->getGround());

    HC_8_r = new SmoothSphereHalfSpaceForce("sphere_8_r", *calcn_r, model->getGround());
    HC_8_r->set_contact_sphere_location(locSphere_8_r);
    HC_8_r->set_contact_sphere_radius(radiusSphere);
    HC_8_r->set_contact_half_space_location(halfSpaceLocation);
    HC_8_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_8_r->set_stiffness(stiffness);
    HC_8_r->set_dissipation(dissipation);
    HC_8_r->set_static_friction(staticFriction);
    HC_8_r->set_dynamic_friction(dynamicFriction);
    HC_8_r->set_viscous_friction(viscousFriction);
    HC_8_r->set_transition_velocity(transitionVelocity);
    HC_8_r->connectSocket_sphere_frame(*calcn_r);
    HC_8_r->connectSocket_half_space_frame(model->getGround());

    HC_9_r = new SmoothSphereHalfSpaceForce("sphere_9_r", *calcn_r, model->getGround());
    HC_9_r->set_contact_sphere_location(locSphere_9_r);
    HC_9_r->set_contact_sphere_radius(radiusSphere);
    HC_9_r->set_contact_half_space_location(halfSpaceLocation);
    HC_9_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_9_r->set_stiffness(stiffness);
    HC_9_r->set_dissipation(dissipation);
    HC_9_r->set_static_friction(staticFriction);
    HC_9_r->set_dynamic_friction(dynamicFriction);
    HC_9_r->set_viscous_friction(viscousFriction);
    HC_9_r->set_transition_velocity(transitionVelocity);
    HC_9_r->connectSocket_sphere_frame(*calcn_r);
    HC_9_r->connectSocket_half_space_frame(model->getGround());

    HC_10_r = new SmoothSphereHalfSpaceForce("sphere_10_r", *calcn_r, model->getGround());
    HC_10_r->set_contact_sphere_location(locSphere_10_r);
    HC_10_r->set_contact_sphere_radius(radiusSphere);
    HC_10_r->set_contact_half_space_location(halfSpaceLocation);
    HC_10_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_10_r->set_stiffness(stiffness);
    HC_10_r->set_dissipation(dissipation);
    HC_10_r->set_static_friction(staticFriction);
    HC_10_r->set_dynamic_friction(dynamicFriction);
    HC_10_r->set_viscous_friction(viscousFriction);
    HC_10_r->set_transition_velocity(transitionVelocity);
    HC_10_r->connectSocket_sphere_frame(*calcn_r);
    HC_10_r->connectSocket_half_space_frame(model->getGround());

    HC_11_r = new SmoothSphereHalfSpaceForce("sphere_11_r", *toes_r, model->getGround());
    HC_11_r->set_contact_sphere_location(locSphere_11_r);
    HC_11_r->set_contact_sphere_radius(radiusSphere);
    HC_11_r->set_contact_half_space_location(halfSpaceLocation);
    HC_11_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_11_r->set_stiffness(stiffness);
    HC_11_r->set_dissipation(dissipation);
    HC_11_r->set_static_friction(staticFriction);
    HC_11_r->set_dynamic_friction(dynamicFriction);
    HC_11_r->set_viscous_friction(viscousFriction);
    HC_11_r->set_transition_velocity(transitionVelocity);
    HC_11_r->connectSocket_sphere_frame(*toes_r);
    HC_11_r->connectSocket_half_space_frame(model->getGround());

    HC_12_r = new SmoothSphereHalfSpaceForce("sphere_12_r", *toes_r, model->getGround());
    HC_12_r->set_contact_sphere_location(locSphere_12_r);
    HC_12_r->set_contact_sphere_radius(radiusSphere);
    HC_12_r->set_contact_half_space_location(halfSpaceLocation);
    HC_12_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_12_r->set_stiffness(stiffness);
    HC_12_r->set_dissipation(dissipation);
    HC_12_r->set_static_friction(staticFriction);
    HC_12_r->set_dynamic_friction(dynamicFriction);
    HC_12_r->set_viscous_friction(viscousFriction);
    HC_12_r->set_transition_velocity(transitionVelocity);
    HC_12_r->connectSocket_sphere_frame(*toes_r);
    HC_12_r->connectSocket_half_space_frame(model->getGround());

    HC_13_r = new SmoothSphereHalfSpaceForce("sphere_13_r", *toes_r, model->getGround());
    HC_13_r->set_contact_sphere_location(locSphere_13_r);
    HC_13_r->set_contact_sphere_radius(radiusSphere);
    HC_13_r->set_contact_half_space_location(halfSpaceLocation);
    HC_13_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_13_r->set_stiffness(stiffness);
    HC_13_r->set_dissipation(dissipation);
    HC_13_r->set_static_friction(staticFriction);
    HC_13_r->set_dynamic_friction(dynamicFriction);
    HC_13_r->set_viscous_friction(viscousFriction);
    HC_13_r->set_transition_velocity(transitionVelocity);
    HC_13_r->connectSocket_sphere_frame(*toes_r);
    HC_13_r->connectSocket_half_space_frame(model->getGround());

    HC_14_r = new SmoothSphereHalfSpaceForce("sphere_14_r", *toes_r, model->getGround());
    HC_14_r->set_contact_sphere_location(locSphere_14_r);
    HC_14_r->set_contact_sphere_radius(radiusSphere);
    HC_14_r->set_contact_half_space_location(halfSpaceLocation);
    HC_14_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_14_r->set_stiffness(stiffness);
    HC_14_r->set_dissipation(dissipation);
    HC_14_r->set_static_friction(staticFriction);
    HC_14_r->set_dynamic_friction(dynamicFriction);
    HC_14_r->set_viscous_friction(viscousFriction);
    HC_14_r->set_transition_velocity(transitionVelocity);
    HC_14_r->connectSocket_sphere_frame(*toes_r);
    HC_14_r->connectSocket_half_space_frame(model->getGround());

    HC_15_r = new SmoothSphereHalfSpaceForce("sphere_15_r", *toes_r, model->getGround());
    HC_15_r->set_contact_sphere_location(locSphere_15_r);
    HC_15_r->set_contact_sphere_radius(radiusSphere);
    HC_15_r->set_contact_half_space_location(halfSpaceLocation);
    HC_15_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_15_r->set_stiffness(stiffness);
    HC_15_r->set_dissipation(dissipation);
    HC_15_r->set_static_friction(staticFriction);
    HC_15_r->set_dynamic_friction(dynamicFriction);
    HC_15_r->set_viscous_friction(viscousFriction);
    HC_15_r->set_transition_velocity(transitionVelocity);
    HC_15_r->connectSocket_sphere_frame(*toes_r);
    HC_15_r->connectSocket_half_space_frame(model->getGround());

    HC_16_r = new SmoothSphereHalfSpaceForce("sphere_16_r", *toes_r, model->getGround());
    HC_16_r->set_contact_sphere_location(locSphere_16_r);
    HC_16_r->set_contact_sphere_radius(radiusSphere);
    HC_16_r->set_contact_half_space_location(halfSpaceLocation);
    HC_16_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_16_r->set_stiffness(stiffness);
    HC_16_r->set_dissipation(dissipation);
    HC_16_r->set_static_friction(staticFriction);
    HC_16_r->set_dynamic_friction(dynamicFriction);
    HC_16_r->set_viscous_friction(viscousFriction);
    HC_16_r->set_transition_velocity(transitionVelocity);
    HC_16_r->connectSocket_sphere_frame(*toes_r);
    HC_16_r->connectSocket_half_space_frame(model->getGround());

    HC_17_r = new SmoothSphereHalfSpaceForce("sphere_17_r", *toes_r, model->getGround());
    HC_17_r->set_contact_sphere_location(locSphere_17_r);
    HC_17_r->set_contact_sphere_radius(radiusSphere);
    HC_17_r->set_contact_half_space_location(halfSpaceLocation);
    HC_17_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_17_r->set_stiffness(stiffness);
    HC_17_r->set_dissipation(dissipation);
    HC_17_r->set_static_friction(staticFriction);
    HC_17_r->set_dynamic_friction(dynamicFriction);
    HC_17_r->set_viscous_friction(viscousFriction);
    HC_17_r->set_transition_velocity(transitionVelocity);
    HC_17_r->connectSocket_sphere_frame(*toes_r);
    HC_17_r->connectSocket_half_space_frame(model->getGround());

    HC_18_r = new SmoothSphereHalfSpaceForce("sphere_18_r", *toes_r, model->getGround());
    HC_18_r->set_contact_sphere_location(locSphere_18_r);
    HC_18_r->set_contact_sphere_radius(radiusSphere);
    HC_18_r->set_contact_half_space_location(halfSpaceLocation);
    HC_18_r->set_contact_half_space_orientation(halfSpaceOrientation);
    HC_18_r->set_stiffness(stiffness);
    HC_18_r->set_dissipation(dissipation);
    HC_18_r->set_static_friction(staticFriction);
    HC_18_r->set_dynamic_friction(dynamicFriction);
    HC_18_r->set_viscous_friction(viscousFriction);
    HC_18_r->set_transition_velocity(transitionVelocity);
    HC_18_r->connectSocket_sphere_frame(*toes_r);
    HC_18_r->connectSocket_half_space_frame(model->getGround());

    /// Add left foot contact spheres to model
    model->addComponent(HC_1_r);
    model->addComponent(HC_2_r);
    model->addComponent(HC_3_r);
    model->addComponent(HC_4_r);
    model->addComponent(HC_5_r);
    model->addComponent(HC_6_r);
    model->addComponent(HC_7_r);
    model->addComponent(HC_8_r);
    model->addComponent(HC_9_r);
    model->addComponent(HC_10_r);
    model->addComponent(HC_11_r);
    model->addComponent(HC_12_r);
    model->addComponent(HC_13_r);
    model->addComponent(HC_14_r);
    model->addComponent(HC_15_r);
    model->addComponent(HC_16_r);
    model->addComponent(HC_17_r);
    model->addComponent(HC_18_r);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);

    // States and controls
    T ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// pro_sup dofs are locked so Qs and Qdots are hard coded (0)
    QsUs[NX] = 1.51;
    QsUs[NX+1] = 0;
    QsUs[NX+2] = 1.51;
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
    /// Add contact forces to appliedBodyForces
    /// Right foot
    Array<osim_double_adouble> Force_values_1_r = HC_1_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_2_r = HC_2_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_3_r = HC_3_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_4_r = HC_4_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_5_r = HC_5_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_6_r = HC_6_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_7_r = HC_7_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_8_r = HC_8_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_9_r = HC_9_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_10_r = HC_10_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_11_r = HC_11_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_12_r = HC_12_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_13_r = HC_13_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_14_r = HC_14_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_15_r = HC_15_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_16_r = HC_16_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_17_r = HC_17_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_18_r = HC_18_r->getRecordValues(*state);
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
    SpatialVec GRF_7_r;
    GRF_7_r[0] = Vec3(Force_values_7_r[3], Force_values_7_r[4], Force_values_7_r[5]);
    GRF_7_r[1] = Vec3(Force_values_7_r[0], Force_values_7_r[1], Force_values_7_r[2]);
    SpatialVec GRF_8_r;
    GRF_8_r[0] = Vec3(Force_values_8_r[3], Force_values_8_r[4], Force_values_8_r[5]);
    GRF_8_r[1] = Vec3(Force_values_8_r[0], Force_values_8_r[1], Force_values_8_r[2]);
    SpatialVec GRF_9_r;
    GRF_9_r[0] = Vec3(Force_values_9_r[3], Force_values_9_r[4], Force_values_9_r[5]);
    GRF_9_r[1] = Vec3(Force_values_9_r[0], Force_values_9_r[1], Force_values_9_r[2]);
    SpatialVec GRF_10_r;
    GRF_10_r[0] = Vec3(Force_values_10_r[3], Force_values_10_r[4], Force_values_10_r[5]);
    GRF_10_r[1] = Vec3(Force_values_10_r[0], Force_values_10_r[1], Force_values_10_r[2]);
    SpatialVec GRF_11_r;
    GRF_11_r[0] = Vec3(Force_values_11_r[3], Force_values_11_r[4], Force_values_11_r[5]);
    GRF_11_r[1] = Vec3(Force_values_11_r[0], Force_values_11_r[1], Force_values_11_r[2]);
    SpatialVec GRF_12_r;
    GRF_12_r[0] = Vec3(Force_values_12_r[3], Force_values_12_r[4], Force_values_12_r[5]);
    GRF_12_r[1] = Vec3(Force_values_12_r[0], Force_values_12_r[1], Force_values_12_r[2]);
    SpatialVec GRF_13_r;
    GRF_13_r[0] = Vec3(Force_values_13_r[3], Force_values_13_r[4], Force_values_13_r[5]);
    GRF_13_r[1] = Vec3(Force_values_13_r[0], Force_values_13_r[1], Force_values_13_r[2]);
    SpatialVec GRF_14_r;
    GRF_14_r[0] = Vec3(Force_values_14_r[3], Force_values_14_r[4], Force_values_14_r[5]);
    GRF_14_r[1] = Vec3(Force_values_14_r[0], Force_values_14_r[1], Force_values_14_r[2]);
    SpatialVec GRF_15_r;
    GRF_15_r[0] = Vec3(Force_values_15_r[3], Force_values_15_r[4], Force_values_15_r[5]);
    GRF_15_r[1] = Vec3(Force_values_15_r[0], Force_values_15_r[1], Force_values_15_r[2]);
    SpatialVec GRF_16_r;
    GRF_16_r[0] = Vec3(Force_values_16_r[3], Force_values_16_r[4], Force_values_16_r[5]);
    GRF_16_r[1] = Vec3(Force_values_16_r[0], Force_values_16_r[1], Force_values_16_r[2]);
    SpatialVec GRF_17_r;
    GRF_17_r[0] = Vec3(Force_values_17_r[3], Force_values_17_r[4], Force_values_17_r[5]);
    GRF_17_r[1] = Vec3(Force_values_17_r[0], Force_values_17_r[1], Force_values_17_r[2]);
    SpatialVec GRF_18_r;
    GRF_18_r[0] = Vec3(Force_values_18_r[3], Force_values_18_r[4], Force_values_18_r[5]);
    GRF_18_r[1] = Vec3(Force_values_18_r[0], Force_values_18_r[1], Force_values_18_r[2]);
    int ncalcn_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
    int ntoes_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
    appliedBodyForces[ncalcn_r] = appliedBodyForces[ncalcn_r] + GRF_1_r + GRF_2_r + GRF_3_r + GRF_4_r + GRF_5_r + GRF_6_r + GRF_7_r + GRF_8_r + GRF_9_r + GRF_10_r;
    appliedBodyForces[ntoes_r] = appliedBodyForces[ntoes_r] + GRF_11_r + GRF_12_r + GRF_13_r + GRF_14_r + GRF_15_r + GRF_16_r + GRF_17_r + GRF_18_r;
    /// Left foot
    Array<osim_double_adouble> Force_values_1_l = HC_1_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_2_l = HC_2_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_3_l = HC_3_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_4_l = HC_4_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_5_l = HC_5_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_6_l = HC_6_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_7_l = HC_7_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_8_l = HC_8_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_9_l = HC_9_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_10_l = HC_10_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_11_l = HC_11_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_12_l = HC_12_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_13_l = HC_13_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_14_l = HC_14_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_15_l = HC_15_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_16_l = HC_16_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_17_l = HC_17_l->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_18_l = HC_18_l->getRecordValues(*state);
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
    SpatialVec GRF_7_l;
    GRF_7_l[0] = Vec3(Force_values_7_l[3], Force_values_7_l[4], Force_values_7_l[5]);
    GRF_7_l[1] = Vec3(Force_values_7_l[0], Force_values_7_l[1], Force_values_7_l[2]);
    SpatialVec GRF_8_l;
    GRF_8_l[0] = Vec3(Force_values_8_l[3], Force_values_8_l[4], Force_values_8_l[5]);
    GRF_8_l[1] = Vec3(Force_values_8_l[0], Force_values_8_l[1], Force_values_8_l[2]);
    SpatialVec GRF_9_l;
    GRF_9_l[0] = Vec3(Force_values_9_l[3], Force_values_9_l[4], Force_values_9_l[5]);
    GRF_9_l[1] = Vec3(Force_values_9_l[0], Force_values_9_l[1], Force_values_9_l[2]);
    SpatialVec GRF_10_l;
    GRF_10_l[0] = Vec3(Force_values_10_l[3], Force_values_10_l[4], Force_values_10_l[5]);
    GRF_10_l[1] = Vec3(Force_values_10_l[0], Force_values_10_l[1], Force_values_10_l[2]);
    SpatialVec GRF_11_l;
    GRF_11_l[0] = Vec3(Force_values_11_l[3], Force_values_11_l[4], Force_values_11_l[5]);
    GRF_11_l[1] = Vec3(Force_values_11_l[0], Force_values_11_l[1], Force_values_11_l[2]);
    SpatialVec GRF_12_l;
    GRF_12_l[0] = Vec3(Force_values_12_l[3], Force_values_12_l[4], Force_values_12_l[5]);
    GRF_12_l[1] = Vec3(Force_values_12_l[0], Force_values_12_l[1], Force_values_12_l[2]);
    SpatialVec GRF_13_l;
    GRF_13_l[0] = Vec3(Force_values_13_l[3], Force_values_13_l[4], Force_values_13_l[5]);
    GRF_13_l[1] = Vec3(Force_values_13_l[0], Force_values_13_l[1], Force_values_13_l[2]);
    SpatialVec GRF_14_l;
    GRF_14_l[0] = Vec3(Force_values_14_l[3], Force_values_14_l[4], Force_values_14_l[5]);
    GRF_14_l[1] = Vec3(Force_values_14_l[0], Force_values_14_l[1], Force_values_14_l[2]);
    SpatialVec GRF_15_l;
    GRF_15_l[0] = Vec3(Force_values_15_l[3], Force_values_15_l[4], Force_values_15_l[5]);
    GRF_15_l[1] = Vec3(Force_values_15_l[0], Force_values_15_l[1], Force_values_15_l[2]);
    SpatialVec GRF_16_l;
    GRF_16_l[0] = Vec3(Force_values_16_l[3], Force_values_16_l[4], Force_values_16_l[5]);
    GRF_16_l[1] = Vec3(Force_values_16_l[0], Force_values_16_l[1], Force_values_16_l[2]);
    SpatialVec GRF_17_l;
    GRF_17_l[0] = Vec3(Force_values_17_l[3], Force_values_17_l[4], Force_values_17_l[5]);
    GRF_17_l[1] = Vec3(Force_values_17_l[0], Force_values_17_l[1], Force_values_17_l[2]);
    SpatialVec GRF_18_l;
    GRF_18_l[0] = Vec3(Force_values_18_l[3], Force_values_18_l[4], Force_values_18_l[5]);
    GRF_18_l[1] = Vec3(Force_values_18_l[0], Force_values_18_l[1], Force_values_18_l[2]);
    int ncalcn_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
    int ntoes_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
    appliedBodyForces[ncalcn_l] = appliedBodyForces[ncalcn_l] + GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l + GRF_5_l + GRF_6_l + GRF_7_l + GRF_8_l + GRF_9_l + GRF_10_l;
    appliedBodyForces[ntoes_l] = appliedBodyForces[ntoes_l] + GRF_11_l + GRF_12_l + GRF_13_l + GRF_14_l + GRF_15_l + GRF_16_l + GRF_17_l + GRF_18_l;
    /// knownUdot
    Vector knownUdot(ndofr);
    knownUdot.setToZero();
    for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
    Vector residualMobilityForces(ndofr);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Extract several joint origins to set constraints in problem
    Vec3 calcn_or_l  = calcn_l->getPositionInGround(*state);
    Vec3 calcn_or_r  = calcn_r->getPositionInGround(*state);
    Vec3 femur_or_l  = femur_l->getPositionInGround(*state);
    Vec3 femur_or_r  = femur_r->getPositionInGround(*state);
    Vec3 hand_or_l   = hand_l->getPositionInGround(*state);
    Vec3 hand_or_r   = hand_r->getPositionInGround(*state);
    Vec3 tibia_or_l  = tibia_l->getPositionInGround(*state);
    Vec3 tibia_or_r  = tibia_r->getPositionInGround(*state);
    Vec3 toes_or_l  = toes_l->getPositionInGround(*state);
    Vec3 toes_or_r  = toes_r->getPositionInGround(*state);

    // Residual forces in OpenSim order
    T res_os[ndofr];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < ndofr; ++i) res_os[i] =
            value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
    // Extract results
    int nc = 3;
    /// Residual forces
    /// We do want to extract the pro_sup torques (last two -> till NU)
    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];
    /// Joint origins
    res[0][NU]     = value<T>(calcn_or_r[0]);   /// calcn_or_r_x
    res[0][NU + 1] = value<T>(calcn_or_r[2]);   /// calcn_or_r_z
    res[0][NU + 2] = value<T>(calcn_or_l[0]);   /// calcn_or_l_x
    res[0][NU + 3] = value<T>(calcn_or_l[2]);   /// calcn_or_l_x
    res[0][NU + 4] = value<T>(femur_or_r[0]);   /// femur_or_r_x
    res[0][NU + 5] = value<T>(femur_or_r[2]);   /// femur_or_r_z
    res[0][NU + 6] = value<T>(femur_or_l[0]);   /// femur_or_l_x
    res[0][NU + 7] = value<T>(femur_or_l[2]);   /// femur_or_l_z
    res[0][NU + 8] = value<T>(hand_or_r[0]);    /// hand_or_r_x
    res[0][NU + 9] = value<T>(hand_or_r[2]);    /// hand_or_r_z
    res[0][NU + 10] = value<T>(hand_or_l[0]);   /// hand_or_l_x
    res[0][NU + 11] = value<T>(hand_or_l[2]);   /// hand_or_l_z
    res[0][NU + 12] = value<T>(tibia_or_r[0]);  /// tibia_or_r_x
    res[0][NU + 13] = value<T>(tibia_or_r[2]);  /// tibia_or_r_z
    res[0][NU + 14] = value<T>(tibia_or_l[0]);  /// tibia_or_l_x
    res[0][NU + 15] = value<T>(tibia_or_l[2]);  /// tibia_or_l_z
    res[0][NU + 16] = value<T>(toes_or_r[0]);  /// tibia_or_r_x
    res[0][NU + 17] = value<T>(toes_or_r[2]);  /// tibia_or_r_z
    res[0][NU + 18] = value<T>(toes_or_l[0]);  /// tibia_or_l_x
    res[0][NU + 19] = value<T>(toes_or_l[2]);  /// tibia_or_l_z

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
    for (int i = 0; i < NR; ++i) Recorder_res[0][i] >>= res[i];

    Recorder::stop_recording();

    return 0;

}
