/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
//#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Common/SimmSpline.h>
//#include <OpenSim/Common/MultiplierFunction.h>
#include <OpenSim/Simulation/Model/HuntCrossleyForce_smooth.h>
#include "SimTKcommon/internal/recorder.h"
//#include <OpenSim/Simulation/Model/PhysicalOffsetFrame.h>
#include <OpenSim/Simulation/Model/Station.h>

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
constexpr int ndof = 23;        // # degrees of freedom
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NR = ndof + 16*3; // # residual torques

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
    //OpenSim::Body* humerus_r;
    //OpenSim::Body* humerus_l;
    //OpenSim::Body* ulna_r;
    //OpenSim::Body* ulna_l;
    //OpenSim::Body* radius_r;
    //OpenSim::Body* radius_l;
    //OpenSim::Body* hand_r;
    //OpenSim::Body* hand_l;
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
    /// Contact elements
    OpenSim::HuntCrossleyForce_smooth* HC_1_r;
    OpenSim::HuntCrossleyForce_smooth* HC_2_r;
    OpenSim::HuntCrossleyForce_smooth* HC_3_r;
    OpenSim::HuntCrossleyForce_smooth* HC_4_r;
    OpenSim::HuntCrossleyForce_smooth* HC_5_r;
    OpenSim::HuntCrossleyForce_smooth* HC_6_r;
    OpenSim::HuntCrossleyForce_smooth* HC_1_l;
    OpenSim::HuntCrossleyForce_smooth* HC_2_l;
    OpenSim::HuntCrossleyForce_smooth* HC_3_l;
    OpenSim::HuntCrossleyForce_smooth* HC_4_l;
    OpenSim::HuntCrossleyForce_smooth* HC_5_l;
    OpenSim::HuntCrossleyForce_smooth* HC_6_l;
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
    pelvis = new OpenSim::Body("pelvis", 10.638762130045265, Vec3(-0.065932147220013546, 0, 0), Inertia(0.092864460131498114, 0.078681852893516396, 0.052304010132429385, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 8.4024269403416003, Vec3(0, -0.17424375354032129, 0), Inertia(0.12707309096764055, 0.033310421904138784, 0.13400089951180619, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 8.4024269403416003, Vec3(0, -0.17424375354032129, 0), Inertia(0.12707309096764055, 0.033310421904138784, 0.13400089951180619, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 3.3491730149565102, Vec3(0, -0.1941032441286128, 0), Inertia(0.049211195113498758, 0.0049797042674373741, 0.049894683934519575, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 3.3491730149565102, Vec3(0, -0.1941032441286128, 0), Inertia(0.049211195113498758, 0.0049797042674373741, 0.049894683934519575, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.09033507794892813, Vec3(0, 0, 0), Inertia(0.0009033507794892812, 0.0009033507794892812, 0.0009033507794892812, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.09033507794892813, Vec3(0, 0, 0), Inertia(0.0009033507794892812, 0.0009033507794892812, 0.0009033507794892812, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 1.1291884743616016, Vec3(0.099874593118053412, 0.029962377935416021, 0), Inertia(0.0012646910912849937, 0.0035230680400081967, 0.0037037381959060534, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 1.1291884743616016, Vec3(0.099874593118053412, 0.029962377935416021, 0), Inertia(0.0012646910912849937, 0.0035230680400081967, 0.0037037381959060534, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.19566577883737829, Vec3(0.034556609218846478, 0.0059924755870832038, 0.018058146706533898), Inertia(9.0335077948928123e-05, 0.00018067015589785625, 0.0009033507794892812, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.19566577883737829, Vec3(0.034556609218846478, 0.0059924755870832038, -0.018058146706533898), Inertia(9.0335077948928123e-05, 0.00018067015589785625, 0.0009033507794892812, 0, 0, 0));
    torso = new OpenSim::Body("torso", 30.927659297062728, Vec3(-0.028028356858651442, 0.29555336623035872, 0), Inertia(1.3319907243569451, 0.6824815139041519, 1.2930563057609572, 0, 0, 0));
    //humerus_l = new OpenSim::Body("humerus_l", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
    //humerus_r = new OpenSim::Body("humerus_r", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
    //ulna_l = new OpenSim::Body("ulna_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //ulna_r = new OpenSim::Body("ulna_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //radius_l = new OpenSim::Body("radius_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //radius_r = new OpenSim::Body("radius_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    //hand_l = new OpenSim::Body("hand_l", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
    //hand_r = new OpenSim::Body("hand_r", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
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
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.065932147220013546, -0.061050242211958476, -0.077868943322081066), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.065932147220013546, -0.061050242211958476, 0.077868943322081066), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *femur_l, Vec3(-0.00461219, -0.40570167, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(-0.00461219, -0.40570167, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.4470508568575442, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.4470508568575442, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
    subtalar_l = new CustomJoint("subtalar_l", *talus_l, Vec3(-0.048708839063674646, -0.041897391813023405, -0.0081726012523284843), Vec3(0), *calcn_l, Vec3(0), Vec3(0), st_subtalar_l);
    subtalar_r = new CustomJoint("subtalar_r", *talus_r, Vec3(-0.048708839063674646, -0.041897391813023405, 0.0081726012523284843), Vec3(0), *calcn_r, Vec3(0), Vec3(0), st_subtalar_r);
    mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.17857577249507947, -0.0019974918623610681, -0.0011144456253175206), Vec3(0), *toes_l, Vec3(0), Vec3(0));
    mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.17857577249507947, -0.0019974918623610681, 0.0011144456253175206), Vec3(0), *toes_r, Vec3(0), Vec3(0));
    back = new CustomJoint("back", *pelvis, Vec3(-0.093909013084234286, 0.075273747961794479, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
    //shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, -0.16500992260263453), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
    //shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, 0.16500992260263453), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
    //elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.014421812854093517, -0.31410344120358441, 0.010527791717515771), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
    //elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.014421812854093517, -0.31410344120358441, -0.010527791717515771), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
    //radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.0078589566368657427, -0.015195696294888168, -0.030472003264362887), Vec3(0), *radius_l, Vec3(0), Vec3(0),st_radioulnar_l);
    //radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.0078589566368657427, -0.015195696294888168, 0.030472003264362887), Vec3(0), *radius_r, Vec3(0), Vec3(0),st_radioulnar_r);
    //radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.010277276874462305, -0.27552611746618899, -0.015900163494535866), Vec3(0), *hand_l, Vec3(0), Vec3(0));
    //radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.010277276874462305, -0.27552611746618899, 0.015900163494535866), Vec3(0), *hand_r, Vec3(0), Vec3(0));
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
    /// Contact elements
    /// Parameters
    osim_double_adouble radiusSphere_s1 = 0.03232;
    osim_double_adouble radiusSphere_s2 = 0.03232;
    osim_double_adouble radiusSphere_s3 = 0.023374;
    osim_double_adouble radiusSphere_s4 = 0.020508;
    osim_double_adouble radiusSphere_s5 = 0.016244;
    osim_double_adouble radiusSphere_s6 = 0.018414;
    osim_double_adouble stiffness = 1000000;
    osim_double_adouble dissipation = 2.0;
    osim_double_adouble staticFriction = 0.8;
    osim_double_adouble dynamicFriction = 0.8;
    osim_double_adouble viscousFriction = 0.5;
    osim_double_adouble transitionVelocity = 0.2;
    Vec3 normal = Vec3(0, 1, 0);
    osim_double_adouble offset = 0;
    Vec3 locSphere_1_r(-0.000460642, -0.01, -0.00564225);
    Vec3 locSphere_2_r(0.0655687, -0.01, 0.0225828);
    Vec3 locSphere_3_r(0.180314, -0.01, 0.0239173);
    Vec3 locSphere_4_r(0.180314, -0.01, -0.0112908);
    Vec3 locSphere_5_r(0.0580873, -0.01, -0.00385841);
    Vec3 locSphere_6_r(1.89942e-06, -0.01, 0.0251717);
    Vec3 locSphere_1_l(locSphere_1_r[0],locSphere_1_r[1],-locSphere_1_r[2]);
    Vec3 locSphere_2_l(locSphere_2_r[0],locSphere_2_r[1],-locSphere_2_r[2]);
    Vec3 locSphere_3_l(locSphere_3_r[0],locSphere_3_r[1],-locSphere_3_r[2]);
    Vec3 locSphere_4_l(locSphere_4_r[0],locSphere_4_r[1],-locSphere_4_r[2]);
    Vec3 locSphere_5_l(locSphere_5_r[0],locSphere_5_r[1],-locSphere_5_r[2]);
    Vec3 locSphere_6_l(locSphere_6_r[0],locSphere_6_r[1],-locSphere_6_r[2]);
    /// Left foot contact shere specifications
    HC_1_l = new HuntCrossleyForce_smooth("sphere_1_l", "calcn_l", locSphere_1_l, radiusSphere_s1,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_2_l = new HuntCrossleyForce_smooth("sphere_2_l", "calcn_l", locSphere_2_l, radiusSphere_s2,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_3_l = new HuntCrossleyForce_smooth("sphere_3_l", "calcn_l", locSphere_3_l, radiusSphere_s3,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_4_l = new HuntCrossleyForce_smooth("sphere_4_l", "calcn_l", locSphere_4_l, radiusSphere_s4,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_5_l = new HuntCrossleyForce_smooth("sphere_5_l", "toes_l", locSphere_5_l, radiusSphere_s5,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_6_l = new HuntCrossleyForce_smooth("sphere_6_l", "toes_l", locSphere_6_l, radiusSphere_s6,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    /// Add left foot contact spheres to model
    model->addComponent(HC_1_l);
    HC_1_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_2_l);
    HC_2_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_3_l);
    HC_3_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_4_l);
    HC_4_l->connectSocket_body_sphere(*calcn_l);
    model->addComponent(HC_5_l);
    HC_5_l->connectSocket_body_sphere(*toes_l);
    model->addComponent(HC_6_l);
    HC_6_l->connectSocket_body_sphere(*toes_l);
    /// Right foot contact shere specifications
    HC_1_r = new HuntCrossleyForce_smooth("sphere_1_r", "calcn_r", locSphere_1_r, radiusSphere_s1,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_2_r = new HuntCrossleyForce_smooth("sphere_2_r", "calcn_r", locSphere_2_r, radiusSphere_s2,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_3_r = new HuntCrossleyForce_smooth("sphere_3_r", "calcn_r", locSphere_3_r, radiusSphere_s3,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_4_r = new HuntCrossleyForce_smooth("sphere_4_r", "calcn_r", locSphere_4_r, radiusSphere_s4,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_5_r = new HuntCrossleyForce_smooth("sphere_5_r", "toes_r", locSphere_5_r, radiusSphere_s5,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    HC_6_r = new HuntCrossleyForce_smooth("sphere_6_r", "toes_r", locSphere_6_r, radiusSphere_s6,
        stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction, transitionVelocity, normal, offset);
    /// Add right foot contact spheres to model
    model->addComponent(HC_1_r);
    HC_1_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_2_r);
    HC_2_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_3_r);
    HC_3_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_4_r);
    HC_4_r->connectSocket_body_sphere(*calcn_r);
    model->addComponent(HC_5_r);
    HC_5_r->connectSocket_body_sphere(*toes_r);
    model->addComponent(HC_6_r);
    HC_6_r->connectSocket_body_sphere(*toes_r);

    /// OpenPose markers
    Neck = new Station(*torso, Vec3(0.0032925873626634683, 0.33154224015512557, 0.0030862732634675538));
    RShoulder = new Station(*torso, Vec3(0.0036590587612567171, 0.34172524342436539, 0.15221053672121387));
    LShoulder = new Station(*torso, Vec3(0.0029260590487432125, 0.32135920356193792, -0.14603793557074141));
    MidHip = new Station(*pelvis, Vec3(-0.068983007808594787, -0.059890551290884941, -0.0026719473665435256));
    RHip = new Station(*femur_r, Vec3(0.01039116651321681, -0.00077433106334712232, 0.012761629729853174));
    LHip = new Station(*femur_l, Vec3(-0.0030674684341940317, -7.894472308711542e-05, -0.023379526614267787));
    RKnee = new Station(*femur_r, Vec3(-0.0057001608129239045, -0.41931702376120189, 0.008422639264188847));
    LKnee = new Station(*femur_l, Vec3(-0.0078508632082388563, -0.41903451190434715, -0.0018849402006098082));
    RAnkle = new Station(*tibia_r, Vec3(0.017612693078604913, -0.41042511166146645, -0.015892845800358768));
    LAnkle = new Station(*tibia_l, Vec3(0.019268688221482477, -0.41172659838098302, 0.0092228736827546509));
    RHeel = new Station(*calcn_r, Vec3(0.0032027230293412146, 0.028165419411338655, -0.0036159145431496897));
    LHeel = new Station(*calcn_l, Vec3(0.0047636780068757156, 0.035314361919311733, 0.0011655514384873999));
    RBigToe = new Station(*toes_r, Vec3(0.041732378213962917, -0.0029887608662524934, -0.023036242341002988));
    LBigToe = new Station(*toes_l, Vec3(0.041717404102602162, 0.0048129334313149197, 0.019597273711784879));
    RSmallToe = new Station(*toes_r, Vec3(-0.0051850627462585175, 0.010042883134173064, 0.057501473092942068));
    LSmallToe = new Station(*toes_l, Vec3(-5.6850859683432731e-05, 0.0036303328102300914, -0.057883223802813522));

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
    T ua[NU]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// Controls
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
    /// Add contact forces to appliedBodyForces
    /// Right foot
    Array<osim_double_adouble> Force_values_1_r = HC_1_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_2_r = HC_2_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_3_r = HC_3_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_4_r = HC_4_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_5_r = HC_5_r->getRecordValues(*state);
    Array<osim_double_adouble> Force_values_6_r = HC_6_r->getRecordValues(*state);
    SpatialVec GRF_1_r;
    GRF_1_r[0] = Vec3(Force_values_1_r[9], Force_values_1_r[10], Force_values_1_r[11]);
    GRF_1_r[1] = Vec3(Force_values_1_r[6], Force_values_1_r[7], Force_values_1_r[8]);
    SpatialVec GRF_2_r;
    GRF_2_r[0] = Vec3(Force_values_2_r[9], Force_values_2_r[10], Force_values_2_r[11]);
    GRF_2_r[1] = Vec3(Force_values_2_r[6], Force_values_2_r[7], Force_values_2_r[8]);
    SpatialVec GRF_3_r;
    GRF_3_r[0] = Vec3(Force_values_3_r[9], Force_values_3_r[10], Force_values_3_r[11]);
    GRF_3_r[1] = Vec3(Force_values_3_r[6], Force_values_3_r[7], Force_values_3_r[8]);
    SpatialVec GRF_4_r;
    GRF_4_r[0] = Vec3(Force_values_4_r[9], Force_values_4_r[10], Force_values_4_r[11]);
    GRF_4_r[1] = Vec3(Force_values_4_r[6], Force_values_4_r[7], Force_values_4_r[8]);
    SpatialVec GRF_5_r;
    GRF_5_r[0] = Vec3(Force_values_5_r[9], Force_values_5_r[10], Force_values_5_r[11]);
    GRF_5_r[1] = Vec3(Force_values_5_r[6], Force_values_5_r[7], Force_values_5_r[8]);
    SpatialVec GRF_6_r;
    GRF_6_r[0] = Vec3(Force_values_6_r[9], Force_values_6_r[10], Force_values_6_r[11]);
    GRF_6_r[1] = Vec3(Force_values_6_r[6], Force_values_6_r[7], Force_values_6_r[8]);
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
    GRF_1_l[0] = Vec3(Force_values_1_l[9], Force_values_1_l[10], Force_values_1_l[11]);
    GRF_1_l[1] = Vec3(Force_values_1_l[6], Force_values_1_l[7], Force_values_1_l[8]);
    SpatialVec GRF_2_l;
    GRF_2_l[0] = Vec3(Force_values_2_l[9], Force_values_2_l[10], Force_values_2_l[11]);
    GRF_2_l[1] = Vec3(Force_values_2_l[6], Force_values_2_l[7], Force_values_2_l[8]);
    SpatialVec GRF_3_l;
    GRF_3_l[0] = Vec3(Force_values_3_l[9], Force_values_3_l[10], Force_values_3_l[11]);
    GRF_3_l[1] = Vec3(Force_values_3_l[6], Force_values_3_l[7], Force_values_3_l[8]);
    SpatialVec GRF_4_l;
    GRF_4_l[0] = Vec3(Force_values_4_l[9], Force_values_4_l[10], Force_values_4_l[11]);
    GRF_4_l[1] = Vec3(Force_values_4_l[6], Force_values_4_l[7], Force_values_4_l[8]);
    SpatialVec GRF_5_l;
    GRF_5_l[0] = Vec3(Force_values_5_l[9], Force_values_5_l[10], Force_values_5_l[11]);
    GRF_5_l[1] = Vec3(Force_values_5_l[6], Force_values_5_l[7], Force_values_5_l[8]);
    SpatialVec GRF_6_l;
    GRF_6_l[0] = Vec3(Force_values_6_l[9], Force_values_6_l[10], Force_values_6_l[11]);
    GRF_6_l[1] = Vec3(Force_values_6_l[6], Force_values_6_l[7], Force_values_6_l[8]);
    int ncalcn_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
    int ntoes_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
    appliedBodyForces[ncalcn_l] = appliedBodyForces[ncalcn_l] + GRF_1_l + GRF_2_l + GRF_3_l + GRF_4_l;
    appliedBodyForces[ntoes_l] = appliedBodyForces[ntoes_l] + GRF_5_l + GRF_6_l;
    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
    Vector residualMobilityForces(ndof);
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

    // Residual forces in OpenSim order
    T res_os[ndof];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < ndof; ++i) res_os[i] =
            value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];

    int nc = 3;
    for (int i = 0; i < nc; ++i) res[0][i + NU + 0*nc] = value<T>(Neck_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 1*nc] = value<T>(RShoulder_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 2*nc] = value<T>(LShoulder_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 3*nc] = value<T>(MidHip_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 4*nc] = value<T>(RHip_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 5*nc] = value<T>(LHip_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 6*nc] = value<T>(RKnee_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 7*nc] = value<T>(LKnee_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 8*nc] = value<T>(RAnkle_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 9*nc] = value<T>(LAnkle_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 10*nc] = value<T>(RHeel_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 11*nc] = value<T>(LHeel_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 12*nc] = value<T>(RBigToe_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 13*nc] = value<T>(LBigToe_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 14*nc] = value<T>(RSmallToe_location[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 15*nc] = value<T>(LSmallToe_location[i]);

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
    }

    Recorder::stop_recording();

    return 0;

}
