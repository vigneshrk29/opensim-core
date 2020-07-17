///*  This code describes the OpenSim model and the skeleton dynamics
//    Author: Antoine Falisse
//    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
//*/
//#include <OpenSim/Simulation/Model/Model.h>
//#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
//#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/BallJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/GimbalJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/FreeJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
//#include <OpenSim/Common/LinearFunction.h>
//#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Common/SimmSpline.h>
//#include <OpenSim/Common/MultiplierFunction.h>
//#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
//
//#include <iostream>
//#include <iterator>
//#include <random>
//#include <cassert>
//#include <algorithm>
//#include <vector>
//#include <fstream>
//
//using namespace SimTK;
//using namespace OpenSim;
//
///*  The function F describes the OpenSim model and, implicitly, the skeleton
//    dynamics. F takes as inputs joint positions and velocities (states x),
//    joint accelerations (controls u), and returns the joint torques as well as
//    several variables for use in the optimal control problems. F is templatized
//    using type T. F(x,u)->(r).
//*/
//
//// Inputs/outputs of function F
///// number of vectors in inputs/outputs of function F
//constexpr int n_in = 2;
//constexpr int n_out = 1;
///// number of elements in input/output vectors of function F
//constexpr int ndof = 6;        // # degrees of freedom (excluding locked)
//constexpr int ndofr = ndof;   // # degrees of freedom (including locked)
//constexpr int NX = ndof*2;      // # states
//constexpr int NU = ndof;        // # controls
//constexpr int NR = ndof;    // # residual torques + # joint origins
//
//// OpenSim and Simbody use different indices for the states/controls when the
//// kinematic chain has joints up and down the origin (e.g., lumbar joint/arms
//// and legs with pelvis as origin).
//// The two following functions allow getting the indices from one reference
//// system to the other. These functions are inspired from
//// createSystemYIndexMap() in Moco.
//// getIndicesOSInSimbody() returns the indices of the OpenSim Qs in the Simbody
//// reference system. Note that we only care about the order here so we divide
//// by 2 because the states include both Qs and Qdots.
//SimTK::Array_<int> getIndicesOSInSimbody(const Model& model) {
//    auto s = model.getWorkingState();
//    const auto svNames = model.getStateVariableNames();
//    SimTK::Array_<int> idxOSInSimbody(s.getNQ());
//    s.updQ() = 0;
//    for (int iy = 0; iy < s.getNQ(); ++iy) {
//        s.updQ()[iy] = SimTK::NaN;
//        const auto svValues = model.getStateVariableValues(s);
//        for (int isv = 0; isv < svNames.size(); ++isv) {
//            if (SimTK::isNaN(svValues[isv])) {
//                s.updQ()[iy] = 0;
//                idxOSInSimbody[iy] = isv/2;
//                break;
//            }
//        }
//    }
//    return idxOSInSimbody;
//}
//// getIndicesSimbodyInOS() returns the indices of the Simbody Qs in the OpenSim
//// reference system.
//SimTK::Array_<int> getIndicesSimbodyInOS(const Model& model) {
//    auto idxOSInSimbody = getIndicesOSInSimbody(model);
//    auto s = model.getWorkingState();
//    SimTK::Array_<int> idxSimbodyInOS(s.getNQ());
//	for (int iy = 0; iy < s.getNQ(); ++iy) {
//		for (int iyy = 0; iyy < s.getNQ(); ++iyy) {
//			if (idxOSInSimbody[iyy] == iy) {
//				idxSimbodyInOS[iy] = iyy;
//				break;
//			}
//		}
//	}
//    return idxSimbodyInOS;
//}
//
////// This function returns the linear acceleration of the imu frame wrt the
////// ground frame expressed in the ground frame.
////// frameName is the name of the body segment to which the imu is attached.
////// imuTranslation_B is the translation offset of the imu frame'origin from the
////// parent (body segment) frame's origin, expressed in the parent frame.
//const SimTK::Vec3 getLinearAccelerationIMUInGround(Model& model, const State& s, const Vector_<SpatialVec>& A_GB , const std::string& frameName, const Vec3& gravity_G, const Vec3& imuTranslation_B)
//{
//        SimTK::Rotation R_GB = model.getBodySet().get(frameName).getMobilizedBody().getBodyTransform(s).R();
//        /// Body linear acceleration in ground
//        SimTK::Vec3 linAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][1];
//         std::cout << "linesAcc: " << linAcc_G << std::endl;
//        /// Body angular acceleration in ground
//        SimTK::Vec3 angAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][0];
//        std::cout << "angAcc: " << angAcc_G << std::endl;
//        /// Body angular velocity in ground
//        SimTK::Vec3 angVel_G = model.getBodySet().get(frameName).getAngularVelocityInGround(s);
//        /// Body angular velocity in body
//        SimTK::Vec3 angVel_B = ~R_GB*angVel_G;
//        /// Body angular acceleration in body
//        SimTK::Vec3 angAcc_B = ~R_GB*angAcc_G;
//        /// Sensor linear acceleration
//        /// See van den Bogert et al. (1995), equation (1), p949.
//        SimTK::Vec3 linAcc_imu_B = ~R_GB * (linAcc_G - gravity_G) + SimTK::cross(angAcc_B, imuTranslation_B) + SimTK::cross(angVel_B, SimTK::cross(angVel_B, imuTranslation_B));
//        SimTK::Vec3 linAcc_imu_G = R_GB * linAcc_imu_B;
//
//        return linAcc_imu_G;
//}
//
//// Function F
//int F_generic(const double** arg, double** res) {
//
//    // OpenSim model: create components
//    /// Model
//    OpenSim::Model* model;
//    /// Bodies
//    OpenSim::Body* pelvis;
//    //OpenSim::Body* femur_r;
//    //OpenSim::Body* femur_l;
//    //OpenSim::Body* tibia_r;
//    //OpenSim::Body* tibia_l;
//    //OpenSim::Body* talus_r;
//    //OpenSim::Body* talus_l;
//    //OpenSim::Body* calcn_r;
//    //OpenSim::Body* calcn_l;
//    //OpenSim::Body* toes_r;
//    //OpenSim::Body* toes_l;
//    //OpenSim::Body* torso;
//    //OpenSim::Body* humerus_r;
//    //OpenSim::Body* humerus_l;
//    //OpenSim::Body* ulna_r;
//    //OpenSim::Body* ulna_l;
//    //OpenSim::Body* radius_r;
//    //OpenSim::Body* radius_l;
//    //OpenSim::Body* hand_r;
//    //OpenSim::Body* hand_l;
//    /// Joints
//    OpenSim::CustomJoint* ground_pelvis;
//    //OpenSim::PinJoint* ground_pelvis;
//    //OpenSim::PlanarJoint* ground_pelvis;
//    //OpenSim::BallJoint* ground_pelvis;
//    //OpenSim::GimbalJoint* ground_pelvis;
//    //OpenSim::FreeJoint* ground_pelvis;
//    //OpenSim::CustomJoint* hip_r;
//    //OpenSim::CustomJoint* hip_l;
//    //OpenSim::CustomJoint* knee_r;
//    //OpenSim::CustomJoint* knee_l;
//    //OpenSim::CustomJoint* ankle_r;
//    //OpenSim::CustomJoint* ankle_l;
//    //OpenSim::PinJoint* subtalar_r;
//    //OpenSim::PinJoint* subtalar_l;
//    //OpenSim::PinJoint* mtp_r;
//    //OpenSim::PinJoint* mtp_l;
//    //OpenSim::CustomJoint* back;
//    //OpenSim::CustomJoint* shoulder_r;
//    //OpenSim::CustomJoint* shoulder_l;
//    //OpenSim::CustomJoint* elbow_r;
//    //OpenSim::CustomJoint* elbow_l;
//    //OpenSim::CustomJoint* radioulnar_r;
//    //OpenSim::CustomJoint* radioulnar_l;
//    //OpenSim::WeldJoint* radius_hand_r;
//    //OpenSim::WeldJoint* radius_hand_l;
//
//    // OpenSim model: initialize components
//    /// Model
//    model = new OpenSim::Model();
//    //model->setUseVisualizer(true);
//    /// Body specifications
//    pelvis = new OpenSim::Body("pelvis", 11.751210011095651, Vec3(0, 0, 0), Inertia(0.099778065737821386, 0.08453958682650041, 0.056197957258948036, 0, 0, 0));
//    //femur_l = new OpenSim::Body("femur_l", 9.2810312301269491, Vec3(0, -0.17281329846170712, 0), Inertia(0.13806543520006995, 0.036191910198076598, 0.14559252763442779, 0, 0, 0));
//    //femur_r = new OpenSim::Body("femur_r", 9.2810312301269491, Vec3(0, -0.17281329846170712, 0), Inertia(0.13806543520006995, 0.036191910198076598, 0.14559252763442779, 0, 0, 0));
//    //tibia_l = new OpenSim::Body("tibia_l", 3.6993810916309018, Vec3(0, -0.20693800230674214, 0), Inertia(0.061783188539925718, 0.0062518702689210552, 0.062641288380758026, 0, 0, 0));
//    //tibia_r = new OpenSim::Body("tibia_r", 3.6993810916309018, Vec3(0, -0.20693800230674214, 0), Inertia(0.061783188539925718, 0.0062518702689210552, 0.062641288380758026, 0, 0, 0));
//    //talus_l = new OpenSim::Body("talus_l", 0.099781013934751236, Vec3(0, 0, 0), Inertia(0.0011448780870697802, 0.0011448780870697802, 0.0011448780870697802, 0, 0, 0));
//    //talus_r = new OpenSim::Body("talus_r", 0.099781013934751236, Vec3(0, 0, 0), Inertia(0.0011448780870697802, 0.0011448780870697802, 0.0011448780870697802, 0, 0, 0));
//    //calcn_l = new OpenSim::Body("calcn_l", 1.2472626741843904, Vec3(0.1071163252191219, 0.032134897565736564, 0), Inertia(0.0016028293218976922, 0.0044650245395721419, 0.0046940001569860989, 0, 0, 0));
//    //calcn_r = new OpenSim::Body("calcn_r", 1.2472626741843904, Vec3(0.1071163252191219, 0.032134897565736564, 0), Inertia(0.0016028293218976922, 0.0044650245395721419, 0.0046940001569860989, 0, 0, 0));
//    //toes_l = new OpenSim::Body("toes_l", 0.21612567618267114, Vec3(0.037062248525816174, 0.006426979513147313, 0.018745356913346334), Inertia(0.000114487808706978, 0.000228975617413956, 0.000114487808706978, 0, 0, 0));
//    //toes_r = new OpenSim::Body("toes_r", 0.21612567618267114, Vec3(0.037062248525816174, 0.006426979513147313, -0.018745356913346334), Inertia(0.000114487808706978, 0.000228975617413956, 0.000114487808706978, 0, 0, 0));
//    //torso = new OpenSim::Body("torso", 26.767853484219973, Vec3(-0.029119398106347269, 0.31060691313437089, 0), Inertia(1.3861651249269105, 0.7102392349150769, 1.345647175191848, 0, 0, 0));
//    //humerus_l = new OpenSim::Body("humerus_l", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
//    //humerus_r = new OpenSim::Body("humerus_r", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
//    //ulna_l = new OpenSim::Body("ulna_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
//    //ulna_r = new OpenSim::Body("ulna_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
//    //radius_l = new OpenSim::Body("radius_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
//    //radius_r = new OpenSim::Body("radius_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
//    //hand_l = new OpenSim::Body("hand_l", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
//    //hand_r = new OpenSim::Body("hand_r", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
//    /// Joints
//    /// Ground-Pelvis transform
//    //SpatialTransform st_ground_pelvis;
//    //st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
//    //st_ground_pelvis[0].setFunction(new LinearFunction());
//    //st_ground_pelvis[0].setAxis(Vec3(1, 0, 0));
//    //st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
//    //st_ground_pelvis[1].setFunction(new LinearFunction());
//    //st_ground_pelvis[1].setAxis(Vec3(0, 1, 0));
//    //st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
//    //st_ground_pelvis[2].setFunction(new LinearFunction());
//    //st_ground_pelvis[2].setAxis(Vec3(0, 0, 1));
//    //st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
//    //st_ground_pelvis[0].setFunction(new LinearFunction());
//    //st_ground_pelvis[0].setAxis(Vec3(0, 0, 1));
//    //st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
//    //st_ground_pelvis[1].setFunction(new LinearFunction());
//    //st_ground_pelvis[1].setAxis(Vec3(1, 0, 0));
//    //st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
//    //st_ground_pelvis[2].setFunction(new LinearFunction());
//    //st_ground_pelvis[2].setAxis(Vec3(0, 1, 0));
//    //st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
//    //st_ground_pelvis[2].setFunction(new LinearFunction());
//    //st_ground_pelvis[2].setAxis(Vec3(0, 0, 1));
//    //st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
//    //st_ground_pelvis[3].setFunction(new LinearFunction());
//    //st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
//    //st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
//    //st_ground_pelvis[4].setFunction(new LinearFunction());
//    //st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));
//    //SpatialTransform st_ground_pelvis;
//    //st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
//    //st_ground_pelvis[0].setFunction(new LinearFunction());
//    //st_ground_pelvis[0].setAxis(Vec3(1, 0, 0));
//    //st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_list", 1, 1));
//    //st_ground_pelvis[1].setFunction(new LinearFunction());
//    //st_ground_pelvis[1].setAxis(Vec3(0, 1, 0));
//    SpatialTransform st_ground_pelvis;
//    st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
//    st_ground_pelvis[0].setFunction(new LinearFunction());
//    st_ground_pelvis[0].setAxis(Vec3(0, 0, 1));
//    st_ground_pelvis[1].setCoordinateNames(OpenSim::Array<std::string>("pelvis_list", 1, 1));
//    st_ground_pelvis[1].setFunction(new LinearFunction());
//    st_ground_pelvis[1].setAxis(Vec3(1, 0, 0));
//    st_ground_pelvis[2].setCoordinateNames(OpenSim::Array<std::string>("pelvis_rotation", 1, 1));
//    st_ground_pelvis[2].setFunction(new LinearFunction());
//    st_ground_pelvis[2].setAxis(Vec3(0, 1, 0));
//    st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
//    st_ground_pelvis[3].setFunction(new LinearFunction());
//    st_ground_pelvis[3].setAxis(Vec3(1, 0, 0));
//    st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
//    st_ground_pelvis[4].setFunction(new LinearFunction());
//    st_ground_pelvis[4].setAxis(Vec3(0, 1, 0));
//    st_ground_pelvis[5].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tz", 1, 1));
//    st_ground_pelvis[5].setFunction(new LinearFunction());
//    st_ground_pelvis[5].setAxis(Vec3(0, 0, 1));
//   /* /// Hip_l transform
//    SpatialTransform st_hip_l;
//    st_hip_l[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_l", 1, 1));
//    st_hip_l[0].setFunction(new LinearFunction());
//    st_hip_l[0].setAxis(Vec3(0, 0, 1));
//    st_hip_l[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_l", 1, 1));
//    st_hip_l[1].setFunction(new LinearFunction());
//    st_hip_l[1].setAxis(Vec3(-1, 0, 0));
//    st_hip_l[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_l", 1, 1));
//    st_hip_l[2].setFunction(new LinearFunction());
//    st_hip_l[2].setAxis(Vec3(0, -1, 0));
//    /// Hip_r transform
//    SpatialTransform st_hip_r;
//    st_hip_r[0].setCoordinateNames(OpenSim::Array<std::string>("hip_flexion_r", 1, 1));
//    st_hip_r[0].setFunction(new LinearFunction());
//    st_hip_r[0].setAxis(Vec3(0, 0, 1));
//    st_hip_r[1].setCoordinateNames(OpenSim::Array<std::string>("hip_adduction_r", 1, 1));
//    st_hip_r[1].setFunction(new LinearFunction());
//    st_hip_r[1].setAxis(Vec3(1, 0, 0));
//    st_hip_r[2].setCoordinateNames(OpenSim::Array<std::string>("hip_rotation_r", 1, 1));
//    st_hip_r[2].setFunction(new LinearFunction());
//    st_hip_r[2].setAxis(Vec3(0, 1, 0));
//    /// Knee_l transform
//    SpatialTransform st_knee_l;
//    st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
//    st_knee_l[0].setFunction(new LinearFunction());
//    st_knee_l[0].setAxis(Vec3(0,0,1));
//    st_knee_l[1].setFunction(new Constant(0));
//    st_knee_l[1].setAxis(Vec3(0,1,0));
//    st_knee_l[2].setFunction(new Constant(0));
//    st_knee_l[2].setAxis(Vec3(1,0,0));
//    st_knee_l[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
//    double knee_X_r_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
//    double knee_X_r_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
//    OpenSim::SimmSpline* knee_X_r = new SimmSpline(12, knee_X_r_x, knee_X_r_y, "function_X");
//    st_knee_l[3].setFunction(new MultiplierFunction(knee_X_r, 1.0165488144806301));
//    st_knee_l[3].setAxis(Vec3(1,0,0));
//    st_knee_l[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
//    double knee_Y_r_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
//    double knee_Y_r_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
//    OpenSim::SimmSpline* knee_Y_r = new SimmSpline(7, knee_Y_r_x, knee_Y_r_y, "function_Y");
//    st_knee_l[4].setFunction(new MultiplierFunction(knee_Y_r, 1.0165488144806301));
//    st_knee_l[4].setAxis(Vec3(0,1,0));
//    st_knee_l[5].setFunction(new Constant(0));
//    st_knee_l[5].setAxis(Vec3(0,0,1));
//    /// Knee_r transform
//    SpatialTransform st_knee_r;
//    st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
//    st_knee_r[0].setFunction(new LinearFunction());
//    st_knee_r[0].setAxis(Vec3(0,0,1));
//    st_knee_r[1].setFunction(new Constant(0));
//    st_knee_r[1].setAxis(Vec3(0,1,0));
//    st_knee_r[2].setFunction(new Constant(0));
//    st_knee_r[2].setAxis(Vec3(1,0,0));
//    st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
//    double knee_X_l_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
//    double knee_X_l_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
//    OpenSim::SimmSpline* knee_X_l = new SimmSpline(12, knee_X_l_x, knee_X_l_y, "function_X");
//    st_knee_r[3].setFunction(new MultiplierFunction(knee_X_l, 1.0165488144806301));
//    st_knee_r[3].setAxis(Vec3(1,0,0));
//    st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
//    double knee_Y_l_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
//    double knee_Y_l_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
//    OpenSim::SimmSpline* knee_Y_l = new SimmSpline(7, knee_Y_l_x, knee_Y_l_y, "function_Y");
//    st_knee_r[4].setFunction(new MultiplierFunction(knee_Y_l, 1.0165488144806301));
//    st_knee_r[4].setAxis(Vec3(0,1,0));
//    st_knee_r[5].setFunction(new Constant(0));
//    st_knee_r[5].setAxis(Vec3(0,0,1));
//    /// Ankle_l transform
//    SpatialTransform st_ankle_l;
//    st_ankle_l[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_l", 1, 1));
//    st_ankle_l[0].setFunction(new LinearFunction());
//    st_ankle_l[0].setAxis(Vec3(0.10501355, 0.17402245, 0.97912632));
//    /// Ankle_r transform
//    SpatialTransform st_ankle_r;
//    st_ankle_r[0].setCoordinateNames(OpenSim::Array<std::string>("ankle_angle_r", 1, 1));
//    st_ankle_r[0].setFunction(new LinearFunction());
//    st_ankle_r[0].setAxis(Vec3(-0.10501355, -0.17402245, 0.97912632));
//    /// Back transform
//    SpatialTransform st_back;
//    st_back[0].setCoordinateNames(OpenSim::Array<std::string>("lumbar_extension", 1, 1));
//    st_back[0].setFunction(new LinearFunction());
//    st_back[0].setAxis(Vec3(0, 0, 1));
//    st_back[1].setCoordinateNames(OpenSim::Array<std::string>("lumbar_bending", 1, 1));
//    st_back[1].setFunction(new LinearFunction());
//    st_back[1].setAxis(Vec3(1, 0, 0));
//    st_back[2].setCoordinateNames(OpenSim::Array<std::string>("lumbar_rotation", 1, 1));
//    st_back[2].setFunction(new LinearFunction());
//    st_back[2].setAxis(Vec3(0, 1, 0));
//    /// Shoulder_l transform
//    SpatialTransform st_sho_l;
//    st_sho_l[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_l", 1, 1));
//    st_sho_l[0].setFunction(new LinearFunction());
//    st_sho_l[0].setAxis(Vec3(0, 0, 1));
//    st_sho_l[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_l", 1, 1));
//    st_sho_l[1].setFunction(new LinearFunction());
//    st_sho_l[1].setAxis(Vec3(-1, 0, 0));
//    st_sho_l[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_l", 1, 1));
//    st_sho_l[2].setFunction(new LinearFunction());
//    st_sho_l[2].setAxis(Vec3(0, -1, 0));
//    /// Shoulder_r transform
//    SpatialTransform st_sho_r;
//    st_sho_r[0].setCoordinateNames(OpenSim::Array<std::string>("arm_flex_r", 1, 1));
//    st_sho_r[0].setFunction(new LinearFunction());
//    st_sho_r[0].setAxis(Vec3(0, 0, 1));
//    st_sho_r[1].setCoordinateNames(OpenSim::Array<std::string>("arm_add_r", 1, 1));
//    st_sho_r[1].setFunction(new LinearFunction());
//    st_sho_r[1].setAxis(Vec3(1, 0, 0));
//    st_sho_r[2].setCoordinateNames(OpenSim::Array<std::string>("arm_rot_r", 1, 1));
//    st_sho_r[2].setFunction(new LinearFunction());
//    st_sho_r[2].setAxis(Vec3(0, 1, 0));
//    /// Elbow_l transform
//    SpatialTransform st_elb_l;
//    st_elb_l[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_l", 1, 1));
//    st_elb_l[0].setFunction(new LinearFunction());
//    st_elb_l[0].setAxis(Vec3(-0.22604696, -0.022269, 0.97386183));
//    /// Elbow_r transform
//    SpatialTransform st_elb_r;
//    st_elb_r[0].setCoordinateNames(OpenSim::Array<std::string>("elbow_flex_r", 1, 1));
//    st_elb_r[0].setFunction(new LinearFunction());
//    st_elb_r[0].setAxis(Vec3(0.22604696, 0.022269, 0.97386183));
//    /// Radioulnar_l transform
//    SpatialTransform st_radioulnar_l;
//    st_radioulnar_l[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_l", 1, 1));
//    st_radioulnar_l[0].setFunction(new LinearFunction());
//    st_radioulnar_l[0].setAxis(Vec3(-0.05639803, -0.99840646, 0.001952));
//    /// Radioulnar_r transform
//    SpatialTransform st_radioulnar_r;
//    st_radioulnar_r[0].setCoordinateNames(OpenSim::Array<std::string>("pro_sup_r", 1, 1));
//    st_radioulnar_r[0].setFunction(new LinearFunction());
//    st_radioulnar_r[0].setAxis(Vec3(0.05639803, 0.99840646, 0.001952));*/
//    /// Joint specifications
//    ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0,5,0), Vec3(0), st_ground_pelvis);
//    //ground_pelvis = new PinJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0,5,0), Vec3(0));
//    /*ground_pelvis = new PlanarJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0,5,0), Vec3(0));*/
//    //ground_pelvis = new BallJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0,5,0), Vec3(0));
//    //ground_pelvis = new GimbalJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0,5,0), Vec3(0));
//    //ground_pelvis = new FreeJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0,5,0), Vec3(0));
//    //hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.069729482228687481, -0.065192627656523949, -0.08235377321209908), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
//    //hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.069729482228687481, -0.065192627656523949, 0.08235377321209908), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
//    //knee_l = new CustomJoint("knee_l", *femur_l, Vec3(0, 0, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
//    //knee_r = new CustomJoint("knee_r", *femur_r, Vec3(0, 0, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
//    //ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.47661136042795454, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
//    //ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.47661136042795454, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
//    //subtalar_l = new PinJoint("subtalar_l", *talus_l, Vec3(-0.052240631809365744, -0.044935298429421636, -0.0084836129573544541), Vec3(1.7681899999999999, -0.906223, 1.8196000000000001),
//    //                            *calcn_l, Vec3(0), Vec3(1.7681899999999999, -0.906223, 1.8196000000000001));
//    //subtalar_r = new PinJoint("subtalar_r", *talus_r, Vec3(-0.052240631809365744, -0.044935298429421636, 0.0084836129573544541), Vec3(-1.7681899999999999, 0.906223, 1.8196000000000001),
//    //                            *calcn_r, Vec3(0), Vec3(-1.7681899999999999, 0.906223, 1.8196000000000001));
//    //mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.19152398949178992, -0.0021423265043824377, -0.0011568563123665165), Vec3(-3.1415899999999999, -0.61990100000000004, 0),
//    //                            *toes_l, Vec3(0), Vec3(-3.1415899999999999, -0.61990100000000004, 0));
//    //mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.19152398949178992, -0.0021423265043824377, 0.0011568563123665165), Vec3(-3.1415899999999999, 0.61990100000000004, 0),
//    //                            *toes_r, Vec3(0), Vec3(-3.1415899999999999, 0.61990100000000004, 0));
//    //back = new CustomJoint("back", *pelvis, Vec3(-0.099317664221058408, 0.080381227745941017, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
//    //shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, -0.16500992260263453), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
//    //shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, 0.16500992260263453), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
//    //elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.014421812854093517, -0.31410344120358441, 0.010527791717515771), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
//    //elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.014421812854093517, -0.31410344120358441, -0.010527791717515771), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
//    //radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.0078589566368657427, -0.015195696294888168, -0.030472003264362887), Vec3(0), *radius_l, Vec3(0), Vec3(0),st_radioulnar_l);
//    //radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.0078589566368657427, -0.015195696294888168, 0.030472003264362887), Vec3(0), *radius_r, Vec3(0), Vec3(0),st_radioulnar_r);
//    //radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.010277276874462305, -0.27552611746618899, -0.015900163494535866), Vec3(0), *hand_l, Vec3(0), Vec3(0));
//    //radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.010277276874462305, -0.27552611746618899, 0.015900163494535866), Vec3(0), *hand_r, Vec3(0), Vec3(0));
//    /// Add bodies and joints to model
//    model->addBody(pelvis);		    model->addJoint(ground_pelvis);
//    //model->addBody(femur_l);		model->addJoint(hip_l);
//    //model->addBody(femur_r);		model->addJoint(hip_r);
//    //model->addBody(tibia_l);		model->addJoint(knee_l);
//    //model->addBody(tibia_r);		model->addJoint(knee_r);
//    //model->addBody(talus_l);		model->addJoint(ankle_l);
//    //model->addBody(talus_r);		model->addJoint(ankle_r);
//    //model->addBody(calcn_l);		model->addJoint(subtalar_l);
//    //model->addBody(calcn_r);		model->addJoint(subtalar_r);
//    //model->addBody(toes_l);		    model->addJoint(mtp_l);
//    //model->addBody(toes_r);		    model->addJoint(mtp_r);
//    //model->addBody(torso);          model->addJoint(back);
//    //model->addBody(humerus_l);      model->addJoint(shoulder_l);
//    //model->addBody(humerus_r);      model->addJoint(shoulder_r);
//    //model->addBody(ulna_l);         model->addJoint(elbow_l);
//    //model->addBody(ulna_r);         model->addJoint(elbow_r);
//    //model->addBody(radius_l);       model->addJoint(radioulnar_l);
//    //model->addBody(radius_r);       model->addJoint(radioulnar_r);
//    //model->addBody(hand_l);         model->addJoint(radius_hand_l);
//    //model->addBody(hand_r);         model->addJoint(radius_hand_r);
//   /* /// Contact elements
//    /// Parameters
//    double radiusSphere_s1 = 0.03232;
//    double radiusSphere_s2 = 0.03232;
//    double radiusSphere_s3 = 0.023374;
//    double radiusSphere_s4 = 0.020508;
//    double radiusSphere_s5 = 0.016244;
//    double radiusSphere_s6 = 0.018414;
//    double stiffness = 1000000;
//    double dissipation = 2.0;
//    double staticFriction = 0.8;
//    double dynamicFriction = 0.8;
//    double viscousFriction = 0.5;
//    double transitionVelocity = 0.2;
//    Vec3 normal = Vec3(0, 1, 0);
//    double offset = 0;
//    Vec3 locSphere_s1_r(-0.00042152, 0, -0.0049972);
//    Vec3 locSphere_s2_r(0.06, 0, 0.020001);
//    Vec3 locSphere_s3_r(0.165, -0.01, 0.021183);
//    Vec3 locSphere_s4_r(0.18, -0.01, -0.01);
//    Vec3 locSphere_s5_r(0.053154, -0.01, -0.0034173);
//    Vec3 locSphere_s6_r(0.01, -0.01, -0.015);
//    Vec3 locSphere_s1_l(locSphere_s1_r[0],locSphere_s1_r[1],-locSphere_s1_r[2]);
//    Vec3 locSphere_s2_l(locSphere_s2_r[0],locSphere_s2_r[1],-locSphere_s2_r[2]);
//    Vec3 locSphere_s3_l(locSphere_s3_r[0],locSphere_s3_r[1],-locSphere_s3_r[2]);
//    Vec3 locSphere_s4_l(locSphere_s4_r[0],locSphere_s4_r[1],-locSphere_s4_r[2]);
//    Vec3 locSphere_s5_l(locSphere_s5_r[0],locSphere_s5_r[1],-locSphere_s5_r[2]);
//    Vec3 locSphere_s6_l(locSphere_s6_r[0],locSphere_s6_r[1],-locSphere_s6_r[2]);
//    /// Left foot contact shere specifications
//
//    Vec3 halfSpaceLocation(0);
//	Vec3 halfSpaceOrientation(0, 0, -0.5 * SimTK::Pi);
//    OpenSim::ContactHalfSpace* contactHalfSpace;
//    contactHalfSpace = new OpenSim::ContactHalfSpace(halfSpaceLocation,halfSpaceOrientation,model->getGround(),"contactHalfSpace");
//    model->addComponent(contactHalfSpace);
//
//    OpenSim::ContactSphere* s1_l;
//    s1_l = new OpenSim::ContactSphere(radiusSphere_s1,locSphere_s1_l,*calcn_l,"s1_l");
//    model->addComponent(s1_l);
//    OpenSim::ContactSphere* s2_l;
//    s2_l = new OpenSim::ContactSphere(radiusSphere_s2,locSphere_s2_l,*calcn_l,"s2_l");
//    model->addComponent(s2_l);
//    OpenSim::ContactSphere* s3_l;
//    s3_l = new OpenSim::ContactSphere(radiusSphere_s3,locSphere_s3_l,*calcn_l,"s3_l");
//    model->addComponent(s3_l);
//    OpenSim::ContactSphere* s4_l;
//    s4_l = new OpenSim::ContactSphere(radiusSphere_s4,locSphere_s4_l,*calcn_l,"s4_l");
//    model->addComponent(s4_l);
//    OpenSim::ContactSphere* s5_l;
//    s5_l = new OpenSim::ContactSphere(radiusSphere_s5,locSphere_s5_l,*toes_l,"s5_l");
//    model->addComponent(s5_l);
//    OpenSim::ContactSphere* s6_l;
//    s6_l = new OpenSim::ContactSphere(radiusSphere_s6,locSphere_s6_l,*toes_l,"s6_l");
//    model->addComponent(s6_l);
//    OpenSim::ContactSphere* s1_r;
//    s1_r = new OpenSim::ContactSphere(radiusSphere_s1,locSphere_s1_r,*calcn_r,"s1_r");
//    model->addComponent(s1_r);
//    OpenSim::ContactSphere* s2_r;
//    s2_r = new OpenSim::ContactSphere(radiusSphere_s2,locSphere_s2_r,*calcn_r,"s2_r");
//    model->addComponent(s2_r);
//    OpenSim::ContactSphere* s3_r;
//    s3_r = new OpenSim::ContactSphere(radiusSphere_s3,locSphere_s3_r,*calcn_r,"s3_r");
//    model->addComponent(s3_r);
//    OpenSim::ContactSphere* s4_r;
//    s4_r = new OpenSim::ContactSphere(radiusSphere_s4,locSphere_s4_r,*calcn_r,"s4_r");
//    model->addComponent(s4_r);
//    OpenSim::ContactSphere* s5_r;
//    s5_r = new OpenSim::ContactSphere(radiusSphere_s5,locSphere_s5_r,*toes_r,"s5_r");
//    model->addComponent(s5_r);
//    OpenSim::ContactSphere* s6_r;
//    s6_r = new OpenSim::ContactSphere(radiusSphere_s6,locSphere_s6_r,*toes_r,"s6_r");
//    model->addComponent(s6_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s1_r;
//    contact_s1_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s1_r",*s1_r,*contactHalfSpace);
//    contact_s1_r->set_stiffness(stiffness);
//    contact_s1_r->set_dissipation(dissipation);
//    contact_s1_r->set_static_friction(staticFriction);
//    contact_s1_r->set_dynamic_friction(dynamicFriction);
//    contact_s1_r->set_viscous_friction(viscousFriction);
//    contact_s1_r->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s1_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s2_r;
//    contact_s2_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s2_r",*s2_r,*contactHalfSpace);
//    contact_s2_r->set_stiffness(stiffness);
//    contact_s2_r->set_dissipation(dissipation);
//    contact_s2_r->set_static_friction(staticFriction);
//    contact_s2_r->set_dynamic_friction(dynamicFriction);
//    contact_s2_r->set_viscous_friction(viscousFriction);
//    contact_s2_r->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s2_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s3_r;
//    contact_s3_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s3_r",*s3_r,*contactHalfSpace);
//    contact_s3_r->set_stiffness(stiffness);
//    contact_s3_r->set_dissipation(dissipation);
//    contact_s3_r->set_static_friction(staticFriction);
//    contact_s3_r->set_dynamic_friction(dynamicFriction);
//    contact_s3_r->set_viscous_friction(viscousFriction);
//    contact_s3_r->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s3_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s4_r;
//    contact_s4_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s4_r",*s4_r,*contactHalfSpace);
//    contact_s4_r->set_stiffness(stiffness);
//    contact_s4_r->set_dissipation(dissipation);
//    contact_s4_r->set_static_friction(staticFriction);
//    contact_s4_r->set_dynamic_friction(dynamicFriction);
//    contact_s4_r->set_viscous_friction(viscousFriction);
//    contact_s4_r->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s4_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s5_r;
//    contact_s5_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s5_r",*s5_r,*contactHalfSpace);
//    contact_s5_r->set_stiffness(stiffness);
//    contact_s5_r->set_dissipation(dissipation);
//    contact_s5_r->set_static_friction(staticFriction);
//    contact_s5_r->set_dynamic_friction(dynamicFriction);
//    contact_s5_r->set_viscous_friction(viscousFriction);
//    contact_s5_r->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s5_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s6_r;
//    contact_s6_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s6_r",*s6_r,*contactHalfSpace);
//    contact_s6_r->set_stiffness(stiffness);
//    contact_s6_r->set_dissipation(dissipation);
//    contact_s6_r->set_static_friction(staticFriction);
//    contact_s6_r->set_dynamic_friction(dynamicFriction);
//    contact_s6_r->set_viscous_friction(viscousFriction);
//    contact_s6_r->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s6_r);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s1_l;
//    contact_s1_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s1_l",*s1_l,*contactHalfSpace);
//    contact_s1_l->set_stiffness(stiffness);
//    contact_s1_l->set_dissipation(dissipation);
//    contact_s1_l->set_static_friction(staticFriction);
//    contact_s1_l->set_dynamic_friction(dynamicFriction);
//    contact_s1_l->set_viscous_friction(viscousFriction);
//    contact_s1_l->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s1_l);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s2_l;
//    contact_s2_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s2_l",*s2_l,*contactHalfSpace);
//    contact_s2_l->set_stiffness(stiffness);
//    contact_s2_l->set_dissipation(dissipation);
//    contact_s2_l->set_static_friction(staticFriction);
//    contact_s2_l->set_dynamic_friction(dynamicFriction);
//    contact_s2_l->set_viscous_friction(viscousFriction);
//    contact_s2_l->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s2_l);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s3_l;
//    contact_s3_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s3_l",*s3_l,*contactHalfSpace);
//    contact_s3_l->set_stiffness(stiffness);
//    contact_s3_l->set_dissipation(dissipation);
//    contact_s3_l->set_static_friction(staticFriction);
//    contact_s3_l->set_dynamic_friction(dynamicFriction);
//    contact_s3_l->set_viscous_friction(viscousFriction);
//    contact_s3_l->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s3_l);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s4_l;
//    contact_s4_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s4_l",*s4_l,*contactHalfSpace);
//    contact_s4_l->set_stiffness(stiffness);
//    contact_s4_l->set_dissipation(dissipation);
//    contact_s4_l->set_static_friction(staticFriction);
//    contact_s4_l->set_dynamic_friction(dynamicFriction);
//    contact_s4_l->set_viscous_friction(viscousFriction);
//    contact_s4_l->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s4_l);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s5_l;
//    contact_s5_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s5_l",*s5_l,*contactHalfSpace);
//    contact_s5_l->set_stiffness(stiffness);
//    contact_s5_l->set_dissipation(dissipation);
//    contact_s5_l->set_static_friction(staticFriction);
//    contact_s5_l->set_dynamic_friction(dynamicFriction);
//    contact_s5_l->set_viscous_friction(viscousFriction);
//    contact_s5_l->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s5_l);
//
//    OpenSim::SmoothSphereHalfSpaceForce* contact_s6_l;
//    contact_s6_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s6_l",*s6_l,*contactHalfSpace);
//    contact_s6_l->set_stiffness(stiffness);
//    contact_s6_l->set_dissipation(dissipation);
//    contact_s6_l->set_static_friction(staticFriction);
//    contact_s6_l->set_dynamic_friction(dynamicFriction);
//    contact_s6_l->set_viscous_friction(viscousFriction);
//    contact_s6_l->set_transition_velocity(transitionVelocity);
//    model->addComponent(contact_s6_l);*/
//
//    // Initialize system and state
//    SimTK::State* state;
//    state = new State(model->initSystem());
//
//    // Read inputs
//    std::vector<double> x(arg[0], arg[0] + NX);
//    std::vector<double> u(arg[1], arg[1] + NU);
//
//    // States and controls
//    double ua[NU]; /// joint accelerations (Qdotdots) - controls
//    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states
//
//    // Assign inputs to model variables
//    /// States
//    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
//    ///// pro_sup dofs are locked so Qs and Qdots are hard coded (0)
//    //QsUs[NX] = SimTK::Pi/2;
//    //QsUs[NX+1] = 0;
//    //QsUs[NX+2] = SimTK::Pi/2;
//    //QsUs[NX+3] = 0;
//    /// Controls
//    double ut[NU];
//    for (int i = 0; i < NU; ++i) ut[i] = u[i];
//    ///// pro_sup dofs are locked so Qdotdots are hard coded (0)
//    ///// Need to have a temporary vector to add 0s to the vector before
//    ///// adjusting for the index difference between OpenSim and Simbody.
//    //ut[NU] = 0;
//    //ut[NU+1] = 0;
//    /// OpenSim and Simbody have different state orders so we need to adjust
//    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
//    for (int i = 0; i < ndofr; ++i) ua[i] = ut[indicesOSInSimbody[i]];
//
//    // Set state variables and realize
//    model->setStateVariableValues(*state, QsUs);
//    model->realizeAcceleration(*state);
//
//    // Compute residual forces
//    /// appliedMobilityForces (# mobilities)
//    Vector appliedMobilityForces(ndofr);
//    appliedMobilityForces.setToZero();
//    /// appliedBodyForces (# bodies + ground)
//    Vector_<SpatialVec> appliedBodyForces;
//    int nbodies = model->getBodySet().getSize() + 1;
//    appliedBodyForces.resize(nbodies);
//    appliedBodyForces.setToZero();
//    /// Set gravity
//    Vec3 gravity(0);
//    gravity[1] = -9.80665;
//    /// Add weights to appliedBodyForces
//    for (int i = 0; i < model->getBodySet().getSize(); ++i) {
//        model->getMatterSubsystem().addInStationForce(*state,
//            model->getBodySet().get(i).getMobilizedBodyIndex(),
//            model->getBodySet().get(i).getMassCenter(),
//            model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
//    }
//    ///// Add contact forces to appliedBodyForces
//    ///// Right foot
//    //Array<double> Force_values_s1_r = contact_s1_r->getRecordValues(*state);
//    //Array<double> Force_values_s2_r = contact_s2_r->getRecordValues(*state);
//    //Array<double> Force_values_s3_r = contact_s3_r->getRecordValues(*state);
//    //Array<double> Force_values_s4_r = contact_s4_r->getRecordValues(*state);
//    //Array<double> Force_values_s5_r = contact_s5_r->getRecordValues(*state);
//    //Array<double> Force_values_s6_r = contact_s6_r->getRecordValues(*state);
//    //SpatialVec GRF_s1_r;
//    //GRF_s1_r[0] = Vec3(Force_values_s1_r[3], Force_values_s1_r[4], Force_values_s1_r[5]);
//    //GRF_s1_r[1] = Vec3(Force_values_s1_r[0], Force_values_s1_r[1], Force_values_s1_r[2]);
//    //SpatialVec GRF_s2_r;
//    //GRF_s2_r[0] = Vec3(Force_values_s2_r[3], Force_values_s2_r[4], Force_values_s2_r[5]);
//    //GRF_s2_r[1] = Vec3(Force_values_s2_r[0], Force_values_s2_r[1], Force_values_s2_r[2]);
//    //SpatialVec GRF_s3_r;
//    //GRF_s3_r[0] = Vec3(Force_values_s3_r[3], Force_values_s3_r[4], Force_values_s3_r[5]);
//    //GRF_s3_r[1] = Vec3(Force_values_s3_r[0], Force_values_s3_r[1], Force_values_s3_r[2]);
//    //SpatialVec GRF_s4_r;
//    //GRF_s4_r[0] = Vec3(Force_values_s4_r[3], Force_values_s4_r[4], Force_values_s4_r[5]);
//    //GRF_s4_r[1] = Vec3(Force_values_s4_r[0], Force_values_s4_r[1], Force_values_s4_r[2]);
//    //SpatialVec GRF_s5_r;
//    //GRF_s5_r[0] = Vec3(Force_values_s5_r[3], Force_values_s5_r[4], Force_values_s5_r[5]);
//    //GRF_s5_r[1] = Vec3(Force_values_s5_r[0], Force_values_s5_r[1], Force_values_s5_r[2]);
//    //SpatialVec GRF_s6_r;
//    //GRF_s6_r[0] = Vec3(Force_values_s6_r[3], Force_values_s6_r[4], Force_values_s6_r[5]);
//    //GRF_s6_r[1] = Vec3(Force_values_s6_r[0], Force_values_s6_r[1], Force_values_s6_r[2]);
//    //int nCalc_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
//    //int nToe_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
//    //appliedBodyForces[nCalc_r] = appliedBodyForces[nCalc_r] + GRF_s1_r + GRF_s2_r + GRF_s3_r + GRF_s4_r;
//    //appliedBodyForces[nToe_r] = appliedBodyForces[nToe_r] + GRF_s5_r + GRF_s6_r;
//
//    ///// Left foot
//    //Array<double> Force_values_s1_l = contact_s1_l->getRecordValues(*state);
//    //Array<double> Force_values_s2_l = contact_s2_l->getRecordValues(*state);
//    //Array<double> Force_values_s3_l = contact_s3_l->getRecordValues(*state);
//    //Array<double> Force_values_s4_l = contact_s4_l->getRecordValues(*state);
//    //Array<double> Force_values_s5_l = contact_s5_l->getRecordValues(*state);
//    //Array<double> Force_values_s6_l = contact_s6_l->getRecordValues(*state);
//    //SpatialVec GRF_s1_l;
//    //GRF_s1_l[0] = Vec3(Force_values_s1_l[3], Force_values_s1_l[4], Force_values_s1_l[5]);
//    //GRF_s1_l[1] = Vec3(Force_values_s1_l[0], Force_values_s1_l[1], Force_values_s1_l[2]);
//    //SpatialVec GRF_s2_l;
//    //GRF_s2_l[0] = Vec3(Force_values_s2_l[3], Force_values_s2_l[4], Force_values_s2_l[5]);
//    //GRF_s2_l[1] = Vec3(Force_values_s2_l[0], Force_values_s2_l[1], Force_values_s2_l[2]);
//    //SpatialVec GRF_s3_l;
//    //GRF_s3_l[0] = Vec3(Force_values_s3_l[3], Force_values_s3_l[4], Force_values_s3_l[5]);
//    //GRF_s3_l[1] = Vec3(Force_values_s3_l[0], Force_values_s3_l[1], Force_values_s3_l[2]);
//    //SpatialVec GRF_s4_l;
//    //GRF_s4_l[0] = Vec3(Force_values_s4_l[3], Force_values_s4_l[4], Force_values_s4_l[5]);
//    //GRF_s4_l[1] = Vec3(Force_values_s4_l[0], Force_values_s4_l[1], Force_values_s4_l[2]);
//    //SpatialVec GRF_s5_l;
//    //GRF_s5_l[0] = Vec3(Force_values_s5_l[3], Force_values_s5_l[4], Force_values_s5_l[5]);
//    //GRF_s5_l[1] = Vec3(Force_values_s5_l[0], Force_values_s5_l[1], Force_values_s5_l[2]);
//    //SpatialVec GRF_s6_l;
//    //GRF_s6_l[0] = Vec3(Force_values_s6_l[3], Force_values_s6_l[4], Force_values_s6_l[5]);
//    //GRF_s6_l[1] = Vec3(Force_values_s6_l[0], Force_values_s6_l[1], Force_values_s6_l[2]);
//    //int nCalc_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
//    //int nToe_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
//    //appliedBodyForces[nCalc_l] = appliedBodyForces[nCalc_l] + GRF_s1_l + GRF_s2_l + GRF_s3_l + GRF_s4_l;
//    //appliedBodyForces[nToe_l] = appliedBodyForces[nToe_l] + GRF_s5_l + GRF_s6_l;
//    /// knownUdot
//    Vector knownUdot(ndofr);
//    knownUdot.setToZero();
//    for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];
//    /// Calculate residual forces
//    Vector residualMobilityForces(ndofr);
//    residualMobilityForces.setToZero();
//    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
//        appliedMobilityForces, appliedBodyForces, knownUdot,
//        residualMobilityForces);
//
//    const SimTK::Vec3 translation_pelvis_imu(-0.17825090079013006, 0.06148338297319611, -0.0039742631657566363);
//    //const SimTK::Vec3 translation_torso_imu(0.11141304895632698, 0.32812980850924067, -0.012040552365984683);
//
//    //const SimTK::Vec3 translation_femur_r_imu(0.043933841399841023, -0.14958344693305592, 0.099866089437989247);
//    //const SimTK::Vec3 translation_femur_l_imu(0.057420458331616797, -0.12466095809783695, -0.1025425014551595);
//
//    //const SimTK::Vec3 translation_tibia_r_imu(0.047080175165178595, -0.11466976609963364, 0.0052928998697051033);
//    //const SimTK::Vec3 translation_tibia_l_imu(0.053945123225004443, -0.12392203671935126, -0.006707066365313652);
//
//    //const SimTK::Vec3 translation_calcn_r_imu(0.14808446197590852, 0.040621301301533644, 0.035088429802814902);
//    //const SimTK::Vec3 translation_calcn_l_imu(0.15253209992514105, 0.045078665693204636, -0.044293377012806667);
//
//    SimTK::Vec3 angVel_pelvis_imu_inG   = model->getBodySet().get("pelvis").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angPos_pelvis_imu_inG   = model->getBodySet().get("pelvis").getPositionInGround(*state);
//    std::cout << "q: " << state->getQ() << std::endl;
//    std::cout << "u: " << state->getU() << std::endl;
//    std::cout << "qdot: " << state->getQDot() << std::endl;
//    std::cout << "udot: " << state->getUDot() << std::endl;
//    std::cout << "theta: " << angPos_pelvis_imu_inG << std::endl;
//    std::cout << "omega: " << angVel_pelvis_imu_inG << std::endl;
//    //SimTK::Vec3 angVel_torso_imu_inG    = model->getBodySet().get("torso").getAngularVelocityInGround(*state);
//    //SimTK::Vec3 angVel_femur_l_imu_inG  = model->getBodySet().get("femur_l").getAngularVelocityInGround(*state);
//    //SimTK::Vec3 angVel_femur_r_imu_inG  = model->getBodySet().get("femur_r").getAngularVelocityInGround(*state);
//    //SimTK::Vec3 angVel_tibia_l_imu_inG  = model->getBodySet().get("tibia_l").getAngularVelocityInGround(*state);
//    //SimTK::Vec3 angVel_tibia_r_imu_inG  = model->getBodySet().get("tibia_r").getAngularVelocityInGround(*state);
//    //SimTK::Vec3 angVel_calcn_l_imu_inG  = model->getBodySet().get("calcn_l").getAngularVelocityInGround(*state);
//    //SimTK::Vec3 angVel_calcn_r_imu_inG  = model->getBodySet().get("calcn_r").getAngularVelocityInGround(*state);
//
//    Vector_<SpatialVec> A_GB(nbodies);
//    model->getMatterSubsystem().calcBodyAccelerationFromUDot(*state, knownUdot, A_GB);
//
//    SimTK::Vec3 linAcc_pelvis_imu_inG   = getLinearAccelerationIMUInGround(*model, *state, A_GB, "pelvis",  gravity, translation_pelvis_imu);
//    //SimTK::Vec3 linAcc_torso_imu_inG    = getLinearAccelerationIMUInGround(*model, *state, A_GB, "torso",   gravity, translation_torso_imu);
//    //SimTK::Vec3 linAcc_femur_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "femur_l", gravity, translation_femur_l_imu);
//    //SimTK::Vec3 linAcc_femur_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "femur_r", gravity, translation_femur_r_imu);
//    //SimTK::Vec3 linAcc_tibia_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "tibia_l", gravity, translation_tibia_l_imu);
//    //SimTK::Vec3 linAcc_tibia_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "tibia_r", gravity, translation_tibia_r_imu);
//    //SimTK::Vec3 linAcc_calcn_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "calcn_l", gravity, translation_calcn_l_imu);
//    //SimTK::Vec3 linAcc_calcn_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "calcn_r", gravity, translation_calcn_r_imu);
//
//    // Residual forces in OpenSim order
//    double res_os[ndofr];
//    /// OpenSim and Simbody have different state orders so we need to adjust
//    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
//    for (int i = 0; i < ndofr; ++i) res_os[i] =
//            (residualMobilityForces[indicesSimbodyInOS[i]]);
//    // Extract results
//    int nc = 3;
//    /// Residual forces
//    /// We do want to extract the pro_sup torques (last two -> till NU)
//    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];
//    //for (int i = 0; i < ndofr; ++i) res[0][i] = (residualMobilityForces[indicesSimbodyInOS[i]]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 0*nc] = (angVel_pelvis_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 1*nc] = (angVel_torso_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 2*nc] = (angVel_femur_l_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 3*nc] = (angVel_femur_r_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 4*nc] = (linAcc_pelvis_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 5*nc] = (linAcc_torso_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 6*nc] = (linAcc_femur_l_imu_inG[i]);
//    //for (int i = 0; i < nc; ++i) res[0][i + ndofr + 7*nc] = (linAcc_femur_r_imu_inG[i]);
//
//    return 0;
//
//}
//
///* In main(), the Recorder is used to save the expression graph of function F.
//This expression graph is saved as a MATLAB function named foo.m. From this
//function, a c-code can be generated via CasADi and then compiled as a dll. This
//dll is then imported in MATLAB as an external function. With this workflow,
//CasADi can use algorithmic differentiation to differentiate the function F.
//*/
//int main() {
//    double x[NX];
//	double u[NU];
//	double tau[NR];
//
//	for (int i = 0; i < NX; ++i) x[i] = 1;
//	for (int i = 0; i < NU; ++i) u[i] = 1;
//
//	const double* Recorder_arg[n_in] = { x,u };
//	double* Recorder_res[n_out] = { tau };
//
//	F_generic(Recorder_arg, Recorder_res);
//
//    for (int i = 0; i < NR; ++i) std::cout << Recorder_res[0][i] << std::endl;
//
//	return 0;
//}





///*  This code describes the OpenSim model and the skeleton dynamics
//    Author: Gil Serrancoli
//    Contributor: Joris Gillis, Antoine Falisse, Chris Dembia
//*/
//#include <OpenSim/Simulation/Model/Model.h>
//#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/Body.h>
//
//#include <iostream>
//#include <iterator>
//#include <random>
//#include <cassert>
//#include <algorithm>
//#include <vector>
//#include <fstream>
//
//using namespace SimTK;
//using namespace OpenSim;
//
///*  The function F describes the OpenSim model and, implicitly, the skeleton
//    dynamics. F takes as inputs joint positions and velocities (states x),
//    joint accelerations (controls u), a platform perturbation value (p), and
//    returns the joint torques. F is templatized using type T. F(x,u,p)->(r).
//*/
//
//// Inputs/outputs of function F
///// number of vectors in inputs/outputs of function F
//constexpr int n_in = 3;
//constexpr int n_out = 1;
///// number of elements in input/output vectors of function F
//constexpr int ndof = 2;     // # degrees of freedom
//constexpr int NX = ndof*2;  // # states
//constexpr int NU = ndof;    // # controls
//constexpr int NR = ndof;    // # residual torques
//
//// Function F
//int F_generic(const double** arg, double** res) {
//
//	// OpenSim model: create components
//    /// Model
//	Model* model;
//    /// Bodies
//	OpenSim::Body* lowerbody;
//	OpenSim::Body* upperbody;
//    /// Joints
//	OpenSim::PinJoint* ankle;
//	OpenSim::PinJoint* hip;
//
//    // OpenSim model: initialize components
//    /// Model
//	model = new OpenSim::Model();
//    /// Body specifications
//	lowerbody = new OpenSim::Body("lowerbody", 26, Vec3(0, 0.55, 0),
//            Inertia(1.4));
//	upperbody = new OpenSim::Body("upperbody", 46, Vec3(0, 0.365, 0),
//            Inertia(2.9));
//    /// Joint specifications
//	ankle = new PinJoint("ankle", model->getGround(), Vec3(0), Vec3(0),
//            *lowerbody, Vec3(0), Vec3(0));
//	hip = new PinJoint("hip", *lowerbody, Vec3(0, 0.853, 0), Vec3(0),
//            *upperbody, Vec3(0, 0, 0), Vec3(0));
//    /// Add bodies and joints to model
//	model->addBody(lowerbody);      model->addBody(upperbody);
//	model->addJoint(ankle);         model->addJoint(hip);
//
//	// Initialize system and state
//	State* state;
//	state = new State(model->initSystem());
//
//	// Read inputs
//	std::vector<double> x(arg[0], arg[0] + NX);
//	std::vector<double> u(arg[1], arg[1] + NU);
//
//    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states
//	double ua[NU]; /// joint accelerations (Qdotdots) - controls
//
//	// Assign inputs to model variables
//    /// States
//	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
//    /// Controls
//	for (int i = 0; i < NU; ++i) ua[i] = u[i];
//
//    // Set state variables and realize
//	model->setStateVariableValues(*state, QsUs);
//	model->realizeVelocity(*state);
//
//	// Compute residual forces
//    /// appliedMobilityForces (# mobilities)
//	Vector appliedMobilityForces(ndof);
//	appliedMobilityForces.setToZero();
//    /// appliedBodyForces (# bodies + ground)
//	Vector_<SpatialVec> appliedBodyForces;
//	int nbodies = model->getBodySet().getSize() + 1;
//	appliedBodyForces.resize(nbodies);
//	appliedBodyForces.setToZero();
//	/// Set gravity
//	Vec3 gravity(0);
//	gravity[1] = -9.81;
//	/// Add to model
//	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
//		model->getMatterSubsystem().addInStationForce(*state,
//            model->getBodySet().get(i).getMobilizedBodyIndex(),
//            model->getBodySet().get(i).getMassCenter(),
//            gravity*model->getBodySet().get(i).getMass(),
//            appliedBodyForces);
//	}
//    /// knownUdot
//	Vector knownUdot(ndof);
//	knownUdot.setToZero();
//    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
//    /// Calculate residual forces
//	Vector residualMobilityForces(ndof);
//	residualMobilityForces.setToZero();
//	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
//        appliedMobilityForces, appliedBodyForces, knownUdot,
//        residualMobilityForces);
//
//    Vector_<SpatialVec> A_GB(nbodies);
//    model->getMatterSubsystem().calcBodyAccelerationFromUDot(*state, knownUdot, A_GB);
//    std::cout << A_GB << std::endl;
//
//	// Extract results
//    /// Residual forces
//	for (int i = 0; i < NR; ++i) res[0][i] = (residualMobilityForces[i]);
//
//	return 0;
//
//}
//
///* In main(), the Recorder is used to save the expression graph of function F.
//This expression graph is saved as a MATLAB function named foo.m. From this
//function, a c-code can be generated via CasADi and then compiled as a dll. This
//dll is then imported in MATLAB as an external function. With this workflow,
//CasADi can use algorithmic differentiation to differentiate the function F.
//*/
//int main() {
//
//	double x[NX];
//	double u[NU];
//	double tau[NR];
//
//	for (int i = 0; i < NX; ++i) x[i] = 1;
//	for (int i = 0; i < NU; ++i) u[i] = 1;
//
//	const double* Recorder_arg[n_in] = { x,u };
//	double* Recorder_res[n_out] = { tau };
//
//	F_generic(Recorder_arg, Recorder_res);
//
//	return 0;
//
//}



//#include <OpenSim/Simulation/Model/Model.h>
//#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
//#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
//#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
//#include <OpenSim/Common/LinearFunction.h>
//#include <OpenSim/Common/Constant.h>
//#include <OpenSim/Common/SimmSpline.h>
//#include <OpenSim/Common/MultiplierFunction.h>
//#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
//#include <OpenSim/Simulation/Model/MovingPathPoint.h>
//#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
//#include <OpenSim/Simulation/InverseDynamicsSolver.h>
//
//#include <iostream>
//#include <iterator>
//#include <random>
//#include <cassert>
//#include <algorithm>
//#include <vector>
//#include <fstream>
//#include <chrono>
//
//using namespace SimTK;
//using namespace OpenSim;
//
//// Declare inputs/outputs of function F
///// number of vectors in inputs/outputs of function F
//constexpr int n_in = 2;
//constexpr int n_out = 1;
///// number of elements in input/output vectors of function F
//constexpr int ndof = 31;        // # degrees of freedom (excluding locked)
//constexpr int ndofr = ndof+2;   // # degrees of freedom (including locked)
//constexpr int NX = 2 * ndof;    // states
//constexpr int NU = ndof;        // controls
//constexpr int NR = ndof+8*3;        // residual forces + segment origin + GRFs + GRMs + PoA
//
//// OpenSim and Simbody use different indices for the states/controls when the
//// kinematic chain has joints up and down the origin (e.g., lumbar joint/arms
//// and legs with pelvis as origin).
//// The two following functions allow getting the indices from one reference
//// system to the other. These functions are inspired from
//// createSystemYIndexMap() in Moco.
//// getIndicesOSInSimbody() returns the indices of the OpenSim Qs in the Simbody
//// reference system. Note that we only care about the order here so we divide
//// by 2 because the states include both Qs and Qdots.
//SimTK::Array_<int> getIndicesOSInSimbody(const Model& model) {
//    auto s = model.getWorkingState();
//    const auto svNames = model.getStateVariableNames();
//    SimTK::Array_<int> idxOSInSimbody(s.getNQ());
//    s.updQ() = 0;
//    for (int iy = 0; iy < s.getNQ(); ++iy) {
//        s.updQ()[iy] = SimTK::NaN;
//        const auto svValues = model.getStateVariableValues(s);
//        for (int isv = 0; isv < svNames.size(); ++isv) {
//            if (SimTK::isNaN(svValues[isv])) {
//                s.updQ()[iy] = 0;
//                idxOSInSimbody[iy] = isv/2;
//                break;
//            }
//        }
//    }
//    return idxOSInSimbody;
//}
//// getIndicesSimbodyInOS() returns the indices of the Simbody Qs in the OpenSim
//// reference system.
//SimTK::Array_<int> getIndicesSimbodyInOS(const Model& model) {
//    auto idxOSInSimbody = getIndicesOSInSimbody(model);
//    auto s = model.getWorkingState();
//    SimTK::Array_<int> idxSimbodyInOS(s.getNQ());
//	for (int iy = 0; iy < s.getNQ(); ++iy) {
//		for (int iyy = 0; iyy < s.getNQ(); ++iyy) {
//			if (idxOSInSimbody[iyy] == iy) {
//				idxSimbodyInOS[iy] = iyy;
//				break;
//			}
//		}
//	}
//    return idxSimbodyInOS;
//}
//
////// This function returns the linear acceleration of the imu frame wrt the
////// ground frame expressed in the ground frame.
////// frameName is the name of the body segment to which the imu is attached.
////// imuTranslation_B is the translation offset of the imu frame'origin from the
////// parent (body segment) frame's origin, expressed in the parent frame.
//const SimTK::Vec3 getLinearAccelerationIMUInGround(Model& model, const State& s, const Vector_<SpatialVec>& A_GB , const std::string& frameName, const Vec3& gravity_G, const Vec3& imuTranslation_B)
//{
//        SimTK::Rotation R_GB = model.getBodySet().get(frameName).getMobilizedBody().getBodyTransform(s).R();
//        /// Body linear acceleration in ground
//        SimTK::Vec3 linAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][1];
//         //std::cout << linAcc_G << std::endl;
//        /// Body angular acceleration in ground
//        SimTK::Vec3 angAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][0];
//        //std::cout << angAcc_G << std::endl;
//        /// Body angular velocity in ground
//        SimTK::Vec3 angVel_G = model.getBodySet().get(frameName).getAngularVelocityInGround(s);
//        /// Body angular velocity in body
//        SimTK::Vec3 angVel_B = ~R_GB*angVel_G;
//        /// Body angular acceleration in body
//        SimTK::Vec3 angAcc_B = ~R_GB*angAcc_G;
//        /// Sensor linear acceleration
//        /// See van den Bogert et al. (1995), equation (1), p949.
//        SimTK::Vec3 linAcc_imu_B = ~R_GB * (linAcc_G - gravity_G) + SimTK::cross(angAcc_B, imuTranslation_B) + SimTK::cross(angVel_B, SimTK::cross(angVel_B, imuTranslation_B));
//        SimTK::Vec3 linAcc_imu_G = R_GB * linAcc_imu_B;
//
//        return linAcc_imu_G;
//}
//
//int F_generic(const double** arg, double** res) {
//	Model* model = new Model("Giacomo_IMU_0_noMuscles_ordered_imu.osim");
//
//	/// Initialize  system and state.
//    SimTK::State* state;
//	state = new State(model->initSystem());
//
//	// Read inputs
//    std::vector<double> x(arg[0], arg[0] + NX);
//    std::vector<double> u(arg[1], arg[1] + NU);
//
//    // States and controls
//    double ua[NU+2]; /// joint accelerations (Qdotdots) - controls
//    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states
//
//    // Assign inputs to model variables
//    /// States
//    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
//    /// pro_sup dofs are locked so Qs and Qdots are hard coded (0)
//    std::string str1 ("/jointset/radioulnar_l/pro_sup_l/value");
//    std::string str2 ("/jointset/radioulnar_r/pro_sup_r/value");
//        auto svNames = model->getStateVariableNames();
//    int idx_radioulnar_l;
//    int idx_radioulnar_r;
//    for (int isv = 0; isv < svNames.size(); ++isv) {
//        if (str1.compare(svNames[isv]) == 0)
//            idx_radioulnar_l = isv;
//        if (str2.compare(svNames[isv]) == 0)
//            idx_radioulnar_r = isv;
//    }
//    QsUs[idx_radioulnar_r] = SimTK::Pi/2;
//    QsUs[idx_radioulnar_r+1] = 0;
//    QsUs[idx_radioulnar_l] = SimTK::Pi/2;
//    QsUs[idx_radioulnar_l+1] = 0;
//    /// Controls
//    double ut[NU+2];
//    for (int i = 0; i < NU; ++i) ut[i] = u[i];
//    /// pro_sup dofs are locked so Qdotdots are hard coded (0)
//    /// Need to have a temporary vector to add 0s to the vector before
//    /// adjusting for the index difference between OpenSim and Simbody.
//    ut[idx_radioulnar_l/2] = 0;
//    ut[idx_radioulnar_r/2] = 0;
//    /// OpenSim and Simbody have different state orders so we need to adjust
//    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
//    for (int i = 0; i < ndofr; ++i) ua[i] = ut[indicesOSInSimbody[i]];
//
//    // Set state variables and realize
//    model->setStateVariableValues(*state, QsUs);
//    model->realizeVelocity(*state);
//
//    //InverseDynamicsSolver* IDsolver;
//    //IDsolver = new InverseDynamicsSolver(*model);
//    //Vector residualMobilityForces = IDsolver->solve(*state);
//
//    // Compute residual forces
//    /// appliedMobilityForces (# mobilities)
//    Vector appliedMobilityForces(ndofr);
//    appliedMobilityForces.setToZero();
//    /// appliedBodyForces (# bodies + ground)
//    Vector_<SpatialVec> appliedBodyForces;
//    int nbodies = model->getBodySet().getSize() + 1;
//    appliedBodyForces.resize(nbodies);
//    appliedBodyForces.setToZero();
//    /// Set gravity
//    const Vec3 gravity(0, -9.80665, 0);
//    /// Add weights to appliedBodyForces
//    for (int i = 0; i < model->getBodySet().getSize(); ++i) {
//        model->getMatterSubsystem().addInStationForce(*state,
//            model->getBodySet().get(i).getMobilizedBodyIndex(),
//            model->getBodySet().get(i).getMassCenter(),
//            model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
//    }
//
//    /// Add contact forces to appliedBodyForces
//    /// Right foot
//    Array<double> Force_values_s1_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s1_r").getRecordValues(*state);
//    Array<double> Force_values_s2_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s2_r").getRecordValues(*state);
//    Array<double> Force_values_s3_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s3_r").getRecordValues(*state);
//    Array<double> Force_values_s4_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s4_r").getRecordValues(*state);
//    Array<double> Force_values_s5_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s5_r").getRecordValues(*state);
//    Array<double> Force_values_s6_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s6_r").getRecordValues(*state);
//    SpatialVec GRF_s1_r;
//    GRF_s1_r[0] = Vec3(Force_values_s1_r[3], Force_values_s1_r[4], Force_values_s1_r[5]);
//    GRF_s1_r[1] = Vec3(Force_values_s1_r[0], Force_values_s1_r[1], Force_values_s1_r[2]);
//    SpatialVec GRF_s2_r;
//    GRF_s2_r[0] = Vec3(Force_values_s2_r[3], Force_values_s2_r[4], Force_values_s2_r[5]);
//    GRF_s2_r[1] = Vec3(Force_values_s2_r[0], Force_values_s2_r[1], Force_values_s2_r[2]);
//    SpatialVec GRF_s3_r;
//    GRF_s3_r[0] = Vec3(Force_values_s3_r[3], Force_values_s3_r[4], Force_values_s3_r[5]);
//    GRF_s3_r[1] = Vec3(Force_values_s3_r[0], Force_values_s3_r[1], Force_values_s3_r[2]);
//    SpatialVec GRF_s4_r;
//    GRF_s4_r[0] = Vec3(Force_values_s4_r[3], Force_values_s4_r[4], Force_values_s4_r[5]);
//    GRF_s4_r[1] = Vec3(Force_values_s4_r[0], Force_values_s4_r[1], Force_values_s4_r[2]);
//    SpatialVec GRF_s5_r;
//    GRF_s5_r[0] = Vec3(Force_values_s5_r[3], Force_values_s5_r[4], Force_values_s5_r[5]);
//    GRF_s5_r[1] = Vec3(Force_values_s5_r[0], Force_values_s5_r[1], Force_values_s5_r[2]);
//    SpatialVec GRF_s6_r;
//    GRF_s6_r[0] = Vec3(Force_values_s6_r[3], Force_values_s6_r[4], Force_values_s6_r[5]);
//    GRF_s6_r[1] = Vec3(Force_values_s6_r[0], Force_values_s6_r[1], Force_values_s6_r[2]);
//    int nCalc_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
//    int nToe_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
//    appliedBodyForces[nCalc_r] = appliedBodyForces[nCalc_r] + GRF_s1_r + GRF_s2_r + GRF_s3_r + GRF_s4_r;
//    appliedBodyForces[nToe_r] = appliedBodyForces[nToe_r] + GRF_s5_r + GRF_s6_r;
//
//    /// Left foot
//    Array<double> Force_values_s1_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s1_l").getRecordValues(*state);
//    Array<double> Force_values_s2_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s2_l").getRecordValues(*state);
//    Array<double> Force_values_s3_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s3_l").getRecordValues(*state);
//    Array<double> Force_values_s4_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s4_l").getRecordValues(*state);
//    Array<double> Force_values_s5_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s5_l").getRecordValues(*state);
//    Array<double> Force_values_s6_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_s6_l").getRecordValues(*state);
//    SpatialVec GRF_s1_l;
//    GRF_s1_l[0] = Vec3(Force_values_s1_l[3], Force_values_s1_l[4], Force_values_s1_l[5]);
//    GRF_s1_l[1] = Vec3(Force_values_s1_l[0], Force_values_s1_l[1], Force_values_s1_l[2]);
//    SpatialVec GRF_s2_l;
//    GRF_s2_l[0] = Vec3(Force_values_s2_l[3], Force_values_s2_l[4], Force_values_s2_l[5]);
//    GRF_s2_l[1] = Vec3(Force_values_s2_l[0], Force_values_s2_l[1], Force_values_s2_l[2]);
//    SpatialVec GRF_s3_l;
//    GRF_s3_l[0] = Vec3(Force_values_s3_l[3], Force_values_s3_l[4], Force_values_s3_l[5]);
//    GRF_s3_l[1] = Vec3(Force_values_s3_l[0], Force_values_s3_l[1], Force_values_s3_l[2]);
//    SpatialVec GRF_s4_l;
//    GRF_s4_l[0] = Vec3(Force_values_s4_l[3], Force_values_s4_l[4], Force_values_s4_l[5]);
//    GRF_s4_l[1] = Vec3(Force_values_s4_l[0], Force_values_s4_l[1], Force_values_s4_l[2]);
//    SpatialVec GRF_s5_l;
//    GRF_s5_l[0] = Vec3(Force_values_s5_l[3], Force_values_s5_l[4], Force_values_s5_l[5]);
//    GRF_s5_l[1] = Vec3(Force_values_s5_l[0], Force_values_s5_l[1], Force_values_s5_l[2]);
//    SpatialVec GRF_s6_l;
//    GRF_s6_l[0] = Vec3(Force_values_s6_l[3], Force_values_s6_l[4], Force_values_s6_l[5]);
//    GRF_s6_l[1] = Vec3(Force_values_s6_l[0], Force_values_s6_l[1], Force_values_s6_l[2]);
//    int nCalc_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
//    int nToe_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
//    appliedBodyForces[nCalc_l] = appliedBodyForces[nCalc_l] + GRF_s1_l + GRF_s2_l + GRF_s3_l + GRF_s4_l;
//    appliedBodyForces[nToe_l] = appliedBodyForces[nToe_l] + GRF_s5_l + GRF_s6_l;
//
//    /// knownUdot
//    Vector knownUdot(ndofr);
//    knownUdot.setToZero();
//    for (int i = 0; i < ndofr; ++i) knownUdot[i] = ua[i];
//    /// Calculate residual forces
//    Vector residualMobilityForces(ndofr);
//    residualMobilityForces.setToZero();
//    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
//        *state, appliedMobilityForces, appliedBodyForces, knownUdot,
//        residualMobilityForces);
//
//    // Virtual sensors
//    // Approach 1: extract from model
//    //SimTK::Vec3 translation_pelvis_imu2(model->getBodySet().get("pelvis").getComponent<OpenSim::PhysicalOffsetFrame>("pelvis_imu").get_translation());
//    //OpenSim::PhysicalOffsetFrame pelvis_imu2(model->getBodySet().get("pelvis").getComponent<OpenSim::PhysicalOffsetFrame>("pelvis_imu"));
//    //std::cout << pelvis_imu2.getOffsetTransform() << std::endl; // same transform
//    // Approach 2: build programmatically
//    //const SimTK::Vec3 translation_pelvis_imu(-0.17825090079013006, 0.06148338297319611, -0.0039742631657566363);
//    //SimTK::Vec3 orientation_pelvis_imu(-1.334791708492248, -1.1869361939759826, 0.23826934527288943);
//    //SimTK::Rotation R_imuB_pelvis;
//    //R_imuB_pelvis.setRotationToBodyFixedXYZ(
//    //    Vec3(std::cos(orientation_pelvis_imu[0]),std::cos(orientation_pelvis_imu[1]),std::cos(orientation_pelvis_imu[2])),
//    //    Vec3(std::sin(orientation_pelvis_imu[0]),std::sin(orientation_pelvis_imu[1]),std::sin(orientation_pelvis_imu[2])));
//    //SimTK::Transform T_pelvis_imu(R_imuB_pelvis, translation_pelvis_imu);
//    //std::cout << T_pelvis_imu << std::endl; // same transform
//    //OpenSim::PhysicalOffsetFrame* pelvis_imu;
//    //pelvis_imu = new PhysicalOffsetFrame("pelvis_imu", model->getBodySet().get("pelvis"), T_pelvis_imu);
//
//    // Now replaced by getLinearAccelerationIMUInGround()
//    //// Compute the linear acceleration of the sensor wrt ground in ground
//    ///// Rotation
//    //SimTK::Rotation R_GB_pelvis = (model->getBodySet().get("pelvis").getMobilizedBody().getBodyTransform(*state).R());
//    ///// Body linear acceleration in ground
//    //SimTK::Vec3 linAcc_pelvis_inG = model->getBodySet().get("pelvis").getLinearAccelerationInGround(*state);
//    ///// Body angular velocity in ground
//    //SimTK::Vec3 angVel_pelvis_inG = model->getBodySet().get("pelvis").getAngularVelocityInGround(*state);
//    ///// Body angular velocity in body
//    //SimTK::Vec3 angVel_pelvis_inB = ~R_GB_pelvis*angVel_pelvis_inG;
//    ///// Body angular acceleration in ground
//    //SimTK::Vec3 angAcc_pelvis_inG = model->getBodySet().get("pelvis").getAngularAccelerationInGround(*state);
//    ///// Body angular acceleration in body
//    //SimTK::Vec3 angAcc_pelvis_inB = ~R_GB_pelvis*angAcc_pelvis_inG;
//    ///// Sensor linear acceleration
//    //SimTK::Vec3 linAcc_pelvis_imu_inB = ~R_GB_pelvis * (linAcc_pelvis_inG - gravity) + SimTK::cross(angAcc_pelvis_inB, translation_pelvis_imu) + SimTK::cross(angVel_pelvis_inB, SimTK::cross(angVel_pelvis_inB, translation_pelvis_imu));
//    //SimTK::Vec3 linAcc_pelvis_imu_inG = R_GB_pelvis * linAcc_pelvis_imu_inB;
//    //std::cout << linAcc_pelvis_imu_inG << std::endl;
//    //// The angular velocity f
//    //std::cout << angAcc_pelvis_inG << std::endl;
//
//
//    const SimTK::Vec3 translation_pelvis_imu(-0.17825090079013006, 0.06148338297319611, -0.0039742631657566363);
//    const SimTK::Vec3 translation_torso_imu(0.11141304895632698, 0.32812980850924067, -0.012040552365984683);
//
//    const SimTK::Vec3 translation_femur_r_imu(0.043933841399841023, -0.14958344693305592, 0.099866089437989247);
//    const SimTK::Vec3 translation_femur_l_imu(0.057420458331616797, -0.12466095809783695, -0.1025425014551595);
//
//    const SimTK::Vec3 translation_tibia_r_imu(0.047080175165178595, -0.11466976609963364, 0.0052928998697051033);
//    const SimTK::Vec3 translation_tibia_l_imu(0.053945123225004443, -0.12392203671935126, -0.006707066365313652);
//
//    const SimTK::Vec3 translation_calcn_r_imu(0.14808446197590852, 0.040621301301533644, 0.035088429802814902);
//    const SimTK::Vec3 translation_calcn_l_imu(0.15253209992514105, 0.045078665693204636, -0.044293377012806667);
//
//    SimTK::Vec3 angVel_pelvis_imu_inG   = model->getBodySet().get("pelvis").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_torso_imu_inG    = model->getBodySet().get("torso").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_femur_l_imu_inG  = model->getBodySet().get("femur_l").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_femur_r_imu_inG  = model->getBodySet().get("femur_r").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_tibia_l_imu_inG  = model->getBodySet().get("tibia_l").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_tibia_r_imu_inG  = model->getBodySet().get("tibia_r").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_calcn_l_imu_inG  = model->getBodySet().get("calcn_l").getAngularVelocityInGround(*state);
//    SimTK::Vec3 angVel_calcn_r_imu_inG  = model->getBodySet().get("calcn_r").getAngularVelocityInGround(*state);
//
//    Vector_<SpatialVec> A_GB(nbodies);
//    model->getMatterSubsystem().calcBodyAccelerationFromUDot(*state, knownUdot, A_GB);
//
//    SimTK::Vec3 linAcc_pelvis_imu_inG   = getLinearAccelerationIMUInGround(*model, *state, A_GB, "pelvis",  gravity, translation_pelvis_imu);
//    SimTK::Vec3 linAcc_torso_imu_inG    = getLinearAccelerationIMUInGround(*model, *state, A_GB, "torso",   gravity, translation_torso_imu);
//    SimTK::Vec3 linAcc_femur_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "femur_l", gravity, translation_femur_l_imu);
//    SimTK::Vec3 linAcc_femur_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "femur_r", gravity, translation_femur_r_imu);
//    SimTK::Vec3 linAcc_tibia_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "tibia_l", gravity, translation_tibia_l_imu);
//    SimTK::Vec3 linAcc_tibia_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "tibia_r", gravity, translation_tibia_r_imu);
//    SimTK::Vec3 linAcc_calcn_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "calcn_l", gravity, translation_calcn_l_imu);
//    SimTK::Vec3 linAcc_calcn_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "calcn_r", gravity, translation_calcn_r_imu);
//
//    //std::cout << model->getGround().getMobilizedBodyIndex() << std::endl;
//    //std::cout << model->getBodySet().get("talus_r").getMobilizedBodyIndex() << std::endl;
//
//    //std::cout << A_GB[model->getGround().getMobilizedBodyIndex()] << std::endl;
//    //std::cout << A_GB[model->getBodySet().get("pelvis").getMobilizedBodyIndex()] << std::endl;
//    //std::cout << model->getBodySet().get("talus_r").getAngularVelocityInGround(*state) << std::endl;
//
//    // Residual forces in OpenSim order
//    double res_os[ndofr];
//    /// OpenSim and Simbody have different state orders so we need to adjust
//    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
//    for (int i = 0; i < ndofr; ++i) res_os[i] =
//            (residualMobilityForces[indicesSimbodyInOS[i]]);
//    // Extract results
//    int nc = 3;
//    /// Residual forces
//    /// We do want to extract the pro_sup torques (last two -> till NU)
//    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];
//    //for (int i = 0; i < ndofr; ++i) res[0][i] = (residualMobilityForces[indicesSimbodyInOS[i]]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 0*nc] = (angVel_pelvis_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 1*nc] = (angVel_torso_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 2*nc] = (angVel_femur_l_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 3*nc] = (angVel_femur_r_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 4*nc] = (linAcc_pelvis_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 5*nc] = (linAcc_torso_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 6*nc] = (linAcc_femur_l_imu_inG[i]);
//    for (int i = 0; i < nc; ++i) res[0][i + NU + 7*nc] = (linAcc_femur_r_imu_inG[i]);
//
//	return 0;
//}
//
//int main() {
//    double x[NX];
//	double u[NU];
//	double tau[NR];
//
//	for (int i = 0; i < NX; ++i) x[i] = -1;
//	for (int i = 0; i < NU; ++i) u[i] = -1;
//
//	const double* Recorder_arg[n_in] = { x,u };
//	double* Recorder_res[n_out] = { tau };
//
//	F_generic(Recorder_arg, Recorder_res);
//
//    for (int i = 0; i < NR; ++i) std::cout << Recorder_res[0][i]  << std::endl;
//        //std::cout << Recorder_res[0][i] << std::endl;
//
//
//	return 0;
//}


/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Common/MultiplierFunction.h>
#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>

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
constexpr int NR = ndof + 8*3;    // # residual torques + # joint origins

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

//// This function returns the linear acceleration of the imu frame wrt the
//// ground frame expressed in the ground frame.
//// frameName is the name of the body segment to which the imu is attached.
//// imuTranslation_B is the translation offset of the imu frame'origin from the
//// parent (body segment) frame's origin, expressed in the parent frame.
const SimTK::Vec3 getLinearAccelerationIMUInGround(Model& model, const State& s, const Vector_<SpatialVec>& A_GB , const std::string& frameName, const Vec3& gravity_G, const Vec3& imuTranslation_B)
{
        SimTK::Rotation R_GB = model.getBodySet().get(frameName).getMobilizedBody().getBodyTransform(s).R();
        /// Body linear acceleration in ground
        SimTK::Vec3 linAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][1];
        /// Body angular acceleration in ground
        SimTK::Vec3 angAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][0];
        /// Body angular velocity in ground
        SimTK::Vec3 angVel_G = model.getBodySet().get(frameName).getAngularVelocityInGround(s);
        /// Body angular velocity in body
        SimTK::Vec3 angVel_B = ~R_GB*angVel_G;
        /// Body angular acceleration in body
        SimTK::Vec3 angAcc_B = ~R_GB*angAcc_G;
        /// Sensor linear acceleration
        /// See van den Bogert et al. (1995), equation (1), p949.
        SimTK::Vec3 linAcc_imu_B = ~R_GB * (linAcc_G - gravity_G) + SimTK::cross(angAcc_B, imuTranslation_B) + SimTK::cross(angVel_B, SimTK::cross(angVel_B, imuTranslation_B));
        SimTK::Vec3 linAcc_imu_G = R_GB * linAcc_imu_B;

        return linAcc_imu_G;
}

// Function F
int F_generic(const double** arg, double** res) {

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
    OpenSim::PinJoint* subtalar_r;
    OpenSim::PinJoint* subtalar_l;
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

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    //model->setUseVisualizer(true);
    /// Body specifications
    pelvis = new OpenSim::Body("pelvis", 11.751210011095651, Vec3(-0.069729482228687481, 0, 0), Inertia(0.099778065737821386, 0.08453958682650041, 0.056197957258948036, 0, 0, 0));
    femur_l = new OpenSim::Body("femur_l", 9.2810312301269491, Vec3(0, -0.17281329846170712, 0), Inertia(0.13806543520006995, 0.036191910198076598, 0.14559252763442779, 0, 0, 0));
    femur_r = new OpenSim::Body("femur_r", 9.2810312301269491, Vec3(0, -0.17281329846170712, 0), Inertia(0.13806543520006995, 0.036191910198076598, 0.14559252763442779, 0, 0, 0));
    tibia_l = new OpenSim::Body("tibia_l", 3.6993810916309018, Vec3(0, -0.20693800230674214, 0), Inertia(0.061783188539925718, 0.0062518702689210552, 0.062641288380758026, 0, 0, 0));
    tibia_r = new OpenSim::Body("tibia_r", 3.6993810916309018, Vec3(0, -0.20693800230674214, 0), Inertia(0.061783188539925718, 0.0062518702689210552, 0.062641288380758026, 0, 0, 0));
    talus_l = new OpenSim::Body("talus_l", 0.099781013934751236, Vec3(0, 0, 0), Inertia(0.0011448780870697802, 0.0011448780870697802, 0.0011448780870697802, 0, 0, 0));
    talus_r = new OpenSim::Body("talus_r", 0.099781013934751236, Vec3(0, 0, 0), Inertia(0.0011448780870697802, 0.0011448780870697802, 0.0011448780870697802, 0, 0, 0));
    calcn_l = new OpenSim::Body("calcn_l", 1.2472626741843904, Vec3(0.1071163252191219, 0.032134897565736564, 0), Inertia(0.0016028293218976922, 0.0044650245395721419, 0.0046940001569860989, 0, 0, 0));
    calcn_r = new OpenSim::Body("calcn_r", 1.2472626741843904, Vec3(0.1071163252191219, 0.032134897565736564, 0), Inertia(0.0016028293218976922, 0.0044650245395721419, 0.0046940001569860989, 0, 0, 0));
    toes_l = new OpenSim::Body("toes_l", 0.21612567618267114, Vec3(0.037062248525816174, 0.006426979513147313, 0.018745356913346334), Inertia(0.000114487808706978, 0.000228975617413956, 0.000114487808706978, 0, 0, 0));
    toes_r = new OpenSim::Body("toes_r", 0.21612567618267114, Vec3(0.037062248525816174, 0.006426979513147313, -0.018745356913346334), Inertia(0.000114487808706978, 0.000228975617413956, 0.000114487808706978, 0, 0, 0));
    torso = new OpenSim::Body("torso", 26.767853484219973, Vec3(-0.029119398106347269, 0.31060691313437089, 0), Inertia(1.3861651249269105, 0.7102392349150769, 1.345647175191848, 0, 0, 0));
    humerus_l = new OpenSim::Body("humerus_l", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
    humerus_r = new OpenSim::Body("humerus_r", 2.028049108223819, Vec3(0, -0.18049429839653774, 0), Inertia(0.014350103502831833, 0.0049503412468751024, 0.016107528701613263, 0, 0, 0));
    ulna_l = new OpenSim::Body("ulna_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    ulna_r = new OpenSim::Body("ulna_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    radius_l = new OpenSim::Body("radius_l", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    radius_r = new OpenSim::Body("radius_r", 0.6061696596536138, Vec3(0, -0.1408058196310753, 0), Inertia(0.0040338514094605477, 0.00084163408880709608, 0.0043756801413223296, 0, 0, 0));
    hand_l = new OpenSim::Body("hand_l", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
    hand_r = new OpenSim::Body("hand_r", 0.45649813875148687, Vec3(0, -0.079553389651757511, 0), Inertia(0.001214785772194061, 0.00074494149931631317, 0.0018249023932063251, 0, 0, 0));
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
    st_knee_l[0].setAxis(Vec3(0,0,1));
    st_knee_l[1].setFunction(new Constant(0));
    st_knee_l[1].setAxis(Vec3(0,1,0));
    st_knee_l[2].setFunction(new Constant(0));
    st_knee_l[2].setAxis(Vec3(1,0,0));
    st_knee_l[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    double knee_X_r_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
    double knee_X_r_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
    OpenSim::SimmSpline* knee_X_r = new SimmSpline(12, knee_X_r_x, knee_X_r_y, "function_X");
    st_knee_l[3].setFunction(new MultiplierFunction(knee_X_r, 1.0165488144806301));
    st_knee_l[3].setAxis(Vec3(1,0,0));
    st_knee_l[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
    double knee_Y_r_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
    double knee_Y_r_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
    OpenSim::SimmSpline* knee_Y_r = new SimmSpline(7, knee_Y_r_x, knee_Y_r_y, "function_Y");
    st_knee_l[4].setFunction(new MultiplierFunction(knee_Y_r, 1.0165488144806301));
    st_knee_l[4].setAxis(Vec3(0,1,0));
    st_knee_l[5].setFunction(new Constant(0));
    st_knee_l[5].setAxis(Vec3(0,0,1));
    /// Knee_r transform
    SpatialTransform st_knee_r;
    st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    st_knee_r[0].setFunction(new LinearFunction());
    st_knee_r[0].setAxis(Vec3(0,0,1));
    st_knee_r[1].setFunction(new Constant(0));
    st_knee_r[1].setAxis(Vec3(0,1,0));
    st_knee_r[2].setFunction(new Constant(0));
    st_knee_r[2].setAxis(Vec3(1,0,0));
    st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    double knee_X_l_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
    double knee_X_l_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
    OpenSim::SimmSpline* knee_X_l = new SimmSpline(12, knee_X_l_x, knee_X_l_y, "function_X");
    st_knee_r[3].setFunction(new MultiplierFunction(knee_X_l, 1.0165488144806301));
    st_knee_r[3].setAxis(Vec3(1,0,0));
    st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
    double knee_Y_l_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
    double knee_Y_l_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
    OpenSim::SimmSpline* knee_Y_l = new SimmSpline(7, knee_Y_l_x, knee_Y_l_y, "function_Y");
    st_knee_r[4].setFunction(new MultiplierFunction(knee_Y_l, 1.0165488144806301));
    st_knee_r[4].setAxis(Vec3(0,1,0));
    st_knee_r[5].setFunction(new Constant(0));
    st_knee_r[5].setAxis(Vec3(0,0,1));
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
    hip_l = new CustomJoint("hip_l", *pelvis, Vec3(-0.069729482228687481, -0.065192627656523949, -0.08235377321209908), Vec3(0), *femur_l, Vec3(0), Vec3(0), st_hip_l);
    hip_r = new CustomJoint("hip_r", *pelvis, Vec3(-0.069729482228687481, -0.065192627656523949, 0.08235377321209908), Vec3(0), *femur_r, Vec3(0), Vec3(0), st_hip_r);
    knee_l = new CustomJoint("knee_l", *femur_l, Vec3(0, 0, 0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
    knee_r = new CustomJoint("knee_r", *femur_r, Vec3(0, 0, 0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
    ankle_l = new CustomJoint("ankle_l", *tibia_l, Vec3(0, -0.47661136042795454, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0), st_ankle_l);
    ankle_r = new CustomJoint("ankle_r", *tibia_r, Vec3(0, -0.47661136042795454, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0), st_ankle_r);
    subtalar_l = new PinJoint("subtalar_l", *talus_l, Vec3(-0.052240631809365744, -0.044935298429421636, -0.0084836129573544541), Vec3(1.7681899999999999, -0.906223, 1.8196000000000001),
                                *calcn_l, Vec3(0), Vec3(1.7681899999999999, -0.906223, 1.8196000000000001));
    subtalar_r = new PinJoint("subtalar_r", *talus_r, Vec3(-0.052240631809365744, -0.044935298429421636, 0.0084836129573544541), Vec3(-1.7681899999999999, 0.906223, 1.8196000000000001),
                                *calcn_r, Vec3(0), Vec3(-1.7681899999999999, 0.906223, 1.8196000000000001));
    mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.19152398949178992, -0.0021423265043824377, -0.0011568563123665165), Vec3(-3.1415899999999999, -0.61990100000000004, 0),
                                *toes_l, Vec3(0), Vec3(-3.1415899999999999, -0.61990100000000004, 0));
    mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.19152398949178992, -0.0021423265043824377, 0.0011568563123665165), Vec3(-3.1415899999999999, 0.61990100000000004, 0),
                                *toes_r, Vec3(0), Vec3(-3.1415899999999999, 0.61990100000000004, 0));
    back = new CustomJoint("back", *pelvis, Vec3(-0.099317664221058408, 0.080381227745941017, 0), Vec3(0), *torso, Vec3(0), Vec3(0), st_back);
    shoulder_l = new CustomJoint("shoulder_l", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, -0.16500992260263453), Vec3(0), *humerus_l, Vec3(0), Vec3(0), st_sho_l);
    shoulder_r = new CustomJoint("shoulder_r", *torso, Vec3(0.0030623900341841875, 0.36059521321693366, 0.16500992260263453), Vec3(0), *humerus_r, Vec3(0), Vec3(0), st_sho_r);
    elbow_l = new CustomJoint("elbow_l", *humerus_l, Vec3(0.014421812854093517, -0.31410344120358441, 0.010527791717515771), Vec3(0), *ulna_l, Vec3(0), Vec3(0), st_elb_l);
    elbow_r = new CustomJoint("elbow_r", *humerus_r, Vec3(0.014421812854093517, -0.31410344120358441, -0.010527791717515771), Vec3(0), *ulna_r, Vec3(0), Vec3(0), st_elb_r);
    radioulnar_l = new CustomJoint("radioulnar_l", *ulna_l, Vec3(-0.0078589566368657427, -0.015195696294888168, -0.030472003264362887), Vec3(0), *radius_l, Vec3(0), Vec3(0),st_radioulnar_l);
    radioulnar_r = new CustomJoint("radioulnar_r", *ulna_r, Vec3(-0.0078589566368657427, -0.015195696294888168, 0.030472003264362887), Vec3(0), *radius_r, Vec3(0), Vec3(0),st_radioulnar_r);
    radius_hand_l = new WeldJoint("radius_hand_l", *radius_l, Vec3(-0.010277276874462305, -0.27552611746618899, -0.015900163494535866), Vec3(0), *hand_l, Vec3(0), Vec3(0));
    radius_hand_r = new WeldJoint("radius_hand_r", *radius_r, Vec3(-0.010277276874462305, -0.27552611746618899, 0.015900163494535866), Vec3(0), *hand_r, Vec3(0), Vec3(0));
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
    double radiusSphere_s1 = 0.03232;
    double radiusSphere_s2 = 0.03232;
    double radiusSphere_s3 = 0.023374;
    double radiusSphere_s4 = 0.020508;
    double radiusSphere_s5 = 0.016244;
    double radiusSphere_s6 = 0.018414;
    double stiffness = 1000000;
    double dissipation = 2.0;
    double staticFriction = 0.8;
    double dynamicFriction = 0.8;
    double viscousFriction = 0.5;
    double transitionVelocity = 0.2;
    Vec3 normal = Vec3(0, 1, 0);
    double offset = 0;
    Vec3 locSphere_s1_r(-0.00042152, 0, -0.0049972);
    Vec3 locSphere_s2_r(0.06, 0, 0.020001);
    Vec3 locSphere_s3_r(0.165, -0.01, 0.021183);
    Vec3 locSphere_s4_r(0.18, -0.01, -0.01);
    Vec3 locSphere_s5_r(0.053154, -0.01, -0.0034173);
    Vec3 locSphere_s6_r(0.01, -0.01, -0.015);
    Vec3 locSphere_s1_l(locSphere_s1_r[0],locSphere_s1_r[1],-locSphere_s1_r[2]);
    Vec3 locSphere_s2_l(locSphere_s2_r[0],locSphere_s2_r[1],-locSphere_s2_r[2]);
    Vec3 locSphere_s3_l(locSphere_s3_r[0],locSphere_s3_r[1],-locSphere_s3_r[2]);
    Vec3 locSphere_s4_l(locSphere_s4_r[0],locSphere_s4_r[1],-locSphere_s4_r[2]);
    Vec3 locSphere_s5_l(locSphere_s5_r[0],locSphere_s5_r[1],-locSphere_s5_r[2]);
    Vec3 locSphere_s6_l(locSphere_s6_r[0],locSphere_s6_r[1],-locSphere_s6_r[2]);
    /// Left foot contact shere specifications

    Vec3 halfSpaceLocation(0);
	Vec3 halfSpaceOrientation(0, 0, -0.5 * SimTK::Pi);
    OpenSim::ContactHalfSpace* contactHalfSpace;
    contactHalfSpace = new OpenSim::ContactHalfSpace(halfSpaceLocation,halfSpaceOrientation,model->getGround(),"contactHalfSpace");
    model->addComponent(contactHalfSpace);

    OpenSim::ContactSphere* s1_l;
    s1_l = new OpenSim::ContactSphere(radiusSphere_s1,locSphere_s1_l,*calcn_l,"s1_l");
    model->addComponent(s1_l);
    OpenSim::ContactSphere* s2_l;
    s2_l = new OpenSim::ContactSphere(radiusSphere_s2,locSphere_s2_l,*calcn_l,"s2_l");
    model->addComponent(s2_l);
    OpenSim::ContactSphere* s3_l;
    s3_l = new OpenSim::ContactSphere(radiusSphere_s3,locSphere_s3_l,*calcn_l,"s3_l");
    model->addComponent(s3_l);
    OpenSim::ContactSphere* s4_l;
    s4_l = new OpenSim::ContactSphere(radiusSphere_s4,locSphere_s4_l,*calcn_l,"s4_l");
    model->addComponent(s4_l);
    OpenSim::ContactSphere* s5_l;
    s5_l = new OpenSim::ContactSphere(radiusSphere_s5,locSphere_s5_l,*toes_l,"s5_l");
    model->addComponent(s5_l);
    OpenSim::ContactSphere* s6_l;
    s6_l = new OpenSim::ContactSphere(radiusSphere_s6,locSphere_s6_l,*toes_l,"s6_l");
    model->addComponent(s6_l);
    OpenSim::ContactSphere* s1_r;
    s1_r = new OpenSim::ContactSphere(radiusSphere_s1,locSphere_s1_r,*calcn_r,"s1_r");
    model->addComponent(s1_r);
    OpenSim::ContactSphere* s2_r;
    s2_r = new OpenSim::ContactSphere(radiusSphere_s2,locSphere_s2_r,*calcn_r,"s2_r");
    model->addComponent(s2_r);
    OpenSim::ContactSphere* s3_r;
    s3_r = new OpenSim::ContactSphere(radiusSphere_s3,locSphere_s3_r,*calcn_r,"s3_r");
    model->addComponent(s3_r);
    OpenSim::ContactSphere* s4_r;
    s4_r = new OpenSim::ContactSphere(radiusSphere_s4,locSphere_s4_r,*calcn_r,"s4_r");
    model->addComponent(s4_r);
    OpenSim::ContactSphere* s5_r;
    s5_r = new OpenSim::ContactSphere(radiusSphere_s5,locSphere_s5_r,*toes_r,"s5_r");
    model->addComponent(s5_r);
    OpenSim::ContactSphere* s6_r;
    s6_r = new OpenSim::ContactSphere(radiusSphere_s6,locSphere_s6_r,*toes_r,"s6_r");
    model->addComponent(s6_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s1_r;
    contact_s1_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s1_r",*s1_r,*contactHalfSpace);
    contact_s1_r->set_stiffness(stiffness);
    contact_s1_r->set_dissipation(dissipation);
    contact_s1_r->set_static_friction(staticFriction);
    contact_s1_r->set_dynamic_friction(dynamicFriction);
    contact_s1_r->set_viscous_friction(viscousFriction);
    contact_s1_r->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s1_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s2_r;
    contact_s2_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s2_r",*s2_r,*contactHalfSpace);
    contact_s2_r->set_stiffness(stiffness);
    contact_s2_r->set_dissipation(dissipation);
    contact_s2_r->set_static_friction(staticFriction);
    contact_s2_r->set_dynamic_friction(dynamicFriction);
    contact_s2_r->set_viscous_friction(viscousFriction);
    contact_s2_r->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s2_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s3_r;
    contact_s3_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s3_r",*s3_r,*contactHalfSpace);
    contact_s3_r->set_stiffness(stiffness);
    contact_s3_r->set_dissipation(dissipation);
    contact_s3_r->set_static_friction(staticFriction);
    contact_s3_r->set_dynamic_friction(dynamicFriction);
    contact_s3_r->set_viscous_friction(viscousFriction);
    contact_s3_r->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s3_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s4_r;
    contact_s4_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s4_r",*s4_r,*contactHalfSpace);
    contact_s4_r->set_stiffness(stiffness);
    contact_s4_r->set_dissipation(dissipation);
    contact_s4_r->set_static_friction(staticFriction);
    contact_s4_r->set_dynamic_friction(dynamicFriction);
    contact_s4_r->set_viscous_friction(viscousFriction);
    contact_s4_r->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s4_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s5_r;
    contact_s5_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s5_r",*s5_r,*contactHalfSpace);
    contact_s5_r->set_stiffness(stiffness);
    contact_s5_r->set_dissipation(dissipation);
    contact_s5_r->set_static_friction(staticFriction);
    contact_s5_r->set_dynamic_friction(dynamicFriction);
    contact_s5_r->set_viscous_friction(viscousFriction);
    contact_s5_r->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s5_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s6_r;
    contact_s6_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_s6_r",*s6_r,*contactHalfSpace);
    contact_s6_r->set_stiffness(stiffness);
    contact_s6_r->set_dissipation(dissipation);
    contact_s6_r->set_static_friction(staticFriction);
    contact_s6_r->set_dynamic_friction(dynamicFriction);
    contact_s6_r->set_viscous_friction(viscousFriction);
    contact_s6_r->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s6_r);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s1_l;
    contact_s1_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s1_l",*s1_l,*contactHalfSpace);
    contact_s1_l->set_stiffness(stiffness);
    contact_s1_l->set_dissipation(dissipation);
    contact_s1_l->set_static_friction(staticFriction);
    contact_s1_l->set_dynamic_friction(dynamicFriction);
    contact_s1_l->set_viscous_friction(viscousFriction);
    contact_s1_l->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s1_l);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s2_l;
    contact_s2_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s2_l",*s2_l,*contactHalfSpace);
    contact_s2_l->set_stiffness(stiffness);
    contact_s2_l->set_dissipation(dissipation);
    contact_s2_l->set_static_friction(staticFriction);
    contact_s2_l->set_dynamic_friction(dynamicFriction);
    contact_s2_l->set_viscous_friction(viscousFriction);
    contact_s2_l->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s2_l);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s3_l;
    contact_s3_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s3_l",*s3_l,*contactHalfSpace);
    contact_s3_l->set_stiffness(stiffness);
    contact_s3_l->set_dissipation(dissipation);
    contact_s3_l->set_static_friction(staticFriction);
    contact_s3_l->set_dynamic_friction(dynamicFriction);
    contact_s3_l->set_viscous_friction(viscousFriction);
    contact_s3_l->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s3_l);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s4_l;
    contact_s4_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s4_l",*s4_l,*contactHalfSpace);
    contact_s4_l->set_stiffness(stiffness);
    contact_s4_l->set_dissipation(dissipation);
    contact_s4_l->set_static_friction(staticFriction);
    contact_s4_l->set_dynamic_friction(dynamicFriction);
    contact_s4_l->set_viscous_friction(viscousFriction);
    contact_s4_l->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s4_l);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s5_l;
    contact_s5_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s5_l",*s5_l,*contactHalfSpace);
    contact_s5_l->set_stiffness(stiffness);
    contact_s5_l->set_dissipation(dissipation);
    contact_s5_l->set_static_friction(staticFriction);
    contact_s5_l->set_dynamic_friction(dynamicFriction);
    contact_s5_l->set_viscous_friction(viscousFriction);
    contact_s5_l->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s5_l);

    OpenSim::SmoothSphereHalfSpaceForce* contact_s6_l;
    contact_s6_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_s6_l",*s6_l,*contactHalfSpace);
    contact_s6_l->set_stiffness(stiffness);
    contact_s6_l->set_dissipation(dissipation);
    contact_s6_l->set_static_friction(staticFriction);
    contact_s6_l->set_dynamic_friction(dynamicFriction);
    contact_s6_l->set_viscous_friction(viscousFriction);
    contact_s6_l->set_transition_velocity(transitionVelocity);
    model->addComponent(contact_s6_l);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<double> x(arg[0], arg[0] + NX);
    std::vector<double> u(arg[1], arg[1] + NU);

    // States and controls
    double ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX+4); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// pro_sup dofs are locked so Qs and Qdots are hard coded (0)
    QsUs[NX] = SimTK::Pi/2;
    QsUs[NX+1] = 0;
    QsUs[NX+2] = SimTK::Pi/2;
    QsUs[NX+3] = 0;
    /// Controls
    double ut[NU+2];
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
    model->realizeAcceleration(*state);

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
    Array<double> Force_values_s1_r = contact_s1_r->getRecordValues(*state);
    Array<double> Force_values_s2_r = contact_s2_r->getRecordValues(*state);
    Array<double> Force_values_s3_r = contact_s3_r->getRecordValues(*state);
    Array<double> Force_values_s4_r = contact_s4_r->getRecordValues(*state);
    Array<double> Force_values_s5_r = contact_s5_r->getRecordValues(*state);
    Array<double> Force_values_s6_r = contact_s6_r->getRecordValues(*state);
    SpatialVec GRF_s1_r;
    GRF_s1_r[0] = Vec3(Force_values_s1_r[3], Force_values_s1_r[4], Force_values_s1_r[5]);
    GRF_s1_r[1] = Vec3(Force_values_s1_r[0], Force_values_s1_r[1], Force_values_s1_r[2]);
    SpatialVec GRF_s2_r;
    GRF_s2_r[0] = Vec3(Force_values_s2_r[3], Force_values_s2_r[4], Force_values_s2_r[5]);
    GRF_s2_r[1] = Vec3(Force_values_s2_r[0], Force_values_s2_r[1], Force_values_s2_r[2]);
    SpatialVec GRF_s3_r;
    GRF_s3_r[0] = Vec3(Force_values_s3_r[3], Force_values_s3_r[4], Force_values_s3_r[5]);
    GRF_s3_r[1] = Vec3(Force_values_s3_r[0], Force_values_s3_r[1], Force_values_s3_r[2]);
    SpatialVec GRF_s4_r;
    GRF_s4_r[0] = Vec3(Force_values_s4_r[3], Force_values_s4_r[4], Force_values_s4_r[5]);
    GRF_s4_r[1] = Vec3(Force_values_s4_r[0], Force_values_s4_r[1], Force_values_s4_r[2]);
    SpatialVec GRF_s5_r;
    GRF_s5_r[0] = Vec3(Force_values_s5_r[3], Force_values_s5_r[4], Force_values_s5_r[5]);
    GRF_s5_r[1] = Vec3(Force_values_s5_r[0], Force_values_s5_r[1], Force_values_s5_r[2]);
    SpatialVec GRF_s6_r;
    GRF_s6_r[0] = Vec3(Force_values_s6_r[3], Force_values_s6_r[4], Force_values_s6_r[5]);
    GRF_s6_r[1] = Vec3(Force_values_s6_r[0], Force_values_s6_r[1], Force_values_s6_r[2]);
    int nCalc_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
    int nToe_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
    appliedBodyForces[nCalc_r] = appliedBodyForces[nCalc_r] + GRF_s1_r + GRF_s2_r + GRF_s3_r + GRF_s4_r;
    appliedBodyForces[nToe_r] = appliedBodyForces[nToe_r] + GRF_s5_r + GRF_s6_r;

    /// Left foot
    Array<double> Force_values_s1_l = contact_s1_l->getRecordValues(*state);
    Array<double> Force_values_s2_l = contact_s2_l->getRecordValues(*state);
    Array<double> Force_values_s3_l = contact_s3_l->getRecordValues(*state);
    Array<double> Force_values_s4_l = contact_s4_l->getRecordValues(*state);
    Array<double> Force_values_s5_l = contact_s5_l->getRecordValues(*state);
    Array<double> Force_values_s6_l = contact_s6_l->getRecordValues(*state);
    SpatialVec GRF_s1_l;
    GRF_s1_l[0] = Vec3(Force_values_s1_l[3], Force_values_s1_l[4], Force_values_s1_l[5]);
    GRF_s1_l[1] = Vec3(Force_values_s1_l[0], Force_values_s1_l[1], Force_values_s1_l[2]);
    SpatialVec GRF_s2_l;
    GRF_s2_l[0] = Vec3(Force_values_s2_l[3], Force_values_s2_l[4], Force_values_s2_l[5]);
    GRF_s2_l[1] = Vec3(Force_values_s2_l[0], Force_values_s2_l[1], Force_values_s2_l[2]);
    SpatialVec GRF_s3_l;
    GRF_s3_l[0] = Vec3(Force_values_s3_l[3], Force_values_s3_l[4], Force_values_s3_l[5]);
    GRF_s3_l[1] = Vec3(Force_values_s3_l[0], Force_values_s3_l[1], Force_values_s3_l[2]);
    SpatialVec GRF_s4_l;
    GRF_s4_l[0] = Vec3(Force_values_s4_l[3], Force_values_s4_l[4], Force_values_s4_l[5]);
    GRF_s4_l[1] = Vec3(Force_values_s4_l[0], Force_values_s4_l[1], Force_values_s4_l[2]);
    SpatialVec GRF_s5_l;
    GRF_s5_l[0] = Vec3(Force_values_s5_l[3], Force_values_s5_l[4], Force_values_s5_l[5]);
    GRF_s5_l[1] = Vec3(Force_values_s5_l[0], Force_values_s5_l[1], Force_values_s5_l[2]);
    SpatialVec GRF_s6_l;
    GRF_s6_l[0] = Vec3(Force_values_s6_l[3], Force_values_s6_l[4], Force_values_s6_l[5]);
    GRF_s6_l[1] = Vec3(Force_values_s6_l[0], Force_values_s6_l[1], Force_values_s6_l[2]);
    int nCalc_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
    int nToe_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
    appliedBodyForces[nCalc_l] = appliedBodyForces[nCalc_l] + GRF_s1_l + GRF_s2_l + GRF_s3_l + GRF_s4_l;
    appliedBodyForces[nToe_l] = appliedBodyForces[nToe_l] + GRF_s5_l + GRF_s6_l;
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

    const SimTK::Vec3 translation_pelvis_imu(-0.17825090079013006, 0.06148338297319611, -0.0039742631657566363);
    const SimTK::Vec3 translation_torso_imu(0.11141304895632698, 0.32812980850924067, -0.012040552365984683);

    const SimTK::Vec3 translation_femur_r_imu(0.043933841399841023, -0.14958344693305592, 0.099866089437989247);
    const SimTK::Vec3 translation_femur_l_imu(0.057420458331616797, -0.12466095809783695, -0.1025425014551595);

    const SimTK::Vec3 translation_tibia_r_imu(0.047080175165178595, -0.11466976609963364, 0.0052928998697051033);
    const SimTK::Vec3 translation_tibia_l_imu(0.053945123225004443, -0.12392203671935126, -0.006707066365313652);

    const SimTK::Vec3 translation_calcn_r_imu(0.14808446197590852, 0.040621301301533644, 0.035088429802814902);
    const SimTK::Vec3 translation_calcn_l_imu(0.15253209992514105, 0.045078665693204636, -0.044293377012806667);

    SimTK::Vec3 angVel_pelvis_imu_inG   = model->getBodySet().get("pelvis").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_torso_imu_inG    = model->getBodySet().get("torso").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_femur_l_imu_inG  = model->getBodySet().get("femur_l").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_femur_r_imu_inG  = model->getBodySet().get("femur_r").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_tibia_l_imu_inG  = model->getBodySet().get("tibia_l").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_tibia_r_imu_inG  = model->getBodySet().get("tibia_r").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_calcn_l_imu_inG  = model->getBodySet().get("calcn_l").getAngularVelocityInGround(*state);
    SimTK::Vec3 angVel_calcn_r_imu_inG  = model->getBodySet().get("calcn_r").getAngularVelocityInGround(*state);

    Vector_<SpatialVec> A_GB(nbodies);
    model->getMatterSubsystem().calcBodyAccelerationFromUDot(*state, knownUdot, A_GB);

    SimTK::Vec3 linAcc_pelvis_imu_inG   = getLinearAccelerationIMUInGround(*model, *state, A_GB, "pelvis",  gravity, translation_pelvis_imu);
    SimTK::Vec3 linAcc_torso_imu_inG    = getLinearAccelerationIMUInGround(*model, *state, A_GB, "torso",   gravity, translation_torso_imu);
    SimTK::Vec3 linAcc_femur_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "femur_l", gravity, translation_femur_l_imu);
    SimTK::Vec3 linAcc_femur_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "femur_r", gravity, translation_femur_r_imu);
    SimTK::Vec3 linAcc_tibia_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "tibia_l", gravity, translation_tibia_l_imu);
    SimTK::Vec3 linAcc_tibia_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "tibia_r", gravity, translation_tibia_r_imu);
    SimTK::Vec3 linAcc_calcn_l_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "calcn_l", gravity, translation_calcn_l_imu);
    SimTK::Vec3 linAcc_calcn_r_imu_inG  = getLinearAccelerationIMUInGround(*model, *state, A_GB, "calcn_r", gravity, translation_calcn_r_imu);

    // Residual forces in OpenSim order
    double res_os[ndofr];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < ndofr; ++i) res_os[i] =
            (residualMobilityForces[indicesSimbodyInOS[i]]);
    // Extract results
    int nc = 3;
    /// Residual forces
    /// We do want to extract the pro_sup torques (last two -> till NU)
    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];
    //for (int i = 0; i < ndofr; ++i) res[0][i] = (residualMobilityForces[indicesSimbodyInOS[i]]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 0*nc] = (angVel_pelvis_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 1*nc] = (angVel_torso_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 2*nc] = (angVel_femur_l_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 3*nc] = (angVel_femur_r_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 4*nc] = (linAcc_pelvis_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 5*nc] = (linAcc_torso_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 6*nc] = (linAcc_femur_l_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 7*nc] = (linAcc_femur_r_imu_inG[i]);

    return 0;

}

/* In main(), the Recorder is used to save the expression graph of function F.
This expression graph is saved as a MATLAB function named foo.m. From this
function, a c-code can be generated via CasADi and then compiled as a dll. This
dll is then imported in MATLAB as an external function. With this workflow,
CasADi can use algorithmic differentiation to differentiate the function F.
*/
int main() {
    double x[NX];
	double u[NU];
	double tau[NR];

	for (int i = 0; i < NX; ++i) x[i] = -0.2;
	for (int i = 0; i < NU; ++i) u[i] = -1;

	const double* Recorder_arg[n_in] = { x,u };
	double* Recorder_res[n_out] = { tau };

	F_generic(Recorder_arg, Recorder_res);

    for (int i = 0; i < NR; ++i) std::cout << Recorder_res[0][i] << std::endl;;

	return 0;
}