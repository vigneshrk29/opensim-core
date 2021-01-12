/*  This code describes the OpenSim model and the skeleton dynamics
    Author: Antoine Falisse
    Contributor: Joris Gillis, Gil Serrancoli, Chris Dembia
*/
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Simulation/SimbodyEngine/ScapulothoracicJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PointConstraint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
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
constexpr int n_in = 4;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 11;        // # coordinates (DoFs) in OpenSim model
constexpr int NX = ndof * 2;      // # states (generalized positions and velocities)
constexpr int NU = ndof;        // # slack controls (generalized accelerations)
constexpr int NLambda = 3;		// # Lagrange multipliers (1 per constraint)
constexpr int NR = ndof + 2 * NLambda + 6 + 2*3 + 3; // # outputs

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
    SimTK::Array_<int> idxOSInSimbody(s.getNU());
    s.updU() = 0;
    for (int iy = 0; iy < s.getNU(); ++iy) {
        s.updU()[iy] = SimTK::NaN;
        const auto svValues = model.getStateVariableValues(s);
        for (int isv = 0; isv < svNames.size(); ++isv) {
            if (SimTK::isNaN(svValues[isv])) {
                s.updU()[iy] = 0;
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
    SimTK::Array_<int> idxSimbodyInOS(s.getNU());
    for (int iy = 0; iy < s.getNU(); ++iy) {
        for (int iyy = 0; iyy < s.getNU(); ++iyy) {
            if (idxOSInSimbody[iyy] == iy) {
                idxSimbodyInOS[iy] = iyy;
                break;
            }
        }
    }
    return idxSimbodyInOS;
}

const SimTK::Vec3 getLinearAccelerationIMUInG(Model& model,
    const State& s, const Vector_<SpatialVec>& A_GB,
    const std::string& frameName, const Vec3& gravity_G,
    const Vec3& imuTranslation_B)
{
    /// Rotation
    const SimTK::Rotation R_GB = model.getBodySet().get(frameName).getMobilizedBody().getBodyTransform(s).R();
    /// Body linear acceleration in ground
    const SimTK::Vec3 linAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][1];
    /// Body angular acceleration in ground
    const SimTK::Vec3 angAcc_G = A_GB[model.getBodySet().get(frameName).getMobilizedBodyIndex()][0];
    /// Body angular velocity in ground
    const SimTK::Vec3 angVel_G = model.getBodySet().get(frameName).getAngularVelocityInGround(s);
    /// Sensor linear acceleration
    /// See van den Bogert et al. (1995), equation (1), p949.
    const SimTK::Vec3 linAcc_imu_G = (linAcc_G - gravity_G) + SimTK::cross(angAcc_G, R_GB*imuTranslation_B) + SimTK::cross(angVel_G, SimTK::cross(angVel_G, R_GB*imuTranslation_B));

    return linAcc_imu_G;
}

// Function F
template<typename T>
int F_generic(const T** arg, T** res) {

    // OpenSim model: create components
    /// Model
    OpenSim::Model* model;
    /// Bodies
    OpenSim::Body* thorax;
    OpenSim::Body* clavicle;
    OpenSim::Body* scapula;
    OpenSim::Body* humerus;
    OpenSim::Body* ulna;
    OpenSim::Body* radius;
    OpenSim::Body* hand;
    /// Joints
    OpenSim::WeldJoint* ground_thorax;
    OpenSim::CustomJoint* sternoclavicular;
    OpenSim::ScapulothoracicJoint* scapulothoracic;
    OpenSim::CustomJoint* GlenoHumeral;
    OpenSim::CustomJoint* elbow;
    OpenSim::CustomJoint* radioulnar;
    OpenSim::WeldJoint* rc;
    /// Constraints
    OpenSim::PointConstraint* ac;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();

    /// Bodies
    thorax = new OpenSim::Body("thorax", 20.866630055554726, Vec3(-0.062718835701941519, -0.15165074549763785, 0), Inertia(1.0142312761589201, 0.43156618119044199, 1.0142312761589201, 0, 0, 0));
    clavicle = new OpenSim::Body("clavicle", 0.23027421524576702, Vec3(-0.014383542433100529, 0.0059574690139403749, 0.061361724799005118), Inertia(0.00023667930695988186, 0.00023667930695988186, 4.3337908599160701e-05, -1.8517650778800334e-05, -6.8236388598144618e-05, 5.2398748794970591e-05));
    scapula = new OpenSim::Body("scapula", 0.35591058418327082, Vec3(-0.047672109433058503, -0.033763145943312339, -0.041755323033245646), Inertia(0.00078311608990092778, 0.00072446748308865937, 0.00078311608990092778, 0.00028312352957027363, 0.00025778440319030612, 0.00015190073029839582));
    humerus = new OpenSim::Body("humerus", 1.7650184258291826, Vec3(0.021950117736946696, -0.11878719323199137, -0.013136206000296318), Inertia(0.010210502436341196, 0.0021229483602644233, 0.010210502436341196, -0.00028909835589940511, -0.00019348702642473181, 0.00010227353907687952));
    ulna = new OpenSim::Body("ulna", 0.86632238509629322, Vec3(0.0097179999999999992, -0.095949999999999994, 0.024289999999999999), Inertia(0.004242651832557889, 0.0009037091378051444, 0.0038742708310241366, 0.00024838284975754575, -5.968556014211773e-05, 0.00085589798654225316));
    radius = new OpenSim::Body("radius", 0.18309319565592519, Vec3(0.03363, -0.18156, 0.015599999999999999), Inertia(0.00034377001547383903, 6.9435899842287734e-05, 3.1555359833508332e-05, 2.3623411460058158e-05, -3.323266907453437e-06, 5.0303601443481501e-05));
    hand = new OpenSim::Body("hand", 0.41148941660685245, Vec3(0.00059999999999999995, -0.090499999999999997, -0.036499999999999998), Inertia(0.00050060626740342222, 0.00014923349508941846, 0.00050060626740342222, 0, 0, 0));

    /// Add bodies to model
    model->addBody(thorax);
    model->addBody(clavicle);
    model->addBody(scapula);
    model->addBody(humerus);
    model->addBody(ulna);
    model->addBody(radius);
    model->addBody(hand);

    /// Joints
    /// Clavicle transform
    SpatialTransform st_clavicle;
    st_clavicle[0].setCoordinateNames(OpenSim::Array<std::string>("clav_prot", 1, 1));
    st_clavicle[0].setFunction(new LinearFunction());
    st_clavicle[0].setAxis(Vec3(0.015299999999999999, 0.98929869999999998, -0.14509996));
    st_clavicle[1].setCoordinateNames(OpenSim::Array<std::string>("clav_elev", 1, 1));
    st_clavicle[1].setFunction(new LinearFunction());
    st_clavicle[1].setAxis(Vec3(-0.99447253999999996, 0, -0.10499695000000001));
    st_clavicle[2].setFunction(new Constant(0));
    st_clavicle[2].setAxis(Vec3(0, 1, 0));
    /// GlenoHumeral transform
    SpatialTransform st_glenoHumeral;
    st_glenoHumeral[0].setCoordinateNames(OpenSim::Array<std::string>("plane_elv", 1, 1));
    st_glenoHumeral[0].setFunction(new LinearFunction());
    st_glenoHumeral[0].setAxis(Vec3(0, 1, 0));
    st_glenoHumeral[1].setCoordinateNames(OpenSim::Array<std::string>("shoulder_elv", 1, 1));
    st_glenoHumeral[1].setFunction(new LinearFunction());
    st_glenoHumeral[1].setAxis(Vec3(-1, 0, 0));
    st_glenoHumeral[2].setCoordinateNames(OpenSim::Array<std::string>("axial_rot", 1, 1));
    st_glenoHumeral[2].setFunction(new LinearFunction());
    st_glenoHumeral[2].setAxis(Vec3(-0.084599999999999995, 0.99470000000000003, -0.058400000000000001));
    /// Elbow transform
    SpatialTransform st_elbow;
    st_elbow[2].setCoordinateNames(OpenSim::Array<std::string>("elbow_flexion", 1, 1));
    st_elbow[2].setFunction(new LinearFunction());
    st_elbow[2].setAxis(Vec3(0.049400010000000001, 0.036600010000000002, 0.99810825000000003));
    /// Radioulnar transform
    SpatialTransform st_radioulnar;
    st_radioulnar[1].setCoordinateNames(OpenSim::Array<std::string>("pro_sup", 1, 1));
    st_radioulnar[1].setFunction(new LinearFunction());
    st_radioulnar[1].setAxis(Vec3(-0.017160990000000001, 0.99266564000000002, -0.11966796));
    /// Joint specifications
    ground_thorax = new WeldJoint("ground_thorax", model->getGround(), Vec3(0), Vec3(0), *thorax, Vec3(0), Vec3(0));
    sternoclavicular = new CustomJoint("sternoclavicular", *thorax, Vec3(0.0067122864657020035, 0.0070722686906692523, 0.027126698859144753), Vec3(0), *clavicle, Vec3(0), Vec3(0), st_clavicle);
    scapulothoracic = new ScapulothoracicJoint("scapulothoracic", *thorax, Vec3(-0.031224595304553068, -0.017655138246985499, 0.074873670018500482), Vec3(0, -1.0700000000000001, 0),
        *scapula, Vec3(-0.052143747840941085, -0.037628215513715965, -0.053471203903169034), Vec3(-0.5181, -1.1415999999999999, -0.28539999999999999),
        Vec3(0.082997500000000002, 0.199991, 0.083001000000000005), Vec2(0), 0);
    GlenoHumeral = new CustomJoint("GlenoHumeral", *scapula, Vec3(-0.0083245245186914282, -0.032770446729409215, 0.0085935842954960595), Vec3(0), *humerus, Vec3(0), Vec3(0), st_glenoHumeral);
    elbow = new CustomJoint("elbow", *humerus, Vec3(0.0074139454488117091, -0.24622301325984491, -0.01267257156393734), Vec3(0), *ulna, Vec3(0), Vec3(0), st_elbow);
    radioulnar = new CustomJoint("radioulnar", *ulna, Vec3(0.00040000000000000002, -0.011502999999999999, 0.019998999999999999), Vec3(0), *radius, Vec3(0), Vec3(0), st_radioulnar);
    rc = new WeldJoint("rc", *radius, Vec3(0.017999999999999999, -0.24199999999999999, 0.025000000000000001), Vec3(0), *hand, Vec3(0), Vec3(0));

    /// Add joints to model
    model->addJoint(ground_thorax);
    model->addJoint(sternoclavicular);
    model->addJoint(scapulothoracic);
    model->addJoint(GlenoHumeral);
    model->addJoint(elbow);
    model->addJoint(radioulnar);
    model->addJoint(rc);

    /// Constraints
    ac = new PointConstraint(*clavicle, Vec3(-0.037887400000000002, 0.018924, 0.135987), *scapula, Vec3(-0.011828699999999999, 0.000106018, -0.0145474));
    model->addConstraint(ac);

    // Initialize system and state
    SimTK::State* state;
    state = new State(model->initSystem());

    // Read inputs
    std::vector<T> x(arg[0], arg[0] + NX);
    std::vector<T> u(arg[1], arg[1] + NU);
    std::vector<T> lambda(arg[2], arg[2] + NLambda);
    std::vector<T> gamma(arg[3], arg[3] + NLambda);

    // States and controls
    T ua[NU]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states
    Vector multipliers(NLambda); /// Lagrange multipliers
    Vector vel_corrs(NLambda); /// Velocity correctors

    // Assign inputs to model variables
    /// States
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// Controls
    T ut[NU];
    for (int i = 0; i < NU; ++i) ut[i] = u[i];
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
    for (int i = 0; i < ndof; ++i) ua[i] = ut[indicesOSInSimbody[i]];
    /// Lagrange multipliers
    for (int i = 0; i < NLambda; ++i) multipliers[i] = lambda[i];
    /// Velocity correctors
    for (int i = 0; i < NLambda; ++i) vel_corrs[i] = gamma[i];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    // Compute residual forces
    /// appliedMobilityForces (# mobilities)
    Vector appliedMobilityForces(ndof);
    appliedMobilityForces.setToZero();
    /// appliedBodyForces (# Mobilized bodies - can be different from # bodies)
    Vector_<SpatialVec> appliedBodyForces;
    appliedBodyForces.resize(model->getMatterSubsystem().getNumBodies());
    appliedBodyForces.setToZero();
    /// Set gravity
    Vec3 gravity(0);
    gravity[1] = -9.80665;
    /// Add weight to appliedBodyForces
    for (MobilizedBodyIndex i(0); i < model->getMatterSubsystem().getNumBodies(); ++i) {
        model->getMatterSubsystem().addInStationForce(*state,
            i,
            model->getMatterSubsystem().getMobilizedBody(i).getBodyMassCenterStation(*state),
            model->getMatterSubsystem().getMobilizedBody(i).getBodyMass(*state) * gravity,
            appliedBodyForces);
    }
    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Compute residual forces.
    /*/// APPROACH 1: Compute constraint forces separately.
    /// Compute constraint forces.
    Vector_<SpatialVec> constraintBodyForces;
    constraintBodyForces.resize(model->getMatterSubsystem().getNumBodies());
    constraintBodyForces.setToZero();
    Vector constraintMobilityForces(ndof);
    constraintMobilityForces.setToZero();
    model->getMatterSubsystem().calcConstraintForcesFromMultipliers(*state,
        -multipliers, constraintBodyForces, constraintMobilityForces);
    /// Add up mobility and body forces.
    appliedBodyForces += constraintBodyForces;
    appliedMobilityForces += constraintMobilityForces;
    /// Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(*state,
        appliedMobilityForces, appliedBodyForces, knownUdot, residualMobilityForces);*/
        /// APPROACH 2: Constraint forces are taken into account in calcResidualForce.
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForce(*state, appliedMobilityForces,
        appliedBodyForces, knownUdot, multipliers, residualMobilityForces);

    // Constraint errors
    /// Position-level constraint errors.
    /// Same as constraint.getPositionErrorsAsVector(*state).
    const auto& qerr = state->getQErr();
    /// Velocity-level constraint errors.
    /// Same as constraint.getVelocityErrorsAsVector(*state).
    const auto& uerr = state->getUErr();
    ///// Acceleration-level constraint errors.
    ///// pvaerr = G udot - b(t,q,u)
    ///// All acceleration-level constraints are included: holonomic second
    ///// derivatives, nonholonomic first derivatives, and acceleration - only constraints.
    //Vector aerr;
    //model->getMatterSubsystem().calcConstraintAccelerationErrors(*state, knownUdot, aerr);

    // Velocity correction based on Posa et al. (2016).
    /// This should be added to the qdots in the collocation equations (see eq.8 in paper).
    Vector qdotCorr;
    model->getMatterSubsystem().multiplyByGTranspose(*state, vel_corrs, qdotCorr);

    /// IMUs
    const SimTK::Vec3 translation_radius_imu(0.018, -0.22, 0.04);
    const SimTK::Vec3 orientation_radius_imu(0, 0, 0);
    // Angular velocity measured with respect to the ground frame and expressed in the GROUND frame.
    SimTK::Vec3 angVel_radius_imu_inG = model->getBodySet().get("radius").getAngularVelocityInGround(*state);

    Vector_<SpatialVec> A_GB(model->getMatterSubsystem().getNumBodies());
    model->getMatterSubsystem().calcBodyAccelerationFromUDot(*state, knownUdot, A_GB);
    // Linear acceleration measured with respect to the ground frame and expressed in the GROUND frame.
    SimTK::Vec3 linAcc_radius_imu_inG = getLinearAccelerationIMUInG(*model, *state, A_GB, "radius", gravity, translation_radius_imu);

    // Orientations (rotation matrices)
    const SimTK::Rotation R_radius_imu_inB(
        SimTK::BodyOrSpaceType::BodyRotationSequence,
        orientation_radius_imu[0], SimTK::XAxis,
        orientation_radius_imu[1], SimTK::YAxis,
        orientation_radius_imu[2], SimTK::ZAxis);
    const SimTK::Rotation R_radius_GB = radius->getTransformInGround(*state).R();
    const SimTK::Rotation R_radius_imu_inG = R_radius_GB * R_radius_imu_inB;
    const SimTK::Vec3 EA_radius_imu_inG = R_radius_imu_inG.convertRotationToBodyFixedXYZ();

    // Outputs
    /// 1. Residual forces.
    /// Adjust for different orders in OpenSim and Simbody.
    T res_os[ndof];
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < ndof; ++i) res_os[i] =
        value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
    for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];
    /// 2. Contraint errors.
    for (int i = 0; i < NLambda; ++i) res[0][i + NU + 0 * NLambda] = value<T>(qerr[i]);
    for (int i = 0; i < NLambda; ++i) res[0][i + NU + 1 * NLambda] = value<T>(uerr[i]);
    //for (int i = 0; i < NLambda; ++i) res[0][i + NU + 2 * NLambda] = value<T>(aerr[i]);
    /// 3. Velocity correctors.
    /// TODO: should sort of automate this.
    res[0][0 + NU + 2 * NLambda] = value<T>(qdotCorr[3]);
    res[0][1 + NU + 2 * NLambda] = value<T>(qdotCorr[4]);
    res[0][2 + NU + 2 * NLambda] = value<T>(qdotCorr[5]);
    res[0][3 + NU + 2 * NLambda] = value<T>(qdotCorr[6]);
    res[0][4 + NU + 2 * NLambda] = value<T>(qdotCorr[7]);
    res[0][5 + NU + 2 * NLambda] = value<T>(qdotCorr[8]);
    /// 4. IMU data
    int nc = 3;
    for (int i = 0; i < nc; ++i) res[0][i + NU + 2 * NLambda + 6 + 0 * nc] = value<T>(angVel_radius_imu_inG[i]);
    for (int i = 0; i < nc; ++i) res[0][i + NU + 2 * NLambda + 6 + 1 * nc] = value<T>(linAcc_radius_imu_inG[i]);
    /// 5. Orientation data (in ground frame) expressed as body fixed XYZ Euler angles
    for (int i = 0; i < nc; ++i) res[0][i + NU + 2 * NLambda + 6 + 2 * nc] = value<T>(EA_radius_imu_inG[i]);

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
    Recorder lambda[NLambda];
    Recorder gamma[NLambda];
    Recorder tau[NR];

    for (int i = 0; i < NX; ++i) x[i] <<= 0;
    for (int i = 0; i < NU; ++i) u[i] <<= 0;

    //for (int i = 1; i < NX; i += 2) x[i] <<= 0;
    //// ground_thorax/ground_thorax_rot_x/value
    //x[0] <<= 3.994 * SimTK::Pi.getValue() / 180;
    //// ground_thorax/ground_thorax_rot_y/value
    //x[2] <<= -3.619 * SimTK::Pi.getValue() / 180;
    //// ground_thorax / ground_thorax_rot_z / value
    //x[4] <<= -10.437 * SimTK::Pi.getValue() / 180;
    //// sternoclavicular/clav_prot/value
    //x[6] <<= -17.183 * SimTK::Pi.getValue() / 180;
    //// sternoclavicular/clav_elev/value
    //x[8] <<= -3.709 * SimTK::Pi.getValue() / 180;
    //// scapulothoracic/scapulothoracic_coord_0/value
    //x[10] <<= -17.851 * SimTK::Pi.getValue() / 180;
    //// scapulothoracic/scapulothoracic_coord_1/value
    //x[12] <<= -2.893 * SimTK::Pi.getValue() / 180;
    //// scapulothoracic/scapulothoracic_coord_2/value
    //x[14] <<= 11.929 * SimTK::Pi.getValue() / 180;
    //// scapulothoracic/scapulothoracic_coord_3/value
    //x[16] <<= 15.675 * SimTK::Pi.getValue() / 180;
    //// GlenoHumeral/plane_elv/value
    //x[18] <<= 0 * SimTK::Pi.getValue() / 180;
    //// GlenoHumeral/shoulder_elv/value
    //x[20] <<= 18.428 * SimTK::Pi.getValue() / 180;
    //// GlenoHumeral/axial_rot/value
    //x[22] <<= 58.948 * SimTK::Pi.getValue() / 180;
    //// elbow/elbow_flexion/value
    //x[24] <<= 0 * SimTK::Pi.getValue() / 180;
    //// radioulnar/pro_sup/value
    //x[26] <<= 34.320 * SimTK::Pi.getValue() / 180;

    for (int i = 0; i < NLambda; ++i) lambda[i] <<= 0;
    for (int i = 0; i < NLambda; ++i) gamma[i] <<= 0;

    const Recorder* Recorder_arg[n_in] = { x,u,lambda,gamma };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) {
        Recorder_res[0][i] >>= res[i];
        //std::cout << res[i] << std::endl;
    }

    Recorder::stop_recording();

    std::vector<std::vector <int> > normal;
    for (int i = 0; i < 2; i++) {
        normal.emplace_back();
        for (int j = 0; j < 3; j++) {
            normal.back().push_back(j);
        }
    }

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << normal[i][j] << " ";
        }
        std::cout << std::endl;
    }

    return 0;

}
