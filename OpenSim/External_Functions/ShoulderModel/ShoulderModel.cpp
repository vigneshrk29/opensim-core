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
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 14;        // # degrees of freedom
constexpr int NX = ndof*2;      // # states
constexpr int NU = ndof;        // # controls
constexpr int NR = ndof;        // # residual torques

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
    OpenSim::CustomJoint* ground_thorax;
    OpenSim::CustomJoint* sternoclavicular;
    OpenSim::ScapulothoracicJoint* scapulothoracic;
    OpenSim::CustomJoint* GlenoHumeral;
    OpenSim::CustomJoint* elbow;
    OpenSim::CustomJoint* radioulnar;
    OpenSim::WeldJoint* rc;

    // OpenSim model: initialize components
    /// Model
    model = new OpenSim::Model();
    /// Body specifications
	thorax = new OpenSim::Body("thorax", 20.45504587680885, Vec3(-0.059099856018640096, -0.14860100703035389, 0), Inertia(1.1470502930519131, 0.488082083684654, 1.1470502930519131, 0, 0, 0));
	clavicle = new OpenSim::Body("clavicle", 0.18977290868104679, Vec3(-0.011100639351336325, 0.0063717609081909501, 0.054170293336786662), Inertia(0.00026767372964065194, 0.00026767372964065194, 4.901323980777646e-05, -2.0942636311609564e-05, -7.7172309095742937e-05, 5.9260645548665655e-05));  
	scapula = new OpenSim::Body("scapula", 0.50159886749985216, Vec3(-0.054689731740849026, -0.035031350908916951, -0.043714611940415045), Inertia(0.00088566933551532175, 0.00081934037957356279, 0.00088566933551532175, 0.00032020007191396071, 0.0002915426512417284, 0.00017179294437008979));
	humerus = new OpenSim::Body("humerus", 1.8801022951697861, Vec3(0.018059994846920106, -0.14010038452973017, -0.012749884499947897), Inertia(0.011547622408340414, 0.0024009598165791765, 0.011547622408340414, -0.00032695733374650814, -0.00021882518867171686, 0.00011566680670100846));
	ulna = new OpenSim::Body("ulna", 0.97977194064209672, Vec3(0.0097179999999999992, -0.095949999999999994, 0.024289999999999999), Inertia(0.004798249809730995, 0.0010220546888268683, 0.0043816273433401652, 0.00028090991404096664, -6.7501703863110161e-05, 0.00096798241127401597)); 
	radius = new OpenSim::Body("radius", 0.20707023010403855, Vec3(0.03363, -0.18156, 0.015599999999999999), Inertia(0.00038878853991280527, 7.8528902760773868e-05, 3.5687703184882663e-05, 2.6717023695786481e-05, -3.7584665053130287e-06, 5.68911274318373e-05));
	hand = new OpenSim::Body("hand", 0.46537615926635378, Vec3(0.00059999999999999995, -0.090499999999999997, -0.036499999999999998), Inertia(0.00056616333890175264, 0.00016877642042726428, 0.00056616333890175264, 0, 0, 0));
	
	/// Joints
    /// Ground-Thorax transform
    SpatialTransform st_ground_thorax;
	st_ground_thorax[0].setCoordinateNames(OpenSim::Array<std::string>("ground_thorax_rot_x", 1, 1));
	st_ground_thorax[0].setFunction(new LinearFunction());
	st_ground_thorax[0].setAxis(Vec3(1, 0, 0));
	st_ground_thorax[1].setCoordinateNames(OpenSim::Array<std::string>("ground_thorax_rot_y", 1, 1));
	st_ground_thorax[1].setFunction(new LinearFunction());
	st_ground_thorax[1].setAxis(Vec3(0, 1, 0));
	st_ground_thorax[2].setCoordinateNames(OpenSim::Array<std::string>("ground_thorax_rot_z", 1, 1));
	st_ground_thorax[2].setFunction(new LinearFunction());
	st_ground_thorax[2].setAxis(Vec3(0, 0, 1));
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
	ground_thorax = new CustomJoint("ground_thorax", model->getGround(), Vec3(0), Vec3(0), *thorax, Vec3(0), Vec3(0), st_ground_thorax);
	sternoclavicular = new CustomJoint("sternoclavicular", *thorax, Vec3(0.006324976528009351, 0.0069300434097705312, 0.025464431687451307), Vec3(0), *clavicle, Vec3(0), Vec3(0), st_clavicle);
	scapulothoracic = new ScapulothoracicJoint("scapulothoracic", *thorax, Vec3(-0.029422884944949112, -0.017300088530082611, 0.070285568667052628), Vec3(0, -0.87, 0),
		*scapula, Vec3(-0.059819622317908099, -0.039041599498770452, -0.055980238178321401), Vec3(-0.5181, -1.1415999999999999, -0.28539999999999999),
		Vec3(0.082997500000000002, 0.199991, 0.083001000000000005), Vec2(0), 0);
	GlenoHumeral = new CustomJoint("GlenoHumeral", *scapula, Vec3(-0.0095499447834722345, -0.034001364113029092, 0.008996821850851211), Vec3(0), *humerus, Vec3(0), Vec3(0), st_glenoHumeral);
	elbow = new CustomJoint("elbow", *humerus, Vec3(0.0061000044831427585, -0.29040116109488806, -0.01229988504701285), Vec3(0), *ulna, Vec3(0), Vec3(0), st_elbow);
	radioulnar = new CustomJoint("radioulnar", *ulna, Vec3(0.00040000000000000002, -0.011502999999999999, 0.019998999999999999), Vec3(0), *radius, Vec3(0), Vec3(0), st_radioulnar);
	rc = new WeldJoint("rc", *radius, Vec3(0.017999999999999999, -0.24199999999999999, 0.025000000000000001), Vec3(0), *hand, Vec3(0), Vec3(0));

    /// Add bodies and joints to model
    model->addBody(thorax);		    model->addJoint(ground_thorax);
    model->addBody(clavicle);		model->addJoint(sternoclavicular);
    model->addBody(scapula);		model->addJoint(scapulothoracic);
    model->addBody(humerus);		model->addJoint(GlenoHumeral);
    model->addBody(ulna);			model->addJoint(elbow);
    model->addBody(radius);			model->addJoint(radioulnar);
    model->addBody(hand);			model->addJoint(rc);

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
    /// appliedBodyForces (# Mobilized bodies - can be different from # bodies)
	Vector_<SpatialVec> appliedBodyForces;
	appliedBodyForces.resize(model->getMatterSubsystem().getNumBodies());
	appliedBodyForces.setToZero();
    /// Set gravity
    Vec3 gravity(0);
    gravity[1] = -9.80665;
	/// Add weight to appliedBodyForces
	for (int i = 0; i < model->getMatterSubsystem().getNumBodies(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
			MobilizedBodyIndex(i),
			model->getMatterSubsystem().getMobilizedBody(MobilizedBodyIndex(i)).getBodyMassCenterStation(*state),
			model->getMatterSubsystem().getMobilizedBody(MobilizedBodyIndex(i)).getBodyMass(*state) * gravity,
			appliedBodyForces);
	}
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
	
	// Residual forces in OpenSim order
	T res_os[ndof];
	/// OpenSim and Simbody have different state orders so we need to adjust
	auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
	for (int i = 0; i < ndof; ++i) res_os[i] =
		value<T>(residualMobilityForces[indicesSimbodyInOS[i]]);
	for (int i = 0; i < NU; ++i) res[0][i] = res_os[i];

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

	//for (int i = 0; i < NX; ++i) x[i] <<= 0;
	for (int i = 0; i < NU; ++i) u[i] <<= 0;

	for (int i = 1; i < NX; i += 2) x[i] <<= 0;
	// ground_thorax/ground_thorax_rot_x/value
	x[0] <<= 3.994 * SimTK::Pi.getValue() / 180;
	// ground_thorax/ground_thorax_rot_y/value
	x[2] <<= -3.619 * SimTK::Pi.getValue() / 180;
	// ground_thorax / ground_thorax_rot_z / value
	x[4] <<= -10.437 * SimTK::Pi.getValue() / 180;
	// sternoclavicular/clav_prot/value
	x[6] <<= -17.183 * SimTK::Pi.getValue() / 180;
	// sternoclavicular/clav_elev/value
	x[8] <<= -3.709 * SimTK::Pi.getValue() / 180;
	// scapulothoracic/scapulothoracic_coord_0/value
	x[10] <<= -17.851 * SimTK::Pi.getValue() / 180;
	// scapulothoracic/scapulothoracic_coord_1/value
	x[12] <<= -2.893 * SimTK::Pi.getValue() / 180;
	// scapulothoracic/scapulothoracic_coord_2/value
	x[14] <<= 11.929 * SimTK::Pi.getValue() / 180;
	// scapulothoracic/scapulothoracic_coord_3/value
	x[16] <<= 15.675 * SimTK::Pi.getValue() / 180;
	// GlenoHumeral/plane_elv/value
	x[18] <<= 0 * SimTK::Pi.getValue() / 180;
	// GlenoHumeral/shoulder_elv/value
	x[20] <<= 18.428 * SimTK::Pi.getValue() / 180;
	// GlenoHumeral/axial_rot/value
	x[22] <<= 58.948 * SimTK::Pi.getValue() / 180;
	// elbow/elbow_flexion/value
	x[24] <<= 0 * SimTK::Pi.getValue() / 180;
	// radioulnar/pro_sup/value
	x[26] <<= 34.320 * SimTK::Pi.getValue() / 180;

    const Recorder* Recorder_arg[n_in] = { x,u };
    Recorder* Recorder_res[n_out] = { tau };

    F_generic<Recorder>(Recorder_arg, Recorder_res);

    double res[NR];
    for (int i = 0; i < NR; ++i) {
        Recorder_res[0][i] >>= res[i];
		std::cout << Recorder_res[0][i] << std::endl;
    }

    Recorder::stop_recording();

    return 0;

}
