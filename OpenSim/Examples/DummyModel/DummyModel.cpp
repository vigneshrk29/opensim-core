#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/SimbodyEngine/PlanarJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/PinJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/WeldJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/SpatialTransform.h>
#include <OpenSim/Simulation/SimbodyEngine/CustomJoint.h>
#include <OpenSim/Simulation/SimbodyEngine/Joint.h>
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/SimmSpline.h>
#include <OpenSim/Simulation/Model/ConditionalPathPoint.h>
#include <OpenSim/Simulation/Model/MovingPathPoint.h>
#include <OpenSim/Simulation/Model/SmoothSphereHalfSpaceForce.h>
#include <OpenSim/Simulation/InverseDynamicsSolver.h>

#include <iostream>
#include <iterator>
#include <random>
#include <cassert>
#include <algorithm>
#include <vector>
#include <fstream>
#include <chrono>

using namespace SimTK;
using namespace OpenSim;

// Declare inputs/outputs of function F
/// number of vectors in inputs/outputs of function F
constexpr int n_in = 2;
constexpr int n_out = 1;
/// number of elements in input/output vectors of function F
constexpr int ndof = 8;         // degrees of freedom
constexpr int NX = 2 * ndof;    // states
constexpr int NU = ndof;        // controls
constexpr int NR = 20;          // residual forces + segment origin + GRFs + GRMs + PoA

int F_generic(const double** arg, double** res) {
	/// bodies
	OpenSim::Model* model;
	OpenSim::Body* Base;
	OpenSim::Body* Segment1;
	OpenSim::Body* Segment2;
	/// joints
	OpenSim::CustomJoint* ground_base;
	OpenSim::CustomJoint* Joint1;
	OpenSim::CustomJoint* Joint2;
	/// contact elements
    OpenSim::ContactSphere* contactSphere;
    OpenSim::ContactHalfSpace* contactHalfSpace;
	OpenSim::SmoothSphereHalfSpaceForce* Endpoint_contact;
	/// states
	SimTK::State* state;

	// OpenSim model
	/// Model
	model = new OpenSim::Model();
	/// Bodies - Definition
	Base = new OpenSim::Body("Base", 3.4465, Vec3(0, -0.062356, 0), Inertia(0.0070757, 0.091049, 0.096688, 0, 0, 0));
	Segment1 = new OpenSim::Body("Segment1", 3.6424, Vec3(0, -0.325, 0), Inertia(0.16228, 0.0037128, 0.16448, 0, 0, 0));
	Segment2 = new OpenSim::Body("Segment2", 3.6424, Vec3(0, -0.325, 0), Inertia(0.16228, 0.0037128, 0.16448, 0, 0, 0));
	/// Joints - Transforms
	// Ground-pelvis
	SpatialTransform st_ground_base;
	st_ground_base[0].setCoordinateNames(OpenSim::Array<std::string>("Base_rx", 1, 1));
	st_ground_base[0].setFunction(new LinearFunction());
	st_ground_base[0].setAxis(Vec3(1, 0, 0));
	st_ground_base[1].setCoordinateNames(OpenSim::Array<std::string>("Base_ry", 1, 1));
	st_ground_base[1].setFunction(new LinearFunction());
	st_ground_base[1].setAxis(Vec3(0, 1, 0));
	st_ground_base[2].setCoordinateNames(OpenSim::Array<std::string>("Base_rz", 1, 1));
	st_ground_base[2].setFunction(new LinearFunction());
	st_ground_base[2].setAxis(Vec3(0, 0, 1));
	st_ground_base[3].setCoordinateNames(OpenSim::Array<std::string>("Base_tx", 1, 1));
	st_ground_base[3].setFunction(new LinearFunction());
	st_ground_base[3].setAxis(Vec3(1, 0, 0));
	st_ground_base[4].setCoordinateNames(OpenSim::Array<std::string>("Base_ty", 1, 1));
	st_ground_base[4].setFunction(new LinearFunction());
	st_ground_base[4].setAxis(Vec3(0, 1, 0));
	st_ground_base[5].setCoordinateNames(OpenSim::Array<std::string>("Base_tz", 1, 1));
	st_ground_base[5].setFunction(new LinearFunction());
	st_ground_base[5].setAxis(Vec3(0, 0, 1));
	// Joint1
	SpatialTransform st_joint1;
	st_joint1[2].setCoordinateNames(OpenSim::Array<std::string>("joint_1", 1, 1));
	st_joint1[2].setFunction(new LinearFunction());
	st_joint1[2].setAxis(Vec3(0, 0, 1));
	// Joint2
	SpatialTransform st_joint2;
	st_joint2[2].setCoordinateNames(OpenSim::Array<std::string>("joint_2", 1, 1));
	st_joint2[2].setFunction(new LinearFunction());
	st_joint2[2].setAxis(Vec3(0, 0, 1));
	/// Joints - Definition
	ground_base = new CustomJoint("ground_base", model->getGround(), Vec3(0), Vec3(0), *Base, Vec3(-0.3, 0, 0), Vec3(0), st_ground_base);
	Joint1 = new CustomJoint("Joint1", *Base, Vec3(0, -0.15, 0), Vec3(0), *Segment1, Vec3(0), Vec3(0), st_joint1);
	Joint2 = new CustomJoint("Joint2", *Segment1, Vec3(0, -0.65, 0), Vec3(0), *Segment2, Vec3(0), Vec3(0), st_joint2);
	/// bodies and joints
	model->addBody(Base);		    model->addJoint(ground_base);
	model->addBody(Segment1);		model->addJoint(Joint1);
	model->addBody(Segment2);		model->addJoint(Joint2);
	/// Contact elements
	/// Parameters
	double radiusSphere = 0.1;
	double stiffness = 250000;
	double dissipation = 2;
	double staticFriction = 0.8;
	double dynamicFriction = 0.8;
	double viscousFriction = 0.5;
	double transitionVelocity = 0.2;
	Vec3 locSphere_contact = Vec3(0, -0.65, 0);
    Vec3 halfSpaceLocation(0);
	Vec3 halfSpaceOrientation(0, 0, -0.5 * SimTK::Pi);
    contactSphere = new OpenSim::ContactSphere(radiusSphere,locSphere_contact,*Segment2,"contacSphere");
    model->addComponent(contactSphere);
    contactHalfSpace = new OpenSim::ContactHalfSpace(halfSpaceLocation,halfSpaceOrientation,model->getGround(),"contactHalfSpace");
    model->addComponent(contactHalfSpace);
    Endpoint_contact = new OpenSim::SmoothSphereHalfSpaceForce("Endoint",*contactSphere,*contactHalfSpace);
    Endpoint_contact->set_stiffness(stiffness);
    Endpoint_contact->set_dissipation(dissipation);
    Endpoint_contact->set_static_friction(staticFriction);
    Endpoint_contact->set_dynamic_friction(dynamicFriction);
    Endpoint_contact->set_viscous_friction(viscousFriction);
    Endpoint_contact->set_transition_velocity(transitionVelocity);
	///Add contact components to model
    Endpoint_contact->connectSocket_half_space(*contactHalfSpace);
    Endpoint_contact->connectSocket_sphere(*contactSphere);
    model->addComponent(Endpoint_contact);
	/// Initialize  system and state.
	state = new State(model->initSystem());
	// Read inputs
    std::vector<double> x(arg[0], arg[0] + NX);
	std::vector<double> u(arg[1], arg[1] + NU);
	// States and controls
	double ua[NU]; /// joint accelerations (Qdotdots), i.e. controls
	Vector QsUs(NX); /// joint positions (Qs) and velocities (Us), i.e. states
	// Assign inputs to model variables
	for (int i = 0; i < NX; ++i) QsUs[i] = x[i];    // states
	for (int i = 0; i < NU; ++i) ua[i] = u[i];      // controls
	model->setStateVariableValues(*state, QsUs);
	model->realizeAcceleration(*state);
	// Residual forces
	/// appliedMobilityForces (# mobilities)
	Vector appliedMobilityForces(ndof);
	appliedMobilityForces.setToZero();
	/// appliedBodyForces (# bodies + ground)
	Vector_<SpatialVec> appliedBodyForces;
	int nbodies = model->getBodySet().getSize() + 1; // including ground
	appliedBodyForces.resize(nbodies);
	appliedBodyForces.setToZero();
	/// Gravity
	Vec3 gravity(0);
	gravity[1] = -9.80665;    // Y direction is up
	/// Weight
	for (int i = 0; i < model->getBodySet().getSize(); ++i) {
		model->getMatterSubsystem().addInStationForce(*state,
			model->getBodySet().get(i).getMobilizedBodyIndex(),
			model->getBodySet().get(i).getMassCenter(),
			model->getBodySet().get(i).getMass()*gravity, appliedBodyForces);
	}
	// Contact forces and moments
    // Approach 1: extract forces and torques and add them to the appliedBodyForces.
    // The torques are expressed in the ground with respect to the segment origin
	Array<double> Force_endpoint = Endpoint_contact->getRecordValues(*state);
	SpatialVec GRF_endpoint;
    GRF_endpoint[0] = Vec3(Force_endpoint[3], Force_endpoint[4], Force_endpoint[5]); // GRMs (x,y,z) applied on the contact body by the plane (ground)
	GRF_endpoint[1] = Vec3(Force_endpoint[0], Force_endpoint[1], Force_endpoint[2]); // GRFs (x,y,z) applied on the contact body by the plane (ground)
	int n_endpoint = model->getBodySet().get("Segment2").getMobilizedBodyIndex();
	appliedBodyForces[n_endpoint] = appliedBodyForces[n_endpoint] + GRF_endpoint;

    // Approach 2: apply the forces as InStationForces at the contact point.
    // The resulting torques will be applied "implicitly".
    Vec3 normal = Vec3(0, 1, 0);
    Vec3 pos_InGround_locSphere = Segment2->findStationLocationInGround(*state, locSphere_contact);
    Vec3 contactPointpos_InGround_locSphere = pos_InGround_locSphere - radiusSphere*normal;
    Vec3 contactPointpos_InGround_locSphere_adj = contactPointpos_InGround_locSphere - 0.5*contactPointpos_InGround_locSphere[1]*normal;
    /*Vec3 contactPointPos_InBody_locSphere_adj = model->getGround().findStationLocationInAnotherFrame(*state, contactPointpos_InGround_locSphere_adj, *Segment2);
    model->getMatterSubsystem().addInStationForce(*state, Segment2->getMobilizedBodyIndex(), contactPointPos_InBody_locSphere_adj, Vec3(Force_endpoint[0], Force_endpoint[1], Force_endpoint[2]), appliedBodyForces);*/

	// knownUdot
	Vector knownUdot(ndof);
	knownUdot.setToZero();
	for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
	// Residual forces
	Vector residualMobilityForces(ndof);
	residualMobilityForces.setToZero();
	model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
		*state, appliedMobilityForces, appliedBodyForces, knownUdot,
		residualMobilityForces);

    // Coordinates of the origin of the segment to which the contact sphere is attached
    Vec3 Segment2_or  = Segment2->getPositionInGround(*state);


    // Virtual sensors
    // Accelerometer data
    /// Transforms: from ground to body
    SimTK::Rotation R_BG_Base = ~(Base->getMobilizedBody().getBodyRotation(*state));
    SimTK::Rotation R_BG_Segment1 = ~(Segment1->getMobilizedBody().getBodyRotation(*state));
    SimTK::Rotation R_BG_Segment2 = ~(Segment2->getMobilizedBody().getBodyRotation(*state));

    SimTK::Rotation R_GB_Base = (Base->getMobilizedBody().getBodyRotation(*state));
    SimTK::Rotation R_GB_Segment1 = (Segment1->getMobilizedBody().getBodyRotation(*state));
    SimTK::Rotation R_GB_Segment2 = (Segment2->getMobilizedBody().getBodyRotation(*state));

    /// Accelerations: of the segment origin relative to the ground in ground
    SimTK::Vec3 linAcc_Base = Base->getLinearAccelerationInGround(*state);
    std::cout << linAcc_Base << std::endl;
    SimTK::Vec3 linVel_Base = Base->getLinearVelocityInGround(*state);
    //std::cout << linVel_Base << std::endl;

    /// Positions: TODO assume it is the body COM
    SimTK::Vec3 com_Base = Base->getMassCenter();
    /// Angular Velocity
    SimTK::Vec3 angVel_Base_inG = Base->getAngularVelocityInGround(*state); // Corresponds to gyroscope data
    SimTK::Vec3 angVel_Base_inB = R_BG_Base*angVel_Base_inG; // TODO
    //std::cout << angVel_Base_inB << std::endl;

    SimTK::Vec3 angVel_Segment1_inG = Segment1->getAngularVelocityInGround(*state); // Corresponds to gyroscope data
    SimTK::Vec3 angVel_Segment1_inB = R_BG_Segment1*angVel_Segment1_inG; // TODO
    //std::cout << angVel_Segment1_inB << std::endl;

    SimTK::Vec3 angVel_Segment2_inG = Segment2->getAngularVelocityInGround(*state); // Corresponds to gyroscope data
    SimTK::Vec3 angVel_Segment2_inB = R_BG_Segment2*angVel_Segment2_inG; // TODO
    //std::cout << angVel_Segment2_inB << std::endl;

    /// Angular acceleration
    SimTK::Vec3 angAcc_Base_inG = Base->getAngularAccelerationInGround(*state);
    SimTK::Vec3 angAcc_Base_inB = R_BG_Base*angAcc_Base_inG; // TODO
    //std::cout << angAcc_Base_inB << std::endl;

    /// Sensor acceleration
    SimTK::Vec3 accSensor_Base_inB = R_BG_Base * (linAcc_Base - gravity) + SimTK::cross(angAcc_Base_inB, com_Base) + SimTK::cross(angVel_Base_inB, SimTK::cross(angVel_Base_inB, com_Base));
    SimTK::Vec3 accSensor_Base_inG = R_GB_Base * accSensor_Base_inB; // Corresponds to accelerometer data TO CONFIRM
    std::cout << accSensor_Base_inB << std::endl;
    std::cout << accSensor_Base_inG << std::endl;

	if (res[0]) {
		for (int i = 0; i < NU; ++i) {
			res[0][i] = (residualMobilityForces[i]); // residual torques
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof] = (Segment2_or[i]); //  Coordinates of the origin of the segment to which the contact sphere is attached (expressed in the ground)
		}
		for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 3] = (GRF_endpoint[1][i]); // GRF (x, y and z)
		}
        for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 6] = (GRF_endpoint[0][i]); // GRM (x, y and z)
		}
        for (int i = 0; i < 3; ++i) {
			res[0][i + ndof + 9] = (contactPointpos_InGround_locSphere_adj[i]); //  Coordinates of the point of application (expressed in the ground)
		}
	}
	return 0;
}

int main() {
    double x[NX];
	double u[NU];
	double tau[NR];

	for (int i = 0; i < NX; ++i) x[i] = 0;
	for (int i = 0; i < NU; ++i) u[i] = 1;

    for (int i = 12; i < NX; ++i) x[i] = 0;

    //x[7] = 2;
    x[8] = 10;
    //x[5] = 1;
    //x[13] = 1;
    //x[15] = 1;

    //u[3] = 1;



    //x[4] = 10*SimTK::Pi/180;
    //x[8] = 1.15;
    //x[12] = -45*SimTK::Pi/180;
    //x[14] = 80*SimTK::Pi/180;

	const double* Recorder_arg[n_in] = { x,u };
	double* Recorder_res[n_out] = { tau };

	F_generic(Recorder_arg, Recorder_res);
/*
    for (int i = 0; i < NR; ++i)
        std::cout << Recorder_res[0][i] << std::endl;*/

	return 0;
}
