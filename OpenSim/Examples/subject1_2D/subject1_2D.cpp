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
#include <OpenSim/Common/MultiplierFunction.h>
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
constexpr int ndof = 12;         // degrees of freedom
constexpr int NX = 2 * ndof;    // states
constexpr int NU = ndof;        // controls
constexpr int NR = 15;          // residual forces + segment origin + GRFs + GRMs + PoA

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

int F_generic(const double** arg, double** res) {
	Model* model = new Model("subject1_scaled_contacts.osim");

 //   // OpenSim model: create components
 //   /// Model
 //   OpenSim::Model* model;
 //   /// Bodies
 //   OpenSim::Body* pelvis;
 //   OpenSim::Body* femur_r;
 //   OpenSim::Body* femur_l;
 //   OpenSim::Body* tibia_r;
 //   OpenSim::Body* tibia_l;
 //   OpenSim::Body* talus_r;
 //   OpenSim::Body* talus_l;
 //   OpenSim::Body* calcn_r;
 //   OpenSim::Body* calcn_l;
 //   OpenSim::Body* toes_r;
 //   OpenSim::Body* toes_l;
 //   OpenSim::Body* torso;
 //   /// Joints
 //   OpenSim::CustomJoint* ground_pelvis;
 //   OpenSim::PinJoint* hip_r;
 //   OpenSim::PinJoint* hip_l;
 //   OpenSim::CustomJoint* knee_r;
 //   OpenSim::CustomJoint* knee_l;
 //   OpenSim::PinJoint* ankle_r;
 //   OpenSim::PinJoint* ankle_l;
 //   OpenSim::WeldJoint* subtalar_r;
 //   OpenSim::WeldJoint* subtalar_l;
 //   OpenSim::PinJoint* mtp_r;
 //   OpenSim::PinJoint* mtp_l;
 //   OpenSim::PinJoint* back;

 //   // OpenSim model: initialize components
 //   /// Model
 //   model = new OpenSim::Model();
 //   /// Body specifications
 //   pelvis = new OpenSim::Body("pelvis", 9.7143336091724048, Vec3(-0.068431260352111167, 0, 0), Inertia(0.084795236055270728, 0.071844990860059146, 0.047759184509729331, 0, 0, 0));
 //   femur_l = new OpenSim::Body("femur_l", 7.6723191502382786, Vec3(0, -0.16734021487797751, 0), Inertia(0.10701920433197534, 0.028053577834595483, 0.11285370912378581, 0, 0, 0));
 //   femur_r = new OpenSim::Body("femur_r", 7.6723191502382786, Vec3(0, -0.16734021487797751, 0), Inertia(0.10701920433197534, 0.028053577834595483, 0.11285370912378581, 0, 0, 0));
 //   tibia_l = new OpenSim::Body("tibia_l", 3.0581550357482121, Vec3(0, -0.18635251395796504, 0), Inertia(0.041418155207686803, 0.004191122848396879, 0.041993407363349118, 0, 0, 0));
 //   tibia_r = new OpenSim::Body("tibia_r", 3.0581550357482121, Vec3(0, -0.18635251395796504, 0), Inertia(0.041418155207686803, 0.004191122848396879, 0.041993407363349118, 0, 0, 0));
 //   talus_l = new OpenSim::Body("talus_l", 0.082485638186061028, Vec3(0), Inertia(0.00076310987643225315, 0.00076310987643225315, 0.00076310987643225315, 0, 0, 0));
 //   talus_r = new OpenSim::Body("talus_r", 0.082485638186061028, Vec3(0), Inertia(0.00076310987643225315, 0.00076310987643225315, 0.00076310987643225315, 0, 0, 0));
 //   calcn_l = new OpenSim::Body("calcn_l", 1.0310704773257626, Vec3(0.096184339663296911, 0.028855301898989071, 0), Inertia(0.0010683538270051544, 0.0029761285180857871, 0.003128750493372238, 0, 0, 0));
 //   calcn_r = new OpenSim::Body("calcn_r", 1.0310704773257626, Vec3(0.096184339663296911, 0.028855301898989071, 0), Inertia(0.0010683538270051544, 0.0029761285180857871, 0.003128750493372238, 0, 0, 0));
 //   toes_l = new OpenSim::Body("toes_l", 0.17866389231100815, Vec3(0.03327978152350073, 0.0057710603797978145, 0.016832259441076958), Inertia(7.631098764322532e-05, 0.00015262197528645064, 7.631098764322532e-05, 0, 0, 0));
 //   toes_r = new OpenSim::Body("toes_r", 0.17866389231100815, Vec3(0.03327978152350073, 0.0057710603797978145, -0.016832259441076958), Inertia(7.631098764322532e-05, 0.00015262197528645064, 7.631098764322532e-05, 0, 0, 0));
 //   torso = new OpenSim::Body("torso", 28.240278003208967, Vec3(-0.028926628525507352, 0.308550704272078427, 0), Inertia(1.1307751170491251, 0.57938324918997219, 1.0977222804639659, 0, 0, 0));
 //   /// Joints
 //   /// Ground_pelvis transform
 //   SpatialTransform st_ground_pelvis;
 //   st_ground_pelvis[0].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tilt", 1, 1));
 //   st_ground_pelvis[0].setFunction(new LinearFunction());
 //   st_ground_pelvis[0].setAxis(Vec3(0,0,1));
 //   st_ground_pelvis[1].setFunction(new Constant(0));
 //   st_ground_pelvis[1].setAxis(Vec3(1,0,0));
 //   st_ground_pelvis[2].setFunction(new Constant(0));
 //   st_ground_pelvis[2].setAxis(Vec3(0,1,0));
 //   st_ground_pelvis[3].setCoordinateNames(OpenSim::Array<std::string>("pelvis_tx", 1, 1));
 //   st_ground_pelvis[3].setFunction(new LinearFunction());
 //   st_ground_pelvis[3].setAxis(Vec3(1,0,0));
 //   st_ground_pelvis[4].setCoordinateNames(OpenSim::Array<std::string>("pelvis_ty", 1, 1));
 //   st_ground_pelvis[4].setFunction(new LinearFunction());
 //   st_ground_pelvis[4].setAxis(Vec3(0,1,0));
 //   st_ground_pelvis[5].setFunction(new LinearFunction());
 //   st_ground_pelvis[5].setAxis(Vec3(0,0,1));
 //   /// Knee_l transform
 //   SpatialTransform st_knee_l;
 //   st_knee_l[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
 //   st_knee_l[0].setFunction(new LinearFunction());
 //   st_knee_l[0].setAxis(Vec3(0,0,1));
 //   st_knee_l[1].setFunction(new Constant(0));
 //   st_knee_l[1].setAxis(Vec3(0,1,0));
 //   st_knee_l[2].setFunction(new Constant(0));
 //   st_knee_l[2].setAxis(Vec3(1,0,0));
 //   st_knee_l[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
 //   double knee_X_r_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
 //   double knee_X_r_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
 //   OpenSim::SimmSpline* knee_X_r = new SimmSpline(12, knee_X_r_x, knee_X_r_y, "function_X");
 //   st_knee_l[3].setFunction(new MultiplierFunction(knee_X_r, 0.9843542051645735));
 //   st_knee_l[3].setAxis(Vec3(1,0,0));
 //   st_knee_l[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_l", 1, 1));
 //   double knee_Y_r_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
 //   double knee_Y_r_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
 //   OpenSim::SimmSpline* knee_Y_r = new SimmSpline(7, knee_Y_r_x, knee_Y_r_y, "function_Y");
 //   st_knee_l[4].setFunction(new MultiplierFunction(knee_Y_r, 0.9843542051645735));
 //   st_knee_l[4].setAxis(Vec3(0,1,0));
 //   st_knee_l[5].setFunction(new Constant(0));
 //   st_knee_l[5].setAxis(Vec3(0,0,1));
 //   /// Knee_r transform
 //   SpatialTransform st_knee_r;
 //   st_knee_r[0].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
 //   st_knee_r[0].setFunction(new LinearFunction());
 //   st_knee_r[0].setAxis(Vec3(0,0,1));
 //   st_knee_r[1].setFunction(new Constant(0));
 //   st_knee_r[1].setAxis(Vec3(0,1,0));
 //   st_knee_r[2].setFunction(new Constant(0));
 //   st_knee_r[2].setAxis(Vec3(1,0,0));
 //   st_knee_r[3].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
 //   double knee_X_l_x[] = { -2.0944, -1.74533, -1.39626, -1.0472, -0.698132, -0.349066, -0.174533, 0.197344, 0.337395, 0.490178, 1.52146, 2.0944 };
 //   double knee_X_l_y[] = { -0.0032, 0.00179, 0.00411, 0.0041, 0.00212, -0.001, -0.0031, -0.005227, -0.005435, -0.005574, -0.005435, -0.00525 };
 //   OpenSim::SimmSpline* knee_X_l = new SimmSpline(12, knee_X_l_x, knee_X_l_y, "function_X");
 //   st_knee_r[3].setFunction(new MultiplierFunction(knee_X_l, 0.9843542051645735));
 //   st_knee_r[3].setAxis(Vec3(1,0,0));
 //   st_knee_r[4].setCoordinateNames(OpenSim::Array<std::string>("knee_angle_r", 1, 1));
 //   double knee_Y_l_x[] = { -2.0944, -1.22173, -0.523599, -0.349066, -0.174533, 0.159149, 2.0944 };
 //   double knee_Y_l_y[] = { -0.4226, -0.4082, -0.399, -0.3976, -0.3966, -0.395264, -0.396 };
 //   OpenSim::SimmSpline* knee_Y_l = new SimmSpline(7, knee_Y_l_x, knee_Y_l_y, "function_Y");
 //   st_knee_r[4].setFunction(new MultiplierFunction(knee_Y_l, 0.9843542051645735));
 //   st_knee_r[4].setAxis(Vec3(0,1,0));
 //   st_knee_r[5].setFunction(new Constant(0));
 //   st_knee_r[5].setAxis(Vec3(0,0,1));
 //   /// Joint specifications
 //   ground_pelvis = new CustomJoint("ground_pelvis", model->getGround(), Vec3(0), Vec3(0), *pelvis, Vec3(0), Vec3(0), st_ground_pelvis);
 //   hip_l = new PinJoint("hip_l", *pelvis, Vec3(-0.068431260352111167, -0.063978872832737607, -0.082439592243860341), Vec3(0), *femur_l, Vec3(0), Vec3(0));
 //   hip_r = new PinJoint("hip_r", *pelvis, Vec3(-0.068431260352111167, -0.063978872832737607, 0.082439592243860341), Vec3(0), *femur_r, Vec3(0), Vec3(0));
 //   knee_l = new CustomJoint("knee_l", *femur_l, Vec3(0), Vec3(0), *tibia_l, Vec3(0), Vec3(0), st_knee_l);
 //   knee_r = new CustomJoint("knee_r", *femur_r, Vec3(0), Vec3(0), *tibia_r, Vec3(0), Vec3(0), st_knee_r);
 //   ankle_l = new PinJoint("ankle_l", *tibia_l, Vec3(0, -0.42919968399531311, 0), Vec3(0), *talus_l, Vec3(0), Vec3(0));
 //   ankle_r = new PinJoint("ankle_r", *tibia_r, Vec3(0, -0.42919968399531311, 0), Vec3(0), *talus_r, Vec3(0), Vec3(0));
 //   subtalar_l = new WeldJoint("subtalar_l", *talus_l, Vec3(-0.046909102453789903, -0.040349330488753055, -0.0076177997013331146), Vec3(1.7681899999999999, -0.906223, 1.8196000000000001), *calcn_l, Vec3(0), Vec3(1.7681899999999999, -0.906223, 1.8196000000000001));
 //   subtalar_r = new WeldJoint("subtalar_r", *talus_r, Vec3(-0.046909102453789903, -0.040349330488753055, 0.0076177997013331146), Vec3(-1.7681899999999999, 0.906223, 1.8196000000000001), *calcn_r, Vec3(0), Vec3(-1.7681899999999999, 0.906223, 1.8196000000000001));
 //   mtp_l = new PinJoint("mtp_l", *calcn_l, Vec3(0.17197759931797485, -0.0019236867932659382, -0.0010387908683636067), Vec3(-3.1415899999999999, -0.61990100000000004, 0), *toes_l, Vec3(0), Vec3(-3.1415899999999999, -0.61990100000000004, 0));
 //   mtp_r = new PinJoint("mtp_r", *calcn_r, Vec3(0.17197759931797485, -0.0019236867932659382, 0.0010387908683636067), Vec3(-3.1415899999999999, 0.61990100000000004, 0), *toes_r, Vec3(0), Vec3(-3.1415899999999999, 0.61990100000000004, 0));
 //   back = new PinJoint("back", *pelvis, Vec3(-0.097468570261069226, 0.078884691919336072, 0), Vec3(0), *torso, Vec3(0), Vec3(0));
 //   /// Add bodies and joints to model
 //   model->addBody(pelvis);		    model->addJoint(ground_pelvis);
 //   model->addBody(femur_l);		model->addJoint(hip_l);
 //   model->addBody(femur_r);		model->addJoint(hip_r);
 //   model->addBody(tibia_l);		model->addJoint(knee_l);
 //   model->addBody(tibia_r);		model->addJoint(knee_r);
 //   model->addBody(talus_l);		model->addJoint(ankle_l);
 //   model->addBody(talus_r);		model->addJoint(ankle_r);
 //   model->addBody(calcn_l);		model->addJoint(subtalar_l);
 //   model->addBody(calcn_r);		model->addJoint(subtalar_r);
 //   model->addBody(toes_l);		    model->addJoint(mtp_l);
 //   model->addBody(toes_r);		    model->addJoint(mtp_r);
 //   model->addBody(torso);          model->addJoint(back);
 //   /// Contact elements
 //   /// Parameters
 //   double radiusSphere_heel = 0.035;
 //   double radiussphere_front = 0.015;
 //   double stiffness = 3000000;
 //   double dissipation = 2.0;
 //   double staticFriction = 0.8;
 //   double dynamicFriction = 0.8;
 //   double viscousFriction = 0.5;
 //   double transitionVelocity = 0.2;
 //   Vec3 locSphere_heel_l = Vec3(0.01, 0, 0);
 //   Vec3 locsphere_front_l = Vec3(0.025, -0.01, -0.005);
 //   Vec3 locSphere_heel_r = Vec3(0.01, 0, 0);
 //   Vec3 locsphere_front_r = Vec3(0.025, -0.01, 0.005);
 //   Vec3 halfSpaceLocation(0);
	//Vec3 halfSpaceOrientation(0, 0, -0.5 * SimTK::Pi);
 //   /// Right foot contact shere specifications
 //   OpenSim::ContactSphere* sphere_heel_l;
 //   OpenSim::ContactSphere* sphere_heel_r;
 //   OpenSim::ContactSphere* sphere_front_l;
 //   OpenSim::ContactSphere* sphere_front_r;
 //   OpenSim::ContactHalfSpace* contactHalfSpace;
 //   OpenSim::SmoothSphereHalfSpaceForce* contact_heel_l;
	//OpenSim::SmoothSphereHalfSpaceForce* contact_heel_r;
 //   OpenSim::SmoothSphereHalfSpaceForce* contact_front_l;
	//OpenSim::SmoothSphereHalfSpaceForce* contact_front_r;

 //   sphere_heel_l = new OpenSim::ContactSphere(radiusSphere_heel,locSphere_heel_l,*calcn_l,"sphere_heel_l");
 //   model->addComponent(sphere_heel_l);
 //   sphere_heel_r = new OpenSim::ContactSphere(radiusSphere_heel,locSphere_heel_r,*calcn_r,"sphere_heel_r");
 //   model->addComponent(sphere_heel_r);
 //   sphere_front_l = new OpenSim::ContactSphere(radiussphere_front,locsphere_front_l,*toes_l,"sphere_front_l");
 //   model->addComponent(sphere_front_l);
 //   sphere_front_r = new OpenSim::ContactSphere(radiussphere_front,locsphere_front_r,*toes_r,"sphere_front_r");
 //   model->addComponent(sphere_front_r);
 //   contactHalfSpace = new OpenSim::ContactHalfSpace(halfSpaceLocation,halfSpaceOrientation,model->getGround(),"contactHalfSpace");
 //   model->addComponent(contactHalfSpace);
 //   contact_heel_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_heel_l",*sphere_heel_l,*contactHalfSpace);
 //   contact_heel_l->set_stiffness(stiffness);
 //   contact_heel_l->set_dissipation(dissipation);
 //   contact_heel_l->set_static_friction(staticFriction);
 //   contact_heel_l->set_dynamic_friction(dynamicFriction);
 //   contact_heel_l->set_viscous_friction(viscousFriction);
 //   contact_heel_l->set_transition_velocity(transitionVelocity);
 //   contact_heel_l->connectSocket_half_space(*contactHalfSpace);
 //   contact_heel_l->connectSocket_sphere(*sphere_heel_l);
 //   model->addComponent(contact_heel_l);

 //   contact_heel_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_heel_r",*sphere_heel_r,*contactHalfSpace);
 //   contact_heel_r->set_stiffness(stiffness);
 //   contact_heel_r->set_dissipation(dissipation);
 //   contact_heel_r->set_static_friction(staticFriction);
 //   contact_heel_r->set_dynamic_friction(dynamicFriction);
 //   contact_heel_r->set_viscous_friction(viscousFriction);
 //   contact_heel_r->set_transition_velocity(transitionVelocity);
 //   contact_heel_r->connectSocket_half_space(*contactHalfSpace);
 //   contact_heel_r->connectSocket_sphere(*sphere_heel_r);
 //   model->addComponent(contact_heel_r);

 //   contact_front_l = new OpenSim::SmoothSphereHalfSpaceForce("contact_front_l",*sphere_front_l,*contactHalfSpace);
 //   contact_front_l->set_stiffness(stiffness);
 //   contact_front_l->set_dissipation(dissipation);
 //   contact_front_l->set_static_friction(staticFriction);
 //   contact_front_l->set_dynamic_friction(dynamicFriction);
 //   contact_front_l->set_viscous_friction(viscousFriction);
 //   contact_front_l->set_transition_velocity(transitionVelocity);
 //   contact_front_l->connectSocket_half_space(*contactHalfSpace);
 //   contact_front_l->connectSocket_sphere(*sphere_front_l);
 //   model->addComponent(contact_front_l);

 //   contact_front_r = new OpenSim::SmoothSphereHalfSpaceForce("contact_front_r",*sphere_front_r,*contactHalfSpace);
 //   contact_front_r->set_stiffness(stiffness);
 //   contact_front_r->set_dissipation(dissipation);
 //   contact_front_r->set_static_friction(staticFriction);
 //   contact_front_r->set_dynamic_friction(dynamicFriction);
 //   contact_front_r->set_viscous_friction(viscousFriction);
 //   contact_front_r->set_transition_velocity(transitionVelocity);
 //   contact_front_r->connectSocket_half_space(*contactHalfSpace);
 //   contact_front_r->connectSocket_sphere(*sphere_front_r);
 //   model->addComponent(contact_front_r);

	/// Initialize  system and state.
    SimTK::State* state;
	state = new State(model->initSystem());
	// Read inputs
    std::vector<double> x(arg[0], arg[0] + NX);
    std::vector<double> u(arg[1], arg[1] + NU);

    // States and controls
    double ua[NU+2]; /// joint accelerations (Qdotdots) - controls
    Vector QsUs(NX); /// joint positions (Qs) and velocities (Us) - states

    // Assign inputs to model variables
    /// States
    QsUs.setToZero();
    for (int i = 0; i < NX; ++i) QsUs[i] = x[i];
    /// Controls
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesOSInSimbody = getIndicesOSInSimbody(*model);
    for (int i = 0; i < NU; ++i) ua[i] = u[indicesOSInSimbody[i]];

    // Set state variables and realize
    model->setStateVariableValues(*state, QsUs);
    model->realizeVelocity(*state);

    //InverseDynamicsSolver* IDsolver;
    //IDsolver = new InverseDynamicsSolver(*model);
    //Vector residualMobilityForces = IDsolver->solve(*state);

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
    Array<double> Force_values_heel_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_heel_r").getRecordValues(*state);
    Array<double> Force_values_front_r = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_front_r").getRecordValues(*state);
    /*Array<double> Force_values_heel_r = contact_heel_r->getRecordValues(*state);
    Array<double> Force_values_front_r = contact_front_r->getRecordValues(*state);*/
    SpatialVec GRF_heel_r;
    GRF_heel_r[0] = Vec3(Force_values_heel_r[3], Force_values_heel_r[4], Force_values_heel_r[5]);
    GRF_heel_r[1] = Vec3(Force_values_heel_r[0], Force_values_heel_r[1], Force_values_heel_r[2]);
    SpatialVec GRF_front_r;
    GRF_front_r[0] = Vec3(Force_values_front_r[3], Force_values_front_r[4], Force_values_front_r[5]);
    GRF_front_r[1] = Vec3(Force_values_front_r[0], Force_values_front_r[1], Force_values_front_r[2]);
    int nCalc_r = model->getBodySet().get("calcn_r").getMobilizedBodyIndex();
    int nToe_r = model->getBodySet().get("toes_r").getMobilizedBodyIndex();
    appliedBodyForces[nCalc_r] = appliedBodyForces[nCalc_r] + GRF_heel_r;
    appliedBodyForces[nToe_r] = appliedBodyForces[nToe_r] + GRF_front_r;

    /// Left foot
    Array<double> Force_values_heel_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_heel_l").getRecordValues(*state);
    Array<double> Force_values_front_l = model->getComponent<OpenSim::SmoothSphereHalfSpaceForce>("forceset/contact_front_l").getRecordValues(*state);
    /*Array<double> Force_values_heel_l = contact_heel_l->getRecordValues(*state);
    Array<double> Force_values_front_l = contact_front_l->getRecordValues(*state);*/
    SpatialVec GRF_heel_l;
    GRF_heel_l[0] = Vec3(Force_values_heel_l[3], Force_values_heel_l[4], Force_values_heel_l[5]);
    GRF_heel_l[1] = Vec3(Force_values_heel_l[0], Force_values_heel_l[1], Force_values_heel_l[2]);
    SpatialVec GRF_front_l;
    GRF_front_l[0] = Vec3(Force_values_front_l[3], Force_values_front_l[4], Force_values_front_l[5]);
    GRF_front_l[1] = Vec3(Force_values_front_l[0], Force_values_front_l[1], Force_values_front_l[2]);
    int nCalc_l = model->getBodySet().get("calcn_l").getMobilizedBodyIndex();
    int nToe_l = model->getBodySet().get("toes_l").getMobilizedBodyIndex();
    appliedBodyForces[nCalc_l] = appliedBodyForces[nCalc_l] + GRF_heel_l;
    appliedBodyForces[nToe_l] = appliedBodyForces[nToe_l] + GRF_front_l;
    /// knownUdot
    Vector knownUdot(ndof);
    knownUdot.setToZero();
    for (int i = 0; i < ndof; ++i) knownUdot[i] = ua[i];
    /// Calculate residual forces
    Vector residualMobilityForces(ndof);
    residualMobilityForces.setToZero();
    model->getMatterSubsystem().calcResidualForceIgnoringConstraints(
        *state, appliedMobilityForces, appliedBodyForces, knownUdot,
        residualMobilityForces);

    // Extract results
    /// Residual forces
    /// OpenSim and Simbody have different state orders so we need to adjust
    auto indicesSimbodyInOS = getIndicesSimbodyInOS(*model);
    for (int i = 0; i < NU; ++i) res[0][i] = (residualMobilityForces[indicesSimbodyInOS[i]]);
    for (int i = 0; i < 3; ++i) res[0][i + NU] = GRF_heel_r[1][i];
	return 0;
}

int main() {
    double x[NX];
	double u[NU];
	double tau[NR];

	for (int i = 0; i < NX; ++i) x[i] = 0;
	for (int i = 0; i < NU; ++i) u[i] = 0;

	const double* Recorder_arg[n_in] = { x,u };
	double* Recorder_res[n_out] = { tau };

	F_generic(Recorder_arg, Recorder_res);

    for (int i = 0; i < NR; ++i)
        std::cout << Recorder_res[0][i] << std::endl;

	return 0;
}
