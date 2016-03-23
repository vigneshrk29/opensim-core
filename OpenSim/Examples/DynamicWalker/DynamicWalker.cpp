/* -------------------------------------------------------------------------- *
 *                    OpenSim:  DynamicWalker.cpp                     *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2012 Stanford University and the Authors                *
 * Author(s): Jeffrey A. Reinbolt, Ayman Habib, Ajay Seth, Jack Middleton,    *
 *            Samuel R. Hamner                                                *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/* 
 *  Below is an example of an OpenSim application that provides its own 
 *  main() routine.  This application is a forward simulation of tug-of-war between two
 *  muscles pulling on a block.
 */

// Author:  Jeff Reinbolt, Ayman Habib, Ajay Seth, Jack Middleton, Samuel Hamner

//==============================================================================
//==============================================================================
#include <OpenSim/OpenSim.h>

using namespace OpenSim;
using namespace SimTK;

int main()
{
    try {
        // Create an OpenSim model and set its name
        Model osimModel;
        // Section: Setup
        // Define key model variables
        double pelvisWidth = 0.20, thighLength = 0.40, shankLength = 0.435;

        osimModel.setName("DynamicWalkerModel");

        // 4.0 API change getGround
        // Get a reference to the ground object
        const Ground& ground = osimModel.getGround();

        // Define the acceleration of gravity
        osimModel.setGravity(Vec3(0, -9.80665, 0));
        // Section: Create the Platform
        double mass = 1;

        // Location of the Center of Mass from the Body Origin expressed in Body Frame
        Vec3 comLocInBody(0.0, 0.0, 0.0);

        // Inertia of the Body expressed in the Body Frame
        Inertia bodyInertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        // Create the body
        OpenSim::Body* platform = new OpenSim::Body("Platform", mass, comLocInBody, bodyInertia);

        // Section: Create the Platform Joint
        // Create the joint connection the platform to the ground
        Vec3 locationInParent(0.0, 0.0, 0.0);
        Vec3 orientationInParent(0.0, 0.0, 0.0);
        Vec3 locationInChild(0.0, 0.0, 0.0);
        Vec3 orientationInChild(0.0, 0.0, 0.0);
        PinJoint *platformToGround = new PinJoint("PlatformToGround", ground, locationInParent, orientationInParent, *platform, locationInChild, orientationInChild);

        // Section: Set the properties of the coordinates that define the joint
        // A pin joint consists of a single coordinate describing a change in orientation about the Z axis
        CoordinateSet &platformJoints = platformToGround->upd_CoordinateSet();
        platformJoints[0].setName("platform_rz");
        double rotRangePlatform[2] = { -Pi / 2.0, 0 };
        platformJoints[0].setRange(rotRangePlatform);
        platformJoints[0].setDefaultValue(convertDegreesToRadians(-10.0));
        platformJoints[0].setDefaultLocked(true);
        // Add and scale model for display in GUI
        // 4.0 API change to attach Geometry to frame
        platform->attachMeshGeometry("box.vtp", Vec3(1, 0.05, 1));

        // Add the platform to the Model
        osimModel.addBody(platform);

        // Section: Create the Pelvis
        mass = 1;

        // Location of the Center of Mass from the Body Origin expressed in Body Frame
        comLocInBody = Vec3(0.0, 0.0, 0.0);

        // Inertia of the Body expressed in the Body Frame
        bodyInertia = Inertia(1.0, 1.0, 1.0, 0.0, 0.0, 0.0);

        // Create the body
        OpenSim::Body* pelvis = new OpenSim::Body("Pelvis", mass, comLocInBody, bodyInertia);

        // Create the joint which connects the Pelvis to the Platform
        locationInParent = Vec3(0.0, 0.0, 0.0);
        orientationInParent = Vec3(0.0, 0.0, 0.0);
        locationInChild = Vec3(0.0, 0.0, 0.0);
        orientationInChild = Vec3(0.0, 0.0, 0.0);
        FreeJoint *pelvisToPlatform = new FreeJoint("PelvisToPlatform", *platform, locationInParent, orientationInParent, *pelvis, locationInChild, orientationInChild, false);

        // A Free joint has six coordinates: (in order) rot_x, rot_y, rot_z, trans_x, trans_y, trans_z
        // Set the properties of the coordinates that define the joint
        CoordinateSet &pelvisJointCoords = pelvisToPlatform->upd_CoordinateSet();
        // Add and scale model for display in GUI
        pelvis->attachMeshGeometry("sphere.vtp", Vec3(pelvisWidth / 2.0, pelvisWidth / 2.0, pelvisWidth));

        // Add the joint
        // 4.0 API now Joint need to be added separately
        osimModel.addJoint(pelvisToPlatform);

        // Add the pelvis to the Model
        osimModel.addBody(pelvis);

        // Section: Add Contact Geometry
        // Add Contact Mesh for Platform
        // The initial orientation is defined by a normal pointing in the positive x direction
        ContactHalfSpace* platformContact = new ContactHalfSpace(Vec3(0.0, 0.0, 0.0), Vec3(0.0, 0.0, -Pi / 2.0), *platform, "PlatformContact");
        osimModel.addContactGeometry(platformContact);

        // Contact Sphere Properties
        double contactSphereRadius = 0.05;

        // Add Contact Sphere for Right Hip
        Vec3 rightHipLocationInPelvis(0.0, 0.0, pelvisWidth / 2.0);
        ContactSphere* rightHipContact = new ContactSphere(contactSphereRadius, rightHipLocationInPelvis, *pelvis, "RHipContact");
        osimModel.addContactGeometry(rightHipContact);

        // Section: Add HuntCrossleyForces
        // Define contact parameters for all the spheres
        double stiffness = 1E7, dissipation = 0.1, staticFriction = 0.6, dynamicFriction = 0.4, viscosity = 0.01;

        // Right Foot Contact Parameters
        OpenSim::HuntCrossleyForce::ContactParameters *rightHipContactParams = new OpenSim::HuntCrossleyForce::ContactParameters(stiffness, dissipation, staticFriction, dynamicFriction, viscosity);
        rightHipContactParams->addGeometry("RHipContact");
        rightHipContactParams->addGeometry("PlatformContact");

        // Right Foot Force
        OpenSim::HuntCrossleyForce* rightHipForce = new OpenSim::HuntCrossleyForce(rightHipContactParams);
        rightHipForce->setName("RightHipForce");

        //Add Force
        osimModel.addForce(rightHipForce);

        // Save the model to a file
        osimModel.print("DynamicWalkerModel.osim");

        Model roundTrip("DynamicWalkerModel.osim");
        roundTrip.initSystem();	// This crashes now due to Exception.

    }
    catch (OpenSim::Exception ex)
    {
        std::cout << ex.getMessage() << std::endl;
        return 1;
    }
    catch (std::exception ex)
    {
        std::cout << ex.what() << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cout << "UNRECOGNIZED EXCEPTION" << std::endl;
        return 1;
    }

    std::cout << "OpenSim example completed successfully.\n";
    std::cin.get();
    return 0;
}
