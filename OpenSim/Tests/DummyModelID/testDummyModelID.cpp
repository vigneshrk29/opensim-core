/* -------------------------------------------------------------------------- *
 *                            OpenSim:  testID.cpp                            *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2017 Stanford University and the Authors                *
 * Author(s): Ayman Habib, Ajay Seth                                          *
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
// INCLUDE
#include <string>
#include <iostream>
#include <OpenSim/version.h>
#include <OpenSim/Common/IO.h>
#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/ForceSet.h>
#include <OpenSim/Tools/InverseDynamicsTool.h>
#include <OpenSim/Auxiliary/auxiliaryTestFunctions.h>

using namespace OpenSim;
using namespace std;

int main()
{
    try {
        InverseDynamicsTool id1("Setup_Standing_GRF_no_GRM_source.xml");
        id1.run();
        //Storage result1("ID_Standing_GRF_no_GRM_source.sto"), standard1("ID_Standing_GRF_no_GRM_GUI.sto");
        //CHECK_STORAGE_AGAINST_STANDARD( result1, standard1,
        //    std::vector<double>(23, 1e-2), __FILE__, __LINE__,
        //    "testNoGRM failed");
        //cout << "testNoGRM passed" << endl;

        InverseDynamicsTool id2("Setup_Standing_GRF_with_GRM_source.xml");
        id2.run();
        //Storage result2("ID_Standing_GRF_with_GRM_source.sto"), standard2("ID_Standing_GRF_with_GRM_GUI.sto");
        //CHECK_STORAGE_AGAINST_STANDARD(result2, standard2,
        //    std::vector<double>(23, 2.0), __FILE__, __LINE__,
        //    "testWithGRM failed");
        //cout << "testWithGRM passed" << endl;
    }
    catch (const Exception& e) {
        e.print(cerr);
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
