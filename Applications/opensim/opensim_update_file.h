#ifndef OPENSIM_UPDATE_FILE_H_
#define OPENSIM_UPDATE_FILE_H_
/* -------------------------------------------------------------------------- *
 *                       OpenSim:  opensim_update_file.h                      *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2016 Stanford University and the Authors                *
 * Author(s): Frank C. Anderson, Ayman Habib, Chris Dembia                    *
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

static const char HELP_UPDATE_FILE[] =
R"(Update an .osim, .xml (e.g., setup) or .sto file to this version's format.

Usage:
  opensim [options]... update-file <input-file> <output-file>
  opensim update-file -h | --help

Options:
  -L <path>, --library <path>  Load a plugin.

Description:
  In an OpenSim XML file, the XML file version appears as
  the "Version" attribute the "OpenSimDocument" element.

Examples:
  opensim update-file lowerlimb_v3.3.osim lowerlimb_updated.osim
  opensim update-file data_v3.3.sto data_updated.sto
)";

int update_file(int argc, const char** argv) {

    using namespace OpenSim;

    std::map<std::string, docopt::value> args = docopt::docopt(
            HELP_UPDATE_FILE, { argv + 1, argv + argc },
            true); // show help if requested

    const std::string inputFile = args["<input-file>"].asString();
    const std::string outputFile = args["<output-file>"].asString();

    // Grab the file extension.
    std::string::size_type extSep = inputFile.rfind(".");
    if (extSep == std::string::npos) {
        throw Exception("Input file '" + inputFile +
                        "' does not have an extension.");
    }
    std::string extension = inputFile.substr(extSep);

    // .osim or .xml file.
    if (extension == ".osim" || extension == ".xml") {
        std::cout << "Loading input file '" << inputFile << "'." << std::endl;
        const auto* obj = Object::makeObjectFromFile(inputFile);
        if (!obj) {
            throw Exception(
                    "Could not make object from file '" + inputFile + "'.\n"
                    "Did you intend to load a plugin (with --library)?");
        }
        std::cout << "Printing updated file to '" << outputFile << "'."
                  << std::endl;
        obj->print(outputFile);
        return EXIT_SUCCESS;
    }

    // .sto file.
    if (extension == ".sto") {
        std::cout << "Loading input file '" << inputFile << "'." << std::endl;
        Storage stg(inputFile);
        std::cout << "Printing updated file to '" << outputFile << "'."
                  << std::endl;
        stg.print(outputFile);
        return EXIT_SUCCESS;
    }

    throw Exception(
            "Input file '" + inputFile + "' has an unrecognized extension.");
}


#endif // OPENSIM_UPDATE_FILE_H_
