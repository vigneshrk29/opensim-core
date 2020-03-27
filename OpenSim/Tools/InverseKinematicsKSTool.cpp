/* -------------------------------------------------------------------------- *
 *                    OpenSim:  InverseKinematicsKSTool.cpp                   *
 * -------------------------------------------------------------------------- *
 * The OpenSim API is a toolkit for musculoskeletal modeling and simulation.  *
 * See http://opensim.stanford.edu and the NOTICE file for more information.  *
 * OpenSim is developed at Stanford University and supported by the US        *
 * National Institutes of Health (U54 GM072970, R24 HD065690) and by DARPA    *
 * through the Warrior Web program.                                           *
 *                                                                            *
 * Copyright (c) 2005-2020 Stanford University and the Authors                *
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

//=============================================================================
// INCLUDES
//=============================================================================
#include "InverseKinematicsKSTool.h"

#include "IKCoordinateTask.h"
#include "IKMarkerTask.h"
#include "IKTaskSet.h"

#include <OpenSim/Analyses/Kinematics.h>
#include <OpenSim/Common/Constant.h>
#include <OpenSim/Common/FunctionSet.h>
#include <OpenSim/Common/GCVSplineSet.h>
#include <OpenSim/Common/IO.h>
#include <OpenSim/Common/Stopwatch.h>
#include <OpenSim/Common/Storage.h>
#include <OpenSim/Common/XMLDocument.h>
#include <OpenSim/Simulation/InverseKinematicsSolver.h>
#include <OpenSim/Simulation/Model/Model.h>

using namespace OpenSim;
using namespace std;
using namespace SimTK;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Destructor.
 */
InverseKinematicsKSTool::~InverseKinematicsKSTool()
{
}
//_____________________________________________________________________________
/**
 * Default constructor.
 */
InverseKinematicsKSTool::InverseKinematicsKSTool() : Tool()
{
    constructProperties();
}
/**
 * Construct from file.
 *
 * The object is constructed from the root element of the XML document.
 * The type of object is the tag name of the XML root element.
 *
 * @param aFileName File name of the document.
 */
InverseKinematicsKSTool::InverseKinematicsKSTool(const string &aFileName, bool aLoadModel) :
    Tool(aFileName, true)
{
    constructProperties();
    updateFromXMLDocument();
}

//_____________________________________________________________________________
/**
 * Construct properties
 */
void InverseKinematicsKSTool::constructProperties()
{
    constructProperty_model_file("");
    constructProperty_constraint_weight(Infinity);
    constructProperty_accuracy(1e-5);
    constructProperty_IKTaskSet(IKTaskSet());
    constructProperty_marker_file("");
    constructProperty_coordinate_file("");
    Array<double> range{Infinity, 2};
    range[0] = -Infinity; // Make range -Infinity to Infinity unless limited by
                          // data
    constructProperty_time_range(range);

    constructProperty_report_errors(true);
    constructProperty_output_motion_file("");
    constructProperty_report_marker_locations(false);
}

// TODO
////=============================================================================
//// OPERATORS
////=============================================================================
////_____________________________________________________________________________
///**
// * Assignment operator.
// *
// * @return Reference to this object.
// */
//InverseKinematicsKSTool& InverseKinematicsKSTool::
//operator=(const InverseKinematicsKSTool &aTool)
//{
//	// BASE CLASS
//	Tool::operator=(aTool);
//
//	// MEMBER VARIABLES
//	_modelFileName = aTool._modelFileName;
//	_constraintWeight = aTool._constraintWeight;
//	_accuracy = aTool._accuracy;
//	_order = aTool._order;
//	_timeScale = aTool._timeScale;
//	_sdProcess = aTool._sdProcess;
//	_sdMeas = aTool._sdMeas;
//	_ikTaskSet = aTool._ikTaskSet;
//	_markerFileName = aTool._markerFileName;
//	_timeRange = aTool._timeRange;
//	_reportErrors = aTool._reportErrors;
//	_coordinateFileName = aTool._coordinateFileName;
//	_reportErrors = aTool._reportErrors;
//	_outputMotionFileName = aTool._outputMotionFileName;
//	_reportMarkerLocations = aTool._reportMarkerLocations;
//
//	return(*this);
//}

//=============================================================================
// GET AND SET
//=============================================================================


//=============================================================================
// RUN
//=============================================================================
//_____________________________________________________________________________
/**
 * Run the inverse kinematics (Kalman smoothing) tool.
 */
bool InverseKinematicsKSTool::run()
{
	bool success = false;
    bool modelFromFile=true;
    Kinematics* kinematicsReporter = nullptr;
    try{
        //Load and create the indicated model
        if (_model.empty()) {
            OPENSIM_THROW_IF_FRMOBJ(get_model_file().empty(), Exception,
                    "No model filename was provided.");
            _model.reset(new Model(get_model_file()));
        }
        else
            modelFromFile = false;

		// although newly loaded model will be finalized
        // there is no guarantee that the _model has not been edited/modified
        _model->finalizeFromProperties();
        _model->printBasicInfo();

		// Do the maneuver to change then restore working directory so that the
        // parsing code behaves properly if called from a different directory.
        string saveWorkingDirectory = IO::getCwd();
        string directoryOfSetupFile = IO::getParentDirectory(getDocumentFileName());
        IO::chDir(directoryOfSetupFile);

		// Define reporter for output
        kinematicsReporter = new Kinematics();
        kinematicsReporter->setRecordAccelerations(false);
        kinematicsReporter->setInDegrees(true);
        _model->addAnalysis(kinematicsReporter);

		cout<<"Running tool "<<getName()<<".\n";

        // Get the trial name to label data written to files
        string trialName = getName();

		_model->setUseVisualizer(true);
		_model->updForceSet().clearAndDestroy(); // Required if using Visualizer

		// Initialize the model's underlying system and get its default state.
        SimTK::State& s = _model->initSystem();

		SimTK::Vector qDefaults = s.getQ(); // TODO

		//Convert old Tasks to references for assembly and tracking
        MarkersReference markersReference;
        SimTK::Array_<CoordinateReference> coordinateReferences;
        // populate the references according to the setting of this Tool
        populateReferences(markersReference, coordinateReferences); // TODO: might replace piece of code below

        // Determine the start time, if the provided time range is not
        // specified then use time from marker reference.
        // Adjust the time range for the tool if the provided range exceeds
        // that of the marker data.
        SimTK::Vec2 markersValidTimeRange = markersReference.getValidTimeRange();
        double start_time = (markersValidTimeRange[0] > get_time_range(0)) ?
            markersValidTimeRange[0] : get_time_range(0);
        double final_time = (markersValidTimeRange[1] < get_time_range(1)) ?
            markersValidTimeRange[1] : get_time_range(1);

        SimTK_ASSERT2_ALWAYS(final_time >= start_time,
            "InverseKinematicsTool final time (%f) is before start time (%f).",
            final_time, start_time);

        const auto& markersTable = markersReference.getMarkerTable();
        const int start_ix = int(
            markersTable.getNearestRowIndexForTime(start_time) );
        const int final_ix = int(
            markersTable.getNearestRowIndexForTime(final_time) );
        const int Nframes = final_ix - start_ix + 1;
        const auto& times = markersTable.getIndependentColumn();

		// create the solver given the input data
        InverseKinematicsSolver ikSolver(*_model, markersReference,
            coordinateReferences, get_constraint_weight());
        ikSolver.setAccuracy(get_accuracy());
        s.updTime() = times[start_ix];
        ikSolver.assemble(s);
        kinematicsReporter->begin(s);

        AnalysisSet& analysisSet = _model->updAnalysisSet();
        analysisSet.begin(s);
        // Get the actual number of markers the Solver is using, which
        // can be fewer than the number of references if there isn't a
        // corresponding model marker for each reference.
        int nm = ikSolver.getNumMarkersInUse();
        SimTK::Array_<double> squaredMarkerErrors(nm, 0.0);
        SimTK::Array_<Vec3> markerLocations(nm, Vec3(0));
        SimTK::Array_<Vec3> measMarkerLocations(nm, Vec3(0));
        std::cout << "number of markers " << nm << std::endl;
		//std::cout << "Marker not used :" << mNotUsed <<std::endl; TODO

		Stopwatch watch;
		double dt = 1.0/markersReference.getSamplingFrequency();
		int Nframes = int((final_time-start_time)/dt)+1;

        // THIS IS WHERE I STOPPED FOR NOW
		// Set up process model
		Matrix F(_order, _order); F=0.0;
		Matrix Q(_order, _order); Q=0.0;
		Vector G(_order); G=0.0;
		for (int j=0; j<_order; ++j) F(j,j) = 1;
		double Ts = dt*_timeScale;
		int fac = 1;
		double diag;
		for (int n=1; n<_order; ++n){
			fac = fac*n;
			diag = pow(Ts,n)/fac;
			for (int j=0; j<_order-n; ++j) F(j,j+n) = diag;
			G(_order-n) = diag;
		}
		fac = fac*_order;
		G(0) = pow(Ts,_order)/fac;
		Q = pow(_sdProcess,2)*G*~G;

		// Find locked joints
		int nq = s.getNQ();
		int nu = s.getNU();
		int ns = nq*_order;

		CoordinateSet& modelCoordinateSet = _model->updCoordinateSet();
		Array<int> iLocked;
		Array<int> iClamped;
		Array<double> lowerClampbounds;
		Array<double> upperClampbounds;
		int nl = 0; // number of locked coordinates
		int k = 0; // Index in the variable x
		// Loop over coordinates
		int n_coordinates = modelCoordinateSet.getSize();
		for (int i = 0; i < n_coordinates; i++) {
			Coordinate& aCoord = modelCoordinateSet[i];
			cout << aCoord.getName() << endl;
			if (aCoord.getLocked(s)) {
				// int iqx = aCoord.getMobilizerQIndex();
				// int iqx = aCoord.getStateVariableSystemIndex(aCoord.getName()); // Did not work anymore.
				iLocked.append(i);
				cout << aCoord.getName() << " (" << i <<") is locked" << endl;
				nl = nl + 1;
			}
			else {
				if (aCoord.getClamped(s)) {
					iClamped.append(k);
					lowerClampbounds.append(aCoord.getRangeMin());
					upperClampbounds.append(aCoord.getRangeMax());
					cout << aCoord.getName() << " (" << k <<") is clamped" << endl;
				}
				k += 1;
			}
		}
		int nClamped = iClamped.size();
		std::cout << "Number locked: " << iLocked.size() << std::endl;
		std::cout << "Number clamped: " << nClamped << std::endl;
		std::cout << "Lower: " << lowerClampbounds << std::endl;
		std::cout << "Upper: " << upperClampbounds << std::endl;

		// Initial state estimate based on IK (global optimization method) solution
		int Nframesinit = 3;
		Matrix qIK(nq,Nframesinit); qIK = 0;
		for (int i = 0; i < Nframesinit; i++) {
			s.updTime() = start_time + i*dt;
			ikSolver.track(s);
			qIK(i) = s.getQ();
		}

		const Vector &q = s.getQ();
		int nqf = nq-nl;		// number of free coordinates
		int nsf = nqf*_order;	// number of free states
		Vector x(nsf); x=0.0;
		int j = 0;
		for(int k=0; k<nq; ++k){
			if (iLocked.findIndex(k) == -1){
				x(j*_order) = qIK(k,0);
				x(j*_order+1) = (qIK(k,2) - qIK(k,0))/(2*Ts);
				x(j*_order+2) = (qIK(k,2) - 2*qIK(k,1) + qIK(k,0))/pow(Ts,2);
				j = j+1;
			}
		}




		s.updTime() = start_time;
		s.updQ() = qIK(0);
		_model->getMultibodySystem().realize(s, SimTK::Stage::Velocity);
		kinematicsReporter.step(s, 0);

		// Initial covariance estimate;
		Matrix P(nsf,nsf); P=0.0;
		for(int n=0; n<nsf; ++n) P(n,n) = 1;

		// Initialization
		Vector x_pred(nsf); x_pred=0.0;
		Matrix P_pred(nsf,nsf); P_pred=0.0;

		SimTK::Array_<Vector> x_pred_array(Nframes,Vector(nsf));
		SimTK::Array_<Matrix> P_pred_array(Nframes,Matrix(nsf,nsf));
		SimTK::Array_<Vector> x_kf_array(Nframes,Vector(nsf));
		SimTK::Array_<Matrix> P_kf_array(Nframes,Matrix(nsf,nsf));

		x_kf_array[0] = x;
		P_kf_array[0] = P;

		////////////
		// FILTER //
		////////////
		double tol_m = 0.000001;
		// Loop over time frames
		for (int n = 1; n < Nframes; n++) {
		// for (int n = 1; n < 2; n++) {

			// PREDICT
			for(int i=0; i<nqf; ++i){
				for(int j=0; j<nqf; ++j)	P_pred(i*_order, j*_order, _order, _order) = F*P(i*_order, j*_order, _order, _order)*~F;
			}

			for(int i=0; i<nqf; ++i){
				x_pred(i*_order,_order) = F*x(i*_order, _order);
				P_pred(i*_order, i*_order, _order, _order) = P_pred(i*_order, i*_order, _order, _order) + Q;
			}


			// Set new state
			s.updTime() = start_time + n*dt;
			Vector qpred(nq); qpred = 0;
			int j=0;
			for(int k=0; k<nq; ++k) {
				if (iLocked.findIndex(k) == -1){
					qpred(k) = x_pred(j*_order);
					j = j+1;
				}
				else { // Update
					qpred(k) = qDefaults(k);
				}
			}
			s.updQ() = qpred;

			ikSolver.setState(s); // Update: Fuction added to AssemblySolver
			// ikSolver.track(s); // failure
			_model->getMultibodySystem().realize(s, SimTK::Stage::Velocity);

			x_pred_array[n] = x_pred;
			P_pred_array[n] = P_pred;

			// UPDATE
			// innovation
			ikSolver.computeCurrentMarkerLocations(markerLocations);

			markersReference.getValues(s, measMarkerLocations);
			// get markers defined by the model
			const MarkerSet& modelMarkerSet = _model->getMarkerSet();

			Array<Vec3> innov;
			Vec3 diff(0);
			int il = 0;
			int index = 0;
			int index_model = 0; //version 2.1
			bool markerNotMeasured = false;
			// Array<double> measlocations(0.0, 3*nm);
			Array<std::string> mNotUsedInStep = mNotUsed;
			for(unsigned int i=0; i < markerNames.size(); i++){ // Loop over all the markers in the marker-file
			// 1. Check if we have this marker in the model, else ignore it
				index_model = modelMarkerSet.getIndex(markerNames[i]);
				// cout << measMarkerLocations[i] << endl;
				if(index_model >= 0){
					// 2. Check if marker used for estimation
					index = mNotUsed.findIndex(markerNames[i]);
					if (index == -1) {
						// 3. Check if marker location is measured in current frame
						// Don't use marker for update if its position is not measured [0 0 0] or if it contains NaN.

						markerNotMeasured = (SimTK::isNaN(measMarkerLocations[i][0]) || SimTK::isNaN(measMarkerLocations[i][1]) || SimTK::isNaN(measMarkerLocations[i][2]))
							||(measMarkerLocations[i][0] < tol_m && measMarkerLocations[i][0] > -tol_m && measMarkerLocations[i][1] < tol_m && measMarkerLocations[i][1] > -tol_m && measMarkerLocations[i][2] < tol_m && measMarkerLocations[i][2] > -tol_m);

						if (markerNotMeasured) {
							// cout << markerNames[i] << " not visible" << endl;
							mNotUsedInStep.append(markerNames[i]);
						}
						else {
							diff = measMarkerLocations[i] - markerLocations[il];
							innov.append(diff);
							//std::cout << markerNames[i] << "  " << measMarkerLocations[i] << "  " << markerLocations[il] << "  " << diff <<std::endl;
						}
						il = il + 1; //version 2.1
					}
					//il = il + 1;
				}
			}

			int nmUsed = innov.getSize();
			// cout << "Time = " << s.getTime() << ": number of used markers = " << nmUsed << endl;

			Vector innovation(3*nmUsed); innovation=0.0;
			for(int j=0; j<nmUsed; ++j) {
				for(int k=0; k<3; ++k)
					innovation[j*3+k] = innov[j][k];
			}

			// Construct Jacobian
			Matrix measModelJ(3*nmUsed, nu); measModelJ=0.0;
			int indexU = 0;
			index = 0;
			il = 0;
			//Loop through all markers in the reference
			for(unsigned int i=0; i < markerNames.size(); i++){
			// Check if we have this marker in the model, else ignore it
				index = modelMarkerSet.getIndex(markerNames[i]);
				if(index >= 0){
					// Check if the marker is used in the current frame
					indexU = mNotUsedInStep.findIndex(markerNames[i]);
					if (indexU == -1) {
						Marker &marker = modelMarkerSet[index];
						// marker.getFrameName(); // frame names correspond to body names
						// cout << marker.getFrameName() << endl;
						// OpenSim::Body& body = marker.getBody(); // had to replace this
						// OpenSim::Frame& frame = _model->getFrameSet().get(marker.getFrameName);
						/*std::string bName = marker.getFrameName();*/
                        std::string bName = marker.getParentFrameName();

						// const OpenSim::Body & body = _model->getBodySet().get(bName);
						const OpenSim::PhysicalFrame & pFrame = _model->getComponent<PhysicalFrame>(bName);
						const SimTK::MobilizedBodyIndex moBodyIndex = pFrame.getMobilizedBodyIndex();
						// cout << "mobilized body index of " << bName << " is " << moBodyIndex << endl;

						Matrix Jm(3, nu); Jm = 0.0;

						SimTK::Vec3 mLoc = marker.get_location();
						_model->getMatterSubsystem().calcStationJacobian(s, moBodyIndex, mLoc, Jm);
						measModelJ(il*3, 0, 3, nu) = Jm;
						il = il + 1;

					}
				}
			}
			//std::cout << "Measured modelJ" << measModelJ << std::endl;



			Matrix P_predr(nsf, nq); P_predr=0.0;
			j=0;
			for(int k=0; k<nq; ++k) {
				if (iLocked.findIndex(k) == -1){
					P_predr(0, k, nsf, 1) = P_pred(0, j*_order, nsf, 1);
					j = j+1;
				}
			}

			Matrix NinvPT(nu, nsf); NinvPT=0.0;
			bool transpose = 0;
			for (int j=0; j<nsf; ++j) {
				_model->getMatterSubsystem().multiplyByNInv(s, transpose, ~P_predr[j], NinvPT(j));
			}


			Matrix PHT(nsf, 3*nmUsed); PHT=0.0;
			PHT = ~NinvPT * ~measModelJ;

			Matrix S(3*nmUsed, 3*nmUsed); S=0.0;
			Matrix PHTr(nq, 3*nmUsed); PHTr=0.0;
			j=0;
			for(int k=0; k<nq; ++k) {
				if (iLocked.findIndex(k) == -1){
					PHTr(k, 0, 1, 3*nmUsed) = PHT(j*_order, 0, 1, 3*nmUsed);
					j = j+1;
				}
			}

			Matrix NinvPHTr(nu, 3*nmUsed); NinvPHTr=0.0;
			for (int j=0; j<3*nmUsed; ++j) {
				_model->getMatterSubsystem().multiplyByNInv(s, transpose, PHTr(j), NinvPHTr(j));
			}
			S = measModelJ * NinvPHTr; // +R!

			SimTK::Array_<double> weights;
			markersReference.getWeights(s, weights);
			index = 0;
			il = 0;
			//Loop through all markers in the reference
			// cout << "Marker names size " << markerNames.size() << endl;
			for(unsigned int j=0; j < markerNames.size(); j++){
			// Check if we have this marker in the model, else ignore it
				index = modelMarkerSet.getIndex(markerNames[j]);
				if(index >= 0){
					// Check if the marker is used in the current frame
					index = mNotUsedInStep.findIndex(markerNames[j]);
					if (index == -1) {
						S(3*il,3*il) = S(3*il,3*il) +  pow(_sdMeas * 1/weights[j],2); // pow(_sdMeas ,2);
						S(3*il+1,3*il+1) = S(3*il+1,3*il+1) + pow(_sdMeas * 1/weights[j],2); // pow(_sdMeas ,2);
						S(3*il+2,3*il+2) = S(3*il+2,3*il+2) + pow(_sdMeas * 1/weights[j],2); // pow(_sdMeas ,2);
						il = il+1;
					}
				}
			}

			// Solve the pseudoinverse problem of ~K = pinv(~S)*~PHT;
			// ( if S has full rank : Matrix K = PHT * S.invert();
			Matrix transposeS = ~S;
			SimTK::FactorQTZ pinvS(transposeS);
			double invCond = pinvS.getRCondEstimate();

			Matrix transposeK(3*nmUsed, nsf); transposeK = 0.0;
			Matrix transposeKm(3*nmUsed, nsf); transposeKm = 0.0;
			Matrix transposePHT = ~PHT;
			pinvS.solve<Real>(transposePHT, transposeKm);
			Matrix K = ~transposeKm;

			x = x_pred + K*innovation;

			Matrix P_predr2(nq, nsf); P_predr2=0.0;
			j=0;
			for(int k=0; k<nq; ++k) {
				if (iLocked.findIndex(k) == -1){
					P_predr2(k, 0, 1, nsf) = P_pred(j*_order, 0, 1, nsf);
					j = j+1;
				}
			}

			Matrix NinvP(nu, nsf); NinvP=0.0;
			for (int j=0; j<nsf; ++j) {
				_model->getMatterSubsystem().multiplyByNInv(s, transpose, P_predr2(j), NinvP(j));
			}

			P = P_pred - K*measModelJ*NinvP;

			// check whether inequality constraints are active
			lowerClampbounds;
            Array<int> iIneqAct;
            Array<double> coordBound;

			for (int clampedNr = 0; clampedNr < nClamped; clampedNr++) {
				int DOFnr = iClamped.get(clampedNr);
				// std::cout << x(DOFnr*_order) << "  " << lowerClampbounds.get(clampedNr) << "  " << upperClampbounds.get(clampedNr) << std::endl;
				if (x(DOFnr*_order) < lowerClampbounds.get(clampedNr) || x(DOFnr*_order) > upperClampbounds.get(clampedNr)) {
					iIneqAct.append(DOFnr);
					if (x(DOFnr*_order) < lowerClampbounds.get(clampedNr)) coordBound.append(lowerClampbounds.get(clampedNr));
					else coordBound.append(upperClampbounds.get(clampedNr));
				}
			}

			int noIneqAct = iIneqAct.getSize(); // number of constraints
			if (noIneqAct > 0) {
				Matrix D(noIneqAct,nsf); D=0.0;
				Vector d(nsf); d=0.0;
				for (int i=0; i<noIneqAct; i++) {
					D(i, iIneqAct.get(i)*_order) = 1.0;
					d(iIneqAct.get(i)*_order) = coordBound[i];
				}

				Matrix DPDt = D * P * ~D;

				SimTK::FactorSVD pinvDPDt(DPDt);
				Matrix invDPDt(noIneqAct,noIneqAct); invDPDt = 0.0;
				pinvDPDt.inverse<Real>(invDPDt);
				x = x - P*~D*invDPDt*D*(x-d);
			}

			x_kf_array[n] = x;
			P_kf_array[n] = P;

			Vector qhat(nq); qhat = 0;
			j=0;
			for(int k=0; k<nq; ++k) {
				if (iLocked.findIndex(k) == -1){
					qhat(k) = x(j*_order);
					j = j+1;
				}
				else { // Update
					qhat(k) = qDefaults(k);
				}
			}


			// Friedl: unnecessary when only KS results are reported.
			 s.updQ() = qhat;
			 _model->getMultibodySystem().realize(s, SimTK::Stage::Velocity);
			 _model->getVisualizer().show(s);


		}


		//////////////
		// SMOOTHER //
		//////////////
		// RAUCH-TUNG-STRIEBEL
		Matrix invF = F.invert();
		Matrix invFQ = invF*Q;
		Vector x_rts(nsf); x_rts = x_kf_array[Nframes-1];
		SimTK::Array_<Vector> x_rts_array(Nframes,Vector(nsf));
		x_rts_array[Nframes-1] = x_rts;
		// Loop over time frames
		for (int n = 1; n < Nframes; n++) {
			int index = Nframes-1-n;

			// Probably not most robust solution!
			Matrix invPpred = P_pred_array[index+1].invert();

			Matrix Ktilde(nsf,nsf); Ktilde = 0;
			for(int i=0; i<nqf; ++i){
				for(int j=0; j<nqf; ++j)	Ktilde(i*_order, j*_order, _order, _order) = invFQ*invPpred(i*_order, j*_order, _order, _order);
			}

			Matrix Ftilde = -Ktilde;
			for(int i=0; i<nqf; ++i){
				Ftilde(i*_order, i*_order, _order, _order) = Ftilde(i*_order, i*_order, _order, _order) + invF;
			}

			x_rts = Ftilde * x_rts + Ktilde * x_pred_array[index+1];
			x_rts_array[index] = x_rts;
		}

		// Report results
		cout << "report results" << endl;
		Storage *modelMarkerLocations = _reportMarkerLocations ? new Storage(Nframes, "ModelMarkerLocations") : NULL;
		Storage *modelMarkerErrors = _reportErrors ? new Storage(Nframes, "ModelMarkerErrors") : NULL;
		// Loop over time frames
		for (int n = 0; n < Nframes; n++) {
			s.updTime() = start_time + n*dt;
			Vector q_ks(nq); q_ks = 0;
			Vector x_ks = x_rts_array[n];

			j=0;
			for(int k=0; k<nq; ++k) {
				if (iLocked.findIndex(k) == -1){
					q_ks(k) = x_ks(j*_order);
					j = j+1;
				}
				else { // Update
					q_ks(k) = qDefaults(k);
				}
			}

			s.updQ() = q_ks;

			_model->getMultibodySystem().realize(s, SimTK::Stage::Velocity);
			_model->getVisualizer().show(s);



			if(_reportErrors){
				Array<double> markerErrors(0.0, 3);
				double totalSquaredMarkerError = 0.0;
				double maxSquaredMarkerError = 0.0;
				int worst = -1;

				ikSolver.setState(s); // Update
				// ikSolver.track(s); // failure
				ikSolver.computeCurrentMarkerLocations(markerLocations);

				markersReference.getValues(s, measMarkerLocations);
				// get markers defined by the model
				const MarkerSet& modelMarkerSet = _model->getMarkerSet();

				int il = 0;
				int index = 0;
				Vec3 diff(0);
				double SquaredError = 0.0;
				nm = 0;
				bool markerNotMeasured = false;
				Array<std::string> mNotUsedInStep = mNotUsed;
				for(unsigned int i=0; i < markerNames.size(); i++){
					// 1. Check if we have this marker in the model, else ignore it
					index = modelMarkerSet.getIndex(markerNames[i]);
					if(index >= 0){
						// 2. Check if marker used for estimation
						index = mNotUsed.findIndex(markerNames[i]);
						if (index == -1) {
							// 3. Check if marker location is measured in current frame
							// Don't use marker for update if its position is not measured [0 0 0] of if it contains NaN.

							markerNotMeasured = (SimTK::isNaN(measMarkerLocations[i][0]) || SimTK::isNaN(measMarkerLocations[i][1]) || SimTK::isNaN(measMarkerLocations[i][2]))
							||(measMarkerLocations[i][0] < tol_m && measMarkerLocations[i][0] > -tol_m && measMarkerLocations[i][1] < tol_m && measMarkerLocations[i][1] > -tol_m && measMarkerLocations[i][2] < tol_m && measMarkerLocations[i][2] > -tol_m);


							if (markerNotMeasured) {
								cout << markerNames[i] << " not visible" << endl;
							}
							else {

								diff = measMarkerLocations[i] - markerLocations[il];
								SquaredError = diff.normSqr();

								totalSquaredMarkerError += SquaredError;
								if(SquaredError > maxSquaredMarkerError){
									maxSquaredMarkerError = SquaredError;
									worst = i;
								}
								nm += 1;
							}
							il = il + 1; // version 2.1
						}
						//il = il + 1;
					}
				}

				double rms = sqrt(totalSquaredMarkerError / nm);
				markerErrors.set(0, totalSquaredMarkerError);
				markerErrors.set(1, rms);
				markerErrors.set(2, sqrt(maxSquaredMarkerError));
				modelMarkerErrors->append(s.getTime(), 3, &markerErrors[0]);

				cout << "Time = " << s.getTime() << ": number of used markers = " << nm << endl;
				cout << "Frame " << n << " (t=" << s.getTime() << "):\t";
				cout << "total squared error = " << totalSquaredMarkerError;
				cout << ", marker error: RMS=" << sqrt(totalSquaredMarkerError/nm);
				cout << ", max=" << sqrt(maxSquaredMarkerError) << " (" << markerNames[worst] << ")" << endl;
			}

			if(_reportMarkerLocations){
				ikSolver.computeCurrentMarkerLocations(markerLocations);
				Array<double> locations(0.0, 3*nm);
				for(int j=0; j<nm; ++j){
					for(int k=0; k<3; ++k)
						locations.set(3*j+k, markerLocations[j][k]);
				}

				modelMarkerLocations->append(s.getTime(), 3*nm, &locations[0]);

			}

			kinematicsReporter.step(s, n);
			analysisSet.step(s, n);

		}

		if (mWeightTooHigh.size()!=0) {
			for (int marker_nr = 0; marker_nr<mWeightTooHigh.size(); marker_nr++) {
				cout << "\nWarning!! Weight for marker " << mWeightTooHigh[marker_nr] << " in TaskSet is too high == unrealistically small measurement error assumed. Choose weights between 0.25 and 4.\n" << endl;
			}
		}

		// Do the maneuver to change then restore working directory
		// so that output files are saved to same folder as setup file.
		if (_outputMotionFileName!= "" && _outputMotionFileName!="Unassigned"){
			kinematicsReporter.getPositionStorage()->print(_outputMotionFileName);
		}

		if (modelMarkerErrors) {
			Array<string> labels("", 4);
			labels[0] = "time";
			labels[1] = "total_squared_error";
			labels[2] = "marker_error_RMS";
			labels[3] = "marker_error_max";

			modelMarkerErrors->setColumnLabels(labels);
			modelMarkerErrors->setName("Model Marker Errors from IK");

			IO::makeDir(getResultsDir());
			string errorFileName = trialName + "_ik_marker_errors";
			Storage::printResult(modelMarkerErrors, errorFileName, getResultsDir(), -1, ".sto");

			delete modelMarkerErrors;
		}

		if(modelMarkerLocations){
			Array<string> labels("", 3*nm+1);
			labels[0] = "time";
			Array<string> XYZ("", 3*nm);
			XYZ[0] = "_tx"; XYZ[1] = "_ty"; XYZ[2] = "_tz";

			for(int j=0; j<nm; ++j){
				for(int k=0; k<3; ++k)
					labels.set(3*j+k+1, markerNames[j]+XYZ[k]);
			}
			modelMarkerLocations->setColumnLabels(labels);
			modelMarkerLocations->setName("Model Marker Locations from IK");

			IO::makeDir(getResultsDir());
			Storage::printResult(modelMarkerLocations, "ks_model_marker_locations", getResultsDir(), -1, ".sto");

			delete modelMarkerLocations;
		}

		IO::chDir(saveWorkingDirectory);

		success = true;

		cout << "InverseKinematicsTool completed " << Nframes << " frames in "
            << watch.getElapsedTimeFormatted() << "\n" <<endl;

	}
	catch (const std::exception& ex) {
		std::cout << "InverseKinematicsKSTool Failed: " << ex.what() << std::endl;
		throw (Exception("InverseKinematicsKSTool Failed, please see messages window for details..."));
	}

	if (modelFromFile) delete _model;
	return success;
}

// Handle conversion from older format
void InverseKinematicsKSTool::updateFromXMLNode(SimTK::Xml::Element& aNode, int versionNumber)
{
	if (versionNumber < XMLDocument::getLatestVersion()) {
		std::string newFileName = getDocumentFileName();
		if (versionNumber < 20300) {
			std::string origFilename = getDocumentFileName();
			newFileName = IO::replaceSubstring(newFileName, ".xml", "_v23.xml");
			cout << "Old version setup file encountered. Converting to new file " << newFileName << endl;
			SimTK::Xml::Document doc = SimTK::Xml::Document(origFilename);
			doc.writeToFile(newFileName);
		}
		if (versionNumber <= 20201) {
			// get filename and use SimTK::Xml to parse it
			SimTK::Xml::Document doc = SimTK::Xml::Document(newFileName);
			Xml::Element root = doc.getRootElement();
			if (root.getElementTag() == "OpenSimDocument") {
				int curVersion = root.getRequiredAttributeValueAs<int>("Version");
				if (curVersion <= 20201) root.setAttributeValue("Version", "20300");
				Xml::element_iterator iter(root.element_begin("IKTool"));
				iter->setElementTag("InverseKinematicsTool");
				Xml::element_iterator toolIter(iter->element_begin("IKTrialSet"));
				// No optimizer_algorithm specification anymore
				Xml::element_iterator optIter(iter->element_begin("optimizer_algorithm"));
				if (optIter != iter->element_end())
					iter->eraseNode(optIter);

				Xml::element_iterator objIter(toolIter->element_begin("objects"));
				Xml::element_iterator trialIter(objIter->element_begin("IKTrial"));
				// Move children of (*trialIter) to root
				Xml::node_iterator p = trialIter->node_begin();
				for (; p != trialIter->node_end(); ++p) {
					iter->insertNodeAfter(iter->node_end(), p->clone());
				}
				// Append constraint_weight of 100 and accuracy of 1e-5
				iter->insertNodeAfter(iter->node_end(), Xml::Comment(_constraintWeightProp.getComment()));
				iter->insertNodeAfter(iter->node_end(), Xml::Element("constraint_weight", "20.0"));
				iter->insertNodeAfter(iter->node_end(), Xml::Comment(_accuracyProp.getComment()));
				iter->insertNodeAfter(iter->node_end(), Xml::Element("accuracy", "1e-4"));
				// erase node for IKTrialSet
				iter->eraseNode(toolIter);
				Xml::Document newDocument;
				Xml::Element docElement = newDocument.getRootElement();
				docElement.setAttributeValue("Version", "20300");
				docElement.setElementTag("OpenSimDocument");
				// Copy all children of root to newRoot
				docElement.insertNodeAfter(docElement.node_end(), iter->clone());
				newDocument.writeToFile(newFileName);
				setDocument(new XMLDocument(newFileName));
				aNode = updDocument()->getRootDataElement();
			}
			else {
				if (root.getElementTag() == "IKTool") {
					root.setElementTag("InverseKinematicsTool");
					Xml::element_iterator toolIter(root.element_begin("IKTrialSet"));
					if (toolIter == root.element_end())
						throw (Exception("Old IKTool setup file doesn't have required IKTrialSet element.. Aborting"));
					// No optimizer_algorithm specification anymore
					Xml::element_iterator optIter(root.element_begin("optimizer_algorithm"));
					if (optIter != root.element_end())
						root.eraseNode(optIter);

					Xml::element_iterator objIter(toolIter->element_begin("objects"));
					Xml::element_iterator trialIter(objIter->element_begin("IKTrial"));
					// Move children of (*trialIter) to root
					Xml::node_iterator p = trialIter->node_begin();
					for (; p != trialIter->node_end(); ++p) {
						root.insertNodeAfter(root.node_end(), p->clone());
					}
					// Append constraint_weight of 100 and accuracy of 1e-5
					root.insertNodeAfter(root.node_end(), Xml::Comment(_constraintWeightProp.getComment()));
					root.insertNodeAfter(root.node_end(), Xml::Element("constraint_weight", "20.0"));
					root.insertNodeAfter(root.node_end(), Xml::Comment(_accuracyProp.getComment()));
					root.insertNodeAfter(root.node_end(), Xml::Element("accuracy", "1e-5"));
					// erase node for IKTrialSet
					root.eraseNode(toolIter);

					// Create an OpenSimDocument node and move root inside it
					Xml::Document newDocument;
					Xml::Element docElement = newDocument.getRootElement();
					docElement.setAttributeValue("Version", "20300");
					docElement.setElementTag("OpenSimDocument");
					// Copy all children of root to newRoot
					docElement.insertNodeAfter(docElement.node_end(), doc.getRootElement().clone());
					newDocument.writeToFile(newFileName);
					setDocument(new XMLDocument(newFileName));
					aNode = updDocument()->getRootDataElement();
				}
				else
					;   // Something wrong! bail out
			}
		}
	}
	Object::updateFromXMLNode(aNode, versionNumber);
}

void InverseKinematicsKSTool::populateReferences(MarkersReference& markersReference,
    SimTK::Array_<CoordinateReference>&coordinateReferences) const
{
    FunctionSet *coordFunctions = NULL;
    // Load the coordinate data
    // bool haveCoordinateFile = false;
    if (get_coordinate_file() != "" && get_coordinate_file() != "Unassigned") {
        Storage coordinateValues(get_coordinate_file());
        // Convert degrees to radian (TODO: this needs to have a check that the storage is, in fact, in degrees!)
        _model->getSimbodyEngine().convertDegreesToRadians(coordinateValues);
        // haveCoordinateFile = true;
        coordFunctions = new GCVSplineSet(5, &coordinateValues);
    }

    Set<MarkerWeight> markerWeights;
    // Loop through old "IKTaskSet" and assign weights to the coordinate and marker references
    // For coordinates, create the functions for coordinate reference values
    int index = 0;
    for (int i = 0; i < get_IKTaskSet().getSize(); i++) {
        if (!get_IKTaskSet()[i].getApply()) continue;
        if (IKCoordinateTask *coordTask = dynamic_cast<IKCoordinateTask *>(&get_IKTaskSet()[i])) {
            CoordinateReference *coordRef = NULL;
            if (coordTask->getValueType() == IKCoordinateTask::FromFile) {
                if (!coordFunctions)
                    throw Exception("InverseKinematicsTool: value for coordinate " + coordTask->getName() + " not found.");

                index = coordFunctions->getIndex(coordTask->getName(), index);
                if (index >= 0) {
                    coordRef = new CoordinateReference(coordTask->getName(), coordFunctions->get(index));
                }
            }
            else if ((coordTask->getValueType() == IKCoordinateTask::ManualValue)) {
                Constant reference(Constant(coordTask->getValue()));
                coordRef = new CoordinateReference(coordTask->getName(), reference);
            }
            else { // assume it should be held at its default value
                double value = _model->getCoordinateSet().get(coordTask->getName()).getDefaultValue();
                Constant reference = Constant(value);
                coordRef = new CoordinateReference(coordTask->getName(), reference);
            }

            if (coordRef == NULL)
                throw Exception("InverseKinematicsTool: value for coordinate " + coordTask->getName() + " not found.");
            else
                coordRef->setWeight(coordTask->getWeight());

            coordinateReferences.push_back(*coordRef);
        }
        else if (IKMarkerTask *markerTask = dynamic_cast<IKMarkerTask *>(&get_IKTaskSet()[i])) {
            if (markerTask->getApply()) {
                // Only track markers that have a task and it is "applied"
                markerWeights.adoptAndAppend(
                    new MarkerWeight(markerTask->getName(), markerTask->getWeight()));
            }
        }
    }

    //Read in the marker data file and set the weights for associated markers.
    //Markers in the model and the marker file but not in the markerWeights are
    //ignored
    markersReference.initializeFromMarkersFile(get_marker_file(), markerWeights);
}

// TODO to be included in populateReferences above
//Array<std::string> mNotUsed;
//Array<std::string> mInTaskSet; //version 2.1
//Array<std::string> mWeightTooHigh; //version 2.1
//
//// Loop through old "IKTaskSet" and assign weights to the coordinate and marker references
//// For coordinates, create the functions for coordinate reference values
//for(int i=0; i < _ikTaskSet.getSize(); i++){
//	if (!_ikTaskSet[i].getApply()) {
//		if(IKMarkerTask *markerTask = dynamic_cast<IKMarkerTask *>(&_ikTaskSet[i])) {
//			mNotUsed.append(markerTask->getName());
//			mInTaskSet.append(markerTask->getName()); //version 2.1
//		}
//	}
//	else if(IKMarkerTask *markerTask = dynamic_cast<IKMarkerTask *>(&_ikTaskSet[i])){
//		double w = markerTask->getWeight();
//		if(w < 0.01) w = 0.01; // prevent that uncertainty becomes too high (or infinity for 0 weight)
//		MarkerWeight *markerWeight = new MarkerWeight(markerTask->getName(), markerTask->getWeight());
//		markerWeights.adoptAndAppend(markerWeight);
//		if (w>4) mWeightTooHigh.append(markerTask->getName()); //version 2.1
//		mInTaskSet.append(markerTask->getName()); //version 2.1
//	}
//}
//
//// Control if all the markers in the model are also in the TaskSet //version 2.1
//const MarkerSet& modelMarkerSet = _model->getMarkerSet();
//for(unsigned int index_model=0; index_model < modelMarkerSet.getSize(); index_model++){
//	string modelMarkerName = modelMarkerSet.get(index_model).getName();
//	bool notInTaskSet = true;
//	for (unsigned int mInTaskSetIndex = 0; mInTaskSetIndex < mInTaskSet.size(); mInTaskSetIndex++) {
//		if (modelMarkerName.compare(mInTaskSet[mInTaskSetIndex]) == 0){
//			notInTaskSet = false;
//			break;
//		}
//	}
//	if (notInTaskSet) {
//		mNotUsed.append(modelMarkerName);
//	}
//}