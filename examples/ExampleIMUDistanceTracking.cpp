/* -------------------------------------------------------------------------- *
 *                     Simbody(tm) Example: IMU Tracking                      *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

/* This example shows a simple example using the Assembler for inverse
kinematics tracking where the input source is orientation observations from
an attached IMU (Inertial Measurement Unit).
*/


#include "Simbody.h"

#include <cstdio>
#include <exception>

using std::cout; using std::endl;

using namespace SimTK;

typedef SimTK::Markers::MarkerIx        MarkerIx;
typedef SimTK::Markers::ObservationIx   MarkerObsIx;

typedef SimTK::OrientationSensors::OSensorIx     IMUIx;
typedef SimTK::OrientationSensors::ObservationIx IMUObsIx;

typedef SimTK::DistanceSensors::DSensorIx     DIx;
typedef SimTK::DistanceSensors::ObservationIx DObsIx;

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity          gravity(forces, matter, -YAxis, 9.81);

    system.setUseUniformBackground(true);

    // Describe mass and visualization properties for a generic body.
    Body::Rigid bodyInfo(MassProperties(1.0, Vec3(0), UnitInertia(1)));
    bodyInfo.addDecoration(Transform(), DecorativeSphere(0.1));


    // An identity Rotation represents an IMU aligned with the body frame.
    const Rotation straightIMU;
    const Rotation tiltedIMU(Pi/4, ZAxis);

    const Rotation p1IMUOri = straightIMU;
    const Rotation p2IMUOri = tiltedIMU;
    const Rotation p3IMUOri = tiltedIMU;

    const Real p3p1Dist = 0.1;

    const Vec3 hdims_head(.3,1,.6);
    const Vec3 marker_head(.3,1,0); // top front

    bodyInfo.addDecoration(DecorativeBrick(hdims_head)
                           .setOpacity(.4).setColor(Cyan));
    bodyInfo.addDecoration(p1IMUOri, 
        DecorativeFrame(1).setColor(Green).setLineThickness(5));
    bodyInfo.addDecoration(marker_head,
        DecorativePoint().setColor(Green));


    const Real arm_length = 3;
    // Create the moving (mobilized) bodies of the pendulum.
    MobilizedBody::Pin pendulum1(matter.Ground(), Vec3(0),
            bodyInfo, Transform(Vec3(0, arm_length, 0)));
    MobilizedBody::Pin pendulum2(pendulum1, Vec3(0),
            bodyInfo, Transform(Vec3(0, arm_length, 0)));
    MobilizedBody::Pin pendulum3(pendulum2, Vec3(0),
        bodyInfo, Transform(Vec3(0, arm_length, 0)));


    // Initialize the system and state.
    State state = system.realizeTopology();


    const Real Accuracy = 1e-5;
    Assembler ik(system);
    ik.setAccuracy(Accuracy);

    Markers*            markers = new Markers();
    OrientationSensors* imus    = new OrientationSensors();
    DistanceSensors* ds    = new DistanceSensors();
    // ik.adoptAssemblyGoal(markers);
    ik.adoptAssemblyGoal(imus);
    ik.adoptAssemblyGoal(ds);
    
    const Real imuWeight = 1;
    const IMUIx    p3IMU    = imus->addOSensor(pendulum3, p3IMUOri, imuWeight);
    const IMUIx    p2IMU   = imus->addOSensor(pendulum2, p2IMUOri, imuWeight);
    const IMUIx    p1IMU   = imus->addOSensor(pendulum1, p1IMUOri, imuWeight);

    const Real dWeight = 0.5;
    const DIx   p3Dist   = ds->addDSensor(pendulum3, pendulum1, p3p1Dist, dWeight);

    // const MarkerIx headMarker = markers->addMarker(pendulum1, marker_head);

    ik.initialize(state);

    const IMUObsIx p3ObsIx = imus->getObservationIxForOSensor(p3IMU);
    const IMUObsIx p2ObsIx = imus->getObservationIxForOSensor(p2IMU);
    const IMUObsIx p1ObsIx  = imus->getObservationIxForOSensor(p1IMU);

    std::cout << "Number of D sensors: " << ds->getNumDSensors() << std::endl;
    std::cout << "Number of DS observations: " << ds->getNumObservations() << "Number of IMU observations: " << imus->getNumObservations() << std::endl;
    std::cout << "Dsensor IX: " << p3Dist << " Dsensor name: " << ds->getDSensorName(p3Dist) << std::endl;
    std::cout << "IMUS Has observation: " << imus->hasObservation(p3IMU) << std::endl;
    std::cout << "DS Has observation: " << ds->hasObservation(p3Dist) << std::endl;
    
    const DObsIx   p3DObsIx = ds->getObservationIxForDSensor(p3Dist);

    // const MarkerObsIx headMarkerObsIx = 
    //     markers->getObservationIxForMarker(headMarker);

    // Try an initial assembly to an arbitrary pose.
    imus->moveOneObservation(p2ObsIx, 
                             Rotation(SimTK::BodyRotationSequence,
                                      Pi/4, ZAxis, 0, YAxis));
    imus->moveOneObservation(p3ObsIx, 
                          Rotation(SimTK::BodyRotationSequence,
                                  -Pi/4, ZAxis, 0, YAxis));
    imus->moveOneObservation(p1ObsIx, 
                             Rotation()); // keep aligned with Ground
    // markers->moveOneObservation(headMarkerObsIx,
    //                             Vec3(0, -arm_length, 0));

    ds->moveOneObservation(p3DObsIx,0.5);

    for (DistanceSensors::DSensorIx mx(0); 
         mx < ds->getNumDSensors(); ++mx)
    {
        printf("mx=%d ox=%d err=%g\n", 
            (int)mx, (int)ds->getObservationIxForDSensor(mx),
            ds->findCurrentDSensorError(mx));
    }

    for (OrientationSensors::OSensorIx mx(0); 
         mx < imus->getNumOSensors(); ++mx)
    {
        printf("mx=%d ox=%d err=%g\n", 
            (int)mx, (int)imus->getObservationIxForOSensor(mx),
            imus->findCurrentOSensorError(mx));
    }

    
    Visualizer viz(system);
    // Show initial configuration
    viz.report(state);
    cout << "Initial state. Type any character to continue:\n";
    getchar();

    printf("Using accuracy=%g\n", ik.getAccuracyInUse());
    ik.assemble(state);

    for (OrientationSensors::OSensorIx mx(0); 
         mx < imus->getNumOSensors(); ++mx)
    {
        printf("mx=%d ox=%d err=%g\n", 
            (int)mx, (int)imus->getObservationIxForOSensor(mx),
            imus->findCurrentOSensorError(mx));
    }

    viz.report(state);
    cout << "ASSEMBLED CONFIGURATION\n";     
    cout << "Type any character to start tracking:\n";
    getchar();
    
    const double startCPU = cpuTime(), startReal = realTime();

    const int NSteps = 200;
    for (int iters=0; iters <= NSteps; ++iters) {
        const Real slow = std::sin(2*Pi*iters/100);
        const Real fast = std::sin(10*Pi*iters/100);
        Rotation p2Obs(SpaceRotationSequence, 
                          (Pi/2)*slow, ZAxis,
                          0*slow, YAxis);
        Rotation p1Obs((Pi/8)*fast, ZAxis); // shake head
        Rotation p3Obs(-(Pi/4)*slow, ZAxis);
        // Real p3dObs((double)iters / 100);
        Real p3dObs((Pi/2) * slow + 3);
        // std::cout << p3dObs << std::endl;
        // Vec3 markerObs(Vec3(0,-12,0) + slow*Vec3(5,0,0));

        // imus->moveOneObservation(p2ObsIx, p2Obs);
        // imus->moveOneObservation(p1ObsIx, p1Obs);
        imus->moveOneObservation(p3ObsIx, p3Obs);
        ds->moveOneObservation(p3DObsIx,p3dObs);
        // markers->moveOneObservation(headMarkerObsIx, markerObs);
                                        
        ik.track();
        ik.updateFromInternalState(state);
        viz.report(state);
    }

    cout << "TRACKED " << NSteps << " steps in " 
         << cpuTime()-startCPU   << " CPU s, " 
         << realTime()-startReal << " REAL s\n";

    cout << "DONE TRACKING ...\n";
    viz.report(state);

    cout << "Type any character to continue:\n";
    getchar();


  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}
