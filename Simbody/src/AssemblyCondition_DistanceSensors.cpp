/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
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

#include "SimTKmath.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/Assembler.h"
#include "simbody/internal/AssemblyCondition.h"
#include "simbody/internal/AssemblyCondition_DistanceSensors.h"
#include <map>
#include <iostream>
using std::cout; using std::endl;

using namespace SimTK;

//------------------------------------------------------------------------------
//                           DISTANCE SENSORS
//------------------------------------------------------------------------------

Rotation DistanceSensors::findCurrentDSensorDistance(DSensorIx mx) const {
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    const DSensor&                dsensor = getDSensor(mx);
    const MobilizedBody&          mobod  = matter.getMobilizedBody(dsensor.bodyB);
    const State&                  state  = getAssembler().getInternalState();
    const Rotation&               R_GB   = mobod.getBodyRotation(state);
    return R_GB * dsensor.distanceInB;
}

// goal = 1/2 sum( wi * ai^2 ) / sum(wi) for WRMS 
// ai == rotation angle between sensor and observation (-Pi:Pi)
int DistanceSensors::calcGoal(const State& state, Real& goal) const {
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();
    goal = 0;
    // Loop over each body that has one or more active dsensors.
    Real wtot = 0;
    PerBodyDSensors::const_iterator bodyp = bodiesWithDSensors.begin();
    for (; bodyp != bodiesWithDSensors.end(); ++bodyp) {
        const MobilizedBodyIndex    mobodIx      = bodyp->first;
        const Array_<DSensorIx>&    bodyDSensors = bodyp->second;
        const MobilizedBody&        mobod = matter.getMobilizedBody(mobodIx);
        const Rotation&             R_GB  = mobod.getBodyRotation(state);
        assert(bodyDSensors.size());
        // Loop over each dsensor on this body.
        for (unsigned m=0; m < bodyDSensors.size(); ++m) {
            const DSensorIx mx = bodyDSensors[m];
            const DSensor&  dsensor = dsensors[mx];
            assert(dsensor.bodyB == mobodIx); // better be on this body!
            const Rotation& R_GO = observations[getObservationIxForDSensor(mx)];
            if (R_GO.isFinite()) { // skip NaNs
                const Rotation R_GS = R_GB * dsensor.distanceInB;
                const Rotation R_SO = ~R_GS*R_GO; // error, in S
                const Vec4 aa_SO = R_SO.convertRotationToAngleAxis();
                goal += dsensor.weight * square(aa_SO[0]);
                wtot += dsensor.weight;
            }
        }
    }

    goal /= (2*wtot);

    return 0;
}
// dgoal/dq = sum( wi * ai * dai/dq ) / sum(wi)
// This calculation is modeled after Peter Eastman's gradient implementation
// in ObservedPointFitter. It treats each dsensor Distance error as a 
// potential energy function whose negative spatial gradient would be a spatial
// force F. We can then use Simbody's spatial force-to-generalized force method
// (using -F instead of F) to obtain the gradient in internal coordinates.
int DistanceSensors::
calcGoalGradient(const State& state, Vector& gradient) const {
    const int np = getNumFreeQs();
    assert(gradient.size() == np);
    const SimbodyMatterSubsystem& matter = getMatterSubsystem();

    Vector_<SpatialVec> dEdR(matter.getNumBodies());
    dEdR = SpatialVec(Vec3(0), Vec3(0));
    // Loop over each body that has one or more active dsensors.
    Real wtot = 0;
    PerBodyDSensors::const_iterator bodyp = bodiesWithDSensors.begin();
    for (; bodyp != bodiesWithDSensors.end(); ++bodyp) {
        const MobilizedBodyIndex    mobodIx     = bodyp->first;
        const Array_<DSensorIx>&     bodyDSensors = bodyp->second;
        const MobilizedBody& mobod = matter.getMobilizedBody(mobodIx);
        const Rotation& R_GB = mobod.getBodyRotation(state);
        assert(bodyDSensors.size());
        // Loop over each dsensor on this body.
        for (unsigned m=0; m < bodyDSensors.size(); ++m) {
            const DSensorIx  mx = bodyDSensors[m];
            const DSensor&   dsensor = dsensors[mx];
            assert(dsensor.bodyB == mobodIx); // better be on this body!
            const Rotation& R_GO = observations[getObservationIxForDSensor(mx)];
            if (R_GO.isFinite()) { // skip NaNs
                const Rotation R_GS = R_GB * dsensor.distanceInB;
                const Rotation R_SO = ~R_GS*R_GO; // error, in S
                const Vec4 aa_SO = R_SO.convertRotationToAngleAxis();
                const Vec3 trq_S = -dsensor.weight * aa_SO[0]
                                    * aa_SO.getSubVec<3>(1);
                const Vec3 trq_G = R_GS * trq_S;
                mobod.applyBodyTorque(state, trq_G, dEdR);
                wtot += dsensor.weight;
            }
        }
    }
    // Convert spatial forces dEdR to generalized forces dEdU.
    Vector dEdU;
    matter.multiplyBySystemJacobianTranspose(state, dEdR, dEdU);

    dEdU /= wtot;

    const int nq = state.getNQ();
    if (np == nq) // gradient is full length
        matter.multiplyByNInv(state, true, dEdU, gradient);
    else { // calculate full gradient; extract the relevant parts
        Vector fullGradient(nq);
        matter.multiplyByNInv(state, true, dEdU, fullGradient);
        for (Assembler::FreeQIndex fx(0); fx < np; ++fx)
            gradient[fx] = fullGradient[getQIndexOfFreeQ(fx)];
    }


    return 0;
}

// TODO: We want the constraint version to minimize the same goal as above. But
// there can never be more than six independent constraints on the pose of
// a rigid body; this method should attempt to produce a minimal set so that
// the optimizer doesn't have to figure it out.
int DistanceSensors::calcErrors(const State& state, Vector& err) const
{   return AssemblyCondition::calcErrors(state,err); } //TODO


int DistanceSensors::calcErrorJacobian(const State& state, Matrix& jacobian) const
{   return AssemblyCondition::calcErrorJacobian(state,jacobian); } //TODO
int DistanceSensors::getNumErrors(const State& state) const
{   return AssemblyCondition::getNumErrors(state); } //TODO

// Run through all the DSensors to find all the bodies that have at least one
// active dsensor. For each of those bodies, we collect all its dsensors so that
// we can process them all at once. Active dsensors are those whose weight is
// greater than zero. Also, if we haven't been given any observation<->dsensor 
// correspondence, we're going to assume they map directly, with each 
// ObservationIx the same as its DSensorIx.
int DistanceSensors::initializeCondition() const {
    // Fill in missing observation information if needed.
    if (observation2dsensor.empty()) {
        const Array_<DSensorIx> zeroLength; // gcc doesn't like this as a temp
        const_cast<DistanceSensors&>(*this)
            .defineObservationOrder(zeroLength);
    }

    bodiesWithDSensors.clear();
    for (DSensorIx mx(0); mx < dsensors.size(); ++mx) {
        const DSensor& dsensor = dsensors[mx];
        if (hasObservation(mx) && dsensor.weight > 0)
            bodiesWithDSensors[dsensor.bodyB].push_back(mx);
    }
    return 0;
}

// Throw away the bodiesWithDSensors map.
void DistanceSensors::uninitializeCondition() const {
    bodiesWithDSensors.clear();
}

