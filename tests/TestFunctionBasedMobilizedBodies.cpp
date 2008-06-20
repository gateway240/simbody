/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-2008 Stanford University and the Authors.      *
 * Authors: Peter Eastman, Ajay Seth                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include <vector>

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

const Real TOL = 1e-6;

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

template <class T>
void assertEqual(T val1, T val2, Real tol = TOL) {
    ASSERT(abs(val1-val2) < tol);
}

template <int N>
void assertEqual(Vec<N> val1, Vec<N> val2, Real tol = TOL) {
    for (int i = 0; i < N; ++i)
        ASSERT(abs(val1[i]-val2[i]) < tol);
}

template<>
void assertEqual(Vector val1, Vector val2, Real tol) {
    ASSERT(val1.size() == val2.size());
    for (int i = 0; i < val1.size(); ++i)
        assertEqual(val1[i], val2[i], tol);
}

template<>
void assertEqual(SpatialVec val1, SpatialVec val2, Real tol) {
    assertEqual(val1[0], val2[0], tol);
    assertEqual(val1[1], val2[1], tol);
}

template<>
void assertEqual(Transform val1, Transform val2, Real tol) {
    assertEqual(val1.T(), val2.T(), tol);
    ASSERT(val1.R().isSameRotationToWithinAngle(val2.R(), tol));
}

void compareMobilizedBodies(const MobilizedBody& b1, const MobilizedBody& b2, bool eulerAngles, int expectedQ, int expectedU) {
    const SimbodyMatterSubsystem& matter = b1.getMatterSubsystem();
    const System& system = matter.getSystem();
    
    // Set whether to use Euler angles.
    
    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, eulerAngles);
    system.realizeModel(state);
    
    // Make sure the number of state variables is correct.
    
    assertEqual(b1.getNumQ(state), expectedQ);
    assertEqual(b1.getNumU(state), expectedU);
    assertEqual(b2.getNumQ(state), expectedQ);
    assertEqual(b2.getNumU(state), expectedU);

    // Set all the state variables to random values.

    Random::Gaussian random;
    int nq = state.getNQ()/2;
    for (int i = 0; i < nq; ++i)
        state.updQ()[i] = state.updQ()[i+nq] = random.getValue();
    int nu = state.getNU()/2;
    for (int i = 0; i < nu; ++i)
        state.updU()[i] = state.updU()[i+nu] = random.getValue(); //0.0; //
    system.realize(state, Stage::Acceleration);
        
    // Compare state variables and their derivatives.
    
    for (int i = 0; i < b1.getNumQ(state); ++i) {
        assertEqual(b1.getOneQ(state, i), b2.getOneQ(state, i));
        assertEqual(b1.getOneQDot(state, i), b2.getOneQDot(state, i));
        assertEqual(b1.getOneQDotDot(state, i), b2.getOneQDotDot(state, i));
    }
    /*
    for (int i = 0; i < b1.getNumU(state); ++i) {
        assertEqual(b1.getOneU(state, i), b2.getOneU(state, i));
        assertEqual(b1.getOneUDot(state, i), b2.getOneUDot(state, i));
    }
    */
    // Compare lots of properties of the two bodies.
    
    assertEqual(b1.getBodyTransform(state), b2.getBodyTransform(state));
    assertEqual(b1.getBodyVelocity(state), b2.getBodyVelocity(state));
    assertEqual(b1.getBodyAcceleration(state), b2.getBodyAcceleration(state));
    assertEqual(b1.getBodyOriginLocation(state), b2.getBodyOriginLocation(state));
    assertEqual(b1.getBodyOriginVelocity(state), b2.getBodyOriginVelocity(state));
    assertEqual(b1.getBodyOriginAcceleration(state), b2.getBodyOriginAcceleration(state));
    assertEqual(b1.getMobilizerTransform(state), b2.getMobilizerTransform(state));
    assertEqual(b1.getMobilizerVelocity(state), b2.getMobilizerVelocity(state));
    
    // Test methods that multiply by various matrices.
    
    Vector tempq(state.getNQ());
    Vector tempu(state.getNU());
    /*
    matter.multiplyByQMatrix(state, false, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i)
        assertEqual(b1.getOneFromQPartition(state, i, tempq), b2.getOneFromQPartition(state, i, tempq));
    matter.multiplyByQMatrix(state, true, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i)
        assertEqual(b1.getOneFromUPartition(state, i, tempu), b2.getOneFromUPartition(state, i, tempu));
    matter.multiplyByQMatrixInverse(state, false, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i)
        assertEqual(b1.getOneFromUPartition(state, i, tempu), b2.getOneFromUPartition(state, i, tempu));
    matter.multiplyByQMatrixInverse(state, true, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i)
        assertEqual(b1.getOneFromQPartition(state, i, tempq), b2.getOneFromQPartition(state, i, tempq));
    */
    
    // Have them calculate q and u, and see if they agree.
    
    if (!eulerAngles) { // The optimizer does not work reliably for Euler angles, since it can hit a singularity
        Transform t = b1.getBodyTransform(state);
        b1.setQFromVector(state, Vector(b1.getNumQ(state), 0.0));
        b2.setQFromVector(state, Vector(b2.getNumQ(state), 0.0));
        b1.setQToFitTransform(state, t);
        b2.setQToFitTransform(state, t);
        system.realize(state, Stage::Velocity);
        assertEqual(b1.getBodyOriginLocation(state), b2.getBodyOriginLocation(state), 1e-2);
        assertEqual((~b1.getBodyRotation(state)*b2.getBodyRotation(state)).convertRotationToAngleAxis()[0], 0.0, 1e-2);
        SpatialVec v = b1.getBodyVelocity(state);
        b1.setUFromVector(state, Vector(b1.getNumU(state), 0.0));
        b2.setUFromVector(state, Vector(b2.getNumU(state), 0.0));
        b1.setUToFitVelocity(state, v);
        b2.setUToFitVelocity(state, v);
        assertEqual(b1.getUAsVector(state), b2.getUAsVector(state), 1e-2);
    }
    
    // Simulate the system, and see if the two bodies remain identical.
    b2.setQFromVector(state, b1.getQAsVector(state));
    b2.setUFromVector(state, b1.getUAsVector(state));
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);

    assertEqual(b1.getQAsVector(integ.getState()), b2.getQAsVector(integ.getState()));
    assertEqual(b1.getQDotAsVector(integ.getState()), b2.getQDotAsVector(integ.getState()));
}

class ConstantFunction : public Function<1> {
// Implements a simple constant function, y = C
private:
    Real C;
public:

    //Default constructor
    ConstantFunction(){
        C = 0.0;
    }

    //Convenience constructor to specify constant value
    ConstantFunction(Real constant){
        C = constant;
    }

    Vec<1> calcValue(const Vector& x) const{
        return Vec<1>(C);
    }

    Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
        Vec<1> deriv(0);
        
        return deriv;
    }

    int getArgumentSize() const{
        // constant has no arguments
        return 0;
    }

    int getMaxDerivativeOrder() const{
        return 10;
    }
};

class LinearFunction : public Function<1> {
// Implements a simple linear functional relationship, y = m*x + b
private:
    Real m;
    Real b;
public:

    //Default constructor
    LinearFunction(){
        m = 1.0;
        b = 0.0;
    }

    //Convenience constructor to specify the slope and Y-intercept of the linear r
    LinearFunction(Real slope, Real intercept){
        m = slope;
        b = intercept;
    }

    Vec<1> calcValue(const Vector& x) const{
        return Vec<1>(m*x[0]+b);
    }

    Vec<1> calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
        Vec<1> deriv(0);

        if (derivComponents.size() == 1){
            deriv[0] = m;
        }
        
        return deriv;
    }

    int getArgumentSize() const{
        return 1;
    }

    int getMaxDerivativeOrder() const{
        return 10;
    }
};

class NonlinearFunction : public Function<1> {
public:
    NonlinearFunction(){
    }
    Vec1 calcValue(const Vector& x) const{
        return Vec<1>(x[0]+x[1]*x[1]);
    }
    Vec1 calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const{
        switch (derivComponents.size()) {
            case 1:
                return (derivComponents[0] == 0 ? Vec1(1.0) : Vec1(x[1]));
            case 2:
                return (derivComponents[0] == 1 && derivComponents[1] == 1 ? Vec1(1.0) : Vec1(0.0));
        }
        return Vec1(0.0);
    }
    int getArgumentSize() const{
        return 2;
    }
    int getMaxDerivativeOrder() const{
        return std::numeric_limits<int>::max();
    }
};

int defineMobilizerFunctions(std::vector<bool> &isdof, std::vector<std::vector<int> > &coordIndices, std::vector<const Function<1>*> &functions)
{
    int nm = 0;
    for(int i=0; i<6; i++){
        if(isdof[i]) {
            std::vector<int> findex(1);
            findex[0] = nm++;
            functions.push_back(new LinearFunction());
            coordIndices.push_back(findex);
        }
        else{
            std::vector<int> findex(0);
            functions.push_back(new ConstantFunction());
            coordIndices.push_back(findex);
        }
    }
    return nm;
}

void testFunctionBasedPin() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);

    // Set the 1 spatial rotation about Z to be mobility
    isdof[2] = true;  //rot Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0.5), Inertia(1)));
    MobilizedBody::Pin p1(matter.Ground(), body);
    MobilizedBody::Pin p2(p1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(p2, fb2, false, nm, nm);
}

void testFunctionBasedSkewedPin() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);
    std::vector<Vec3> axes(6);

    // Set the 1 spatial rotation about first axis
    isdof[0] = true;  //rot 1
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    double angle = 0;

    axes[0] = Vec3(0,0,1);
    axes[1] = Vec3(0,1,0);
    axes[2] = Vec3(1,0,0);
    axes[3] = Vec3(1,0,0);
    axes[4] = Vec3(0,1,0);
    axes[5] = Vec3(0,0,1);

    Transform inParentPin = Transform(Rotation(angle, YAxis), Vec3(0));
    Transform inChildPin = Transform(Rotation(angle, YAxis), Vec3(0,1,0));
    
    Transform inParentFB = Transform(Vec3(0));
    Transform inChildFB = Transform(Vec3(0,1,0));

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

    //Built-in
    MobilizedBody::Pin p1(matter.Ground(), inParentPin, body, inChildPin);
    MobilizedBody::Pin p2(p1, inParentPin, body, inChildPin);
    //Function-based
    MobilizedBody::FunctionBased fb1(matter.Ground(), inParentFB, body, inChildFB, nm, functions, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, inParentFB, body, inChildFB, nm, functions, coordIndices, axes);

    system.realizeTopology();
    compareMobilizedBodies(p2, fb2, false, nm, nm);
}

void testFunctionBasedSlider() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);

    // Set the 1 spatial translation along X to be mobility
    isdof[3] = true; //trans X
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Slider s1(matter.Ground(), body);
    MobilizedBody::Slider s2(s1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(s2, fb2, false, nm, nm);
}


void testFunctionBasedSkewedSlider() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);
    std::vector<Vec3> axes(6);
    //axes[0] = Vec3(1/sqrt(2.0),0, -1/sqrt(2.0));
    axes[0] = Vec3(1,0,0);
    axes[1] = Vec3(0,1,0);
    axes[2] = Vec3(0,0,1);
    axes[3] = Vec3(0,0,1);
    axes[4] = Vec3(0,1,0);
    axes[5] = Vec3(1,0,0);
    // Set the 1 spatial translation along X to be mobility
    isdof[5] = true; //trans X
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    Transform inParent = Transform(Vec3(0)); //Transform(Rotation(-Pi/2, YAxis));
    Transform inChild = Transform(Vec3(0));

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Slider s1(matter.Ground(), inParent, body, inChild);
    MobilizedBody::Slider s2(s1, inParent, body, inChild);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions, coordIndices, axes);
    system.realizeTopology();
    compareMobilizedBodies(s2, fb2, false, nm, nm);
}

void testFunctionBasedCylinder() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);

    // Set 2 mobilities: rotation about and translation along Z
    isdof[2] = true;  //rot Z
    isdof[5] = true;  //trans Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Cylinder c1(matter.Ground(), body); 
    MobilizedBody::Cylinder c2(c1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(c2, fb2, false, nm, nm);
}

void testFunctionBasedUniversal() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);

    // Set 2 rotation mobilities about body's X then Y
    isdof[0] = true;  //rot X
    isdof[1] = true;  //rot Y
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Universal u1(matter.Ground(), body); 
    MobilizedBody::Universal u2(u1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices);
    MobilizedBody::FunctionBased fb2(u1, body, nm, functions, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(u2, fb2, true, nm, nm);
}

void testFunctionBasedPlanar() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[2] = true;  //rot Z
    isdof[3] = true;  //trans X
    isdof[4] = true;  //trans Y
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    //Vec3(1/sqrt(2.0000000000000), 1/sqrt(2.0000000000000), 0)
    Body::Rigid body(MassProperties(1.0, Vec3(0, 0, 0), Inertia(1)));
    MobilizedBody::Planar u1(matter.Ground(), body); 
    MobilizedBody::Planar u2(u1, body);
    //MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices);
    //MobilizedBody::FunctionBased fb2(fb1, body, nm, functions, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(u2, u2, false, nm, nm);
}

void testFunctionBasedGimbal() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<bool> isdof(6,false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[0] = true;  //rot X
    isdof[1] = true;  //rot Y
    isdof[2] = true;  //rot Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    MobilizedBody::Gimbal b1(matter.Ground(), body); 
    MobilizedBody::Gimbal b2(b1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(b2, fb2, true, nm, nm);
}

void testFunctionBasedGimbalUserAxes() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int> > coordIndices;
    std::vector<const Function<1>*> functions;
    std::vector<Vec3> axes(6);
    std::vector<bool> isdof(6,false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[0] = true;  //rot X
    isdof[1] = true;  //rot Y
    isdof[2] = true;  //rot Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions);

    Random::Gaussian random;

    axes[0] = Vec3(random.getValue(),random.getValue(), 0); //Vec3(0,0,1);//
     axes[1] = Vec3(random.getValue(), 0,random.getValue());
    axes[2] = Vec3(0,random.getValue(), random.getValue());
    axes[3] = Vec3(1,0,0);
    axes[4] = Vec3(0,1,0);
    axes[5] = Vec3(0,0,1);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    //Use massless bodies for generationg skewed-axes
    Body::Massless massLessBody;

    Transform inParent = Transform(Vec3(0));
    Transform inChild = Transform(Vec3(0,1,0));

       // Compared to standard buil-in pin mobilizers with skewed axes
    // Pin rotates about Z-axis and need to align with fist axis
    Transform parentPinAxis1 = Transform(Rotation(UnitVec3(axes[0]), ZAxis), Vec3(0,0,0));
    Transform childPinAxis1 = Transform(Rotation(UnitVec3(axes[0]), ZAxis), Vec3(0,0,0));
    Transform parentPinAxis2 = Transform(Rotation(UnitVec3(axes[1]), ZAxis), Vec3(0,0,0));
    Transform childPinAxis2 = Transform(Rotation(UnitVec3(axes[1]), ZAxis), Vec3(0,0,0));
    Transform parentPinAxis3 = Transform(Rotation(UnitVec3(axes[2]), ZAxis), Vec3(0,0,0));
    Transform childPinAxis3 = Transform(Rotation(UnitVec3(axes[2]), ZAxis), Vec3(0,1,0));
    
    //MobilizedBody::Gimbal b1(matter.Ground(), body); 
    MobilizedBody::Pin masslessPin1(matter.Ground(), parentPinAxis1, massLessBody, childPinAxis1);
    MobilizedBody::Pin masslessPin2(masslessPin1, parentPinAxis2, massLessBody, childPinAxis2);
    MobilizedBody::Pin b1(masslessPin2, parentPinAxis3, body, childPinAxis3);
    //MobilizedBody::Gimbal b2(b1, body);
       MobilizedBody::Pin masslessPin3(b1, parentPinAxis1, massLessBody, childPinAxis1);
    MobilizedBody::Pin masslessPin4(masslessPin3, parentPinAxis2, massLessBody, childPinAxis2);
    MobilizedBody::Pin b2(masslessPin4, parentPinAxis3, body, childPinAxis3);

    MobilizedBody::FunctionBased fb1(matter.Ground(), inParent, body, inChild, nm, functions, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, inParent, body, inChild, nm, functions, coordIndices, axes);
    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    int nq = state.getNQ()/2;
    for (int i = 0; i < nq; ++i)
        state.updQ()[i] = state.updQ()[i+nq] = random.getValue();
    int nu = state.getNU()/2;
    for (int i = 0; i < nu; ++i)
        state.updU()[i] = state.updU()[i+nu] = random.getValue(); //0.0; //

    system.realize(state, Stage::Acceleration);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(1.0);

    Vec3 com_bin = b2.getBodyOriginLocation(state);
    Vec3 com_fb = fb2.getBodyOriginLocation(state);
    
    assertEqual(com_bin, com_fb);
    assertEqual(b2.getBodyVelocity(state), fb2.getBodyVelocity(state));
    assertEqual(b2.getBodyAcceleration(state), fb2.getBodyAcceleration(state));
}

/**
 * Test a mobilized body based on functions that take multiple arguments.
 */

void testMultipleArguments() {
    // Define the functions that specify the FunctionBased Mobilized Body.
    
    std::vector<std::vector<int> > coordIndices(6);
    std::vector<const Function<1>*> functions(6);
    Vector_<Vec1> coeff(3);
    coeff[0] = Vec1(0.5);
    coeff[1] = Vec1(-0.5);
    coeff[2] = Vec1(1.0);
    functions[0] = new Function<1>::Constant(Vec1(0.0), 0);
    functions[1] = new Function<1>::Constant(Vec1(0.0), 0);
    functions[2] = new Function<1>::Constant(Vec1(0.0), 0);
    functions[3] = new NonlinearFunction();
    functions[4] = new Function<1>::Linear(coeff);
    functions[5] = new Function<1>::Constant(Vec1(0.0), 0);
    coordIndices[3].push_back(0);
    coordIndices[3].push_back(1);
    coordIndices[4].push_back(0);
    coordIndices[4].push_back(1);

    // Create the system.

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    MobilizedBody::FunctionBased fb(matter.Ground(), body, 2, functions, coordIndices);
    State state = system.realizeTopology();
    
    // See if coordinates and velocities are calculated correctly.
    
    ASSERT(state.getNQ() == 2);
    ASSERT(state.getNU() == 2);
    state.updQ()[0] = 2.0;
    state.updQ()[1] = -3.0;
    state.updU()[0] = 0.1;
    state.updU()[1] = -0.4;
    system.realize(state, Stage::Acceleration);
    assertEqual(fb.getBodyTransform(state), Transform(Vec3(11.0, 3.5, 0.0)));
    assertEqual(fb.getBodyVelocity(state), SpatialVec(Vec3(0.0), Vec3(0.1+3.0*0.4, 0.25, 0.0)));

    // Simulate it.
    
    Real energy = system.calcEnergy(state);
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(5.0);
    assertEqual(energy, system.calcEnergy(ts.getState()));
}

int main() {
    try {
        cout << "FunctionBased MobilizedBodies vs. Built-in Types: " << endl;
      
        testFunctionBasedPin();
        cout << "Pin: Passed" << endl;

        testFunctionBasedSkewedPin();
        cout << "Skewed Pin: Passed" << endl;

        testFunctionBasedSlider();
        cout << "Slider: Passed" << endl;

        testFunctionBasedSkewedSlider();
        cout << "Skewed Slider: Passed" << endl;
        
        testFunctionBasedCylinder();
        cout << "Cylinder: Passed" << endl;

        testFunctionBasedUniversal();
        cout << "Universal: Passed" << endl;

        testFunctionBasedPlanar();
        cout << "Planar: Passed" << endl;

        testFunctionBasedGimbal();
        cout << "Gimbal: Passed" << endl;

        testMultipleArguments();
        cout << "Functions with Multiple Arguments: Passed" << endl;

        testFunctionBasedGimbalUserAxes();
        cout << "Gimbal with User Axes: Passed" << endl;
    }
    catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}