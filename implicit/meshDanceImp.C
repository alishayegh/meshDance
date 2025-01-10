/**
* @brief Description
* A code to read a scalar field on a particular mesh, move the mesh nodes while preserving topology, and calculate the field on the new mesh.
* 
* - Theory
* Writing the integral transport equation for a variable $T$ leads to
* 
* \begin{align}
* \frac{d}{dt}\int_{AR}\rho T \mathrm{d}V = 0 + \int_{AR}\rho T n_i (u_i - {u_b}_i)\mathrm{d}S
* \end{align}
* 
* where AR is an arbitrary region, $u_i$ is the material velocity, {u_b}_i is the boundary velocity for the arbitrary region, and $0$ on the RHS stands for
* $\int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V$, which is $0$, because we are solving for a fixed time, i.e., $t$ does not exist in our problem.
* 
* - Dependencies
* foam-extend-4.0
* 
* - Roadmap
* For a given mesh associated with some fields, increase mesh quality,
* read all fields into a \c for loop,
* advect them implicitly (2nd-order accurate), and write the results.
*
* - Signature
* Maalik, Maxwell Corner, Cambridge, 040924
* ali@tensorfields.com
*/

#include "argList.H"
#include "foamTime.H"
#include "fvMesh.H"
#include "IOobject.H"
#include "GeometricField.H"
#include "vector.H"
#include "fvPatchField.H"
#include "volMesh.H"
#include "pointField.H"
#include "scalar.H"
#include "faceList.H"
#include "labelList.H"
#include "Xfer.H"
#include "pointPatchField.H"
#include "pointMesh.H"
#include "dimensionedVector.H"
#include "dimensionSets.H"
#include "fvc.H"
#include "fvm.H"
//#include "boundaryMesh.H"
//#include "fvBoundaryMesh.H"
//#include "meshReader.H"
/** forAllIter definition */
//#include "UList.H"
#include "IOobjectList.H"

using namespace Foam;

typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;

typedef GeometricField<vector, pointPatchField, pointMesh> pointVectorField;

typedef GeometricField<vector, fvPatchField, volMesh> volVectorField;

//typedef HashTable
//        <
//         const GeometricField <scalar, fvPatchField, volMesh>*
//        > ScalarFields;

/// A point motion engine
#include "newPoints.H"

/// Another point motion engine
//#include "readNewPoints.H"
pointField deltaPoints(const fvMesh& startMesh, const fvMesh& targetMesh)
{
    /// Calculate displacement increment
    scalar nTimeSteps = 
      (
          startMesh.time().endTime().value() 
        - startMesh.time().startTime().value()
      ) / startMesh.time().deltaT().value();

    return (targetMesh.points() - startMesh.points()) / nTimeSteps;
}

#include "CreateFields.H"
#include "AdvectFields.H"

//#include "CreateAndAdvectFields.H"

/**
 * Algorithm:
 *
 * Read fields every time step, advect them in the
 * memory; write to the disk every time step,
 * update fields by reading them again 
 * from the disk.
 *
 * For fields to get updated you need to write
 * every field to the disk and read it again,
 * because `fields` is a HashTable of const 
 * `field`s.
 * 
 * .--> Read fields from the disk into `fields`
 * |    container
 * |    |
 * |    V
 * |    for `field` in `fields`
 * |    |
 * |    .-> const-cast every single field into 
 * |        `field` container
 * |        |
 * |        V
 * |        Advect `field`
 * |        |
 * |        V
 * |    .-- Write `field`
 * |    |
 * |    V
 * .-<- Next timestep
 *
 * Possible alternative algorithms: 
 *
 * Read fields only once, and advect them in the
 * memory; update fields in the memory after
 * advection; write based on time controls.
 *
 *      Read fields from the disk into `fields`
 *      container
 *      |
 *      V
 * .--> for `field` in `fields`
 * |    |
 * |    .-> const-cast every single field into 
 * |        `field` container
 * |        |
 * |        V
 * |        Advect `field`
 * |        |
 * |        V
 * |    .-- Update `fields`
 * |    |
 * |    V
 * |    Write `fields`
 * |    |
 * |    V
 * .-<- Next timestep
 */

int main(int argc, char* argv[])
{
    argList::validArgs.append("Smooth mesh directory name (in constant)");

    argList args(argc, argv);

    const word& smoothMeshDir = args.additionalArgs()[0];

    #include "createTime.H"
    #include "createMesh.H"

    /// Calculate mesh motion increments
    pointField deltaPoint = deltaPoints
    (
        mesh,
        fvMesh
        (
            IOobject
            (
                smoothMeshDir,
                runTime.constant(),
                runTime,
                IOobject::MUST_READ
            )
        )
    );

    /// Initialize the new points
    pointField newPoints = mesh.points();

    //IOobjectList objects(mesh, runTime.timeName());
        
    /// Create fields, register them with 'mesh'
    //CreateFields<scalar>(objects, mesh);

    //Info<< "Fields created." 
    //    << endl;

    /// Create IOobjects from existing field files
    IOobjectList objects(mesh, runTime.timeName());

    CreateFields<scalar>(objects, mesh);
    //CreateAndAdvectFields<scalar>
    //(
    //    objects,
    //    mesh
    //);

    //Info<< "mesh.toc() after create fields:"
    //    << mesh.toc()
    //    << endl;

    while (/*runTime.run()*/runTime.loop())
    {
        // Force V allocation for V0 to be stored
        //const scalar totalVol = gSum(mesh.V());

        //mesh.movePoints(avgPoints(mesh, runTime));
        //mesh.movePoints(newPoints(mesh));

        // Move mesh points
        newPoints += deltaPoint;

        // Move mesh
        mesh.movePoints(newPoints);
        
        Info<< "Time = "
            << runTime.timeName()
            << endl;

        //ScalarFields sFields = 

        //HashTable
        //<
        //    GeometricField<scalar, fvPatchField, volMesh>
        //> sFields = ReadFields<scalar>(runTime, mesh);

        //runTime++; // Only necessary with runTime.run()

        // Advect fields
        //CreateAndAdvectFields<scalar>
        //(
        //       objects,
        //       mesh
        //

        //Info<< "mesh.toc() before advect fields:"
        //    << mesh.toc()
        //    << endl;

	AdvectFields<scalar>(mesh);

        //Info<< "mesh.toc() after advect fields:"
        //    << mesh.toc()
        //    << endl;

        //Info<< "Fields advected." 
        //    << endl;

        //mesh.write(); 
        // Time dirs written only when mesh changes (what if feild changes?)

        runTime.write();// Time dirs written without fields
        //mesh.write();// Time dirs written without fields

        Info<< "Fields written to "
            << runTime.timeName()
            << endl;

        //T.write();
    
        //mesh["T"]->write(); // Magic :]
    }
}
