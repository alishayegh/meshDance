/*
- Description
A code to read a scalar field on a particular mesh, move the mesh nodes while preserving topology, and calculate the field on the new mesh.

- Theory
Writing the integral transport equation for a variable $T$ leads to

\begin{align}
\frac{d}{dt}\int_{AR}\rho T \mathrm{d}V = 0 + \int_{AR}\rho T n_i (u_i - {u_b}_i)\mathrm{d}S
\end{align}

where AR is an arbitrary region, $u_i$ is the material velocity, {u_b}_i is the boundary velocity for the arbitrary region, and $0$ on the RHS stands for
$\int_{AR} \rho \frac{\mathrm{d}T}{\mathrm{d}t}\mathrm{d}V$, which is $0$, because we are solving for a fixed time, i.e., $t$ does not exist in our problem.

- Dependencies
foam-extend-4.0

- Signature
Maalik, Maxwell Corner, Cambridge, 040924
ali@tensorfields.com
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

using namespace Foam;

typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;

typedef GeometricField<vector, pointPatchField, pointMesh> pointVectorField;

typedef GeometricField<vector, fvPatchField, volMesh> volVectorField;

// Calc/Assume new mesh points (only internal points, not boundary points)
pointField newPoints(const fvMesh& mesh)
{
    // Find max coordinates
    scalar xMax = -GREAT;
    scalar yMax = -GREAT;
    scalar zMax = -GREAT;

    scalar xMin = GREAT;
    scalar yMin = GREAT;
    scalar zMin = GREAT;

    forAll(mesh.points(), pointI)
    {
        xMax = max(mesh.points()[pointI][0], xMax);
        yMax = max(mesh.points()[pointI][1], yMax);
        zMax = max(mesh.points()[pointI][2], zMax);

        xMin = min(mesh.points()[pointI][0], xMin);
        yMin = min(mesh.points()[pointI][1], yMin);
        zMin = min(mesh.points()[pointI][2], zMin);
    }

    scalar scale(5e-3);

    tmp<pointField> pfPtr(new pointField(mesh.points()));

    forAll(pfPtr(), pointI)
    {
    	// Avoid moving boundary points
		if (pointI > mesh.points().size())
		{
		    FatalErrorIn("point number") << exit(FatalError);
		}

        if 
    	(
             pfPtr()[pointI][0] != xMax || pfPtr()[pointI][0] != xMin
          || pfPtr()[pointI][1] != yMax || pfPtr()[pointI][1] != yMin
          || pfPtr()[pointI][2] != zMax || pfPtr()[pointI][2] != zMin
    	)
        pfPtr()[pointI][0] = mesh.points()[pointI][0] * (1 - (xMax - mesh.points()[pointI][0])/(VSMALL + xMax) * scale);
        pfPtr()[pointI][1] = mesh.points()[pointI][1] * (1 - (yMax - mesh.points()[pointI][1])/(VSMALL + yMax) * scale);
        pfPtr()[pointI][2] = mesh.points()[pointI][2] * (1 - (zMax - mesh.points()[pointI][2])/(VSMALL + zMax) * scale);
    }

	return pfPtr();
}

int main(int argc, char* argv[])
{
	// Create time
    argList args(argc, argv);

	Time runTime
	(
	    "controlDict",
		args.rootPath(),
		args.caseName(),
		"system",
		"constant",
		!args.optionFound("noFunctionObjects")
	);

    //
	// Read material velocity field and mesh velocity field
	// and the scalar field
	//

	// Read the mesh
	IOobject meshIO
	(
		fvMesh::defaultRegion,
	    runTime.timeName(),
		runTime,
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	);

	//fvMesh mesh(meshIO);
    autoPtr<fvMesh> meshPtr
	(
	    new fvMesh(meshIO)
	);

	// Read material velocity
	volVectorField U
	(
	    IOobject
		(
		    "U",
			runTime.timeName(),
			meshPtr(),
			IOobject::READ_IF_PRESENT
		),
		meshPtr(),
		dimensionedVector("0", dimLength/dimTime, vector::zero)
	);

	// Read the existing scalar field
	Info << "Read the existing scalar field" << endl;

	autoPtr<volScalarField> TPtr//OldPtr
	(
	    new volScalarField
	    (
	        IOobject
	    	(
	    	    "T",
	    		runTime.timeName(),
	    		meshPtr(),
	    		IOobject::MUST_READ,
	    		IOobject::AUTO_WRITE
	    	),
	    	meshPtr()
	    )
	);

	// Create boundaries to be added to the newMesh
	Info << "Init list of patch pointers" << endl;

	//scalar maxIter(10);

	while (runTime.run())
	{
	    runTime++;

		meshPtr().movePoints(newPoints(meshPtr()));
		
	    // Advect field
	    Info << "Calculate the new values for the scalar field" << endl;

	    //TPtr.reset(new volScalarField(advectField(TOldPtr(), meshPtr(), newMeshPtr(), Ub(newMeshPtr(), meshPtr(), runTime), U, runTime)));

        fvScalarMatrix TEq
	    (
	        fvm::ddt(TPtr()) == fvc::ddt(TPtr()) /*- fvm::div(phi, T)*/ + fvm::div(meshPtr().phi(), TPtr())
	    );

	    TEq.solve();

	    TPtr().correctBoundaryConditions();

	    //runTime++;

	    Info << "Write mesh to " << runTime.timeName() << endl;

		if(meshPtr.valid())
		{
	        meshPtr().write(); // <- does this write T, as T is registered with newMesh? no!
		}

		
		if (TPtr.valid())
		{
			TPtr().write();
		}
	}
}
