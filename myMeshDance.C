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
Maalik, Maxwell corner, Cambridge, 040924
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
#include "boundaryMesh.H"

using namespace Foam;

typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;

typedef GeometricField<vector, pointPatchField, pointMesh> pointVectorField;

typedef GeometricField<vector, fvPatchField, volMesh> volVectorField;

// Move mesh (with corresponding field T) to newMesh and calc Tnew
void moveMesh
(
    volScalarField& T,
    volScalarField& Tnew,
	fvMesh& mesh,
	fvMesh& newMesh,
	pointVectorField& Ub,
	volVectorField& U,
	Time& runTime
)
{
	forAll(newMesh.cells(), cellI)
	{
		//debug
	    //Info << "Cell " << newMesh.cells()[cellI] << " has faces "<< endl;
		
		Tnew.internalField()[cellI] = T.internalField()[cellI] * mesh.V()[cellI];

	    forAll(newMesh.cells()[cellI], faceI)
		{
		
			//Info << newMesh.faces()[newMesh.cells()[cellI][faceI]] << endl;//.normal(newMesh.points()) &
			
			Tnew.internalField()[cellI] +=
			(
				// T on face
                (fvc::interpolate(T)())[mesh.cells()[cellI][faceI]]//;
                //(fvc::interpolate(T)()).internalField()[mesh.cells()[cellI][faceI]]//;
                //(fvc::interpolate(T)()).internalField()[newMesh.cells()[cellI][faceI]]//;
                *
			    // face area
				mesh.faces()[mesh.cells()[cellI][faceI]].normal(mesh.points())
				//newMesh.faces()[newMesh.cells()[cellI][faceI]].normal(newMesh.points())
			    &
				(
			        //newMesh.cells()[cellI][faceI] <- Does not return a face, but faceLabel

					// boundary velocity
                    -mesh.faces()[newMesh.cells()[cellI][faceI]].average(mesh.points(), Ub) 
				    //(fvc::interpolate(Ub)())[mesh.cells()[cellI][faceI]] // <- Error: No matching for interpolate pointVectorField
				    //(fvc::interpolate(Ub)()).internalField()[mesh.cells()[cellI][faceI]]
			        // field velocit at face centre is
			        //- Uf[newMesh.cells()[cellI][faceI]]

					// material velocity
				    //+(fvc::interpolate(U)()).internalField()[mesh.cells()[cellI][faceI]]

			        // or
                    //- newMesh.faces()[newMesh.cells()[cellI][faceI]].average(newMesh.points(), U) 
		        )
			    *runTime.deltaT()
			).value();
		
		}

	    Tnew.internalField()[cellI] /= (VSMALL + newMesh.V()[cellI]); 
	}
}

int main(int argc, char* argv[])
{
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

	fvMesh mesh(meshIO);

	// Read material velocity
	//typedef GeometricField<vector, fvPatchField, volMesh> volVectorField;
    //      GeometricField<vector, fvPatchField, volMesh> volVectorField;
	volVectorField U
	(
	    IOobject
		(
		    "U",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ//,
			//IOobject::AUTO_WRITE
		),
		mesh
	);

	// Calc/Assume new mesh points (only internal points, not boundary points)
	pointField newPoints = mesh.points();

	// Find max coordinates
	scalar xMax = -GREAT;
	scalar yMax = -GREAT;
	scalar zMax = -GREAT;

	scalar xMin = GREAT;
	scalar yMin = GREAT;
	scalar zMin = GREAT;

	forAll(newPoints, pointI)
	{
	    xMax = max(mesh.points()[pointI][0], xMax);
	    yMax = max(mesh.points()[pointI][1], yMax);
	    zMax = max(mesh.points()[pointI][2], zMax);

	    xMin = min(mesh.points()[pointI][0], xMin);
	    yMin = min(mesh.points()[pointI][1], yMin);
	    zMin = min(mesh.points()[pointI][2], zMin);
	}

	scalar scale(50e-2);

	forAll(newPoints, pointI)
	{
		// Avoid moving boundary points
	    if 
		(
             newPoints[pointI][0] != xMax || newPoints[pointI][0] != xMin
          || newPoints[pointI][1] != yMax || newPoints[pointI][1] != yMin
          || newPoints[pointI][2] != zMax || newPoints[pointI][2] != zMin
		)
	    newPoints[pointI][0] = mesh.points()[pointI][0] * (1 - (xMax - mesh.points()[pointI][0])/(VSMALL + xMax) * runTime.deltaT().value());
	    newPoints[pointI][1] = mesh.points()[pointI][1] * (1 - (yMax - mesh.points()[pointI][1])/(VSMALL + yMax) * runTime.deltaT().value());
	    newPoints[pointI][2] = mesh.points()[pointI][2] * (1 - (zMax - mesh.points()[pointI][2])/(VSMALL + zMax) * runTime.deltaT().value());
	}

	// Create the new mesh
	fvMesh newMesh
	(
		meshIO,
		xferCopy(newPoints),
		xferCopy(mesh.faces()),
		xferCopy(mesh.faceOwner()),
		xferCopy(mesh.faceNeighbour())
	);

	boundaryMesh bMesh;//(mesh.boundaryMesh().size()); // <- no such
	                                                   //    a constructor
	bMesh.read(mesh);

	List<polyPatch*> patches(mesh.boundaryMesh().size());

	forAll(patches, patchI)
	{
		autoPtr<polyPatch> newPatchPtr
		(
		    new polyPatch
			(
			    bMesh.patches()[patchI].name(),
				bMesh.patches()[patchI].size(),
				bMesh.patches()[patchI].start(),
				bMesh.patches()[patchI].index(),
				mesh.boundaryMesh()
			)
		);

	    patches[patchI] = newPatchPtr.ptr();
	}

	newMesh.addPatches(patches);

	// Read/Calculate/Assume mesh velocity field
	//typedef GeometricField<vector, pointPatchField, pointMesh> pointVectorField;

	pointMesh pMesh(mesh);

    //wordList types
	//(
	//    pMesh.boundary().size(),
	//	calculatedFvPatchVectorField::typeName
	//);

	pointVectorField Ub
	(
	    IOobject
		(
		    "Ub",
			runTime.timeName(),
			mesh//,
			//IOobject::NO_READ,
			//IOobject::AUTO_WRITE
		),
		pMesh,
		dimensionedVector("Ub", dimLength/dimTime, vector::zero)//,
		//types
	);

	Info << "Read/Calculate/Assume mesh velocity field" << endl;
	forAll(Ub, UbI)
	{
	    //Ub[UbI][0] (newMesh.points()[UbI][0] - mesh.points()[UbI][0]) / (VSMALL + runTime.deltaT());
	    //Ub[UbI][0] = double((newMesh.points()[UbI][0] - mesh.points()[UbI][0]) / (VSMALL + runTime.deltaT()));
	    //Ub[UbI][0] = dimensionedScalar("Ub0", dimLength, (newMesh.points()[UbI][0] - mesh.points()[UbI][0])) / (VSMALL + runTime.deltaT());
	    Ub.internalField()[UbI][0] = ((newMesh.points()[UbI][0] - mesh.points()[UbI][0]) / (VSMALL + runTime.deltaT().value()));//.value();
	    Ub.internalField()[UbI][1] = ((newMesh.points()[UbI][1] - mesh.points()[UbI][1]) / (VSMALL + runTime.deltaT().value()));//.value();
	    Ub.internalField()[UbI][2] = ((newMesh.points()[UbI][2] - mesh.points()[UbI][2]) / (VSMALL + runTime.deltaT().value()));//.value();
	}
	
	// Return Ub at the face centre
	// face.average(mesh.points(), Ub);

	// Read the existing scalar field
	Info << "Read the existing scalar field" << endl;
	//typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;

	volScalarField T
	(
	    IOobject
		(
		    "T",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);

	volScalarField Tnew
	(
	    IOobject
		(
		    //"Tnew",
		    "T",
			runTime.timeName(),
			newMesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		newMesh,
		dimensionedScalar("Tnew", dimTemperature, 0)
	);

	Info << "Calculate the new values for the scalar field" << endl;

	moveMesh(T, Tnew, mesh, newMesh, Ub, U, runTime);

	runTime++;

	Info << "Correct boundary conditions" << endl;
	Tnew.correctBoundaryConditions();

	Info << "Write results" << endl;
	//T.write();
	//Tnew.write();
	//Tf.write();
	newMesh.write(); // <- does this write T, as T is registered with newMesh? no!
}
