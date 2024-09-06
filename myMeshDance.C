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
void advectField
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

// Calc/Assume new mesh points (only internal points, not boundary points)
pointField& calcNewPoints(fvMesh& mesh, Time& runTime, pointField& newPoints)
{
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

    scalar scale(1e-2);

    forAll(newPoints, pointI)
    {
    	// Avoid moving boundary points
        if 
    	(
             newPoints[pointI][0] != xMax || newPoints[pointI][0] != xMin
          || newPoints[pointI][1] != yMax || newPoints[pointI][1] != yMin
          || newPoints[pointI][2] != zMax || newPoints[pointI][2] != zMin
    	)
        newPoints[pointI][0] = mesh.points()[pointI][0] * (1 - (xMax - mesh.points()[pointI][0])/(VSMALL + xMax) * scale);
        newPoints[pointI][1] = mesh.points()[pointI][1] * (1 - (yMax - mesh.points()[pointI][1])/(VSMALL + yMax) * scale);
        newPoints[pointI][2] = mesh.points()[pointI][2] * (1 - (zMax - mesh.points()[pointI][2])/(VSMALL + zMax) * scale);
    }

	return newPoints;
}

pointVectorField& calcUb
(
    pointVectorField& Ub,
	fvMesh& newMesh,
	fvMesh& mesh,
	Time& runTime
)
{
	Info << "Read/Calculate/Assume mesh velocity field" << endl;
	forAll(Ub, UbI)
	{
	    Ub.internalField()[UbI][0] = ((newMesh.points()[UbI][0] - mesh.points()[UbI][0]) / (VSMALL + runTime.deltaT().value()));//.value();
	    Ub.internalField()[UbI][1] = ((newMesh.points()[UbI][1] - mesh.points()[UbI][1]) / (VSMALL + runTime.deltaT().value()));//.value();
	    Ub.internalField()[UbI][2] = ((newMesh.points()[UbI][2] - mesh.points()[UbI][2]) / (VSMALL + runTime.deltaT().value()));//.value();
	}

	return Ub;
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
			IOobject::MUST_READ//,
			//IOobject::AUTO_WRITE
		),
		meshPtr()
	);

	// Calculate cell face velocity field
	pointMesh pMesh(meshPtr());

	pointVectorField Ub
	(
	    IOobject
		(
		    "Ub",
			runTime.timeName(),
			meshPtr()//,
			//IOobject::NO_READ,
			//IOobject::AUTO_WRITE
		),
		pMesh,
		dimensionedVector("Ub", dimLength/dimTime, vector::zero)//,
	);

	// Read the existing scalar field
	Info << "Read the existing scalar field" << endl;

	volScalarField T
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
	);

    //autoPtr<fvMesh> newMeshPtr(&mesh);

	/*
	Info << "mesh faces are " << endl;
	forAll(meshPtr().faces(), faceI)
	{
	    Info << faceI << endl;
	}
	*/
	Info << "Read boundaryMesh" << endl;
	boundaryMesh bMesh;//(mesh.boundaryMesh().size()); // <- no such

	List<polyPatch*> patches(meshPtr().boundaryMesh().size());

	scalar maxIter(10);

	for (int i=0; i < maxIter; ++i)
	{
	    // Create the new mesh

		Info << "Read points of the old mesh" << endl;
        pointField newPoints = meshPtr().points();

		Info << "Create the new mesh" << endl;
		autoPtr<fvMesh> newMeshPtr
		(
		    new fvMesh
	        (
	        	meshIO,
	        	xferCopy(calcNewPoints(meshPtr(), runTime, newPoints)),
	        	xferCopy(meshPtr().faces()),
	        	xferCopy(meshPtr().faceOwner()),
	        	xferCopy(meshPtr().faceNeighbour())
	        )
		);

		/*
		forAll(newMeshPtr().faces(), faceI)
		{
		    Info << faceI << endl;
		}
		*/

	    // Create boundaries to be added to the newMesh
		//Info << meshPtr().boundaryMesh() << endl;
		//Info << meshPtr() << endl;
	    Info << "Init list of patch pointers" << endl;
        //Info << "meshPtr().boundaryMesh().size() = " << meshPtr().boundaryMesh().size() << endl;

	    bMesh.read(meshPtr());

		Info << "Create the list of patch pointers" << endl;
	    forAll(patches, patchI)
	    {
	    	//Info << "bMesh.patches()[patchI].name(), = " << bMesh.patches()[patchI].name()  << endl;
	    	//Info << "bMesh.patches()[patchI].size(), = " << bMesh.patches()[patchI].size()  << endl;
	    	//Info << "bMesh.patches()[patchI].start(),= " << bMesh.patches()[patchI].start() << endl;
	    	//Info << "bMesh.patches()[patchI].index(),= " << bMesh.patches()[patchI].index() << endl;
	    	autoPtr<polyPatch> newPatchPtr
	    	(
	    	    new polyPatch
	    		(
	    		    bMesh.patches()[patchI].name(),
	    			bMesh.patches()[patchI].size(),
	    			bMesh.patches()[patchI].start() 
				  + meshPtr().nInternalFaces(),
	    			bMesh.patches()[patchI].index(),
	    			meshPtr().boundaryMesh()
	    		)
	    	);

	        patches[patchI] = newPatchPtr.ptr();
	    }

		//if (i != 0)
		//{
        Info << "Add patches" << endl;
		//Info << newMeshPtr() << endl;
	    //newMeshPtr().removePatches();
		//Info << "boundary size before adding patches = " << newMeshPtr().boundaryMesh().size() << endl;
	    newMeshPtr().addPatches(patches);
		//}
	    
	    // Create the scalarField on the new mesh
		Info << "Create the new scalar field" << endl;
	    volScalarField Tnew
	    (
	        IOobject
	    	(
	    	    //"Tnew",
	    	    "T",
	    		runTime.timeName(),
	    		newMeshPtr(),
	    		IOobject::NO_READ,
	    		IOobject::AUTO_WRITE
	    	),
	    	newMeshPtr(),
	    	dimensionedScalar("Tnew", dimTemperature, 0)
	    );

	    // Advect field
	    Info << "Calculate the new values for the scalar field" << endl;

	    advectField(T, Tnew, meshPtr(), newMeshPtr(), calcUb(Ub, newMeshPtr(), meshPtr(), runTime), U, runTime);

		//Info << "before" <<  meshPtr() << endl;
		//meshPtr = newMeshPtr; // <- segfaults at the second iter
		//meshPtr.reset(newMeshPtr); // <- Error: no matching function
		meshPtr.reset(newMeshPtr.ptr());

		//meshPtr.reset(newMeshPtr.ptr());
		//Info << "after" << meshPtr() << endl;
		//meshPtr.reset(newMeshPtr);

	    Info << "Correct boundary conditions" << endl;
	    Tnew.correctBoundaryConditions();

		//if(i == maxIter)
		//{
	        runTime++;

	        Info << "Write results" << endl;
	        meshPtr().write(); // <- does this write T, as T is registered with newMesh? no!
	                 // Sorry! yes! Previously Tnew was registered with mesh, not newMesh!
		//}
	}
}
