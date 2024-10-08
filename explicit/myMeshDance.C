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
//#include "boundaryMesh.H"
//#include "fvBoundaryMesh.H"
//#include "meshReader.H"

using namespace Foam;

typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;

typedef GeometricField<vector, pointPatchField, pointMesh> pointVectorField;

typedef GeometricField<vector, fvPatchField, volMesh> volVectorField;

// Move mesh (with corresponding field T) to newMesh and calc Tnew
volScalarField advectField
(
    const volScalarField& T,
	const fvMesh& mesh,
	const fvMesh& newMesh,
	const pointVectorField& Ub,
	const volVectorField& U,
	const Time& runTime
)
{
	/*
	tmp<IOobject> tIO (new IOobject 
	                       (
	                           T.name(), 
	                       	   runTime.timeName(), 
	                       	   newMesh, 
	                       	   IOobject::NO_READ, 
	                       	   IOobject::AUTO_WRITE
					       )
	                   );
	*/

	IOobject IO 
	(
	    T.name(), 
		runTime.timeName(), 
		newMesh, 
		IOobject::NO_READ, 
		IOobject::AUTO_WRITE
	);

    tmp<volScalarField> Tnew (new volScalarField
	                              (
								      //tIO(), 
								      IO, 
									  newMesh,
									  dimensionedScalar("0", T.dimensions(), 0)
								  )
						     );
    const cellList& cells = newMesh.cells();

	forAll(cells, cellI)
	{
		Tnew().internalField()[cellI] = T.internalField()[cellI] * mesh.V()[cellI];

	    const cell& cell = cells[cellI];

	    forAll(cell, faceI)
		{
			const face& cellFace = mesh.faces()[cell[faceI]];
			Tnew().internalField()[cellI] +=
			(
			    // T on face
                (fvc::interpolate(T)())[cell[faceI]] // <- Unconditionally unstable
			    //0;

			    // face area
			  * cellFace.normal(mesh.points())

			  & (
			      	// boundary velocity
                    cellFace.average(mesh.points(), Ub) 
				    // material velocity
                  //- cellFace.average(mesh.points(), U) 
		        )
			  * runTime.deltaT()
			).value();

			//        /*newMesh.cells()[cellI][faceI] <- Does not return a face, but faceLabel*/

			//	    (fvc::interpolate(Ub)())[mesh.cells()[cellI][faceI]] // <- Error: No matching for interpolate pointVectorField
		}

		//Tnew() += 
		//((fvc::ddt(T)())[cellI] * /* mesh.V()[cellI]*/  runTime.deltaT());//value();

	    Tnew().internalField()[cellI] /= (VSMALL + newMesh.V()[cellI]); 
	}

	Tnew().correctBoundaryConditions();
	
	return Tnew();
}

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

pointVectorField Ub
(
	const fvMesh& newMesh,
	const fvMesh& mesh,
	const Time& runTime
)
{
	// Calculate cell face velocity field
	pointMesh pMesh(mesh);

	tmp<pointVectorField> tUb
	(
	    new pointVectorField
	    (
	        IOobject
	    	(
	    	    "Ub",
	    		runTime.timeName(),
	    		mesh
			),
	    	pMesh,
	    	dimensionedVector("Ub", dimLength/dimTime, vector::zero)//,
	    )
	);

	const label nMesh = mesh.points().size();
	const label nNewMesh = newMesh.points().size();

	Info << "Read/Calculate/Assume mesh velocity field" << endl;
	forAll(tUb(), UbI)
	{
	    if(UbI > nMesh || UbI > nNewMesh)
		{
		    FatalErrorIn("Ub size") << exit(FatalError);
		}

	    tUb().internalField()[UbI][0] = ((newMesh.points()[UbI][0] - mesh.points()[UbI][0]) / (VSMALL + runTime.deltaT().value()));//.value();
	    tUb().internalField()[UbI][1] = ((newMesh.points()[UbI][1] - mesh.points()[UbI][1]) / (VSMALL + runTime.deltaT().value()));//.value();
	    tUb().internalField()[UbI][2] = ((newMesh.points()[UbI][2] - mesh.points()[UbI][2]) / (VSMALL + runTime.deltaT().value()));//.value();
	}
	//Info << "Ub.size() = " << Ub.size() << endl;

	return tUb();
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
			meshPtr()//,
			//IOobject::MUST_READ
		),
		meshPtr(),
		dimensionedVector("0", dimLength/dimTime, vector::zero)
	);

	// Read the existing scalar field
	Info << "Read the existing scalar field" << endl;

	autoPtr<volScalarField> TOldPtr
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

	//autoPtr<boundaryMesh> bMeshPtr(new boundaryMesh);//(mesh.boundaryMesh().size()); // <- no such

	// Create boundaries to be added to the newMesh
	Info << "Init list of patch pointers" << endl;

	List<polyPatch*> patches(meshPtr().boundaryMesh().size());

	scalar maxIter(10);

	for (int i=0; i < maxIter; ++i)
	{
	    // Create the new mesh

		//if(meshPtr.valid())
		//{
			
		    Info << "Create the new mesh" << endl;
		    autoPtr<fvMesh> newMeshPtr
		    (
		        new fvMesh
	            (
	            	meshIO,
	            	xferCopy(newPoints(meshPtr())),
		    		//faces,
		    		//owners,
		    		//neighbours
                    xferCopy(meshPtr().faces()),
                    xferCopy(meshPtr().faceOwner()),
                    xferCopy(meshPtr().faceNeighbour())
	            )
		    );

			if(newMeshPtr.valid())
			{
			    Info << "newMeshPtr assigned successfully" << endl;
			}
			else
			{
		        FatalErrorIn("newMeshPtr is not set")
			        << exit(FatalError);
			}
		//}
		//else
		//{
		    //autoPtr<fvMesh> newMeshPtr;

		    //FatalErrorIn("meshPtr is not valid")
			//    << exit(FatalError);
		//}

		//if(newMeshPtr.valid())
		//{
		    //newMeshPtr -> boundaryMesh().transfer(meshPtr().boundaryMesh()); // <- Error: discards qualifiers

	        Info << "Create the list of patch pointers" << endl;
		
	        forAll(patches, patchI)
	        {
	            patches[patchI] = meshPtr().boundaryMesh()[patchI].clone(newMeshPtr->boundaryMesh()).ptr();

	        	//&meshPtr().boundaryMesh()[patchI]; // <- Error: const to non-const conversion
	        }

            Info << "Add patches" << endl;
	        //newMeshPtr().addFvPatches(meshReader::polyBoundaryPatches(meshPtr())); // <- the function can be used for objects of type meshReader only
	        newMeshPtr->addFvPatches(patches);
		//}
		//else
		//{
		//    FatalErrorIn("newMeshPtr is not valid")
		//	    << exit(FatalError);
		//}
	    
        autoPtr<volScalarField> TPtr;

		if (TOldPtr.valid() /*&& newMeshPtr.valid()*/)
		{
            autoPtr<volScalarField> TPtr 
		    (
		        new volScalarField
    	        (
    	            //tIO(), 
	                IOobject
		            (
		                TOldPtr().name(),
		            	runTime.timeName(),
		            	newMeshPtr(),
		            	IOobject::NO_READ,
		            	IOobject::AUTO_WRITE
		            ),
    	       	    newMeshPtr(),
    	       	    dimensionedScalar("0", TOldPtr().dimensions(), 0)
    	        )
    	    );
		}
		else
		{
		    FatalErrorIn("newMeshPtr or TOldPtr are not valid")
			    << exit(FatalError);
		}

	    // Advect field
	    Info << "Calculate the new values for the scalar field" << endl;

	    TPtr.reset(new volScalarField(advectField(TOldPtr(), meshPtr(), newMeshPtr(), Ub(newMeshPtr(), meshPtr(), runTime), U, runTime)));

	    runTime++;

	    Info << "Write mesh to " << runTime.timeName() << endl;

		if(newMeshPtr.valid())
		{
	        newMeshPtr().write(); // <- does this write T, as T is registered with newMesh? no!
			TPtr().write();
		}
	    //Info << "Write field to " << runTime.timeName() << endl;
		//TPtr().write();
	                 // Sorry! yes! Previously Tnew was registered with mesh, not newMesh!

		// Update mesh and field pointers
		TOldPtr = TPtr;
		//TOldPtr.reset(TPtr.ptr());
		if(!TOldPtr.valid())
		{
		    FatalErrorIn("TOldPtr not valid")
			    << exit(FatalError);
		}

		//meshPtr.reset(newMeshPtr.ptr());  
		//meshPtr = newMeshPtr; // <- Causes newMeshPtr become invalid
		//meshPtr = newMeshPtr; // <- segfaults at the second iter
		//meshPtr.reset(newMeshPtr); // <- Error: no matching function
		meshPtr.reset(newMeshPtr.ptr()); 
		if(!meshPtr.valid())
		{
		    FatalErrorIn("meshPtr not valid")
			    << exit(FatalError);
		}
	}
}
