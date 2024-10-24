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

using namespace Foam;

typedef GeometricField<scalar, fvPatchField, volMesh> volScalarField;

typedef GeometricField<vector, pointPatchField, pointMesh> pointVectorField;

typedef GeometricField<vector, fvPatchField, volMesh> volVectorField;

#include "newPoints.H"
#include "AdvectFields.H"

int main(int argc, char* argv[])
{
	#include "createTime.H"
	#include "createMesh.H"
	#include "createFields.H"

	while (runTime.run())
	{
		//mesh.movePoints(avgPoints(mesh, runTime));
		//mesh.movePoints(newPoints(mesh));
		
	    // Advect field
	    Info << "runTime before AdvectFields = " << runTime.timeName();
		AdvectFields<scalar>(mesh/*Ptr()*/);
		//AdvectFields<vector>(meshPtr/*()*/);
		//AdvectFields<tensor>(meshPtr/*()*/);

	    runTime++;

	    Info << "Write mesh and updated field to " << runTime.timeName() << endl;

		runTime.write();
	}
}
