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

typedef HashTable
        <
            const GeometricField
        	<scalar, fvPatchField, volMesh>*
        > ScalarFields;

#include "newPoints.H"

/// Read fields
#include "ReadFields.H"
//{
    /** Read fields from the first time step */
    //HashTable
    //<
    //    const GeometricField
	//	<Type, fvPatchField, volMesh>*
    //>
    //fields
    //(
    //    mesh.thisDb().objectRegistry::template lookupClass 
    //    <
    //        GeometricField
	//		<Type, fvPatchField, volMesh>
    //    > ()
    //);

    //typename 
    //HashTable
    //<
    //    const GeometricField<Type, fvPatchField, volMesh>*
    //>::iterator fieldIter;

    /** Read the first field */
    //GeometricField<Type, fvPatchField, volMesh>& field =
    //const_cast
    //<
    //    GeometricField<Type, fvPatchField, volMesh>&
    //>(*fields.begin()());

    ///** Read mesh from the first field */
    //fvMesh& mesh = field.parent();

    /** Read time from mesh */
    //const Time& runTime = mesh.time(); // From objectRegistry.H

    //for
    //(
    //    fieldIter = fields.begin();
    //    fieldIter != fields.end();
    //    ++fieldIter
    //)
    //{
    //    GeometricField<Type, fvPatchField, volMesh>& field =
    //    const_cast
    //    <
    //        GeometricField<Type, fvPatchField, volMesh>&
    //    >(*fieldIter());

    //    Info << "Field read: " << field.name() << endl;

    //    IOobject fieldIOobject
    //    (
    //        field.name(),
    //        runTime.timeName(),
    //        mesh,
    //        IOobject::MUST_READ,
    //        IOobject::AUTO_WRITE
    //    );
    //}
//}

/// Advect fields
#include "AdvectFields.H"

int main(int argc, char* argv[])
{
    #include "createTime.H"
    #include "createMesh.H"

	/// Read fields from the disk
    //#include "createFields.H"
	//HashTable
	//<
    //    const GeometricField<Type, fvPatchField, volMesh>*
	//>& fields 
    //=
	//mesh.thisDb().objectRegistry::template lookupClass
	//<
    //    GeometricField<Type, fvPatchField, volMesh>
	//> ();

	/**
	 * Algorithm:
     *
	 * Read fields every time step, advect them in the
	 * memory; write to the disk every time step,
	 * update fields by reading the again 
	 * from the disk.
	 *
	 * For fields to get updated you need to write
	 * every field to the disk and read it again,
	 * because fields is a HashTable of const 
	 * `field`s.
	 * 
	 * .--> Read fields from the disk into `fields`
	 * |    container
	 * |    |
	 * |    V
	 * |    for `field` in `fields`
     * |    |
	 * |    .-> const-cast every single field into 
	 * |  	    `field` container
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
	 * |  	    `field` container
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

    while (runTime.run())
    {
        runTime++;
        //mesh.movePoints(avgPoints(mesh, runTime));
        //mesh.movePoints(newPoints(mesh));
        
        ScalarFields fields = 
		ReadFields<scalar>(mesh);

        // Advect fields
        AdvectFields<scalar>
        (
            fields, 
			//ReadFields<scalar>(mesh),
			mesh
        );

        //runTime.write(); // Moved to AdvectFields
        //T.write();
    }
}
