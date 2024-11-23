template<class Type>
void AdvectFields
(
	
	HashTable
	<
        const GeometricField<Type, fvPatchField, volMesh>*
	>& fields,
	const fvMesh& mesh
)
{
	/** Read fields from the first time step */
	//HashTable
	//<
    //    const GeometricField<Type, fvPatchField, volMesh>*
	//>
	//fields
	//(
	//    mesh.thisDb().objectRegistry::template lookupClass 
	//    <
    //        GeometricField<Type, fvPatchField, volMesh>
	//    > ()
	//);

	typename 
	HashTable
	<
        const GeometricField<Type, fvPatchField, volMesh>*
	>::iterator fieldIter;

	///** Read the first field */
	//GeometricField<Type, fvPatchField, volMesh>& field =
	//const_cast
	//<
	//    GeometricField<Type, fvPatchField, volMesh>&
	//>(*fields.begin()());

	///** Read mesh from the first field */
	//fvMesh& mesh = field.parent();

	/** Read time from mesh */
	const Time& runTime = mesh.time(); // From objectRegistry.H

	for
	(
	    fieldIter = fields.begin();
		fieldIter != fields.end();
		++fieldIter
	)
	{
	    GeometricField<Type, fvPatchField, volMesh>& field =
	    const_cast
	    <
	        GeometricField<Type, fvPatchField, volMesh>&
	    >(*fieldIter());

	    Info << "Calculate the new values for the scalar field" << endl;

        //fvScalarMatrix TEq
        fvScalarMatrix fieldEq
	    (
	        //fvm::ddt(TPtr()) == /*fvc::ddt(TPtr()) - fvm::div(phi, TPtr()) +*/ fvm::div(mesh.phi(), TPtr())
	        fvm::ddt(field) == /*fvc::ddt(TPtr()) - fvm::div(phi, TPtr()) +*/ fvm::div(mesh.phi(), field)
	    );

	    //TEq.solve();
	    fieldEq.solve();

	    //TPtr().correctBoundaryConditions();
	    field.correctBoundaryConditions();

        Info<< "Write "
	        << field.name()
	    	<< " to " 
	        << runTime.timeName()
	    	<< endl;

	    field.write();
	}
}
