template<class Type>
void AdvectFields
(
    //const autoPtr<fvMesh> meshPtr
    const fvMesh& mesh
)
{
	/** Read time from mesh */
	const Time& runTime = mesh.db().time(); // RHS in const

	Info << "runTime in AdvectFields = " << runTime.timeName() << endl;
	//Time& runTime = mesh.time(); //?

	HashTable
	<
        const GeometricField<Type, fvPatchField, volMesh>*
	>
	fields
	(
	    mesh.thisDb().objectRegistry::template lookupClass 
	    <
            GeometricField<Type, fvPatchField, volMesh>
	    > ()
	);

	Info << "fields: " << fields << endl;

	typename 
	HashTable
	<
        const GeometricField<Type, fvPatchField, volMesh>*
	>::iterator fieldIter;
	
	/*
	forAllIter
	(
        HashTable<const GeometricField <Type, fvPatchField, volMesh>* >,
		fields,
		fieldI
	)
	*/
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

		Info << "Field read: " << field.name() << endl;

		IOobject fieldIOobject
        (
	        field.name(),
	        runTime.timeName(),
	        mesh,
	        IOobject::MUST_READ,
	        IOobject::AUTO_WRITE
	    );

		/*
	    autoPtr
		<
		    GeometricField<Type, fvPatchField, volMesh>
		> TPtr//OldPtr
	    (
	        new GeometricField<Type, fvPatchField, volMesh>
	        (
	            fieldIOobject,
	        	mesh,
				dimensioned<Type>
				(
				    "zero",
					field.dimensions(),
					pTraits<Type>::zero
				)
	        )
	    );
		*/

	    Info << "Calculate the new values for the scalar field" << endl;

	    //TPtr.reset(new volScalarField(advectField(TOldPtr(), mesh, newMeshPtr(), Ub(newMeshPtr(), mesh, runTime), U, runTime)));

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
	}
}
