template<class Type>
void AdvectFields(const fvMesh& mesh)
{
    typedef GeometricField<Type, fvPatchField, volMesh> fieldType;

    HashTable
    <
        const GeometricField<Type, fvPatchField, volMesh>*
    > fields = mesh.thisDb().lookupClass
    <
        GeometricField<Type, fvPatchField, volMesh>
    >();

    typename
    HashTable
    <
        const GeometricField<Type, fvPatchField, volMesh>*
    >::iterator fieldIter;

    for
    (
        fieldIter = fields.begin();
        fieldIter !=fields.end();
        ++fieldIter
    )
    {

        GeometricField
	<
	    Type, fvPatchField, volMesh
	>& field = const_cast
	<
            GeometricField
	    <
	        Type, fvPatchField, volMesh
	    >&
	>(*fieldIter());

        fvScalarMatrix fieldEq
        (
            fvm::div(mesh.phi(), field)
         == fvm::ddt(field)
        );

        Info<< "fieldEq "
            << "created."
            << endl;

        //TEq.solve();
        fieldEq.solve();

        Info<< "Field " 
        << field.name()
        //<< fieldPtr->name()
        << " advected."
        << endl;
        
        //TPtr().correctBoundaryConditions();
        field.correctBoundaryConditions();
        //fieldPtr->correctBoundaryConditions();
        
        Info<< "Boundary condition of " 
        << field.name()
        //<< fieldPtr->name()
        << " corrected."
        << endl;

        if 
        (
            mesh.thisDb().foundObject<fieldType>
            (word(field.name() + "_0"))
        )
        {
            const_cast<objectRegistry&>(mesh.thisDb()).checkOut
            (
                const_cast<fieldType&>
                (
                    mesh.thisDb().lookupObject<fieldType>
                    (word(field.name() + "_0"))
                )
            );
        }
        
        if 
        (
            mesh.thisDb().foundObject<fieldType>
            (word(field.name() + "_0_0"))
        )
        {
            const_cast<objectRegistry&>(mesh.thisDb()).checkOut
            (
                const_cast<fieldType&>
                (
                    mesh.thisDb().lookupObject<fieldType>
                    (word(field.name() + "_0_0"))
                )
            );
        }
    }
}
