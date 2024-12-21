template<class Type>
void AdvectFields(const fvMesh& mesh)
{
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

        const_cast<objectRegistry&>(mesh.thisDb()).checkOut(field);

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

        const_cast<objectRegistry&>(mesh.thisDb()).checkIn(field);
        
        //Info<< "fieldPtr->sizr():"
        //<< fieldPtr->size()
        //<< endl;
        
        //if (mesh.time().write())
        //{
        //    //fieldPtr->write();
        //    field.write();
        //    
        //    //Info<< fieldPtr->name()
        //    Info<< field.name()
        //        << " written to " 
        //        << mesh.time().timeName()
        //        << endl;
        //}
    }
}
