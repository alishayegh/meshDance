#include "IOobjectList.H"

template<class Type>
void CreateAndAdvectFields
(
    const IOobjectList& objects,
    const fvMesh& mesh
)
{
    word fieldClassName
    (
        GeometricField<Type, fvPatchField, volMesh>::typeName
    );

	   Info<< "Find fields"
	       << endl;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    /// \note
    /// List of fields needs you to create all fields once, or one field N 
				/// times 
    /// and then overwrite them with
    /// their actual properties, i.e., their respective IO, through 
				/// `for loop`. The 
    /// alternatives are HashTable and PtrList.

    /// Create an empty list of fields
    //List
    //<
    //    GeometricField<Type, fvPatchField, volMesh> 
				//    // Error: The constructor GeometricField<>() does not exist
    //    GeometricField<Type, fvPatchField, volMesh>
    //      (*fieldStart(), mesh)                       
				//    // Not tested; assuming that fieldStart is defined.
    //>

    /// Create an empty HashTable of fields
    //HashTable <GeometricField<Type, fvPatchField, volMesh> > 

    /// Create an empty HashTable of fields
    /// PtrList also can take the advantage of forAll
    /// But has the same redundancy as List, i.e., you need to construct a 
				/// fieldStart as a place holder when creating the list, and
    /// then overwrite, or you get `error: hanging pointer, cannot 
				/// dereference`
    /// PtrList can be used by `set` member function to avoid `[]` operator 
				/// when assigning an element of the list to avoid
    /// hanging pointr error (not checked).
    
    //PtrList< GeometricField<Type, fvPatchField, volMesh> > 

    //List< GeometricField<Type, fvPatchField, volMesh>* >  
				// Exchanged for List<autoPtr<...>>

    	//Info<< "Create list of field pointers"
	    //<< endl;

    //List< autoPtr<GeometricField<Type, fvPatchField, volMesh> > >  
				// Exchanged for List<autoPtr<...>>

    //List< GeometricField<Type, fvPatchField, volMesh>* >  
    //fieldPtrList(fields.size());
    //fieldList(fields.size());

    /// Create fields
    //label fieldI = 0;

    for
    (
        IOobjectList::iterator fieldIter = fields.begin();
        fieldIter !=fields.end();
        ++fieldIter
    )
    {
		      Info<< "fieldIOobject created"
		          << endl;
			
        IOobject fieldIOobject
        (
            fieldIter()->name(),
            mesh.time().timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        //if (fieldIOobject.headerOk())
        {
            Info<< "    Reading "
                << fieldIter()->name()
                << endl;

            //fieldList[fieldI] = // if it is a List
            //fieldList.insert    // if it is a HashTable
            //(
            //    fieldIter()->name(),
            //    GeometricField<Type, fvPatchField, volMesh> //field
            //    (
            //        *fieldIter(),
            //        mesh
            //    )
            //);

            //fieldList[fieldI] = // 
												// if it is a PtrList 
												// (hanging pointer error); alternatively, fieldList.set(...)

												GeometricField<Type, fvPatchField, volMesh>* fieldPtr
												(
												    new GeometricField<Type, fvPatchField, volMesh>
                (
                    *fieldIter(),
                    mesh
                )
												);

            //fieldPtrList[fieldI] = new 
												//    GeometricField<Type, fvPatchField, volMesh> 
            //    (
            //        *fieldIter(),
            //        mesh
            //    );

                //autoPtr< GeometricField<Type, fvPatchField, volMesh> >
                //(
                //    new GeometricField<Type, fvPatchField, volMesh> 
                //    (
                //        *fieldIter(),
                //        mesh
                //    )
                //);
            
            //++fieldI;

    /// Advect fields
    //forAll(fieldList, fieldI)
    //forAll(fieldPtrList, fieldI)
								    /// ddt(*fieldPtr) is not understandable for the compiler, so I
								    /// create a copy
								    //GeometricField<Type, fvPatchField, volMesh> 
								    //fieldCopy = *fieldPtrList[fieldI];

            fvScalarMatrix fieldEq
            (
                fvm::div(mesh.phi(), *(fieldPtr))
             == fvm::ddt(*(fieldPtr))
            );

            Info<< "fieldEq "
                << "created."
                << endl;

            //TEq.solve();
            fieldEq.solve();

            Info<< "Field " 
                //<< field().name()
                << fieldPtr->name()
                << " advected."
                << endl;

            //TPtr().correctBoundaryConditions();
            //field.correctBoundaryConditions();
            fieldPtr->correctBoundaryConditions();

            Info<< "Boundary condition of " 
                //<< field().name()
                << fieldPtr->name()
                << " corrected."
                << endl;

												Info<< "fieldPtr->sizr():"
												    << fieldPtr->size()
																<< endl;

            fieldPtr->write();

            Info<< fieldPtr->name()
                << " written to " 
                << mesh.time().timeName()
                << endl;
				    }
    }
}
