/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenFOAM Foundation
\*---------------------------------------------------------------------------*/

#ifndef atmTurbulentMoistureFluxHumidityFvPatchScalarField_H
#define atmTurbulentMoistureFluxHumidityFvPatchScalarField_H

#include "fixedGradientFvPatchFields.H"
#include "PatchFunction1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
      Class atmTurbulentMoistureFluxHumidityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class atmTurbulentMoistureFluxHumidityFvPatchScalarField
:
    public fixedGradientFvPatchScalarField
{
    // Private Data

        //- Name of the effective mass diffusivity field
        word gammaEffName_;

        //- Effective latent heat scaling factor (rho * Lv)
        scalar Lv0_;

        //- Humidity flux [kg/kg . m/s] or [kg/m2/s]
        autoPtr<PatchFunction1<scalar>> HumidityFlux_;


public:

    //- Runtime type information
    TypeName("atmTurbulentMoistureFluxHumidity");


    // Constructors

        //- Construct from patch and internal field
        atmTurbulentMoistureFluxHumidityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        atmTurbulentMoistureFluxHumidityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        atmTurbulentMoistureFluxHumidityFvPatchScalarField
        (
            const atmTurbulentMoistureFluxHumidityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        atmTurbulentMoistureFluxHumidityFvPatchScalarField
        (
            const atmTurbulentMoistureFluxHumidityFvPatchScalarField&
        );

        //- Construct as copy setting internal field reference
        atmTurbulentMoistureFluxHumidityFvPatchScalarField
        (
            const atmTurbulentMoistureFluxHumidityFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new atmTurbulentMoistureFluxHumidityFvPatchScalarField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new atmTurbulentMoistureFluxHumidityFvPatchScalarField(*this, iF)
            );
        }


    // Member Functions

        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper&);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
