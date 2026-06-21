/*---------------------------------------------------------------------------*\
atmTurbulentMoistureFluxHumidity: Implementation of a moisture flux
                                  boundary condition in OpenFOAM.
-------------------------------------------------------------------------------
    Copyright (C) 2026 Tariq Ridwan
-------------------------------------------------------------------------------
License
    This file is part of buoyantBoussinesqHumidityPimpleFoam.

    buoyantBoussinesqHumidityPimpleFoam is free software: you can redistribute
    it and/or modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    buoyantBoussinesqHumidityPimpleFoam is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM. If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField

Group
    grpAtmWallFunctions

Description
    This boundary condition provides a fixed moisture constraint on specific
    humidity (i.e. \c sh) to specify specific humidity gradient through an
    input moisture source which can be specified as flux [kg/kg·m/s].

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            atmTurbulentMoistureFluxHumidity;
        MoistureSource  flux;
        gammaEff        gammaEff; // from shEqn.H line 8 of the
                                  // buoyantBoussinesqHumidityPimpleFoam solver
        value           uniform 0.00904; // initial guess of specific humidity
                                         // typically from initial condition
        MoistureFlux    uniform 1.4435E-04; // Kinematic moisture flux in kg/kg·m/s
        
        // Optional (inherited) entries
        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property       | Description                         | Type | Reqd  | Dflt
      MoistureSource | Moisture source type                | word |  yes  | -
      gammaEff       | Name of turbulent mass disffusivity field [kg/m/s] <!--
                                                       --> | word | yes   | -
      MoistureFlux   | Moisture source value [kg/kg·m/s (flux)] <!--
                                     --> | PatchFunction1<scalar> | yes | -
    \endtable

    The inherited entries are elaborated in:
     - \link fixedGradientFvPatchScalarField.H \endlink
     - \link PatchFunction1.H \endlink

SourceFiles
    atmTurbulentMoistureFluxHumidityFvPatchScalarField.C
    
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

        //- Humidity flux [kg/kg·m/s]
        autoPtr<PatchFunction1<scalar>> MoistureFlux_;


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
