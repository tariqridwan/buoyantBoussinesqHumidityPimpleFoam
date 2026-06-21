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
    
\*---------------------------------------------------------------------------*/

#include "atmTurbulentMoistureFluxHumidityFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::
atmTurbulentMoistureFluxHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    gammaEffName_("undefinedGammaEff"),
    Lv0_(1.0),
    MoistureFlux_(nullptr)
{}


Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::
atmTurbulentMoistureFluxHumidityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF), // Bypass dict constructor to handle startup safely
    gammaEffName_(dict.get<word>("gammaEff")),
    Lv0_(dict.getOrDefault<scalar>("Lv0", 1.0)),
    MoistureFlux_(PatchFunction1<scalar>::New(p.patch(), "MoistureFlux", dict))
{
    // Ensure the user supplied 'flux' mode if they included the keyword
    if (dict.found("MoistureSource"))
    {
        word sourceMode = dict.get<word>("MoistureSource");
        if (sourceMode != "flux")
        {
            FatalIOErrorInFunction(dict)
                << "MoistureSource must be 'flux'. 'power' is not supported."
                << exit(FatalIOError);
        }
    }

    // Safety initialization rule for time 0
    if (!this->readGradientEntry(dict) || !this->readValueEntry(dict))
    {
        extrapolateInternal();
        gradient() = Zero;
    }
}


Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::
atmTurbulentMoistureFluxHumidityFvPatchScalarField
(
    const atmTurbulentMoistureFluxHumidityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper),
    gammaEffName_(ptf.gammaEffName_),
    Lv0_(ptf.Lv0_),
    MoistureFlux_(ptf.MoistureFlux_->clone(p.patch()).ptr())
{}


Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::
atmTurbulentMoistureFluxHumidityFvPatchScalarField
(
    const atmTurbulentMoistureFluxHumidityFvPatchScalarField& ptf
)
:
    fixedGradientFvPatchScalarField(ptf),
    gammaEffName_(ptf.gammaEffName_),
    Lv0_(ptf.Lv0_),
    MoistureFlux_(ptf.MoistureFlux_->clone(this->patch().patch()).ptr())
{}


Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::
atmTurbulentMoistureFluxHumidityFvPatchScalarField
(
    const atmTurbulentMoistureFluxHumidityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(ptf, iF),
    gammaEffName_(ptf.gammaEffName_),
    Lv0_(ptf.Lv0_),
    MoistureFlux_(ptf.MoistureFlux_->clone(this->patch().patch()).ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    if (MoistureFlux_)
    {
        MoistureFlux_->autoMap(m);
    }
}


void Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchScalarField::rmap(ptf, addr);

    const auto& tiptf =
        refCast<const atmTurbulentMoistureFluxHumidityFvPatchScalarField>(ptf);

    if (MoistureFlux_)
    {
        MoistureFlux_->rmap(tiptf.MoistureFlux_(), addr);
    }
}


void Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Lookup the boundary field slice directly from the database matching OpenFOAM design rules
    const scalarField& gammaEffp =
        patch().lookupPatchField<volScalarField>(gammaEffName_);

    // Get current simulation time
    const scalar t = this->db().time().timeOutputValue();
    
    // Evaluate current flux value (e.g. reading from CSV row)
    scalarField HumFlux(MoistureFlux_->value(t));

    // Compute the boundary gradient: ∇sh = MoistureFlux / (Lv0 * gammaEff)
    gradient() = HumFlux / (Lv0_*gammaEffp + SMALL);

    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    
    os.writeEntry("MoistureSource", "flux");
    os.writeEntry("gammaEff", gammaEffName_);
    os.writeEntry("Lv0", Lv0_);

    if (MoistureFlux_)
    {
        MoistureFlux_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        atmTurbulentMoistureFluxHumidityFvPatchScalarField
    );
}
