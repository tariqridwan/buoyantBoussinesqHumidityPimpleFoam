/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenFOAM Foundation
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
    HumidityFlux_(nullptr)
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
    HumidityFlux_(PatchFunction1<scalar>::New(p.patch(), "HumidityFlux", dict))
{
    // Ensures the user supplied 'flux' mode if they included the keyword
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
    HumidityFlux_(ptf.HumidityFlux_->clone(p.patch()).ptr())
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
    HumidityFlux_(ptf.HumidityFlux_->clone(this->patch().patch()).ptr())
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
    HumidityFlux_(ptf.HumidityFlux_->clone(this->patch().patch()).ptr())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmTurbulentMoistureFluxHumidityFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchScalarField::autoMap(m);
    if (HumidityFlux_)
    {
        HumidityFlux_->autoMap(m);
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

    if (HumidityFlux_)
    {
        HumidityFlux_->rmap(tiptf.HumidityFlux_(), addr);
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
    scalarField HumFlux(HumidityFlux_->value(t));

    // Compute the boundary gradient: ∇sh = HumidityFlux / (Lv0 * gammaEff)
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

    if (HumidityFlux_)
    {
        HumidityFlux_->writeData(os);
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
