Info<< "Reading thermophysical properties\n" << endl;

Info<< "Reading field T\n" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field sh\n" << endl;
volScalarField sh
(
    IOobject
    (
        "sh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field p_rgh\n" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

// 1. Explicitly create the transportProperties dictionary
Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// 2. Read standard Boussinesq properties (Required for T and rhok)
dimensionedScalar beta
(
    "beta",
    dimless/dimTemperature,
    transportProperties
);

dimensionedScalar TRef
(
    "TRef",
    dimTemperature,
    transportProperties
);

dimensionedScalar Pr
(
    "Pr",
    dimless,
    transportProperties
);

dimensionedScalar Prt
(
    "Prt",
    dimless,
    transportProperties
);

// 3. Read the new humidity properties
dimensionedScalar Sch
(
    "Sch",
    dimless,
    transportProperties
);

dimensionedScalar Scht
(
    "Scht",
    dimless,
    transportProperties
);

dimensionedScalar betaSh
(
    "betaSh",
    dimless,
    transportProperties
);

dimensionedScalar shRef
(
    "shRef",
    dimless,
    transportProperties
);

// 4. Create the Laminar Transport Model
Info<< "Creating laminar transport model\n" << endl;
singlePhaseTransportModel laminarTransport(U, phi);

// 5. Create the Turbulence Model
Info<< "Creating turbulence model\n" << endl;
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);

// Kinematic density for buoyancy force
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    // 1.0 - beta*(T - TRef)
    1.0 - beta*(T - TRef) - betaSh*(sh - shRef)
);

// kinematic turbulent thermal thermal conductivity m2/s
Info<< "Reading field alphat\n" << endl;
volScalarField alphat
(
    IOobject
    (
        "alphat",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"


volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rhok*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    pimple.dict(),
    pRefCell,
    pRefValue
);

if (p_rgh.needReference())
{
    p += dimensionedScalar
    (
        "p",
        p.dimensions(),
        pRefValue - getRefCellValue(p, pRefCell)
    );
}

mesh.setFluxRequired(p_rgh.name());

#include "createMRF.H"
#include "createIncompressibleRadiationModel.H"
#include "createFvOptions.H"
