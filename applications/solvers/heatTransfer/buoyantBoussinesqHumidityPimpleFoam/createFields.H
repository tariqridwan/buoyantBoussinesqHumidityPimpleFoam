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

Info<< "Reading field sh\n" << endl; // to read the specific humidity field, added for the buoyantBoussinesqHumidityPimpleFoam solver
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
dimensionedScalar beta // Thermal expansion coefficient
(
    "beta",
    dimless/dimTemperature,
    transportProperties
);

dimensionedScalar TRef // Reference potential temperature
(
    "TRef",
    dimTemperature,
    transportProperties
);

dimensionedScalar Pr // Laminar Prandtl number
(
    "Pr",
    dimless,
    transportProperties
);

dimensionedScalar Prt // Turbulent Prandtl number
(
    "Prt",
    dimless,
    transportProperties
);

// 3. Read the new specific humidity properties
dimensionedScalar Sch // Laminar Schmidt number
(
    "Sch",
    dimless,
    transportProperties
);

dimensionedScalar Scht // Turbulent Schmidt number
(
    "Scht",
    dimless,
    transportProperties
);

dimensionedScalar betaSh // Mass expansion coefficient
(
    "betaSh",
    dimless,
    transportProperties
);

dimensionedScalar shRef // Reference specific humidity
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
    1.0 - beta*(T - TRef) - betaSh*(sh - shRef) // modified to consider density change by specific humidity "sh", for the buoyantBoussinesqHumidityPimpleFoam solver
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
