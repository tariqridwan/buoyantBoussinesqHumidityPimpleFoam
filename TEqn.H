{
    alphat = turbulence->nut()/Prt;
    alphat.correctBoundaryConditions();

    volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(alphaEff, T)
     ==
        radiation->ST(rhoCpRef, T)
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    radiation->correct();

    fvOptions.correct(T);

    //rhok = 1.0 - beta*(T - TRef);
    rhok = 1.0 - beta*(T - TRef) - betaSh*(sh - shRef); // modified to consider density change by specific humidity "sh", for the buoyantBoussinesqHumidityPimpleFoam solver, used in the last part of Equation (21) in the doc file. Equation (20) is yet to be tested.
}
