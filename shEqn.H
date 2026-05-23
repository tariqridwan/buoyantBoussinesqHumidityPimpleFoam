    // This file is unique for this solver, to consider specific humidity
{
    //alphat = turbulence->nut()/Prt;
    //alphat.correctBoundaryConditions();
    volScalarField gammat = turbulence->nut()/Scht; // similar to alphat, here Scht is the turbulent Schmidt Number
    //gammat.correctBoundaryConditions(); // not needed as wall function for gammat is not considered. For alphat, this line was necessary to correct the alphat only between surface and first cell-center, using the alphat wall function, e.g., atmAlphatkWallFunction.

    //volScalarField alphaEff("alphaEff", turbulence->nu()/Pr + alphat);
    volScalarField gammaEff("gammaEff", turbulence->nu()/Sch + gammat); // similar to alphaEff, here Sch is the Schmidt Number

    //fvScalarMatrix TEqn
    //(
    //    fvm::ddt(T)
    //  + fvm::div(phi, T)
    //  - fvm::laplacian(alphaEff, T)
    // ==
    //    radiation->ST(rhoCpRef, T)
    //  + fvOptions(T)
    //);
    fvScalarMatrix shEqn
    (
        fvm::ddt(sh)
      + fvm::div(phi, sh)
      - fvm::laplacian(gammaEff, sh)
     ==
        fvOptions(sh)
    );

    //TEqn.relax();
    shEqn.relax();

    //fvOptions.constrain(TEqn);
    fvOptions.constrain(shEqn);

    //TEqn.solve();
    shEqn.solve();
    
    //radiation->correct(); // not necessary for specific humidity

    //fvOptions.correct(T);
    fvOptions.correct(sh);

    //rhok = 1.0 - beta*(T - TRef); // this line should be "1.0 - beta*(T - TRef) - betaq*(sh - shRef)" now in TEqn.H file
}
