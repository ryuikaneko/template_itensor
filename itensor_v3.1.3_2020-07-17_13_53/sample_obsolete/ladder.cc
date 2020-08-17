// ladder DMRG:
// https://itensor.org/docs/cppv3/formulas/ladder/ladder.cc
// https://itensor.org/docs.cgi?vers=cppv3&page=formulas/ladder

// onsite and 2-point correlarions:
// http://www.itensor.org/docs.cgi?vers=cppv3&page=formulas/measure_mps
// http://www.itensor.org/docs.cgi?vers=cppv3&page=formulas/correlator_mps

// chain energy: -ln2+0.25 = -0.44314718056
// 2-leg ladder: -0.578043140180 (PhysRevB.89.094424)

#include "itensor/all.h"

using namespace itensor;

#include <iostream>
using namespace std;

int
//main()
main(int argc, char** argv)
    {
    for (int i = 0; i < argc; ++i)
        cout << argv[i] << "\n";

    printfln("argc %d\n",argc);

    auto Nx = 5;
    auto N = 2*Nx;
    auto Jx = 1.0;// leg
    auto Jy = 1.0;// rung
    auto J4 = 1.23;

    if(argc>4){
        Nx = atoi(argv[1]);
        N = 2*Nx;
        Jx = atof(argv[2]);
        Jy = atof(argv[3]);
        J4 = atof(argv[4]);
    }

    printfln("# Nx N Jx Jy J4: %4d %4d %.10f %.10f %.10f",Nx,N,Jx,Jy,J4);
    println("\n");

    // QNs are conserved by default. Use
    // the arg {"ConserveQNs=",false} to not
    // conserve QNs
//    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto sites = SpinHalf(N);

    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N-3; j += 2)
        {
        ampo +=   Jx,"Sz",j,"Sz",j+2;
        ampo += Jx/2,"S+",j,"S-",j+2;
        ampo += Jx/2,"S-",j,"S+",j+2;

        ampo +=   Jx,"Sz",j+1,"Sz",j+3;
        ampo += Jx/2,"S+",j+1,"S-",j+3;
        ampo += Jx/2,"S-",j+1,"S+",j+3;

        ampo +=   J4,"Sz",j,"Sz",j+2,"Sz",j+1,"Sz",j+3;
        ampo += J4/2,"Sz",j,"Sz",j+2,"S+",j+1,"S-",j+3;
        ampo += J4/2,"Sz",j,"Sz",j+2,"S-",j+1,"S+",j+3;
        ampo += J4/2,"S+",j,"S-",j+2,"Sz",j+1,"Sz",j+3;
        ampo += J4/4,"S+",j,"S-",j+2,"S+",j+1,"S-",j+3;
        ampo += J4/4,"S+",j,"S-",j+2,"S-",j+1,"S+",j+3;
        ampo += J4/2,"S-",j,"S+",j+2,"Sz",j+1,"Sz",j+3;
        ampo += J4/4,"S-",j,"S+",j+2,"S+",j+1,"S-",j+3;
        ampo += J4/4,"S-",j,"S+",j+2,"S-",j+1,"S+",j+3;
        }
    for(int j = 1; j <= N-1; j += 2)
        {
        ampo +=   Jy,"Sz",j,"Sz",j+1;
        ampo += Jy/2,"S+",j,"S-",j+1;
        ampo += Jy/2,"S-",j,"S+",j+1;
        }
    auto H = toMPO(ampo);

    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i)
        {
        if(i%2 == 1) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    auto psi0 = MPS(state);
//    auto psi0 = randomMPS(sites);

    auto sweeps = Sweeps(7);
    sweeps.maxdim() = 1200,1600,2000,2400,2800,3200,3600;
    sweeps.mindim() = 1,1,1,1,1,1,1;
    sweeps.cutoff() = 1E-8,1E-9,1E-10,1E-11,1E-11,1E-12,1E-12;
    sweeps.niter() = 4,3,3,2,2,2,2;
    sweeps.noise() = 1E-5,1E-5,1E-8,1E-9,1E-10,1E-10,1E-10;

    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",true});

    println("\n");

    printfln("## Nx N Jx Jy J4 energy %4d %4d %.10f %.10f %.10f %.10f",Nx,N,Jx,Jy,J4,energy/N);
    println("\n");

    // save a wave function
//    writeToFile(format("dat_sites_N%d_Jx%.10f_Jy%.10f_Jf%.10f",N,Jx,Jy,J4),sites);
//    writeToFile(format("dat_psi_N%d_Jx%.10f_Jy%.10f_Jf%.10f",N,Jx,Jy,J4),psi);



    // Measure Sz on every bond
    println("## mag j t_or_b num_from_left sx sy sz");
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        auto ket = psi.A(j);
        auto bra = dag(prime(ket,"Site"));
        auto Sxjop = 0.5*(op(sites,"S+",j)+op(sites,"S-",j));
        auto Syjop = -0.5*Cplx_i*(op(sites,"S+",j)-op(sites,"S-",j));
        auto Szjop = op(sites,"Sz",j);
        auto sxj = elt(bra*Sxjop*ket);
        auto syj = elt((bra*Syjop*ket).takeReal());
        auto szj = elt(bra*Szjop*ket);
        printfln("# mag %4d %4d %4d %.10f %.10f %.10f",j,(j-1)%2,(j-1)/2,sxj,syj,szj);
        }
    println("\n");

    // Measure S.S on every bond
    println("## dimer_top b1 b2 num_from_left SS SzSz SpSm SmSp");
    for(int j = 1; j < Nx; ++j)
        {
        auto b1 = 2*(j-1) + 1;
        auto b2 = b1 + 2;
//        auto b2 = 2*j + 1;
        auto Sp1op = op(sites,"S+",b1);
        auto Sm1op = op(sites,"S-",b1);
        auto Sz1op = op(sites,"Sz",b1);
        auto Sp2op = op(sites,"S+",b2);
        auto Sm2op = op(sites,"S-",b2);
        auto Sz2op = op(sites,"Sz",b2);
        psi.position(b1);
        auto psidag = dag(psi);
        psidag.prime("Link");
        auto l1minus1 = leftLinkIndex(psi,b1);
        auto Cpm = prime(psi(b1),l1minus1)*Sp1op;
        auto Cmp = prime(psi(b1),l1minus1)*Sm1op;
        auto Cz = prime(psi(b1),l1minus1)*Sz1op;
        Cpm *= prime(psidag(b1),"Site");
        Cmp *= prime(psidag(b1),"Site");
        Cz *= prime(psidag(b1),"Site");
        for(int k = b1+1; k < b2; ++k)
            {
            Cpm *= psi(k);
            Cmp *= psi(k);
            Cz *= psi(k);
            Cpm *= psidag(k);
            Cmp *= psidag(k);
            Cz *= psidag(k);
            }
        auto l2 = rightLinkIndex(psi,b2);
        Cpm *= prime(psi(b2),l2)*Sm2op;
        Cmp *= prime(psi(b2),l2)*Sp2op;
        Cz *= prime(psi(b2),l2)*Sz2op;
        Cpm *= prime(psidag(b2),"Site");
        Cmp *= prime(psidag(b2),"Site");
        Cz *= prime(psidag(b2),"Site");
        auto pm = 0.5*elt(Cpm);
        auto mp = 0.5*elt(Cmp);
        auto zz = elt(Cz);
        auto SdS = zz+pm+mp;
        printfln("# dimer_top %d %d %d %.10f %.10f %.10f %.10f",b1,b2,j-1,SdS,zz,pm,mp);
        }
    println("\n");

    // Measure S.S on every bond
    println("## dimer_bottom b1 b2 num_from_left SS SzSz SpSm SmSp");
    for(int j = 1; j < Nx; ++j)
        {
        auto b1 = 2*(j-1) + 2;
        auto b2 = b1 + 2;
//        auto b2 = 2*j + 2;
        auto Sp1op = op(sites,"S+",b1);
        auto Sm1op = op(sites,"S-",b1);
        auto Sz1op = op(sites,"Sz",b1);
        auto Sp2op = op(sites,"S+",b2);
        auto Sm2op = op(sites,"S-",b2);
        auto Sz2op = op(sites,"Sz",b2);
        psi.position(b1);
        auto psidag = dag(psi);
        psidag.prime("Link");
        auto l1minus1 = leftLinkIndex(psi,b1);
        auto Cpm = prime(psi(b1),l1minus1)*Sp1op;
        auto Cmp = prime(psi(b1),l1minus1)*Sm1op;
        auto Cz = prime(psi(b1),l1minus1)*Sz1op;
        Cpm *= prime(psidag(b1),"Site");
        Cmp *= prime(psidag(b1),"Site");
        Cz *= prime(psidag(b1),"Site");
        for(int k = b1+1; k < b2; ++k)
            {
            Cpm *= psi(k);
            Cmp *= psi(k);
            Cz *= psi(k);
            Cpm *= psidag(k);
            Cmp *= psidag(k);
            Cz *= psidag(k);
            }
        auto l2 = rightLinkIndex(psi,b2);
        Cpm *= prime(psi(b2),l2)*Sm2op;
        Cmp *= prime(psi(b2),l2)*Sp2op;
        Cz *= prime(psi(b2),l2)*Sz2op;
        Cpm *= prime(psidag(b2),"Site");
        Cmp *= prime(psidag(b2),"Site");
        Cz *= prime(psidag(b2),"Site");
        auto pm = 0.5*elt(Cpm);
        auto mp = 0.5*elt(Cmp);
        auto zz = elt(Cz);
        auto SdS = zz+pm+mp;
        printfln("# dimer_bottom %d %d %d %.10f %.10f %.10f %.10f",b1,b2,j-1,SdS,zz,pm,mp);
        }
    println("\n");

    // Measure S.S on every bond
    println("## dimer_rung b1 b2 num_from_left SS SzSz SpSm SmSp");
    for(int j = 1; j <= Nx; ++j)
        {
//        auto b1 = 2*(j-1) + 1;
//        auto b2 = b1 + 1;
//        auto b2 = 2*(j-1) + 2;
        auto b = 2*(j-1) + 1;
        psi.position(b);
        auto bondket = psi(b)*psi(b+1);
        auto bondbra = dag(prime(bondket,"Site"));
        auto zzop = op(sites,"Sz",b)*op(sites,"Sz",b+1);
        auto pmop = 0.5*op(sites,"S+",b)*op(sites,"S-",b+1);
        auto mpop = 0.5*op(sites,"S-",b)*op(sites,"S+",b+1);
        auto zz = elt(bondbra*zzop*bondket);
        auto pm = elt(bondbra*pmop*bondket);
        auto mp = elt(bondbra*mpop*bondket);
        auto SdS = zz+pm+mp;
//        printfln("# dimer_rung %d %d %d %.10f",b1,b2,j-1,SdS);
//        printfln("# dimer_rung %d %d %d %.10f",b,b+1,j-1,SdS);
        printfln("# dimer_rung %d %d %d %.10f %.10f %.10f %.10f",b,b+1,j-1,SdS,pm,mp,zz);
        }
    println("\n");


    return 0;
    }
