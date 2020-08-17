#include "itensor/all.h"

using namespace itensor;

int main(int argc, char* argv[])
    {
    if(argc < 2) { printfln("Usage: %s inputfile_dmrg_table",argv[0]); return 0; }
    auto input = InputGroup(argv[1],"input");

    auto Nleg = input.getInt("Nleg");
    auto Jleg = input.getReal("Jleg");
    auto Jrung = input.getReal("Jrung");
    auto Dleg = input.getReal("Dleg");
    auto Drung = input.getReal("Drung");
    auto J4 = input.getReal("J4");
    auto nsweeps = input.getInt("nsweeps");
    auto quiet = input.getYesNo("quiet",true);
    auto table = InputGroup(input,"sweeps");
    auto sweeps = Sweeps(nsweeps,table);
    auto N = 2*Nleg;
    println(sweeps);

    println("\n");
    printfln("# Nleg %d",Nleg);
    printfln("# N %d",N);
    printfln("# Jleg %.16f",Jleg);
    printfln("# Jrung %.16f",Jrung);
    printfln("# Dleg %.16f",Dleg);
    printfln("# Drung %.16f",Drung);
    printfln("# J4 %.16f",J4);
    println("\n");

//    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    auto sites = SpinHalf(N,{"ConserveQNs=",true});
    auto ampo = AutoMPO(sites);
    // rung
    for(int j = 1; j <= N-1; j += 2)
        {
        ampo +=   Jrung*Drung,"Sz",j,"Sz",j+1;
        ampo += Jrung/2,"S+",j,"S-",j+1;
        ampo += Jrung/2,"S-",j,"S+",j+1;
        }
    // leg
    for(int j = 1; j <= N-2; ++j)
        {
        ampo +=   Jleg*Dleg,"Sz",j,"Sz",j+2;
        ampo += Jleg/2,"S+",j,"S-",j+2;
        ampo += Jleg/2,"S-",j,"S+",j+2;
        }
    // J4
    for(int j = 1; j <= N-3; j += 2)
        {
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
    auto H = toMPO(ampo);

    auto state = InitState(sites);
    for(int i = 1; i <= N; ++i)
        {
//        if(i%2 == 1) state.set(i,"Up");
//        else         state.set(i,"Dn");
        if(i%4 == 1 || i%4 == 0) state.set(i,"Up");
        else         state.set(i,"Dn");
        }
    auto psi0 = MPS(state);
//    auto psi0 = randomMPS(sites);

    println("\n");
    println("## show initial mag\n");
    for(int j = 1; j <= N; ++j)
        {
        psi0.position(j);
        auto wf = psi0.A(j);
        auto mz = (dag(prime(wf,"Site")) * sites.op("Sz",j) * wf).real();
        printfln("# initial mag %4d %.16f",j,mz);
        }
    println("\n");

    auto [energy,psi] = dmrg(H,psi0,sweeps,{"Quiet",quiet});

    println("\n");
    auto avemz = 0.0;
    auto avestagmz = 0.0;
    auto sign = 1.0;
    for(int j = 1; j <= N; ++j)
        {
        psi.position(j);
        auto wf = psi.A(j);
        auto mz = (dag(prime(wf,"Site")) * sites.op("Sz",j) * wf).real();
        avemz += mz;
        if(j%4 == 1 || j%4 == 0) sign = 1.0;
        else sign = -1.0;
        avestagmz += mz*sign;
        printfln("# mag %4d %.16f",j,mz);
        }
    println("\n");

    // Measure S.S on every bond
    println("## dimer_top b1 b2 num_from_left SS SzSz SpSm SmSp");
    for(int j = 1; j < Nleg; ++j)
        {
        auto b1 = 2*(j-1) + 1;
        auto b2 = b1 + 2;
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
    for(int j = 1; j < Nleg; ++j)
        {
        auto b1 = 2*(j-1) + 2;
        auto b2 = b1 + 2;
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
    for(int j = 1; j <= Nleg; ++j)
        {
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
        printfln("# dimer_rung %d %d %d %.16f %.16f %.16f %.16f",b,b+1,j-1,SdS,zz,pm,mp);
        }
    println("\n");

    println("\n");
    printfln("# N,ene,ene/N,mag/N,stagmag/N %4d %.16f %.16f %.16f %.16f",N,energy,energy/N,avemz/N,avestagmz/N);
    println("\n");

    return 0;
    }
