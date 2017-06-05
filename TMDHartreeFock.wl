(* ::Package:: *)

BeginPackage["TMDHartreeFock`"];


FermiDirac::usage = "FermiDirac[T] returns a Fermi-Dirac distribution function for temperature T.";
CompiledFermiDirac::usage = "FermiDirac[T] returns a compiled version of Fermi-Dirac distribution function for temperature T.";

NewLattice::usage = "NewLattice[a1, a2, N1, N2] returns an Association of properties of a lattice on a torus of size N1, N2 with lattice constants a1, a2.";

InitializeTMD::usage = "Define some functions.
InitializeTMD[N1_Integer, N2_Integer, OptionsPattern[{IsingSpinOrbitCoupling->0, 
                                                       ChargeTransferEnergy->0,
                                                       ChemicalPotential->0,
                                                       RashbaSpinOrbitCoupling->0,
                                                       OnSiteInteraction->0,
                                                       NearestNeighborInteraction->0,
                                                       Temperature->0}]]";
TMDLattice::usage="";
TMDHamNonInt::usage = "";
TMDHamMF::usage = "";
ComputeMF::usage = "";
CollectMF::usage = "";

CollectMF2::usage = "";


CompiledTMDHamNonInt::usage = "";
CompiledTMDHamMF::usage = "";
CompiledComputeMF::usage = "";
CompiledCollectMF::usage = "";



Begin["`Private`"];


FermiDirac[0,OptionsPattern[{Epsilon->$MachineEpsilon}]] = With[{\[Epsilon]=OptionValue[Epsilon]},Function[{x},If[Abs[x]<\[Epsilon], 1/2, UnitStep[-x]]]];
FermiDirac[temperature_] := With[{\[Beta]=1/temperature}, 1/(Exp[\[Beta] #]+1)&];


(* Compiled Version *)
CompiledFermiDirac[0] = Compile[{{x,_Real}},If[Abs[x]<1.0*10^-12, 0.5, UnitStep[-x]]];
CompiledFermiDirac[temperature_] := Module[{\[Beta]=1/temperature}, Compile[{{e,_Real}},1.0/(Exp[\[Beta] e]+1.0)]];



NewLattice[a1_List, a2_List, N1_Integer, N2_Integer] :=
  Module[{G1, G2, kIndexSpan, kVectorSpan},
    {G1,G2} = 2\[Pi] Transpose[Inverse[{a1,a2}]];
    kIndexSpan = Flatten[Table[{i1, i2},{i1, 0, N1-1},{i2, 0, N2-1}], 1];
    kVectorSpan = Map[(#.{G1/N1, G2/N2})&, kIndexSpan];
    <|
      "a1"->a1, "a2"->a2,
      "G1"->G1, "G2"->G2,
      "LatticeVectors"->{a1, a2},
      "ReciprocalLatticeVectors"->{G1, G2},
      "N1"->N1, "N2"->N2,
      "kIndexSpan"->kIndexSpan,
      "kVectorSpan"->kVectorSpan
    |>
  ]


(*InitializeTMD[N1_Integer, N2_Integer, \[Beta]Ising_, mz_, \[Mu]_, \[Alpha]_, U_, V_, T_] := *)
InitializeTMD[N1_Integer, N2_Integer, OptionsPattern[{IsingSpinOrbitCoupling->0, 
                                                       ChargeTransferEnergy->0,
                                                       ChemicalPotential->0,
                                                       RashbaSpinOrbitCoupling->0,
                                                       OnSiteInteraction->0,
                                                       NearestNeighborInteraction->0,
                                                       Temperature->0,
                                                       FindHartreeShift->True,
                                                       FindPairing->True}]] := 
Module[{\[Beta]Ising=OptionValue[IsingSpinOrbitCoupling],
        mz=OptionValue[ChargeTransferEnergy],
        \[Mu]=OptionValue[ChemicalPotential],
        \[Alpha]=OptionValue[RashbaSpinOrbitCoupling],
        U=OptionValue[OnSiteInteraction],
        V=OptionValue[NearestNeighborInteraction],
        T=OptionValue[Temperature],
        Up=1, Dn=2, A=1, B=2,
        \[Sigma], \[Sigma]\[Sigma],
        a1, a2, a3, b1, b2, b3, tmdLattice,
        fermi,compiledFermi,
        numNambu=2, numSpin=2, numOrbital=2,
        vecSpinSpin},

(* Use \[Sigma] and \[Sigma]\[Sigma] For Pauli and Double Pauli *)
\[Sigma]["+"] = {{0,1},{0,0}};
\[Sigma]["-"] = {{0,0},{1,0}};
\[Sigma][i_] := \[Sigma][i] = PauliMatrix[i];
Do[\[Sigma]\[Sigma][i,j] = KroneckerProduct[\[Sigma][i],\[Sigma][j]], 
  {i, {"+","-",0,1,2,3,4}},
  {j, {"+","-",0,1,2,3,4}}
];

(* `fermi` is the Fermi-Dirac distribution function for temperature T *)
fermi = FermiDirac[T];
compiledFermi = CompiledFermiDirac[T];

(* Lattice constants *)
a1 = { Sqrt[3]/2, 1/2};
a2 = {-Sqrt[3]/2, 1/2};
a3 = { 0, -1};
b1 = { Sqrt[3], 0};
b2 = {-Sqrt[3]/2,  3/2};
b3 = {-Sqrt[3]/2, -3/2};
tmdLattice = NewLattice[b1, b2, N1, N2];
TMDLattice[]=tmdLattice;

TMDHamNonInt = Function[{kx, ky},
  Module[{k={kx, ky}, T1, T1c, T2, \[Lambda]p, \[Lambda]m, \[Lambda]pc, \[Lambda]mc},
    T1 = -(Exp[I k.a1] + Exp[I k.a2] + Exp[I k.a3]);
    T1c = Conjugate[T1];
    T2 = 2\[Beta]Ising * (Sin[b1.k] + Sin[b2.k] + Sin[b3.k]);
    \[Lambda]p = \[Alpha] * Sum[(r[[2]]+I r[[1]])Exp[I k.r],{r,{ a1, a2, a3}}];
    \[Lambda]m = \[Alpha] * Sum[(r[[2]]+I r[[1]])Exp[I k.r],{r,{-a1,-a2,-a3}}];
    \[Lambda]pc = Conjugate[\[Lambda]p];
    \[Lambda]mc = Conjugate[\[Lambda]m];
    {{-\[Mu]+mz-T2,        T1,         0,        \[Lambda]p},
     {      T1c, -\[Mu]-mz+T2,        \[Lambda]m,         0},
     {        0,       \[Lambda]mc, -\[Mu]+mz+T2,        T1},
     {      \[Lambda]pc,         0,       T1c, -\[Mu]-mz-T2}}
  ]
];


CompiledTMDHamNonInt = Compile[{{kx,_Real}, {ky,_Real}},
  Module[{k={kx, ky}, T1, T1c, T2, \[Lambda]p, \[Lambda]m, \[Lambda]pc, \[Lambda]mc},
    T1 = -(Exp[I k.a1] + Exp[I k.a2] + Exp[I k.a3]);
    T1c = Conjugate[T1];
    T2 = 2\[Beta]Ising * (Sin[b1.k] + Sin[b2.k] + Sin[b3.k]);
    \[Lambda]p = \[Alpha] * Sum[(r[[2]]+I r[[1]])Exp[I k.r],{r,{ a1, a2, a3}}];
    \[Lambda]m = \[Alpha] * Sum[(r[[2]]+I r[[1]])Exp[I k.r],{r,{-a1,-a2,-a3}}];
    \[Lambda]pc = Conjugate[\[Lambda]p];
    \[Lambda]mc = Conjugate[\[Lambda]m];
    {{-\[Mu]+mz-T2,        T1,         0,        \[Lambda]p},
     {      T1c, -\[Mu]-mz+T2,        \[Lambda]m,         0},
     {        0,       \[Lambda]mc, -\[Mu]+mz+T2,        T1},
     {      \[Lambda]pc,         0,       T1c, -\[Mu]-mz-T2}}
  ]
];

If[OptionValue[FindHartreeShift],
  TMDKineticWithHartreeShift[\[CapitalGamma]_] := Function[{kx,ky},
    TMDHamNonInt[ kx,  ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + DiagonalMatrix[\[CapitalGamma]]
  ],
  TMDKineticWithHartreeShift[\[CapitalGamma]_] := Function[{kx,ky},
    TMDHamNonInt[ kx,  ky]
  ]
];

If[OptionValue[FindPairing],
  TMDPairingMatrix[\[CapitalDelta]os_, \[CapitalDelta]nn_] := Function[{kx,ky},
    Module[{k={kx,ky},\[CapitalDelta]AB,\[CapitalDelta]BA,\[CapitalDelta]0, zero=0.0+0.0I},
      Do[\[CapitalDelta]AB[\[Sigma]1, \[Sigma]2] =  \[CapitalDelta]nn[[1, \[Sigma]1, \[Sigma]2]] Exp[ I k.a1] + \[CapitalDelta]nn[[2, \[Sigma]1, \[Sigma]2]] Exp[ I k.a2] + \[CapitalDelta]nn[[3, \[Sigma]1, \[Sigma]2]] Exp[ I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
      Do[\[CapitalDelta]BA[\[Sigma]1, \[Sigma]2] = -\[CapitalDelta]nn[[1, \[Sigma]2, \[Sigma]1]] Exp[-I k.a1] - \[CapitalDelta]nn[[2, \[Sigma]2, \[Sigma]1]] Exp[-I k.a2] - \[CapitalDelta]nn[[3, \[Sigma]2, \[Sigma]1]] Exp[-I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
      \[CapitalDelta]0[sub_] := \[CapitalDelta]os[[sub]];
      {{        zero,  \[CapitalDelta]AB[Up,Up],       \[CapitalDelta]0[A],  \[CapitalDelta]AB[Up,Dn]},
       {  \[CapitalDelta]BA[Up,Up],        zero,  \[CapitalDelta]BA[Up,Dn],       \[CapitalDelta]0[B]},
       {      -\[CapitalDelta]0[A],  \[CapitalDelta]AB[Dn,Up],        zero,  \[CapitalDelta]AB[Dn,Dn]},
       {  \[CapitalDelta]BA[Dn,Up],      -\[CapitalDelta]0[B],  \[CapitalDelta]BA[Dn,Dn],        zero}}
    ]
  ],
  TMDPairingMatrix[\[CapitalDelta]os_, \[CapitalDelta]nn_] := Function[{kx,ky},ConstantArray[0.0+0.0I, {4,4}]]
];

(* TMDHamMF[\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn] returns a mean field Hamiltonian with mean field parameters \[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn.
\[CapitalGamma] is a real number, \[CapitalDelta]os is a 2 component array of complex numbers, and \[CapitalDelta]nn is a 3x2x2 array of complex numbers. *)
TMDHamMF = Function[{\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn},
  With[{hamNonIntShift=TMDKineticWithHartreeShift[\[CapitalGamma]], hamPairing=TMDPairingMatrix[\[CapitalDelta]os,\[CapitalDelta]nn]},
    Function[{kx, ky},
      Module[{k={kx,ky}, hKinetic1, hKinetic2, hShift, hPairing, \[CapitalDelta]AB, \[CapitalDelta]BA, \[CapitalDelta]0},
        (*
        hKinetic1 = TMDHamNonInt[ kx,  ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + DiagonalMatrix[\[CapitalGamma]];
        hKinetic2 = TMDHamNonInt[-kx, -ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + DiagonalMatrix[\[CapitalGamma]];
      
        Do[\[CapitalDelta]AB[\[Sigma]1, \[Sigma]2] =  \[CapitalDelta]nn[[1, \[Sigma]1, \[Sigma]2]] Exp[ I k.a1] + \[CapitalDelta]nn[[2, \[Sigma]1, \[Sigma]2]] Exp[ I k.a2] + \[CapitalDelta]nn[[3, \[Sigma]1, \[Sigma]2]] Exp[ I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
        Do[\[CapitalDelta]BA[\[Sigma]1, \[Sigma]2] = -\[CapitalDelta]nn[[1, \[Sigma]2, \[Sigma]1]] Exp[-I k.a1] - \[CapitalDelta]nn[[2, \[Sigma]2, \[Sigma]1]] Exp[-I k.a2] - \[CapitalDelta]nn[[3, \[Sigma]2, \[Sigma]1]] Exp[-I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
        \[CapitalDelta]0[sub_] := \[CapitalDelta]os[[sub]];
        hPairing = {{           0,  \[CapitalDelta]AB[Up,Up],       \[CapitalDelta]0[A],  \[CapitalDelta]AB[Up,Dn]},
                    {  \[CapitalDelta]BA[Up,Up],           0,  \[CapitalDelta]BA[Up,Dn],       \[CapitalDelta]0[B]},
                    {      -\[CapitalDelta]0[A],  \[CapitalDelta]AB[Dn,Up],           0,  \[CapitalDelta]AB[Dn,Dn]},
                    {  \[CapitalDelta]BA[Dn,Up],     -\[CapitalDelta]0[B],  \[CapitalDelta]BA[Dn,Dn],            0}}; (*TODO(kyungminlee): Verify *)
        *)
        hKinetic1=hamNonIntShift[kx,ky];
        hKinetic2=hamNonIntShift[-kx,-ky];
        hPairing=hamPairing[kx,ky];
        ArrayFlatten[{{hKinetic1,hPairing},{ConjugateTranspose[hPairing],-Transpose[hKinetic2]}}]
      ]
    ]
  ]
];


CompiledTMDHamMF = Function[{\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn},
  Compile[{{kx,_Real}, {ky,_Real}},
    Module[{
        k={kx,ky},
        hKinetic1 = CompiledTMDHamNonInt[ kx,  ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + DiagonalMatrix[\[CapitalGamma]],
        hKinetic2 = CompiledTMDHamNonInt[-kx, -ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + DiagonalMatrix[\[CapitalGamma]],
        hPairing=Table[0.0+0.0I,{i,1,4},{j,1,4}],
        hOut=Table[0.0+0.0I,{i,1,8},{j,1,8}],
        \[CapitalDelta]AB={{0.0+0.0 I, 0.0+0.0 I},{0.0+0.0 I, 0.0+0.0 I}},
        \[CapitalDelta]BA={{0.0+0.0 I, 0.0+0.0 I},{0.0+0.0 I, 0.0+0.0 I}}
      },
      Do[\[CapitalDelta]AB[[\[Sigma]1, \[Sigma]2]] =  \[CapitalDelta]nn[[1, \[Sigma]1, \[Sigma]2]] Exp[ I k.a1] + \[CapitalDelta]nn[[2, \[Sigma]1, \[Sigma]2]] Exp[ I k.a2] + \[CapitalDelta]nn[[3, \[Sigma]1, \[Sigma]2]] Exp[ I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
      Do[\[CapitalDelta]BA[[\[Sigma]1, \[Sigma]2]] = -\[CapitalDelta]nn[[1, \[Sigma]2, \[Sigma]1]] Exp[-I k.a1] - \[CapitalDelta]nn[[2, \[Sigma]2, \[Sigma]1]] Exp[-I k.a2] - \[CapitalDelta]nn[[3, \[Sigma]2, \[Sigma]1]] Exp[-I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
      hPairing = {{     0.0+0.0 I,  \[CapitalDelta]AB[[Up,Up]],      \[CapitalDelta]os[[A]],  \[CapitalDelta]AB[[Up,Dn]]},
                  {  \[CapitalDelta]BA[[Up,Up]],     0.0+0.0 I,  \[CapitalDelta]BA[[Up,Dn]],      \[CapitalDelta]os[[B]]},
                  {     -\[CapitalDelta]os[[A]],  \[CapitalDelta]AB[[Dn,Up]],     0.0+0.0 I,  \[CapitalDelta]AB[[Dn,Dn]]},
                  {  \[CapitalDelta]BA[[Dn,Up]],    -\[CapitalDelta]os[[B]],  \[CapitalDelta]BA[[Dn,Dn]],    0.0+0.0 I}};
      hOut[[1;;4, 1;;4]]=hKinetic1;
      hOut[[1;;4, 5;;8]]=hPairing;
      hOut[[5;;8, 1;;4]]=ConjugateTranspose[hPairing];
      hOut[[5;;8, 5;;8]]=-Transpose[hKinetic2];
      hOut
    ]
  ]
];

(* vecSpinSpin=Table[{r,\[Sigma]1,\[Sigma]2},{r, {a1, a2, a3}}, {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}]; *)

(* Assume eigenvalues sorted in descending order :
   eigenvectors matrix (returned from `Eigensystem`)  = [ U  V  ;  V^*  U^*] *)
   (*Return "average" density AND mean fields *)
ComputeMF = Function[{kx, ky, eigenvalues, eigenvectors},
  Module[{k={kx,ky}, \[Psi], e, f, u, v, \[Rho], t, \[CapitalGamma]part, \[CapitalDelta]part, \[CapitalDelta]nnpart},
    \[Psi] = ArrayReshape[eigenvectors, {numNambu*numSpin*numOrbital,  numNambu,  numSpin,  numOrbital}];
    e = ArrayReshape[eigenvalues,  {numNambu*numSpin*numOrbital}];
    u = \[Psi][[All, 1, All, All]];
    v = \[Psi][[All, 2, All, All]];
    f = fermi /@ e;
    (* TODO(kmlee): check Sign of I k.r in Exp *)
    \[Rho][s1_,l1_,s2_,l2_] := \[Rho][s1,l1,s2,l2] = Total[ u[[All, s1, l1]] * f * Conjugate[ u[[All, s2, l2]] ] ];
    t[s1_,l1_,s2_,l2_] := t[s1,l1,s2,l2] = Total[ u[[All, s1, l1]] * f * Conjugate[ v[[All, s2, l2]] ] ];

    \[Rho][s1_,l1_,s2_,l2_,r_] := \[Rho][s1,l1,s2,l2] * Exp[-I k.r];
    t[s1_,l1_,s2_,l2_,r_] := t[s1,l1,s2,l2] * Exp[-I k.r];

	\[CapitalGamma]part = {U*\[Rho][Dn,A,Dn,A] + 3*V*(\[Rho][Up,B,Up,B] + \[Rho][Dn,B,Dn,B]),
             U*\[Rho][Dn,B,Dn,B] + 3*V*(\[Rho][Up,A,Up,A] + \[Rho][Dn,A,Dn,A]),
             U*\[Rho][Up,A,Up,A] + 3*V*(\[Rho][Up,B,Up,B] + \[Rho][Dn,B,Dn,B]),
             U*\[Rho][Up,B,Up,B] + 3*V*(\[Rho][Up,A,Up,A] + \[Rho][Dn,A,Dn,A])};
    \[CapitalDelta]part = (U/2) * {t[Up,A,Dn,A] - t[Dn,A,Up,A], t[Up,B,Dn,B] - t[Dn,B,Up,B]};
    \[CapitalDelta]nnpart = (V/2) * Table[t[\[Sigma]1, A, \[Sigma]2, B, r] - t[\[Sigma]2, B, \[Sigma]1, A, -r], {r, {a1, a2, a3}}, {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
    {0.25*(\[Rho][Up,A,Up,A]+\[Rho][Up,B,Up,B]+\[Rho][Dn,A,Dn,A]+\[Rho][Dn,B,Dn,B]),\[CapitalGamma]part, \[CapitalDelta]part, \[CapitalDelta]nnpart}
  ]
];


CompiledComputeMF = Compile[{{kx,_Real}, {ky,_Real}, {eigenvalues,_Real, 1}, {eigenvectors,_Complex, 4}},
(*CompiledComputeMF = Function[{kx, ky, eigenvalues, eigenvectors},*)
  Module[{
      k={kx,ky},
      n=Length[eigenvalues],
      \[Rho]=Table[0.0+0.0 I, {s1,1,2}, {f1,1,2}, {s2,1,2}, {f2,1,2}],
      t=Table[0.0+0.0 I, {s1,1,2}, {f1,1,2}, {s2,1,2}, {f2,1,2}],
      \[Psi]=eigenvectors,
      f=Table[0.0+0.0 I, {i,1,Length[eigenvalues]}],
      \[CapitalGamma]part=Table[0.0+0.0 I, {i,1,4}],
      \[CapitalDelta]part=Table[0.0+0.0 I, {i,1,2}],
      \[CapitalDelta]nnpart=Table[0.0+0.0 I, {i1,1,3},{i2,1,2},{i3,1,2}]
    },
	Do[f[[i]] = compiledFermi[ eigenvalues[[i]] ], {i, 1, n}];
    Do[\[Rho][[s1,f1,s2,f2]] = Sum[\[Psi][[i, 1, s1, f1]] * f[[i]] * Conjugate[ \[Psi][[i, 1, s2, f2]] ], {i,1,n}], {s1,1,2}, {f1,1,2}, {s2,1,2}, {f2,1,2}];
    Do[t[[s1,f1,s2,f2]] = Sum[\[Psi][[i, 1, s1, f1]] * f[[i]] * Conjugate[ \[Psi][[i, 2, s2, f2]] ], {i,1,n}], {s1,1,2}, {f1,1,2}, {s2,1,2}, {f2,1,2}];
    
	\[CapitalGamma]part[[1]] = U*\[Rho][[Dn,A,Dn,A]] + 3*V*(\[Rho][[Up,B,Up,B]] + \[Rho][[Dn,B,Dn,B]]);
    \[CapitalGamma]part[[2]] = U*\[Rho][[Dn,B,Dn,B]] + 3*V*(\[Rho][[Up,A,Up,A]] + \[Rho][[Dn,A,Dn,A]]);
    \[CapitalGamma]part[[3]] = U*\[Rho][[Up,A,Up,A]] + 3*V*(\[Rho][[Up,B,Up,B]] + \[Rho][[Dn,B,Dn,B]]);
    \[CapitalGamma]part[[4]] = U*\[Rho][[Up,B,Up,B]] + 3*V*(\[Rho][[Up,A,Up,A]] + \[Rho][[Dn,A,Dn,A]]);
    
    \[CapitalDelta]part[[1]] = (U/2) * (t[[Up,A,Dn,A]] - t[[Dn,A,Up,A]]);
    \[CapitalDelta]part[[2]] = (U/2) * (t[[Up,B,Dn,B]] - t[[Dn,B,Up,B]]);
    
    Do[
      \[CapitalDelta]nnpart[[1,\[Sigma]1,\[Sigma]2]] = (V/2) * (t[[\[Sigma]1, A, \[Sigma]2, B]] * Exp[-I k.a1] - t[[\[Sigma]2, B, \[Sigma]1, A]]* Exp[I k.a1]);
      \[CapitalDelta]nnpart[[2,\[Sigma]1,\[Sigma]2]] = (V/2) * (t[[\[Sigma]1, A, \[Sigma]2, B]] * Exp[-I k.a2] - t[[\[Sigma]2, B, \[Sigma]1, A]]* Exp[I k.a2]);
      \[CapitalDelta]nnpart[[3,\[Sigma]1,\[Sigma]2]] = (V/2) * (t[[\[Sigma]1, A, \[Sigma]2, B]] * Exp[-I k.a3] - t[[\[Sigma]2, B, \[Sigma]1, A]]* Exp[I k.a3]);
    ,{\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
    
    Join[\[CapitalGamma]part,\[CapitalDelta]part,Flatten[\[CapitalDelta]nnpart]]
    
  ]
];



(* Given \[CapitalGamma], \[CapitalDelta],
   scan over BZ,
   compute contributionos of new \[CapitalGamma] and new \[CapitalDelta] at each (kx, ky),
   and return their means (over BZ) *)
CollectMF = Function[{\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn},
  Module[{hamMF= TMDHamMF[\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn]},
    1/(N1*N2) Sum[
      Module[{hk=hamMF@@k, eigenvalues, eigenvectors},
        {eigenvalues, eigenvectors} = Eigensystem[hk];
        ComputeMF[Sequence@@k, eigenvalues, eigenvectors]
      ],
      {k, N[tmdLattice["kVectorSpan"]]}
    ]
  ]
];



(* Given \[CapitalGamma], \[CapitalDelta],
   scan over BZ,
   compute contributionos of new \[CapitalGamma] and new \[CapitalDelta] at each (kx, ky),
   and return their means (over BZ) *)
CompiledCollectMF = Function[{\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn},
  Module[{hamMF=CompiledTMDHamMF[\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn],out},
    out=1/(N1*N2) Sum[
      Module[{hk=hamMF@@k, eigenvalues, eigenvectors},
        {eigenvalues, eigenvectors} = Eigensystem[hk];
        CompiledComputeMF[Sequence@@k,
          Re[eigenvalues],
          ArrayReshape[eigenvectors,{numNambu*numSpin*numOrbital,  numNambu,  numSpin,  numOrbital}]]
      ],
      {k, N[tmdLattice["kVectorSpan"]]}
    ];
    {out[[1;;4]],
     out[[5;;6]],
     ArrayReshape[out[[7;;]],{3,2,2}]}
  ]
];

] (* InitializeTMD *)



End[];
EndPackage[];


(* Example *)
(*
onsiteInteraction = -3.0;
nearestNeighborInteraction = -0.2;
InitializeTMD[24, 24, OnSiteInteraction->onsiteInteraction, NearestNeighborInteraction->nearestNeighborInteraction, ChemicalPotential->0.5, Temperature->0.1];
{\[CapitalGamma]c2, \[CapitalDelta]osc, \[CapitalDelta]nnc} = {0.5*(onsiteInteraction + 6*nearestNeighborInteraction), RandomReal[{-0.5, 0.5}, 2], RandomReal[{-0.5, 0.5}, {3,2,2}]};
allResults = NestList[(CollectMF2@@#)&,{\[CapitalGamma]c2, \[CapitalDelta]osc, \[CapitalDelta]nnc}, 50];
results=allResults[[;;]];
ArrayReshape[ Flatten[results], {Length[results], 1+2+3*2*2}]//Chop//TableForm
*)
