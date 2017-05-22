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

TMDHamNonInt::usage = "";

TMDHamMF::usage = "";
ComputeMF::usage = "";
CollectMF::usage = "";

TMDHamMF2::usage = "";
ComputeMF2::usage = "";
CollectMF2::usage = "";


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
                                                       Temperature->0}]] := 
Module[{\[Beta]Ising=OptionValue[IsingSpinOrbitCoupling],
        mz=OptionValue[ChargeTransferEnergy],
        \[Mu]=OptionValue[ChemicalPotential],
        \[Alpha]=OptionValue[RashbaSpinOrbitCoupling],
        U=OptionValue[OnSiteInteraction],
        V=OptionValue[NearestNeighborInteraction],
        T=OptionValue[Temperature],
        \[Sigma], \[Sigma]\[Sigma],
        Up=1, Dn=2, A=1, B=2,
        a1, a2, a3, b1, b2, b3, tmdLattice,
        fermi,
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

(* Lattice constants *)
a1 = { Sqrt[3]/2, 1/2};
a2 = {-Sqrt[3]/2, 1/2};
a3 = { 0, -1};
b1 = { Sqrt[3], 0};
b2 = {-Sqrt[3]/2,  3/2};
b3 = {-Sqrt[3]/2, -3/2};
tmdLattice = NewLattice[b1, b2, N1, N2];

TMDHamNonInt = Function[{kx, ky},
  Module[{k={kx, ky}, T1, T1c, T2, \[Lambda]},
    T1 = -(Exp[I k.a1] + Exp[I k.a2] + Exp[I k.a3]);
    T1c = Conjugate[T1];
    T2 = 2\[Beta]Ising * (Sin[b1.k] + Sin[b2.k] + Sin[b3.k]);
    \[Lambda] = 2\[Alpha] * (Sin[a1.k] + Sin[a2.k] + Sin[a3.k]);
    {{-\[Mu]+mz-T2,        T1,         0,         \[Lambda]},
     {      T1c, -\[Mu]-mz+T2,        -\[Lambda],         0},
     {        0,        -\[Lambda], -\[Mu]+mz+T2,        T1},
     {        \[Lambda],         0,       T1c, -\[Mu]-mz-T2}}
  ]
];

TMDHamMF = Function[{\[CapitalGamma],\[CapitalDelta]},
  Function[{kx, ky},
    Module[{hKinetic1, hKinetic2, hShift, hPairing},
      hKinetic1 = TMDHamNonInt[ kx,  ky] - U/2 IdentityMatrix[numSpin * numOrbital] + \[CapitalGamma] * IdentityMatrix[numSpin * numOrbital];
      hKinetic2 = TMDHamNonInt[-kx, -ky] - U/2 IdentityMatrix[numSpin * numOrbital] + \[CapitalGamma] * IdentityMatrix[numSpin * numOrbital];
      hPairing = {{ 0, 0, \[CapitalDelta], 0},
                  { 0, 0, 0, \[CapitalDelta]},
                  {-\[CapitalDelta], 0, 0, 0},
                  { 0,-\[CapitalDelta], 0, 0}};
      KroneckerProduct[{{1,0},{0,0}}, hKinetic1]
        + KroneckerProduct[{{0,0},{0,1}}, -Transpose[hKinetic2]]
        + KroneckerProduct[{{0,1},{0,0}}, hPairing]
        + KroneckerProduct[{{0,0},{1,0}}, ConjugateTranspose[hPairing]]
    ]
  ]
];

(* Assume eigenvalues sorted in descending order :
   eigenvectors matrix (returned from Eigensystem)  = [ U  V  ;  V^*  U^*] *)
ComputeMF = Function[{kx, ky, eigenvalues, eigenvectors},
  Module[{\[Psi], e, f, u, v, \[Rho], t, \[CapitalGamma]part, \[CapitalDelta]part},
    \[Psi] = ArrayReshape[eigenvectors, {numNambu,  numSpin*numOrbital,  numNambu,  numSpin,  numOrbital}];
    e = ArrayReshape[eigenvalues, {numNambu,  numSpin*numOrbital}];
    u = \[Psi][[1, All, 1, All, All]];
    v = \[Psi][[1, All, 2, All, All]];
    f = Map[fermi, e[[1, All]]]; (* Only positive eigenvalues *)
    \[Rho][i_,j_,k_,l_] := (* \[Rho][i, j, k, l] =*) Total[u[[All, i, j]]*f*Conjugate[u[[All, k, l]]] + Conjugate[v[[All, i, j]]] * (1 - f) * v[[All, k, l]]];
    t[i_,j_,k_,l_] := (* t[i, j, k, l] =*) Total[u[[All, i, j]]*f*Conjugate[v[[All, k, l]]] + Conjugate[v[[All, i, j]]] * (1 - f) * u[[All, k, l]]];
    \[CapitalGamma]part = U (\[Rho][Up,A,Up,A] + \[Rho][Up,B,Up,B] + \[Rho][Dn,A,Dn,A] + \[Rho][Dn,B,Dn,B])/4; (* Note: U can be replaced by a function of (kx, ky) *)
    \[CapitalDelta]part = U (t[Up,A,Dn,A] + t[Up,B,Dn,B])/2;
    {\[CapitalGamma]part, \[CapitalDelta]part}
  ]
];


(* Given \[CapitalGamma], \[CapitalDelta],
   scan over BZ,
   compute contributionos of new \[CapitalGamma] and new \[CapitalDelta] at each (kx, ky),
   and return their means (over BZ) *)
CollectMF = Function[{\[CapitalGamma],\[CapitalDelta]},
  Module[{hamMF= TMDHamMF[\[CapitalGamma],\[CapitalDelta]]},
    1/(N1*N2) Sum[
      Module[{hk = hamMF@@k, eigenvalues, eigenvectors, idxPerm},
        {eigenvalues, eigenvectors} = Eigensystem[hk];
        
        (* essential step for extracting U and V correctly *)
        idxPerm = Ordering[eigenvalues, numNambu * numSpin * numOrbital, Greater];
        eigenvalues = eigenvalues[[idxPerm]];
        eigenvectors = eigenvectors[[idxPerm]];

        (*Print["eigenvalues=", eigenvalues];Print["eigenvectors=", eigenvectors];*)
        ComputeMF[Sequence@@k, eigenvalues, eigenvectors]
      ],
      {k, tmdLattice["kVectorSpan"]}
    ]
  ]
];

(* TMDHamMF2[\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn] returns a mean field Hamiltonian with mean field parameters \[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn.
\[CapitalGamma] is a real number, \[CapitalDelta]os is a 2 component array of complex numbers, and \[CapitalDelta]nn is a 3x2x2 array of complex numbers. *)
TMDHamMF2 = Function[{\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn},
  Function[{kx, ky},
    Module[{k={kx,ky}, hKinetic1, hKinetic2, hShift, hPairing, \[CapitalDelta]AB, \[CapitalDelta]BA, \[CapitalDelta]0},
      hKinetic1 = TMDHamNonInt[ kx,  ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + \[CapitalGamma] * IdentityMatrix[numSpin * numOrbital];
      hKinetic2 = TMDHamNonInt[-kx, -ky] - ((U+6V)/2) * IdentityMatrix[numSpin * numOrbital] + \[CapitalGamma] * IdentityMatrix[numSpin * numOrbital];
      Do[\[CapitalDelta]AB[\[Sigma]1, \[Sigma]2] =  \[CapitalDelta]nn[[1, \[Sigma]1, \[Sigma]2]] Exp[ I k.a1] + \[CapitalDelta]nn[[2, \[Sigma]1, \[Sigma]2]] Exp[ I k.a2] + \[CapitalDelta]nn[[3, \[Sigma]1, \[Sigma]2]] Exp[ I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
      Do[\[CapitalDelta]BA[\[Sigma]1, \[Sigma]2] = -\[CapitalDelta]nn[[1, \[Sigma]2, \[Sigma]1]] Exp[-I k.a1] - \[CapitalDelta]nn[[2, \[Sigma]2, \[Sigma]1]] Exp[-I k.a2] - \[CapitalDelta]nn[[3, \[Sigma]2, \[Sigma]1]] Exp[-I k.a3], {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
      \[CapitalDelta]0[sub_] := \[CapitalDelta]os[[sub]];
      hPairing = {{           0,  \[CapitalDelta]AB[Up,Up],       \[CapitalDelta]0[A],  \[CapitalDelta]AB[Up,Dn]},
                  {  \[CapitalDelta]BA[Up,Up],           0,  \[CapitalDelta]BA[Up,Dn],       \[CapitalDelta]0[B]},
                  {      -\[CapitalDelta]0[A],  \[CapitalDelta]AB[Dn,Up],           0,  \[CapitalDelta]AB[Dn,Dn]},
                  {  \[CapitalDelta]BA[Dn,Up],     -\[CapitalDelta]0[B],  \[CapitalDelta]BA[Dn,Dn],          0}}; (*TODO(kyungminlee): Verify *)
      ArrayFlatten[{{hKinetic1,hPairing},{ConjugateTranspose[hPairing],-Transpose[hKinetic2]}}]
    ]
  ]
];

(* vecSpinSpin=Table[{r,\[Sigma]1,\[Sigma]2},{r, {a1, a2, a3}}, {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}]; *)

(* Assume eigenvalues sorted in descending order :
   eigenvectors matrix (returned from `Eigensystem`)  = [ U  V  ;  V^*  U^*] *)
ComputeMF2 = Function[{kx, ky, eigenvalues, eigenvectors},
  Module[{k={kx,ky}, \[Psi], e, f, u, v, \[Rho], t, \[CapitalGamma]part, \[CapitalDelta]part, \[CapitalDelta]nnpart},
    (*
    \[Psi] = ArrayReshape[eigenvectors, {numNambu,  numSpin*numOrbital,  numNambu,  numSpin,  numOrbital}];
    e = ArrayReshape[eigenvalues,  {numNambu,  numSpin*numOrbital}];
    u = \[Psi][[1, All, 1, All, All]];
    v = \[Psi][[1, All, 2, All, All]];
    f = Map[fermi, e[[1, All]]]; (* Only positive eigenvalues *)
    \[Rho][s1_,l1_,s2_,l2_] :=                  Total[u[[All, s1, l1]]*f*Conjugate[u[[All, s2, l2]]]] + Total[Conjugate[v[[All, s1, l1]]] * (1 - f) * v[[All, s2, l2]]];
    t[s1_,l1_,s2_,l2_] := t[s1,l1,s2,l2] = Total[u[[All, s1, l1]]*f*Conjugate[v[[All, s2, l2]]]] + Total[Conjugate[v[[All, s1, l1]]] * (1 - f) * u[[All, s2, l2]]];
    *)
    \[Psi] = ArrayReshape[eigenvectors, {numNambu*numSpin*numOrbital,  numNambu,  numSpin,  numOrbital}];
    e = ArrayReshape[eigenvalues,  {numNambu*numSpin*numOrbital}];
    u = \[Psi][[All, 1, All, All]];
    v = \[Psi][[All, 2, All, All]];
    f = fermi /@ e;
    (* TODO(kmlee): check Sign of I k.r in Exp *)
    \[Rho][s1_,l1_,s2_,l2_,r_] := Total[ u[[All, s1, l1]] * (1 - f) * Conjugate[ u[[All, s2, l2]] ] ] * Exp[I k.r];
    t[s1_,l1_,s2_,l2_,r_] := Total[ Conjugate[ v[[All, s1, l1]] ] * (1 - f) * u[[All, s2, l2]] ] * Exp[I k.r];
                                       
    (* Print[StringForm["\[Rho]=``",(\[Rho][Up,A,Up,A,{0,0}] + \[Rho][Up,B,Up,B,{0,0}] + \[Rho][Dn,A,Dn,A,{0,0}] + \[Rho][Dn,B,Dn,B,{0,0}])]]; *)
    \[CapitalGamma]part = ((U + 6 V)/4) * (\[Rho][Up,A,Up,A,{0,0}] + \[Rho][Up,B,Up,B,{0,0}] + \[Rho][Dn,A,Dn,A,{0,0}] + \[Rho][Dn,B,Dn,B,{0,0}]); (* Note: U can be replaced by a function of (kx, ky) *)
    \[CapitalDelta]part = (U/2) * {t[Up,A,Dn,A,{0,0}] - t[Dn,A,Up,A,{0,0}], t[Up,B,Dn,B,{0,0}] - t[Dn,B,Up,B,{0,0}]};
    \[CapitalDelta]nnpart = (V/2) * Table[t[\[Sigma]1, A, \[Sigma]2, B, r] - t[\[Sigma]2, B, \[Sigma]1, A, -r], {r, {a1, a2, a3}}, {\[Sigma]1, {Up, Dn}}, {\[Sigma]2, {Up, Dn}}];
    {\[CapitalGamma]part, \[CapitalDelta]part, \[CapitalDelta]nnpart}
  ]
];


(* Given \[CapitalGamma], \[CapitalDelta],
   scan over BZ,
   compute contributionos of new \[CapitalGamma] and new \[CapitalDelta] at each (kx, ky),
   and return their means (over BZ) *)
CollectMF2 = Function[{\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn},
  Module[{hamMF= TMDHamMF2[\[CapitalGamma], \[CapitalDelta]os, \[CapitalDelta]nn]},
    1/(N1*N2) Sum[
      Module[{hk=hamMF@@k, eigenvalues, eigenvectors},
        {eigenvalues, eigenvectors} = Eigensystem[hk];
        ComputeMF2[Sequence@@k, eigenvalues, eigenvectors]
      ],
      {k, N[tmdLattice["kVectorSpan"]]}
    ]
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
