(* ::Package:: *)

FermiDirac[0] = UnitStep[-#]&;
FermiDirac[temperature_] := Module[{\[Beta]=1/temperature}, 1/(Exp[\[Beta] #]+1)&];


(*Construct a new 2D Bravais lattice of size N1, N2 with lattice constants a1, a2 on a torus*)
NewLattice[a1_List,a2_List,N1_Integer,N2_Integer] :=
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


(*InitializeTMD[N1_Integer, N2_Integer, \[Beta]Ising_, mz_, \[Mu]_, \[Alpha]_, U_, T_] := *)
InitializeTMD[N1_Integer, N2_Integer, OptionsPattern[{\[Beta]Ising->0, mz->0, \[Mu]->0, \[Alpha]->0, U->0, T->0}]] := 

Module[{\[Beta]Ising=OptionValue[\[Beta]Ising],
        mz=OptionValue[mz],
        \[Mu]=OptionValue[\[Mu]],
        \[Alpha]=OptionValue[\[Alpha]],
        U=OptionValue[U],
        T=OptionValue[T],
        \[Sigma], \[Sigma]\[Sigma],
        Up=1, Dn=2, A=1, B=2,
        a1, a2, a3, b1, b2, b3, tmdLattice,
        fermi,
        numNambu=2, numSpin=2, numOrbital=2},

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
    T1 = -(Exp[-I k.a1] + Exp[-I k.a2] + Exp[-I k.a3]);
    T1c = Conjugate[T1];
    T2 = 2\[Beta]Ising * (Sin[b1.k]+Sin[b2.k]+Sin[b3.k]);
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
      hKinetic1 = TMDHamNonInt[kx, ky]-U/2 IdentityMatrix[numSpin * numOrbital] + \[CapitalGamma] * IdentityMatrix[numSpin * numOrbital];
      hKinetic2 = TMDHamNonInt[-kx, -ky]-U/2 IdentityMatrix[numSpin * numOrbital] + \[CapitalGamma] * IdentityMatrix[numSpin * numOrbital];
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
    (*
    Do[If[Abs[t[i,j,k,l]]>10^-8,Print[StringForm["t(``,``,``,``)=``",i,j,k,l,t[i,j,k,l]]]],{i,1,2},{j,1,2},{k,1,2},{l,1,2}];
    Print["t1=", t[Up,A,Dn,A]];Print["t2=", t[Up,B,Dn,B]];
    *)
    \[CapitalGamma]part = U ( \[Rho][Up,A,Up,A]+\[Rho][Up,B,Up,B]+\[Rho][Dn,A,Dn,A]+\[Rho][Dn,B,Dn,B])/4; (* Note: U can be replaced by a function of (kx, ky) *)
    \[CapitalDelta]part = U (t[Up,A,Dn,A] + t[Up,B,Dn,B])/2;
    {\[CapitalGamma]part, \[CapitalDelta]part}
  ]
];

(* Given \[CapitalGamma], \[CapitalDelta],
   scan over BZ,
   compute contributionos of new \[CapitalGamma] and new \[CapitalDelta] at each (kx, ky),
   and return their means (over BZ) *)
CollectMF = Function[{\[CapitalGamma],\[CapitalDelta]},
  Module[{hamMF= TMDHamMF[\[CapitalGamma],\[CapitalDelta]];},
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

] (* InitializeTMD *)





(* Example *)

currentValue={0.5, 0.1};
result = Table[
	InitializeTMD[16, 16, \[Alpha]-> 0.0, U -> -4.0, T->temperature];
	{temperature, currentValue=Nest[CollectMF@@#&, currentValue, 200]},
	{temperature, {0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0}}
	]
ListPlot[ Map[{#[[1]],Abs[#[[2]][[2]]]}&, result] ]



