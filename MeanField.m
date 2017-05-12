(* ::Package:: *)

(* ::Code::Initialization:: *)
Needs["ErrorBarPlots`"]


(* ::Code::Initialization:: *)
Fermi[En_] := HeavisideTheta[-En]


Nsites=20;
kSpan=Flatten[Table[i*4\[Pi]/3{1/Sqrt[3],1}+j*{8\[Pi]/(3Sqrt[3]),0},{i,0,0.999,1/Nsites},{j,0,0.999,1/Nsites}],1];


TMDHamNonInt[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[Alpha]_:0] :=
  With[{
      Tk = -(1+2Cos[Sqrt[3]kx/2]Exp[-I 3ky/2]),
      Tkstar = -(1+2Cos[Sqrt[3]kx/2]Exp[I 3ky/2]),
      \[Lambda]kx = \[Alpha](1-Cos[Sqrt[3]kx/2]Exp[3I ky/2]),
      \[Lambda]ky = I \[Alpha] Sqrt[3]Sin[Sqrt[3]kx/2]Exp[3I ky/2],
      \[Lambda]kxstar = \[Alpha](1-Cos[Sqrt[3]kx/2]Exp[-3I ky/2]),
      \[Lambda]kystar = -I \[Alpha] Sqrt[3]Sin[Sqrt[3]kx/2]Exp[-3I ky/2],
      \[Beta]SOC = (2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]kx]-2Sin[Sqrt[3]kx/2]Cos[3ky/2])
    },
    KroneckerProduct[Tk*{{0,1},{0,0}},IdentityMatrix[2]]
    +KroneckerProduct[Tkstar*{{0,0},{1,0}},IdentityMatrix[2]]
    +KroneckerProduct[mz*{{1,0},{0,-1}},IdentityMatrix[2]]
    +KroneckerProduct[\[Beta]SOC*{{1,0},{0,-1}},{{1,0},{0,-1}}]
    +KroneckerProduct[I \[Lambda]kx*{{0,1},{0,0}},{{0,1},{1,0}}]
    +KroneckerProduct[I \[Lambda]ky*{{0,1},{0,0}},{{0,-I},{I,0}}]
    +KroneckerProduct[-I \[Lambda]kxstar*{{0,0},{1,0}},{{0,1},{1,0}}]
    +KroneckerProduct[-I \[Lambda]kystar*{{0,0},{1,0}},{{0,-I},{I,0}}]
    -\[Mu]*IdentityMatrix[4]
  ]


TMDHamMF[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[CapitalDelta]_:0,\[CapitalDelta]x_:0,U_:0,nA_:0,nB_:0] :=
  With[{
      Tk = -(1+2Cos[Sqrt[3]kx/2]Exp[-I 3ky/2]),
      Tkstar = -(1+2Cos[Sqrt[3]kx/2]Exp[I 3ky/2]),
      \[Beta]SOC = (2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]kx]-2Sin[Sqrt[3]kx/2]Cos[3ky/2]),
      \[CapitalDelta]star = Conjugate[\[CapitalDelta]],
      \[CapitalDelta]xstar = Conjugate[\[CapitalDelta]x],
      mztilde = mz-U/4*(nA-nB),
      \[Mu]tilde = \[Mu]-U/4*(nA+nB-2)
      (*,\[CapitalDelta]=0,\[CapitalDelta]star=0,U=0*)(*,mz=0.2*)
    },
    Ek=Sqrt[(\[Beta]SOC+mztilde)^2+Tk Tkstar];
    SparseArray[{{Ek-\[Mu]tilde, \[CapitalDelta], 0, \[CapitalDelta]x},
                 {\[CapitalDelta]star, -(Ek-\[Mu]tilde), \[CapitalDelta]x, 0},
                 {0, \[CapitalDelta]xstar, -(Ek+\[Mu]tilde), \[CapitalDelta]},
                 {\[CapitalDelta]xstar, 0, \[CapitalDelta]star, Ek+\[Mu]tilde}}]
  ]


TMDHamMFDn[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[CapitalDelta]_:0,\[CapitalDelta]x_:0,U_:0,nA_:0,nB_:0] :=
  With[{
      Tk=-(1+2Cos[Sqrt[3]kx/2]Exp[-I 3ky/2]),
      Tkstar=-(1+2Cos[Sqrt[3]kx/2]Exp[I 3ky/2]),
      \[Beta]SOC=(2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]kx]-2Sin[Sqrt[3]kx/2]Cos[3ky/2]),
      \[CapitalDelta]star=Conjugate[\[CapitalDelta]],\[CapitalDelta]xstar=Conjugate[\[CapitalDelta]x],
      mztilde=mz-U/4*(nA-nB),
      \[Mu]tilde=\[Mu]-U/4*(nA+nB-2)
      (*,\[CapitalDelta]=0,\[CapitalDelta]star=0,U=0*)(*,mz=0.2*)
    },
    Ek=Sqrt[(-\[Beta]SOC+mztilde)^2+Tk Tkstar];
    SparseArray[{{Ek-\[Mu]tilde, \[CapitalDelta], 0, \[CapitalDelta]x},
                 {\[CapitalDelta]star, -(Ek-\[Mu]tilde), \[CapitalDelta]x, 0},
                 {0, \[CapitalDelta]xstar, -(Ek+\[Mu]tilde), \[CapitalDelta]},
                 {\[CapitalDelta]xstar, 0, \[CapitalDelta]star, Ek+\[Mu]tilde}}]
  ]


find\[CapitalDelta][\[Beta]Ising_,mz_,\[Mu]_,\[CapitalDelta]_,\[CapitalDelta]x_,U_,nA_,nB_] := -U/Length[kSpan]*Sum[k=Simplify[k];


With[{
    Ham=TMDHamMF[k[[1]],k[[2]],\[Beta]Ising,mz,\[Mu],\[CapitalDelta],\[CapitalDelta]x,U,nA,nB],
    sin=Simplify[Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])]/\[Sqrt](((2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]k[[1]]]-2Sin[Sqrt[3]k[[1]]/2]Cos[3k[[2]]/2])+mz)^2+Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])]^2)]
  },
  eigs=Eigensystem[Ham];
  Sum[Fermi[eigs[[1,i]]]*{eigs[[2,i]].SparseArray[{{2,1}->1,{4,3}->1},{4,4}].eigs[[2,i]],(sin)/2*eigs[[2,i]].SparseArray[{{3,2}->1,{4,1}->1},{4,4}].eigs[[2,i]]},
    {i,Length[eigs[[1]]]}]]
,{k,kSpan}]


findn[\[Beta]Ising_,mz_,\[Mu]_,\[CapitalDelta]_,\[CapitalDelta]x_,U_,nA_,nB_]:=2/(Length[kSpan])*Sum[k=Simplify[k];


With[{
    Ham=TMDHamMF[k[[1]],k[[2]],\[Beta]Ising,mz,\[Mu],\[CapitalDelta],\[CapitalDelta]x,U,nA,nB],
    phase=If[Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])]>0,-(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])/Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])],1],sin=Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])]/\[Sqrt](((2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]k[[1]]]-2Sin[Sqrt[3]k[[1]]/2]Cos[3k[[2]]/2])+mz)^2+Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])]^2),cosby2=(1+((2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]k[[1]]]-2Sin[Sqrt[3]k[[1]]/2]Cos[3k[[2]]/2])+mz)/\[Sqrt](((2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]k[[1]]]-2Sin[Sqrt[3]k[[1]]/2]Cos[3k[[2]]/2])+mz)^2+Abs[(1+2Cos[Sqrt[3]k[[1]]/2]Exp[-I 3k[[2]]/2])]^2))/2
  },
  eigs=Eigensystem[Ham];
  Sum[Check[Fermi[eigs[[1,i]]]*{eigs[[2,i]].SparseArray[{{1,1}->cosby2,{2,2}->-0cosby2,{3,3}->1-cosby2,{4,4}->0(cosby2-1),{1,3}->sin/2*phase,{2,4}->-0sin/2*Conjugate[phase],{3,1}->sin/2*Conjugate[phase],{4,2}->-0sin/2*phase},{4,4}].eigs[[2,i]],eigs[[2,i]].SparseArray[{{1,1}->1-cosby2,{2,2}->0(cosby2-1),{3,3}->cosby2,{4,4}->-0cosby2,{1,3}->-sin/2*phase,{2,4}->0sin/2*Conjugate[phase],{3,1}->-sin/2*Conjugate[phase],{4,2}->0sin/2*phase},{4,4}].eigs[[2,i]]},Print[{k,nA,nB}]],
  {i,Length[eigs[[1]]]}]]
,{k,kSpan}]


del={1.0,1.0};nanb={0.9,1.1};


Optimal\[CapitalDelta][\[Beta]Ising_,mz_,\[Mu]_,U_,nA_,nB_] := 
  FixedPoint[Re[find\[CapitalDelta][\[Beta]Ising,mz,\[Mu],#[[1]],#[[2]],U,nA,nB]]&,{del[[1]]+0.01,del[[2]]+0.01},20,SameTest->(Norm[Re[#1-#2]]<1*^-20&)]


OptimalnAnB[\[Beta]Ising_,mz_,\[Mu]_,U_] :=
  {FixedPoint[(del=Optimal\[CapitalDelta][\[Beta]Ising,mz,\[Mu],U,#[[1]],#[[2]]];nanb=Re[findn[\[Beta]Ising,mz,\[Mu],del[[1]],del[[2]],U,#[[1]],#[[2]]]])&,{nanb[[1]]-0.001,nanb[[2]]+0.001},20,SameTest->(Norm[Re[#1-#2]]<1*^-3&)],del}//Flatten


UbetaPhaseDiagram={};


Do[AppendTo[UbetaPhaseDiagram,{U,\[Beta]Ising,Check[OptimalnAnB[\[Beta]Ising,0.1,0.0,U],Continue[]]}],{U,0.01,1.01,0.1},{\[Beta]Ising,0.01,0.31,0.05}]
