(* ::Package:: *)

(* ::Text:: *)
(*50*50 lattice : \[CapitalGamma]KK'\[CapitalGamma] cut*)


(* ::Input:: *)
(*Nsites=50;kCut=Join[Table[i*2*)
(*\[Pi]/3{1/Sqrt[3],1},{i,0,1,2/Nsites}],Table[2\[Pi]/3{1/Sqrt[3],1}-i{4\[Pi]/(3Sqrt[3]),0},{i,2/Nsites,1-2/Nsites,2/Nsites}],Table[i*2\[Pi]/3{-1/Sqrt[3],1},{i,1,0,-2/Nsites}]];*)


(* ::Text:: *)
(*4*4 Non-Interacting k-space Hamiltonian *)


(* ::Input:: *)
(*TMDHamNonInt[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[Alpha]_:0]:=With[{Tk=-(1+2Cos[Sqrt[3]kx/2]Exp[-I 3ky/2]),Tkstar=-(1+2Cos[Sqrt[3]kx/2]Exp[I 3ky/2]),\[Beta]SOC=(2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]kx]-2Sin[Sqrt[3]kx/2]Cos[3ky/2])},KroneckerProduct[Tk*{{0,1},{0,0}},IdentityMatrix[2]]+KroneckerProduct[Tkstar*{{0,0},{1,0}},IdentityMatrix[2]]+KroneckerProduct[mz*{{1,0},{0,-1}},IdentityMatrix[2]]+KroneckerProduct[\[Beta]SOC*{{1,0},{0,-1}},{{1,0},{0,-1}}]-\[Mu]*IdentityMatrix[4]]*)


(* ::Text:: *)
(*8*8 Mean Field Hamiltonian in Nambu space*)


(* ::Input:: *)
(*TMDHamMF[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[Alpha]_:0,\[CapitalDelta]_:0,\[CapitalDelta]x_:0,U_:0,nA_:0,nB_:0]:=ArrayFlatten[{{TMDHamNonInt[kx,ky,\[Beta]Ising,mz-U/4*(nA-nB),\[Mu]-U/4*(nA+nB-2),\[Alpha]],SparseArray[{{1,2},{2,1},{3,4},{4,3}}->{-\[CapitalDelta],\[CapitalDelta],-\[CapitalDelta],\[CapitalDelta]}]},{SparseArray[{{1,2},{2,1},{3,4},{4,3}}->-Conjugate[{-\[CapitalDelta],\[CapitalDelta],-\[CapitalDelta],\[CapitalDelta]}]],-TMDHamNonInt[-kx,-ky,\[Beta]Ising,mz-U/4*(nA-nB),\[Mu]-U/4*(nA+nB-2),\[Alpha]]}}]*)


(* ::Input:: *)
(*With[{\[Beta]Ising=0.5,mz=0.5,mu=0.0,\[CapitalDelta]=0.2,\[CapitalDelta]x=0,U=0,nA=1.0857984307519946`,nB=0.9187476777471787`},ListPlot[Join[Table[Table[Tooltip[(Eigenvalues[TMDHamMF[i[[1]],i[[2]],\[Beta]Ising,mz,mu,0,\[CapitalDelta],\[CapitalDelta]x,U,nA,nB]]//Sort)[[n]]],{i,kCut}],{n,8}]],PlotStyle->{Red,Red,Red,Red,Blue,Blue,Blue,Blue},(*Joined\[Rule]True,*)ImageSize->Large,GridLines->Automatic,PlotRange->All(*{-0.5,0.5}*)]]*)


(* ::Text:: *)
(*Alternatively if the non-interacting Hamiltonian is first diagonalised and then a pair potential added to this Hamiltonian*)


(* ::Input:: *)
(*TMDHamMFUp[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[CapitalDelta]_:0,\[CapitalDelta]x_:0,U_:0,nA_:0,nB_:0]:=With[{Tk=-(1+2Cos[Sqrt[3]kx/2]Exp[-I 3ky/2]),Tkstar=-(1+2Cos[Sqrt[3]kx/2]Exp[I 3ky/2]),\[Beta]SOC=(2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]kx]-2Sin[Sqrt[3]kx/2]Cos[3ky/2]),\[CapitalDelta]star=Conjugate[\[CapitalDelta]],\[CapitalDelta]xstar=Conjugate[\[CapitalDelta]x],mztilde=mz-U/4*(nA-nB),\[Mu]tilde=\[Mu]-U/4*(nA+nB-2)(*,\[CapitalDelta]=0,\[CapitalDelta]star=0,U=0*)(*,mz=0.2*)},*)
(*Ek=Sqrt[(\[Beta]SOC+mztilde)^2+Tk Tkstar];SparseArray[({*)
(* {Ek-\[Mu]tilde, \[CapitalDelta], 0, \[CapitalDelta]x},*)
(* {\[CapitalDelta]star, -(Ek-\[Mu]tilde), \[CapitalDelta]x, 0},*)
(* {0, \[CapitalDelta]xstar, -(Ek+\[Mu]tilde), \[CapitalDelta]},*)
(* {\[CapitalDelta]xstar, 0, \[CapitalDelta]star, Ek+\[Mu]tilde}*)
(*})]]*)


(* ::Input:: *)
(*TMDHamMFDn[kx_,ky_,\[Beta]Ising_,mz_,\[Mu]_:0,\[CapitalDelta]_:0,\[CapitalDelta]x_:0,U_:0,nA_:0,nB_:0]:=With[{Tk=-(1+2Cos[Sqrt[3]kx/2]Exp[-I 3ky/2]),Tkstar=-(1+2Cos[Sqrt[3]kx/2]Exp[I 3ky/2]),\[Beta]SOC=(2\[Beta]Ising/(3Sqrt[3]))*(Sin[Sqrt[3]kx]-2Sin[Sqrt[3]kx/2]Cos[3ky/2]),\[CapitalDelta]star=Conjugate[\[CapitalDelta]],\[CapitalDelta]xstar=Conjugate[\[CapitalDelta]x],mztilde=mz-U/4*(nA-nB),\[Mu]tilde=\[Mu]-U/4*(nA+nB-2)(*,\[CapitalDelta]=0,\[CapitalDelta]star=0,U=0*)(*,mz=0.2*)},*)
(*Ek=Sqrt[(-\[Beta]SOC+mztilde)^2+Tk Tkstar];SparseArray[({*)
(* {Ek-\[Mu]tilde, \[CapitalDelta], 0, \[CapitalDelta]x},*)
(* {\[CapitalDelta]star, -(Ek-\[Mu]tilde), \[CapitalDelta]x, 0},*)
(* {0, \[CapitalDelta]xstar, -(Ek+\[Mu]tilde), \[CapitalDelta]},*)
(* {\[CapitalDelta]xstar, 0, \[CapitalDelta]star, Ek+\[Mu]tilde}*)
(*})]]*)


(* ::Input:: *)
(*With[{\[Beta]Ising=0.5,mz=0.5,mu=0.0,\[CapitalDelta]=0.2,\[CapitalDelta]x=0,U=0,nA=1.0857984307519946`,nB=0.9187476777471787`},*)
(*(*Column[{*)Show[{ListPlot[Join[Table[Table[Tooltip[(Eigenvalues[TMDHamMFUp[i[[1]],i[[2]],\[Beta]Ising,mz,mu,\[CapitalDelta],\[CapitalDelta]x,U,nA,nB]]//Sort)[[n]]],{i,kCut}],{n,4}]],PlotStyle->Red,(*Joined\[Rule]True,*)ImageSize->Large,GridLines->Automatic,PlotRange->All(*{-0.5,0.5}*)],ListPlot[Join[Table[Table[Tooltip[(Eigenvalues[TMDHamMFDn[i[[1]],i[[2]],\[Beta]Ising,mz,mu,\[CapitalDelta],\[CapitalDelta]x,U,nA,nB]]//Sort)[[n]]],{i,kCut}],{n,4}]],PlotStyle->Blue,(*Joined\[Rule]True,*)ImageSize->Large,GridLines->Automatic,PlotRange->All(*{-0.5,0.5}*)]}]]*)
