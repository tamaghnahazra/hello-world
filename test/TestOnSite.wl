(* ::Package:: *)

(* ::Text:: *)
(*Test based on Tsuchiya et al.*)


TablePrint[t_Association] := Module[{keys = Sort[Keys[t]]},
  Do[Print[k, " : ", t[k]], {k, keys}]
]

<<TMDHartreeFock`


Us=Range[-8.0, 8.0, 0.1];
Vs={0.0};


spinString[1]="UP";
spinString[2]="DN";



results=Dataset@Flatten@
Module[{result,diff,\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc,\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp},
	\[CapitalGamma]c = 0.5*(Us[[1]]+6Vs[[1]])*{1,1,1,1}+RandomReal[{-0.1,0.1},4];
	\[CapitalDelta]osc = RandomReal[{-0.5, 0.5}, 2];
	\[CapitalDelta]nnc = RandomReal[{-0.5, 0.5}, {3,2,2}];
	Table[
		InitializeTMD[24, 24,
			OnSiteInteraction->U,
			NearestNeighborInteraction->V,
			ChemicalPotential->0.0,
			IsingSpinOrbitCoupling->0.1,
			Temperature->0];
        Print["Starting calculation with U = ", U, ", V = ", V];
        WriteString["stdout", "Stabilizing..."];
		Do[
			{\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp}={\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc};
			{\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc}=CompiledCollectMF@@{\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp};
			,{i,1,20}
		];
        Print["Done."];

        WriteString["stdout", "Finding solution..."];
		Do[
			{\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp}={\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc};
			{\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc}=CompiledCollectMF@@{\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp};
			diff=Max[Abs[Flatten[{\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc}]-Flatten[{\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp}]]];
			If[diff<10^-8, WriteString["stdout", "Within ", i, " steps,"]; Break[]]
			,{i,1,3000}
		];
        Print["Finished."];
		diff=Flatten[{\[CapitalGamma]c, \[CapitalDelta]osc, \[CapitalDelta]nnc}]-Flatten[{\[CapitalGamma]p, \[CapitalDelta]osp, \[CapitalDelta]nnp}];
		result=Association[{
			"PARAMETER:ON-SITE INTERACTION"->U,
			"PARAMETER:NEAREST-NEIGHBOR INTERACTION"->V,
			"MEAN FIELD:DENSITY A UP"->\[CapitalGamma]c[[1]],
			"MEAN FIELD:DENSITY B UP"->\[CapitalGamma]c[[2]],
			"MEAN FIELD:DENSITY A DN"->\[CapitalGamma]c[[3]],
			"MEAN FIELD:DENSITY B DN"->\[CapitalGamma]c[[4]],
			"MEAN FIELD:ON-SITE PAIRING A"->\[CapitalDelta]osc[[1]],
			"MEAN FIELD:ON-SITE PAIRING B"->\[CapitalDelta]osc[[2]]
		}
		~Join~Flatten[Table[("MEAN FIELD:NEAREST-NEIGHBOR PAIRING BOND-"<>ToString[i]<>" A-"<>spinString[a]<>" B-"<>spinString[b])->\[CapitalDelta]nnc[[i,a,b]],{i,1,3},{a,1,2},{b,1,2}]]
		~Join~{"SELF-CONSISTENT LOOP ERROR"->Max[Abs[Flatten[diff]]]}];
        TablePrint[result];
		result
		,{U,Us},{V,Vs}
	]
];
Export["result_onsite.wl",results]
