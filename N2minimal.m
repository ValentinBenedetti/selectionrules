(* ::Package:: *)
(* N=2 minimal models  (N=2)_k = [su(2)_k + so(2)_1]/u(1)   (notation of N2Gab.pdf, sec. 2.1)
   Fields (l,m,s): l=0..k, m mod 2(k+2), s mod 4, l+m+s even,
   field identification (l,m,s) ~ (k-l, m+k+2, s+2).
   Usage (same philosophy as minimal[p,q] in NonUnitary.nb):
     minimalN2[k]   -> defines rs (labels without duplicates, identity first),
                       hh (weights mod 1), qq (U(1) charges), smat (numerical S matrix ordered as rs), tmat
   then e.g.  verlinde[smat]; dimensions[smat]; loopednu[smat]  *)

cN2[k_]:=3k/(k+2);

hN2[l_,m_,s_,k_]:=(l(l+2)-m^2)/(4(k+2))+s^2/8;   (* mod 1 *)

qN2[l_,m_,s_,k_]:=s/2-m/(k+2);                   (* mod 2 *)

(* reduce m to -(k+1)..k+2 and s to -1..2 *)
redmN2[m_,k_]:=Mod[m+k+1,2(k+2)]-(k+1);
redsN2[s_]:=Mod[s+1,4]-1;

(* canonical key of a field: the smaller of the two reduced representatives of its identification orbit *)
labelN2[{l_,m_,s_},k_]:=First[Sort[{{l,redmN2[m,k],redsN2[s]},{k-l,redmN2[m+k+2,k],redsN2[s+2]}}]];

(* rs = list of primary labels {l,m,s}, each field appearing only once, identity {0,0,0} first *)
fieldsN2[k_]:=Module[{all},all=Select[Flatten[Table[{l,m,s},{l,0,k},{s,-1,2},{m,-(k+1),k+2}],2],EvenQ[Total[#]]&];rs=DeleteDuplicatesBy[all,labelN2[#,k]&];rs=Join[{{0,0,0}},DeleteCases[rs,{0,0,0}]];];

scalingN2[k_,rs_]:=Module[{},hh=Table[Mod[hN2[rs[[i,1]],rs[[i,2]],rs[[i,3]],k],1],{i,1,Length[rs]}];qq=Table[qN2[rs[[i,1]],rs[[i,2]],rs[[i,3]],k],{i,1,Length[rs]}];];

(* eq. (2.6) of N2Gab.pdf; normalisation 1/(k+2) makes smat unitary on the 2(k+1)(k+2) fields
   (the prefactor 1/Sqrt[2(k+2)] printed in the paper is not unitary) *)
modularN2[k_,rs_]:=Module[{},smat=Table[N[1/(k+2) Sin[Pi(rs[[i,1]]+1)(rs[[j,1]]+1)/(k+2)] Exp[I Pi rs[[i,2]] rs[[j,2]]/(k+2)] Exp[-I Pi rs[[i,3]] rs[[j,3]]/2]],{i,1,Length[rs]},{j,1,Length[rs]}];tmat=DiagonalMatrix[Table[N[Exp[2 Pi I (hN2[rs[[i,1]],rs[[i,2]],rs[[i,3]],k]-cN2[k]/24)]],{i,1,Length[rs]}]];];

minimalN2[k_]:=Module[{},fieldsN2[k];scalingN2[k,rs];modularN2[k,rs];]

(* consistency checks: all should return 0 *)
checkN2[k_]:=Module[{},minimalN2[k];{Length[rs]-2(k+1)(k+2),Chop[Norm[smat.ConjugateTranspose[smat]-IdentityMatrix[Length[rs]]]],Chop[Norm[MatrixPower[smat.tmat,3]-smat.smat]]}]
