(****************************************************************************)
(*                                                                          *)
(* Application: PSA interface under Tcl/Tk                                  *)
(*                                                                          *)
(* Module:      Integrations.cmd                                            *)
(*                                                                          *)
(* Use:         This file contains all integrations needed for the          *)
(*              phase space analysis                                        *)
(*                                                                          *)
(* Author:      Peter J. Daly and Andreas 'Mounty' Bergmann                 *)
(*                                                                          *)
(* Date:        January, 1995                                               *)
(*                                                                          *)
(* Copyright:   ESRF, BP220, F-38043 GRENOBLE cedex, France.                *)
(*                                                                          *)
(****************************************************************************)


	(****************)
	(* Integrations *)
	(****************)


If[FreeQ[NewSourceX[x,xp,dl],dl],
	IntegralXXP[x_,xp_]=NewSourceX[x,xp,0],
	IntegralXXP[x_,xp_]=
        Integrate[NewSourceX[x,xp,dl],{dl,-Infinity,Infinity},GenerateConditions -> False]
  ];
        IntegralXXP[x_,xp_]=Simplify[IntegralXXP[x,xp]];


If[FreeQ[NewSourceY[y,yp,dl],dl],
	IntegralYYP[y_,yp_]=NewSourceY[y,yp,0],
	IntegralYYP[y_,yp_]=
        Integrate[NewSourceY[y,yp,dl],{dl,-Infinity,Infinity},GenerateConditions -> False]
  ];
        IntegralYYP[y_,yp_]=Simplify[IntegralYYP[y,yp]];


IntegralXYXPYPdl[x_,y_,xp_,yp_,dl_]=
NewSourceX[x,xp,dl]*NewSourceY[y,yp,dl]*SourceLambda[dl]*SourceIntensity;


If[FreeQ[IntegralXYXPYPdl[x,y,xp,yp,dl],dl],
	IntegralXYXPYP[x_,y_,xp_,yp_]=IntegralXYXPYPdl[x,y,xp,yp,0],
	IntegralXYXPYP[x_,y_,xp_,yp_]=
        Integrate[IntegralXYXPYPdl[x,y,xp,yp,dl],{dl,-Infinity,Infinity},GenerateConditions -> False]
  ]; 
        IntegralXYXPYP[x_,y_,xp_,yp_]=Simplify[IntegralXYXPYP[x,y,xp,yp]];

(****************************************************************************)
(*                                                                          *)
(*  Calculations for the FWHM-values (geometric size) of the beam           *)
(*                                                                          *)
(****************************************************************************)

If[FreeQ[IntegralXYXPYP[x,y,xp,yp],yp],
	IntegralXYXP[x_,y_,xp_]=IntegralXYXPYP[x,y,xp,0],
	IntegralXYXP[x_,y_,xp_]=
        Integrate[IntegralXYXPYP[x,y,xp,yp],{yp,-Infinity,Infinity},GenerateConditions -> False] 
  ];
        IntegralXYXP[x_,y_,xp_]=Simplify[IntegralXYXP[x,y,xp]];


If[FreeQ[IntegralXYXP[x,y,xp],xp],
	IntegralXY[x_,y_]=IntegralXYXP[x,y,0],
	IntegralXY[x_,y_]=
        Integrate[IntegralXYXP[x,y,xp],{xp,-Infinity,Infinity},GenerateConditions -> False] 
  ];
        IntegralXY[x_,y_]=Simplify[IntegralXY[x,y]];

	MaxFlux=N[IntegralXY[0,0]];



FunctionX[x_]:=IntegralXY[x,0];
FunctionY[y_]:=IntegralXY[0,y];
ValueAX=FunctionX[0];
ValueAY=FunctionY[0];
ValueExponentX=FunctionX[10^-2];
ValueExponentY=FunctionY[10^-2];
FWHMx=N[Sqrt[4*Log[2]*10^-4/-Log[ValueExponentX/ValueAX]],6]
FWHMy=N[Sqrt[4*Log[2]*10^-4/-Log[ValueExponentY/ValueAY]],6]

(****************************************************************************)
(*                                                                          *)
(*    Calculations for the FWHMp-values (angle size) of the beam            *)
(*                                                                          *)
(****************************************************************************)

If[FreeQ[IntegralXYXPYP[x,y,xp,yp],y],
	IntegralXXPYP[x_,xp_,yp_]=IntegralXYXPYP[x,0,xp,yp],
	IntegralXXPYP[x_,xp_,yp_]=
        Integrate[IntegralXYXPYP[x,y,xp,yp],{y,-Infinity,Infinity},GenerateConditions -> False]
  ];
        IntegralXXPYP[x_,xp_,yp_]=Simplify[IntegralXXPYP[x,xp,yp]];


If[FreeQ[IntegralXXPYP[x,xp,yp],x],
	IntegralXPYP[xp_,yp_]=IntegralXYXPYP[0,xp,yp],
	IntegralXPYP[xp_,yp_]=
        Integrate[IntegralXXPYP[x,xp,yp],{x,-Infinity,Infinity},GenerateConditions -> False]
  ];
        IntegralXPYP[xp_,yp_]=Simplify[IntegralXPYP[xp,yp]];

        MaxFluxP=N[IntegralXPYP[0,0]];



FunctionXP[xp_]:=IntegralXPYP[xp,0];
FunctionYP[yp_]:=IntegralXPYP[0,yp];
ValueAXP=FunctionXP[0];
ValueAYP=FunctionYP[0];
ValueExponentXP=FunctionXP[10^-6];
ValueExponentYP=FunctionYP[10^-6];
FWHMxp=N[Sqrt[4*Log[2]/-Log[ValueExponentXP/ValueAXP]/10^12],8]
FWHMyp=N[Sqrt[4*Log[2]/-Log[ValueExponentYP/ValueAYP]/10^12],8]
FWHMxpp=N[FWHMxp*1000000,8]
FWHMypp=N[FWHMyp*1000000,8]


During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
Out[274]= 0.825823
Out[275]= 0.210807
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
During evaluation of In[214]:= $MaxPrecision::prec: In increasing internal precision while attempting to evaluate 2 π, the limit $MaxPrecision = 14.` was reached. Increasing the value of $MaxPrecision may help resolve the uncertainty. >>
Out[298]= 5.9572*10^21
Out[299]= 5.92287*10^21
Out[300]= 0.0000561856
Out[301]= 0.0000204078
Out[302]= 56.1856
Out[303]= 20.4078