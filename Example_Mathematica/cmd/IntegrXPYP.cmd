
IntegralXYXPYPdl[x_,y_,xp_,yp_,dl_]=
NewSourceX[x,xp,dl]*NewSourceY[y,yp,dl]*SourceLambda[dl]*SourceIntensity;


If[FreeQ[IntegralXYXPYPdl[x,y,xp,yp,dl],dl],
	IntegralXYXPYP[x_,y_,xp_,yp_]=IntegralXYXPYPdl[x,y,xp,yp,0],
	IntegralXYXPYP[x_,y_,xp_,yp_]=
        Integrate[IntegralXYXPYPdl[x,y,xp,yp,dl],{dl,-Infinity,Infinity},GenerateConditions -> False]
  ]; 
        IntegralXYXPYP[x_,y_,xp_,yp_]=Simplify[IntegralXYXPYP[x,y,xp,yp]];


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
Fxp=N[Sqrt[4*Log[2]/-Log[ValueExponentXP/ValueAXP]/10^12]]
Fyp=N[Sqrt[4*Log[2]/-Log[ValueExponentYP/ValueAYP]/10^12]]
