
IntegralXYXPYPdl[x_,y_,xp_,yp_,dl_]=
NewSourceX[x,xp,dl]*NewSourceY[y,yp,dl]*SourceLambda[dl]*SourceIntensity;


If[FreeQ[IntegralXYXPYPdl[x,y,xp,yp,dl],dl],
	IntegralXYXPYP[x_,y_,xp_,yp_]=IntegralXYXPYPdl[x,y,xp,yp,0],
	IntegralXYXPYP[x_,y_,xp_,yp_]=
        Integrate[IntegralXYXPYPdl[x,y,xp,yp,dl],{dl,-Infinity,Infinity},GenerateConditions -> False]
  ]; 
        IntegralXYXPYP[x_,y_,xp_,yp_]=Simplify[IntegralXYXPYP[x,y,xp,yp]];


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
        FluxIntegralXY[x_,y_]=IntegralXY[x,y];



FunctionX[x_]:=FluxIntegralXY[x,0];
FunctionY[y_]:=FluxIntegralXY[0,y];
ValueAX=FunctionX[0];
ValueAY=FunctionY[0];
ValueExponentX=FunctionX[1];
ValueExponentY=FunctionY[1];
Fx=N[Sqrt[4*Log[2]/-Log[ValueExponentX/ValueAX]]]
Fy=N[Sqrt[4*Log[2]/-Log[ValueExponentY/ValueAY]]]

