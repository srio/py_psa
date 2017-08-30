IntegralXYXPYPdl[x_,y_,xp_,yp_,dl_]=
NewSourceX[x,xp,dl]*NewSourceY[y,yp,dl]*SourceLambda[dl]*SourceIntensity;


If[FreeQ[IntegralXYXPYPdl[x,y,xp,yp,dl],xp],
	IntegralXYYPdl[x_,y_,yp_,dl_]=IntegralXYXPYPdl[x,y,0,yp,dl],
	IntegralXYYPdl[x_,y_,yp_,dl_]=
        Integrate[IntegralXYXPYPdl[x,y,xp,yp,dl],{xp,-Infinity,Infinity},GenerateConditions -> False]
  ]; 
        IntegralXYYPdl[x_,y_,yp_,dl_]=Simplify[IntegralXYYPdl[x,y,yp,dl]];


If[FreeQ[IntegralXYYPdl[x,y,yp,dl],yp],
	IntegralXYdl[x_,y_,dl_]=IntegralXYYPdl[x,y,0,dl],
	IntegralXYdl[x_,y_,dl_]=
        Integrate[IntegralXYYPdl[x,y,yp,dl],{yp,-Infinity,Infinity},GenerateConditions -> False] 
  ];
        IntegralXYdl[x_,y_,dl_]=Simplify[IntegralXYdl[x,y,dl]];


If[FreeQ[IntegralXYdl[x,y,dl],x],
	IntegralYdl[y_,dl_]=IntegralXYdl[0,y,dl],
	IntegralYdl[y_,dl_]=
        Integrate[IntegralXYdl[x,y,dl],{x,-Infinity,Infinity},GenerateConditions -> False] 
  ];
        IntegralYdl[y_,dl_]=Simplify[IntegralYdl[y,dl]];


If[FreeQ[IntegralYdl[y,dl],y],
	Integraldl[dl_]=IntegralYdl[0,dl],
	Integraldl[dl_]=
        Integrate[IntegralYdl[y,dl],{y,-Infinity,Infinity},GenerateConditions -> False]
  ];
     	Integraldl[dl_]=Simplify[Integraldl[dl]];

ValueAL=Integraldl[0];
ValueExponentL=Integraldl[10^-3];
FWHMl=N[Sqrt[4*Log[2]/-Log[ValueExponentL/ValueAL]/10^6],6];
MaxFluxL=ValueAL;

If[FreeQ[Integraldl[dl],dl],
	FluxPhi=Integraldl[0],
	FluxPhi=
        Integrate[Integraldl[dl],{dl,DlLow,DlUp},GenerateConditions -> False]
  ];
  
FluxPhi = CoefAtten*CoefMonoX*CoefMonoY*FluxPhi
