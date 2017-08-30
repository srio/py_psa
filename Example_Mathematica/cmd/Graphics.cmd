(****************************************************************************)
(*                                                                          *)
(* Application: PSA interface under Tcl/Tk                                  *)
(*                                                                          *)
(* Module:      Graphics.cmd                                                *)
(*                                                                          *)
(* Use:         Graphics module for the Math.cmd file, which is send to     *)
(*              Mathematica                                                 *)
(*                                                                          *)
(* Author:      Peter J. Daly and Andreas 'Mounty' Bergmann                 *)
(*                                                                          *)
(* Date:        January, 1995                                               *)
(*                                                                          *)
(* Copyright:   ESRF, BP220, F-38043 GRENOBLE cedex, France.                *)
(*                                                                          *)
(****************************************************************************)


(****************************************************************************)
(*                                                                          *)
(*  Calculations for the ranges in x,xp and y,yp directions                 *)
(*                                                                          *)
(****************************************************************************)


Coef=2.576; (* 99% of the Gaussian *)


RangeX=N[Coef*FWHMx/(2*Sqrt[2*Log[2]])];

RangeXP=N[Coef*FWHMxp/(2*Sqrt[2*Log[2]])];

RangeY=N[Coef*FWHMy/(2*Sqrt[2*Log[2]])];

RangeYP=N[Coef*FWHMyp/(2*Sqrt[2*Log[2]])];

RangeL=N[Coef*FWHMl/(2*Sqrt[2*Log[2]])];



(****************************************************************************)
(*                                                                          *)
(*  Definings for the graphics of the result-function                       *)
(*                                                                          *)
(****************************************************************************)

(* <<Graphics`Graphics3D` *)

Graph3DXY:=Plot3D[IntegralXY[x,y]/MaxFlux,
  {x,-RangeX,RangeX},
  {y,-RangeY,RangeY},
  Axes->True,
  AxesLabel->{"X","Y","  Distribution"},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

Graph3DXYtxt:=Table[{x,y,IntegralXY[x,y]/MaxFlux},
  {x,-RangeX,RangeX,RangeX/20},
  {y,-RangeY,RangeY,RangeY/20}
];

ContourGraphXY:=ContourPlot[IntegralXY[x,y],
  {x,-RangeX,RangeX},
  {y,-RangeY,RangeY},
  Contours->20,
  ColorFunction->Hue,
  ContourLines->False,
  Axes->True,
  AxesLabel->{"X","Y"},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

DensityGraphXY:=DensityPlot[IntegralXY[x,y],
  {x,-RangeX,RangeX},
  {y,-RangeY,RangeY},
  PlotRange->{0,MaxFlux},
  Mesh->False,
  ColorFunction->Hue,
  PlotPoints->20,
  Axes->True,
  AxesLabel->{"X","Y"},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

GraphXYtxt:=Table[{x,y,IntegralXY[x,y]},
  {x,-RangeX,RangeX,RangeX/20},
  {y,-RangeY,RangeY,RangeY/20}
];

GraphX:=Plot[IntegralXY[x,0]/MaxFlux,
  {x,-RangeX,RangeX},
  AxesLabel->{"X","Distribution"},
  PlotRange->{0,1},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

GraphXtxt:=Table[{x,IntegralXY[x,0]/MaxFlux},
  {x,-RangeX,RangeX,RangeX/20}
];

GraphY:=Plot[IntegralXY[0,y]/MaxFlux,
  {y,-RangeY,RangeY},
  AxesLabel->{"Y","Distribution"},
  PlotRange->{0,1},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

GraphYtxt:=Table[{y,IntegralXY[0,y]/MaxFlux},
  {y,-RangeY,RangeY,RangeY/20}
];

GraphL:=Plot[Integraldl[dl]/MaxFluxL,
  {dl,-RangeL,RangeL},
  AxesLabel->{"Lambda","Distribution"},
  PlotRange->{0,1},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

GraphLtxt:=Table[{dl,Integraldl[dl]/MaxFluxL},
  {dl,-RangeL,RangeL,RangeL/20}
];

GraphXP:=Plot[IntegralXPYP[xp,0],
  {xp,-RangeXP,RangeXP},
  AxesLabel->{"XP","Z"},
  PlotRange->{0,MaxFluxP},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

GraphXPtxt:=Table[{xp,IntegralXPYP[xp,0]},
  {xp,-RangeXP,RangeXP,RangeXP/20}
];

GraphYP:=Plot[IntegralXPYP[0,yp],
  {yp,-RangeYP,RangeYP},
  AxesLabel->{"YP","Z"},
  PlotRange->{0,MaxFluxP},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

GraphYPtxt:=Table[{yp,IntegralXPYP[0,yp]},
  {yp,-RangeYP,RangeYP,RangeYP/20}
];

ContourGraphXXP:=ContourPlot[IntegralXXP[x,xp],
  {x,-RangeX,RangeX},
  {xp,-RangeXP,RangeXP},
  Contours->20,
  ColorFunction->Hue,
  ContourLines->False,
  Axes->True,
  AxesLabel->{"X","XP"},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

ContourGraphXXPtxt:=Table[{x,xp,IntegralXXP[x,xp]},
  {x,-RangeX,RangeX,RangeX/20},
  {xp,-RangeXP,RangeXP,RangeXP/20}
];

ContourGraphYYP:=ContourPlot[IntegralYYP[y,yp],
  {y,-RangeY,RangeY},
  {yp,-RangeYP,RangeYP},
  Contours->20,
  ColorFunction->Hue,
  ContourLines->False,
  Axes->True,
  AxesLabel->{"Y","YP"},
  PlotLabel->$DisplayTitle,
  DisplayFunction->Identity
];

ContourGraphYYPtxt:=Table[{y,yp,IntegralYYP[y,yp]},
  {y,-RangeY,RangeY,RangeY/20},
  {yp,-RangeYP,RangeYP,RangeYP/20}
];
