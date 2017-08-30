$MaxPrecision = 14;
        (*********************************************)
        (* This is the command-file for Mathematica. *)
        (* It was created of the program PSA,        *)
        (* programmed by                             *)
        (* Peter J. DALY and Andreas Mounty BERGMANN *)
        (* under Tcl/Tk.                             *)
        (* Modified August 1998 by G. Gatta          *)
        (*                                           *)
        (* This command-file will be sent to         *)
        (* Mathematica to start the calculation      *)
        (* for the simulated beamline.               *)
        (*                                           *)
        (*                                           *)
        (* This file was created by:                 *)
        (* bmiller                                   *)
        (*                                           *)
        (*                                           *)
        (*                                           *)
        (*********************************************)
                                            
                                            
(*********************************************)
(* Definiton of the parameters of the source *)
(*********************************************)
                                            
SigmaXSource    = 1.000000*10^-1;
SigmaYSource    = 1.000000*10^-2;
SigmaXPSource   = 3.000000*10^-5;
SigmaYPSource   = 1.000000*10^-4;
GammaSource     = 0;
SigmaSLambda    = 1.000000*10^-3;
SourceIntensity = 1.000000*10^25;
Energy          = 1.100000*10^4;
                                            
FilterCoef      = 1;
CoefMonoX       = 1;
CoefMonoY       = 1;
BMonoX          = 1;
BMonoY          = 1;
                                            
(*********************************************)
(* Definition of the Flight-Paths            *)
(*********************************************)
                                            
z0 = 10000;
z1 = 5000;
z2 = 3000;
                                            
(*************************************************************)
(* Definition of the parameters of the used optical elements *)
(*************************************************************)
                                            
                                            
(* Parameters for the 0. Element: Toroidal Mirror *)
                                            
SigmaX0 = 4.000000*10^-7;
SigmaY0 = 4.000000*10^-7;
DeltaX0 = 1.000000*10^-12;
DeltaY0 = 1.000000*10^-11;
IncAng0 = 1.000000*10^-1*Degree;
Lambda0 = 1.127140*10^-7;
FmX0 = 1.000000*10^4;
FmY0 = 8.000000*10^3;
SX = 1;
SY = 1;
                                            
(* Parameters for the 1. Element: Plane Monochromator *)
                                            
Alpha1  = 10^-15*Degree;
Wd1     = 1.000000*10^-6;
WidthX1 = 2.000000*10^1;
WidthY1 = 2.000000*10^1;
RMono1  = 9.600000*10^-1;
Rint1   = 1.000000*10^-5;
thetaB1 = 1.035423*10^1*Degree;
bMono1  = Sin[thetaB1+Alpha1]/Sin[thetaB1-Alpha1];
                                            
                                            
                                            
Get["/users/bmiller/psa_gatta2/cmd/MatrixDefinition.cmd"];
                                            
Get["/mntdirect/_users/bmiller/psa_gatta2/bin/math_PreCalc.cmd"];
                                            
(*********************************************)
(* Transformation of the source              *)
(*********************************************)
                                            
MatrixTempX=
MatrixFlight[z0].
MatrixMirror[FmX0,SigmaX0, IncAng0,Lambda0,DeltaX0,SX].
MatrixFlight[z1].
MatrixMonoFirst[CoefX1].
MatrixFlight[z2].
{x,xp,dl};
                                            
                                            
MatrixTempY=
MatrixFlight[z0].
MatrixMirror[FmY0,SigmaY0, IncAng0,Lambda0,DeltaY0,SY].
MatrixFlight[z1].
MatrixMonoPlane[bMono1,thetaB1, CoefY1].
MatrixFlight[z2].
{y,yp,dl};
                                            
MatrixTempX=Rationalize[Simplify[MatrixTempX],0];
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceX[x_,xp_,dl_]=SourceXXP[MatrixTempX[[1]],MatrixTempX[[2]]];
NewSourceY[y_,yp_,dl_]=SourceYYP[MatrixTempY[[1]],MatrixTempY[[2]]];
                                            
NewSourceX[x_,xp_,dl_]=Simplify[NewSourceX[x,xp,dl]];
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
(*******************************************************)
(* Now the source has been transformed into NewSource  *)
(*******************************************************)
                                            
                                            
                                            
MatrixTempY=
MatrixMonoPlane[bMono1,thetaB1, CoefY1].
MatrixFlight[z2].
{y,yp,dl};
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceAngleMonoPlane[ MatrixTempY[[2]], MatrixTempY[[3]], thetaB1, Wd1, RMono1, Rint1 ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceWaveMonoPlane[ MatrixTempY[[3]], BMonoY*SigmaYPSource, thetaB1, Wd1 ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
(*****************************)
(* Starting the Integrations *)
(*****************************)
                                           
DlLow = -9.218515*10^-4;
DlUp = 9.218515*10^-4;
CoefAtten = 1;
Get["/users/bmiller/psa_gatta2/cmd/Flux.cmd"];
Get["/users/bmiller/psa_gatta2/cmd/Integrations.cmd"];
                                           
                                           
(****************************)
(* Commands for the display *)
(****************************)
                                           
                                           
$DisplayTitle="Intensity at position 18000";
$DisplayWidth=600;
$DisplayHeight=300;
  
Get["/users/bmiller/psa_gatta2/cmd/Graphics.cmd"];
  
  
