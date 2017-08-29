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
z1 = 20000;
                                            
(*************************************************************)
(* Definition of the parameters of the used optical elements *)
(*************************************************************)
                                            
                                            
(* Parameters for the 0. Element: Bent Monochromator *)
                                            
Alpha0  = 10^-15*Degree;
Wd0     = 1.000000*10^-6;
RowlandRadius0 = 5.000000*10^3;
WidthX0 = 1.000000*10^1;
WidthY0 = 1.000000*10^1;
RMono0  = 9.000000*10^-1;
Rint0   = 1.300000*10^-5;
thetaB0 = 1.035423*10^1*Degree;
Fc0     = 4.493334*10^2;
DfS0    = 1.000000*10^4;
bMono0  = Sin[thetaB0+Alpha0]/Sin[thetaB0-Alpha0];
                                            
                                            
                                            
Get["/users/bmiller/psa_gatta2/cmd/MatrixDefinition.cmd"];
                                            
Get["/mntdirect/_users/bmiller/psa_gatta2/bin/math_PreCalc.cmd"];
                                            
(*********************************************)
(* Transformation of the source              *)
(*********************************************)
                                            
MatrixTempX=
MatrixFlight[z0].
MatrixMonoFirst[CoefX0].
MatrixFlight[z1].
{x,xp,dl};
                                            
                                            
MatrixTempY=
MatrixFlight[z0].
MatrixMonoBent[bMono0,thetaB0, Fc0,CoefY0].
MatrixFlight[z1].
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
MatrixMonoBent[bMono0,thetaB0, Fc0,CoefY0].
MatrixFlight[z1].
{y,yp,dl};
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceAngleMonoBent[ MatrixTempY[[1]], MatrixTempY[[2]], MatrixTempY[[3]], Alpha0, thetaB0, Wd0, RowlandRadius0, RMono0, Rint0 ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceWaveMonoBent[ MatrixTempY[[3]], thetaB0, DfS0, Wd0, SigmaYSource ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
(*****************************)
(* Starting the Integrations *)
(*****************************)
                                           
DlLow = -2.225721*10^2;
DlUp = 2.225721*10^2;
CoefAtten = 1;
Get["/users/bmiller/psa_gatta2/cmd/Flux.cmd"];
Get["/users/bmiller/psa_gatta2/cmd/Integrations.cmd"];
                                           
                                           
(****************************)
(* Commands for the display *)
(****************************)
                                           
                                           
$DisplayTitle="Intensity at position 30000";
$DisplayWidth=600;
$DisplayHeight=300;
  
Get["/users/bmiller/psa_gatta2/cmd/Graphics.cmd"];
  
  
