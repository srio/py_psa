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
                                            
z0 = 28000;
z1 = 2000;
z2 = 2000;
                                            
(*************************************************************)
(* Definition of the parameters of the used optical elements *)
(*************************************************************)
                                            
                                            
(* Parameters for the 0. Element: Slit *)
                                            
AperturSlitX0 = 1.000000*10^-1;
AperturSlitY0 = 1.000000*10^-1;
                                            
(* Parameters for the 1. Element: Mosaic Monochromator *)
                                            
eta1  = 5.240000*10^-3;
WidthX1 = 1.000000*10^1;
WidthY1 = 1.100000*10^-2;
RMono1  = 1.000000*10^-2;
Rint1   = 1.100000*10^-2;
thetaB1 = 9.136369*10^-2*Degree;
                                            
                                            
                                            
Get["/users/bmiller/psa_gatta2/cmd/MatrixDefinition.cmd"];
                                            
Get["/mntdirect/_users/bmiller/psa_gatta2/bin/math_PreCalc.cmd"];
                                            
(*********************************************)
(* Transformation of the source              *)
(*********************************************)
                                            
MatrixTempX=
MatrixFlight[z0].
MatrixFlight[z1].
MatrixMonoFirst[CoefX1].
MatrixFlight[z2].
{x,xp,dl};
                                            
                                            
MatrixTempY=
MatrixFlight[z0].
MatrixFlight[z1].
MatrixMonoMosaic[thetaB1,CoefY1].
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
                                            
MatrixTempX=
MatrixFlight[z1].
MatrixMonoFirst[CoefX1].
MatrixFlight[z2].
{x,xp,dl};
                                            
                                            
MatrixTempY=
MatrixFlight[z1].
MatrixMonoMosaic[thetaB1,CoefY1].
MatrixFlight[z2].
{y,yp,dl};
                                            
MatrixTempX=Rationalize[Simplify[MatrixTempX],0];
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceX[x_,xp_,dl_]=NewSourceX[x,xp,dl]*AcceptanceSlit[ MatrixTempX[[1]], AperturSlitX0, 0 ];
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceSlit[ MatrixTempY[[1]], AperturSlitY0, 0 ];
                                            
NewSourceX[x_,xp_,dl_]=Simplify[NewSourceX[x,xp,dl]];
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
                                            
                                            
MatrixTempY=
MatrixMonoMosaic[thetaB1,CoefY1].
MatrixFlight[z2].
{y,yp,dl};
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceAngleMonoMosaic[ MatrixTempY[[2]], MatrixTempY[[3]], thetaB1, eta1, RMono1, Rint1 ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceWaveMonoMosaic[ MatrixTempY[[3]], BMonoY*SigmaYPSource, thetaB1, eta1 ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
(*****************************)
(* Starting the Integrations *)
(*****************************)
                                           
DlLow = -3.534557*10^2;
DlUp = 3.534557*10^2;
CoefAtten = 1;
Get["/users/bmiller/psa_gatta2/cmd/Flux.cmd"];
Get["/users/bmiller/psa_gatta2/cmd/Integrations.cmd"];
                                           
                                           
(****************************)
(* Commands for the display *)
(****************************)
                                           
                                           
$DisplayTitle="Intensity at position 32000";
$DisplayWidth=600;
$DisplayHeight=300;
  
Get["/users/bmiller/psa_gatta2/cmd/Graphics.cmd"];
  
  
