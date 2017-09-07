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
                                            
z0 = 32000;
z1 = 3000;
                                            
(*************************************************************)
(* Definition of the parameters of the used optical elements *)
(*************************************************************)
                                            
                                            
(* Parameters for the 0. Element: MultiLayer *)
                                            
WsMu0       = 1.653743*10^-4;
RML0        = 1.000000*10^0;
RowlandMul0 = 5.000000*10^3;
thetaMu0    = 6.720000*10^-1*Degree;
Fmu0        = 1.556381*10^3;
DfS0        = 3.200000*10^4;
                                            
                                            
                                            
Get["/users/bmiller/psa_gatta3/cmd/MatrixDefinition.cmd"];
                                            
                                            
(*********************************************)
(* Transformation of the source              *)
(*********************************************)
                                            
MatrixTempX=
MatrixFlight[z0].
MatrixFlight[z1].
{x,xp,dl};
                                            
                                            
MatrixTempY=
MatrixFlight[z0].
MatrixMultilayer[Fmu0].
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
MatrixMultilayer[Fmu0].
MatrixFlight[z1].
{y,yp,dl};
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceAngleMulti[ MatrixTempY[[1]], MatrixTempY[[2]], MatrixTempY[[3]], thetaMu0, RML0, RowlandMul0, WsMu0 ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
                                            
MatrixTempY=Rationalize[Simplify[MatrixTempY],0];
                                            
NewSourceY[y_,yp_,dl_]=NewSourceY[y,yp,dl]*AcceptanceWaveMulti[ MatrixTempY[[3]], thetaMu0, DfS0, WsMu0, SigmaYSource ];
                                            
NewSourceY[y_,yp_,dl_]=Simplify[NewSourceY[y,yp,dl]];
                                            
(*****************************)
(* Starting the Integrations *)
(*****************************)
                                           
DlLow = -1.326872*10^-4;
DlUp = 1.326872*10^-4;
CoefAtten = 1;
Get["/users/bmiller/psa_gatta3/cmd/Flux.cmd"];
Get["/users/bmiller/psa_gatta3/cmd/Integrations.cmd"];
                                           
                                           
(****************************)
(* Commands for the display *)
(****************************)
                                           
                                           
$DisplayTitle="Intensity at position 35000";
$DisplayWidth=600;
$DisplayHeight=300;
  
Get["/users/bmiller/psa_gatta3/cmd/Graphics.cmd"];
  
  
