(****************************************************************************)
(*                                                                          *)
(* Application: PSA interface under Tcl/Tk                                  *)
(*                                                                          *)
(* Module:      MatrixDefinition.cmd                                        *)
(*                                                                          *)
(* Use:         see below                                                   *)
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
(*  This file contains all definitions of the beam-source, matricies        *)
(*  and acceptance-windows of the optical elements, which are used          *)
(*  in the simulated beamline.                                              *)
(*                                                                          *)
(****************************************************************************)



(*** Definition of the Beam-Source ***)

SourceXXP[x_,xp_]:=Exp[-(x^2/(SigmaXSource)^2+xp^2/(SigmaXPSource)^2)/2];

SourceYYP[y_,yp_]:=
     Exp[-(y^2/(SigmaYSource)^2+(yp-GammaSource*y)^2/SigmaYPSource^2)/2];

SourceLambda[dl_]:=Exp[-dl^2/(2*SigmaSLambda^2)];




(*** Definition of the Matrices ***)

     (* Matrix for the flight-path between the optical elements *)
MatrixFlight[L_]:={{1,-L,0},{0,1,0},{0,0,1}};


     (* Matrix for the first Monochromator *)
MatrixMonoFirst[c_]:=
     {{c,0,0},{0,1,0},{0,0,1}};


     (* Matrix for the plane Monochromator *)
MatrixMonoPlane[b_,thetaBragg_,c_]:=
     {{b*c,0,0},{0,1/b,(1-1/b)*Tan[thetaBragg]},{0,0,1}};


     (* Matrix for the mosaic Monochromator *)
MatrixMonoMosaic[thetaBragg_,c_]:=
     {{c,0,0},{0,-1,2*Tan[thetaBragg]},{0,0,1}};


     (* Matrix for the bent Monochromator *)
MatrixMonoBent[b_,thetaBragg_,Fc_,c_]:=
     {{b*c,0,0},{1/Fc,1/b,(1-1/b)*Tan[thetaBragg]},{0,0,1}};


     (* Matrix for bent and toroidal mirrors *)
MatrixMirror[Fm_,Sig_,IncAng_,Lambda_,Delta_,S_]:=
	Power[E,-(4*Pi*Sin[IncAng]*Sig/Lambda)^2]*
	{{1,0,0},{(1+S*Fm*Delta)/Fm,1,0},{0,0,1}};


     (* Matrix for the plane mirror *)
MatrixMirrorPlane[Sig_,IncAng_,Lambda_,Delta_]:=
	Power[E,-(4*Pi*Sin[IncAng]*Sig/Lambda)^2]*
	{{1,0,0},{Delta,1,0},{0,0,1}};


     (* Matrix for the compound refractive lens with parabolic holes *)
MatrixCompLensPara[F_]:={{1,0,0},{1/F,1,0},{0,0,1}};


     (* Matrix for the compound refractive lens with circular holes *)
MatrixCompLensCirc[F_,Coef_]:={{1,0,0},{Coef/F,1,0},{0,0,1}};


     (* Matrix for the multilayer *)
MatrixMultilayer[Fmu_]:={{1,0,0},{1/Fmu,1,0},{0,0,1}};


(*** Definitions of the acceptance-windows of the optical parts ***)


     (* Definition of the acceptance of the slits *)
AcceptanceSlit[t_,aperture_,calctype_]:=
If[calctype==0,
	Sqrt[6/Pi]/Sqrt[6*Log[2]/Pi]*Exp[-t^2/(2*aperture^2/12)],
	1/Sqrt[6*Log[2]/Pi]*Exp[-t^2/(2*aperture^2/(2*Pi))]
  ];


     (* Definition of the acceptance of the pinhole *)
AcceptancePin[t_,diameter_]:=Sqrt[8/Pi]*Exp[-t^2/(2*diameter^2/16)];


 (* Definition of the angle acceptance of the plane Monochromator *)
AcceptanceAngleMonoPlane[tp_,deltalambda_,theta_,wd_,Rmono_,Rint_]:=
Rmono*Rint*Sqrt[6/Pi]/wd*
Exp[-(tp-deltalambda*Tan[theta])^2/(2*wd^2/12)];
(*Rmono*4*Sqrt[2/(3*Pi)]*
Exp[-(tp-deltalambda*Tan[theta])^2/(2*wd^2/12)];*)


 (* Definition of the wave acceptance of the plane Monochromator *)
AcceptanceWaveMonoPlane[deltalambda_,SigmaTP_,theta_,wd_]:=
Sqrt[6/Pi]*
Exp[-deltalambda^2/(2*((SigmaTP^2+wd^2/12)*Cot[theta]^2))];
 
  
 (* Definition of the angle acceptance of the mosaic Monochromator *)
AcceptanceAngleMonoMosaic[tp_,deltalambda_,theta_,eta_,Rmono_,Rint_]:=
Rint*Sqrt[6/Pi]/eta*
Exp[-(tp-deltalambda*Tan[theta])^2/(2*eta^2)]; 
(* eta or eta*Sqrt[12] ??? *)

 (* Definition of the wave acceptance of the mosaic Monochromator *)
AcceptanceWaveMonoMosaic[deltalambda_,SigmaTP_,theta_,eta_]:=
Sqrt[6/Pi]*
Exp[-deltalambda^2/(2*((SigmaTP^2+eta^2)*Cot[theta]^2))];
 
  
 (* Definition of the angle acceptance of the bent Monochromator *)
AcceptanceAngleMonoBent[t_,tp_,deltalambda_,alpha_,theta_,wd_,r_,Rmono_,Rint_]:=
Rmono*Rint*Sqrt[6/Pi]/wd*
Exp[-(tp-t/(r*Sin[theta+alpha])-deltalambda*Tan[theta])^2/(2*wd^2/12)];


 (* Definition of the wavelength acceptance of the bent Monochromator *)
AcceptanceWaveMonoBent[dl_,theta_,DistanceFromSource_,wd_,SigmaTSource_]:=
Sqrt[6/Pi]*
Exp[-dl^2/(2*(Cot[theta]^2*((SigmaTSource/DistanceFromSource)^2+wd^2/12)))];


    (* Definition of the angle acceptance of the compound refractive lens *)
    (* with a parabolic shape of holes *)
AcceptanceCompLensPara[t_,Radius_,Distance_,Sigma_]:=
Exp[-(t^2+Radius*Distance)/(2*Sigma^2)];


    (* Definition of the angle acceptance of the compound refractive lens *)
    (* with circular holes *)
AcceptanceCompLensCirc[t_,Radius_,Distance_,Sigma_,FWHM_]:=
Exp[-(t^2/(2*Sigma^2))-(t^2*FWHM^2/(8*(Radius*Sigma)^2))-
(t^2*FWHM^4/(16*Sigma^2*Radius^4))-(Radius*Distance/(2*Sigma^2))];


    (* Definition of the angle acceptance of the multilayer *)
AcceptanceAngleMulti[t_,tp_,dl_,theta_,Rml_,Rowland_,ws_]:=
Rml*8/3*Sqrt[Log[2]/Pi]*
Exp[-(-tp-t/(Rowland*Sin[theta])-dl*Tan[theta])^2*8*Log[2]/(2*ws^2)];


    (* Definition of the wavelength acceptance of the multilayer *)
AcceptanceWaveMulti[dl_,theta_,DistanceFromSource_,ws_,SigmaTSource_]:=
Sqrt[6/Pi]*
Exp[-dl^2/(2*(SigmaTSource^2/DistanceFromSource^2+ws^2/(8*Log[2]))*Cot[theta]^2)];




