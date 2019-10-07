#pragma rtGlobals=3		// Use modern global access method and strict wave access.
//*******************************************************************************************************************
//FUNCTIONS AVAILABLE IN THIS PROCEDURE FILE

//			auto_solveJ_eigenvalue(x1,x2)   (x1= J value minima ; x2 = J value maxima)

//			JfoldersPrep(x2)			(x2 = J value upto which Hamiltonian will be made)

//			JfoldersPrep_N (x2,A,B)    (x2 = J value upto which Hamiltonian will be made, A,B= {1=H, 2=D})

//			MatrixPrep_auto(J,ma,mb, r_w, pot_wave)    (this takes repective inputs from JfoldersPrep(x2)	procedure
//			 										and makes the Hamilatonian )

//*********************************************************************************************************************
function ini_hm()
DFREF cdf1 = GetDataFolderDFR()
 variable checkf=datafolderexists ("root:Params:RoVibHamiltonian")
 if (checkf == 0)
 	setdatafolder  root:
 	newdatafolder /O /s Params
 	newdatafolder /O Customspline
 	newdatafolder /O /s RoVibHamiltonian
 	
 	variable /G j_min=0 ; variable /G j_max=0 , massA=0,massB=0,n_electronA=0,n_electronB=0;
 //	variable hartree_asym,AsympCorr
 	variable /G pmin=0,pmax=0,atom1=1,atom2=1;
 	variable /G llim=0, ulim=0, int_pnts=25,int_thresh=14
 	string /G nwave,distance_wave
 	string /G iwv1,iwv2,iwv3,iwv4
 	variable /G int_min,int_max,step_h=0.005
 	variable /G J_ini=0, Bint_min=0.20,Bint_max=3.000
 	variable /G J_final=2
 	string /G root:Params:RoVibHamiltonian:b_norm_distance
 	string potential_rwave, potential_energywave
 	variable	 /G  hartree_asym=-1
 	variable /G   AsympCorr=117.1522
 	endif
 	setdatafolder root: 
 	
 	string cmd="RoVib_Diatomic()"
 	execute cmd
 	
 	//-------------------
 	 cmd=" Prep_HMatrixOther()"
 	execute cmd
 	
 	setdatafolder cdf1
end

Window RoVib_Diatomic() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1717,54,1909,988) as "RoVib_Diatomic"
	ShowInfo/W=$WinName(0,64)
	SetDrawLayer UserBack
	DrawLine 4,70,181,70
	DrawLine 2,148,175,148
	DrawLine 2,294,179,294
	DrawLine 5,549,177,549
	SetDrawEnv dash= 1
	DrawLine 104,455,195,455
	SetDrawEnv dash= 1
	DrawLine 2,487,180,487
	SetDrawEnv dash= 2
	DrawLine 6,578,117,578
	DrawLine 5,690,177,690
	SetDrawEnv dash= 2
	DrawLine 107,596,107,578
	SetDrawEnv dash= 2
	DrawLine 119,624,186,624
	DrawLine -1,861,110,861
	DrawLine 109,837,109,861
	DrawLine 184,836,111,836
	SetVariable llim5,pos={100.00,628.00},size={59.00,18.00},title="--"
	SetVariable llim5,limits={1,25,0.005},value= root:Params:RoVibHamiltonian:ulim,live= 1
	Button JFoldersPrep,pos={27.00,49.00},size={125.00,20.00},proc=JfoldersPrep_N2,title="JFoldersPrep"
	Button JFoldersPrep,fColor=(65280,43520,0)
	Button AutoSolve_EigenVal,pos={5.00,96.00},size={118.00,20.00},proc=auto_solveJ_eigenval,title="AutoSolve_EigenVal"
	Button AutoSolve_EigenVal,fColor=(0,43520,65280)
	SetVariable J_max,pos={124.00,75.00},size={60.00,18.00},title="Max"
	SetVariable J_max,limits={1,20,1},value= root:Params:RoVibHamiltonian:j_max,live= 1
	SetVariable j_min,pos={59.00,75.00},size={65.00,18.00},title="J_min"
	SetVariable j_min,limits={0,20,1},value= root:Params:RoVibHamiltonian:j_min
	TitleBox Solution,pos={4.00,74.00},size={54.00,23.00},title="Solution",fSize=12
	TitleBox Solution,fStyle=1
	TitleBox Prep,pos={3.00,3.00},size={74.00,23.00},title="Preparing H",fStyle=1
	SetVariable Atom1,pos={6.00,28.00},size={86.00,18.00},title="Atom_A"
	SetVariable Atom1,limits={1,2,1},value= root:Params:RoVibHamiltonian:atom1,styledText= 1,live= 1
	SetVariable pmax,pos={84.00,4.00},size={75.00,18.00},title="P_max"
	SetVariable pmax,limits={0,20,1},value= root:Params:RoVibHamiltonian:pmax,styledText= 1,live= 1
	SetVariable Atom2,pos={93.00,28.00},size={84.00,18.00},title="Atom_B"
	SetVariable Atom2,limits={1,2,1},value= root:Params:RoVibHamiltonian:atom2,live= 1
	TitleBox Normalization,pos={2.00,149.00},size={86.00,23.00},title="Normalization"
	TitleBox Normalization,fStyle=1
	PopupMenu popup_wavelist,pos={92.00,152.00},size={67.00,19.00},proc=NormWaveSelection,title="wfn"
	PopupMenu popup_wavelist,mode=3,popvalue="v0J0",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu popup_wavelist1,pos={4.00,176.00},size={105.00,19.00},proc=Norm_distance_Selection,title="distance"
	PopupMenu popup_wavelist1,mode=2,popvalue="r_wave",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	SetVariable llim,pos={10.00,199.00},size={110.00,18.00},title="Low_limit"
	SetVariable llim,limits={0,25,0.005},value= root:Params:RoVibHamiltonian:llim,live= 1
	SetVariable llim1,pos={2.00,217.00},size={116.00,18.00},title="Upper_limit"
	SetVariable llim1,limits={1,25,0.005},value= root:Params:RoVibHamiltonian:ulim,live= 1
	SetVariable Npnts,pos={9.00,235.00},size={110.00,18.00},title="Int_points"
	SetVariable Npnts,limits={25,1000,25},value= root:Params:RoVibHamiltonian:int_pnts,styledText= 1,live= 1
	Button button_Norm,pos={15.00,272.00},size={90.00,20.00},proc=norm_wfn,title="Wfn_Normalize"
	Button button_Norm,fColor=(16384,65280,41216)
	SetVariable setvar_thsr,pos={22.00,253.00},size={96.00,18.00},title="Thresh: 1e-"
	SetVariable setvar_thsr,limits={12,16,1},value= root:Params:RoVibHamiltonian:int_thresh,live= 1
	TitleBox title_expVal,pos={3.00,297.00},size={105.00,23.00},title="ExpectationValue"
	TitleBox title_expVal,fStyle=1
	PopupMenu popup_wavelist2,pos={29.00,321.00},size={143.00,19.00},proc=ExpVal_wf1_Selection,title="Norm_wfn1"
	PopupMenu popup_wavelist2,mode=150,popvalue="v2J0_norm",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu popup_wavelist3,pos={30.00,345.00},size={143.00,19.00},proc=ExpVal_wf2_Selection,title="Norm_wfn2"
	PopupMenu popup_wavelist3,mode=150,popvalue="v2J0_norm",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu popup_wavelist5,pos={35.00,386.00},size={168.00,19.00},proc=ExpVal_parameter_Selection,title="Parameter"
	PopupMenu popup_wavelist5,mode=74,popvalue="isotropy_234_984",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu popup_wavelist6,pos={6.00,406.00},size={164.00,19.00},proc=ExpVal_parameter_distance_sel,title="distance_expVal"
	PopupMenu popup_wavelist6,mode=68,popvalue="g_distance",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	Button button_expVal,pos={13.00,526.00},size={115.00,20.00},proc=expectation_val,title="Calc_expectn_value"
	Button button_expVal,fColor=(44032,29440,58880)
	SetVariable setvar_pdmin1,pos={61.00,427.00},size={82.00,18.00},title="rangm"
	SetVariable setvar_pdmin1,limits={0,20,1e-008},value= root:Params:RoVibHamiltonian:pd_vmax,noedit= 1
	SetVariable setvar_pdmin,pos={6.00,427.00},size={89.00,18.00},title="Range"
	SetVariable setvar_pdmin,limits={0,20,1e-008},value= root:Params:RoVibHamiltonian:pd_vmin,noedit= 1
	SetVariable llim2,pos={5.00,468.00},size={65.00,18.00},title="a="
	SetVariable llim2,limits={0.1,5,0.005},value= root:Params:RoVibHamiltonian:int_min,live= 1
	SetVariable llim3,pos={79.00,468.00},size={65.00,18.00},title="b="
	SetVariable llim3,limits={0.5,25,0.005},value= root:Params:RoVibHamiltonian:int_max,live= 1
	SetVariable step,pos={13.00,489.00},size={100.00,18.00},title="step_h"
	SetVariable step,limits={0.001,2,0.001},value= root:Params:RoVibHamiltonian:step_h
	TitleBox title0,pos={4.00,446.00},size={99.00,23.00},title="Integration range"
	PopupMenu popup_wavelist4,pos={48.00,366.00},size={105.00,19.00},proc=ExpVal_wf_distance,title="distance"
	PopupMenu popup_wavelist4,mode=13,popvalue="r_wave",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	SetVariable Npnts1,pos={4.00,507.00},size={110.00,18.00},title="Int_points"
	SetVariable Npnts1,limits={25,1000,25},value= root:Params:RoVibHamiltonian:int_pnts,styledText= 1,live= 1
	TitleBox Bnormalize,pos={6.00,693.00},size={181.00,21.00},title="Batch_ExpVal_for Specific gamma"
	TitleBox Bnormalize,fSize=11,fStyle=1
	PopupMenu popup_wavelist7,pos={8.00,607.00},size={105.00,19.00},proc=Batch_norm_disatance,title="distance"
	PopupMenu popup_wavelist7,mode=2,popvalue="r_wave",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	SetVariable llim4,pos={3.00,629.00},size={96.00,18.00},title="Low_lim"
	SetVariable llim4,limits={0,25,0.005},value= root:Params:RoVibHamiltonian:llim,live= 1
	SetVariable Npnts2,pos={3.00,649.00},size={98.00,18.00},title="Int_points"
	SetVariable Npnts2,limits={25,1000,25},value= root:Params:RoVibHamiltonian:int_pnts,styledText= 1,live= 1
	SetVariable setvar_thsr1,pos={102.00,648.00},size={96.00,18.00},title="Thresh: 1e-"
	SetVariable setvar_thsr1,limits={12,16,1},value= root:Params:RoVibHamiltonian:int_thresh,live= 1
	Button button_BatchNormalize,pos={13.00,669.00},size={115.00,20.00},proc=norm_wfn_batch,title="Batch_Normalize"
	Button button_BatchNormalize,fColor=(65535,32768,32768)
	Button Get_wfns,pos={79.00,552.00},size={70.00,20.00},proc=get_wfns_v0,title="Wfns_v0J(i)"
	Button Get_wfns,fColor=(49163,65535,32768)
	PopupMenu popup_wavelist8,pos={8.00,714.00},size={168.00,19.00},proc=BExpVal_parameter_Selection,title="Parameter"
	PopupMenu popup_wavelist8,mode=76,popvalue="isotropy_265_986",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu popup_wavelist9,pos={11.00,735.00},size={164.00,19.00},proc=BExpVal_parameter_Selection_d,title="distance_expVal"
	PopupMenu popup_wavelist9,mode=68,popvalue="g_distance",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	SetVariable J_ini,pos={6.00,757.00},size={80.00,18.00},title="J_ini"
	SetVariable J_ini,limits={0,45,1},value= root:Params:RoVibHamiltonian:J_ini
	SetVariable J_final,pos={93.00,757.00},size={90.00,18.00},title="J_final"
	SetVariable J_final,limits={1,45,1},value= root:Params:RoVibHamiltonian:J_final
	Button Batch_ExpVal_Spec_Gamma,pos={4.00,818.00},size={99.00,20.00},proc=batch_integral_me_gamma,title="ME_Spec.Gamma"
	Button Batch_ExpVal_Spec_Gamma,fColor=(32768,40777,65535)
	PopupMenu Wfn_distance_wave,pos={4.00,776.00},size={164.00,19.00},proc=BExpVal_wfn_distance_Selection,title="Wfn_distance_wave"
	PopupMenu Wfn_distance_wave,mode=13,popvalue="r_wave",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	TitleBox batch_norm,pos={7.00,581.00},size={93.00,23.00},title="Batch normalize"
	TitleBox batch_norm,fStyle=0
	TitleBox ext_wfns,pos={3.00,552.00},size={73.00,23.00},title="Extract_wfns"
	TitleBox ext_wfns,fStyle=0
	SetVariable llim6,pos={107.00,816.00},size={66.00,18.00},title="b="
	SetVariable llim6,limits={0.5,25,0.005},value= root:Params:RoVibHamiltonian:Bint_max,live= 1
	SetVariable llim7,pos={106.00,796.00},size={66.00,18.00},title="a="
	SetVariable llim7,limits={0.1,5,0.005},value= root:Params:RoVibHamiltonian:Bint_min,live= 1
	TitleBox title1,pos={3.00,795.00},size={99.00,23.00},title="Integration range"
	Button CleanCSpline,pos={123.00,177.00},size={65.00,40.00},proc=cleanCustomSpline,title="Clean_CS"
	Button CleanCSpline,fColor=(65535,49151,49151)
	Button Clean_HM,pos={122.00,229.00},size={65.00,45.00},proc=clean_H_Matrix,title="Clean_HM"
	Button Clean_HM,fColor=(49151,49152,65535)
	Button Print_eval,pos={3.00,117.00},size={62.00,20.00},proc=print_eval_J0,title="Print_eval"
	Button Print_eval,fColor=(65535,32764,16385)
	SetVariable AsymptCorrn,pos={60.00,132.00},size={130.00,16.00},title="AsymptCorrn"
	SetVariable AsymptCorrn,fSize=11,format="%7.7f"
	SetVariable AsymptCorrn,limits={-inf,inf,0},value= root:Params:RoVibHamiltonian:AsympCorr
	SetVariable H_v,pos={68.00,115.00},size={108.00,16.00},title="H_v",fSize=11
	SetVariable H_v,format="%7.9f"
	SetVariable H_v,limits={0,inf,0.05},value= root:Params:RoVibHamiltonian:hartree_asym
	Button Opt_ExpVal,pos={4.00,840.00},size={99.00,20.00},proc=batch_integral_me_meanpol,title="ME_Spec.Alpha"
	Button Opt_ExpVal,fColor=(65535,43688,32768)
	Button Get_wfns1,pos={150.00,552.00},size={37.00,20.00},proc=get_wfns_v1,title="v1J(i)"
	Button Get_wfns1,fColor=(32768,40777,65535)
	Button Batch_ExpVal_Spec_Gamma1,pos={4.00,865.00},size={65.00,20.00},proc=Set_integral_me_gamma_v0,title="Gamma:v0"
	Button Batch_ExpVal_Spec_Gamma1,fColor=(65535,49151,62258)
	Button Opt_ExpVal1,pos={77.00,872.00},size={100.00,20.00},proc=Set_integral_me_meanPolv0v1,title="Set:ME:alphav0v1"
	Button Opt_ExpVal1,fColor=(52428,52425,1)
	Button Get_wfns_v2,pos={113.00,575.00},size={37.00,20.00},proc=get_wfns_v2,title="v2J(i)"
	Button Get_wfns_v2,fColor=(32768,65535,65535)
	Button Opt_ExpVal_alphav2,pos={74.00,895.00},size={104.00,20.00},proc=Set_integral_me_meanPolv1v2,title="Set:ME:alphav1v2"
	Button Opt_ExpVal_alphav2,fColor=(65535,43690,0)
	TitleBox Bnormalize1,pos={113.00,839.00},size={70.00,21.00},title="Full set MEs"
	TitleBox Bnormalize1,fSize=11,fStyle=1
	Button Batch_ExpVal_Spec_Gamma2,pos={5.00,886.00},size={65.00,20.00},proc=Set_integral_me_gamma_v1,title="Gamma:v1"
	Button Batch_ExpVal_Spec_Gamma2,fColor=(65535,0,0)
	Button Batch_ExpVal_Spec_Gamma3,pos={6.00,907.00},size={65.00,20.00},proc=Set_integral_me_gamma_v2,title="Gamma:v2"
	Button Batch_ExpVal_Spec_Gamma3,fColor=(0,43690,65535)
	Button Get_wfns_v3,pos={152.00,576.00},size={37.00,20.00},proc=get_wfns_v3,title="v3J(i)"
	Button Get_wfns_v3,fColor=(65535,43690,0)
	Button Get_wfns_v4,pos={114.00,597.00},size={37.00,20.00},proc=get_wfns_v4,title="v4J(i)"
	Button Get_wfns_v4,fColor=(19675,39321,1)
	Button Get_wfns_v5,pos={152.00,597.00},size={37.00,20.00},proc=get_wfns_v5,title="v5J(i)"
	Button Get_wfns_v5,fColor=(1,12815,52428)
EndMacro






//*********************************************************************************************************************
Function NormWaveSelection(norm_sel,popNumA,popStr1) : PopupMenuControl
String norm_sel
Variable popNumA
String popStr1
string /G root:Params:RoVibHamiltonian:nwave =popStr1
End

Function Norm_distance_Selection(norm_sel2,popNumB,popStr2) : PopupMenuControl
String norm_sel2
Variable popNumB
String popStr2
string /G root:Params:RoVibHamiltonian:distance = popStr2
End
//*********************************************************************************************************************
Function ExpVal_wf1_Selection(ev_sel1,popNum_ev1,popStr_ev1) : PopupMenuControl	//wavefunction 1
String ev_sel1
Variable popNum_ev1
String popStr_ev1
string /G root:Params:RoVibHamiltonian:iwv1 = popStr_ev1
End
//*********************************************************************************************************************

Function ExpVal_wf2_Selection(ev_sel2,popNum_ev2,popStr_ev2) : PopupMenuControl	//wavefunction 2
String ev_sel2
Variable popNum_ev2
String popStr_ev2
string /G root:Params:RoVibHamiltonian:iwv2 = popStr_ev2
End
//*********************************************************************************************************************

Function ExpVal_wf_distance(ev_sel2d,popNum_ev2d,popStr_ev2d) : PopupMenuControl  //wavefunction distance
String ev_sel2d
Variable popNum_ev2d
String popStr_ev2d
string /G root:Params:RoVibHamiltonian:int_distance = popStr_ev2d
End

//*********************************************************************************************************************

Function ExpVal_parameter_Selection(ev_sel3,popNum_ev3,popStr_ev3) : PopupMenuControl	//parameter (eg; gamma)
String ev_sel3
Variable popNum_ev3
String popStr_ev3
string /G root:Params:RoVibHamiltonian:iwv3 = popStr_ev3
End

//*********************************************************************************************************************

Function ExpVal_parameter_distance_sel(ev_sel4,popNum_ev4,popStr_ev4) : PopupMenuControl	//parameter distance
String ev_sel4
Variable popNum_ev4
String popStr_ev4
string /G root:Params:RoVibHamiltonian:iwv4 = popStr_ev4
string cdf_ev=GetDataFolder(1)
string s1,s2
sprintf s1,"%s%s",cdf_ev,popStr_ev4
wave w1a=$s1
StatsQuantiles /Q w1a
string stat ; sprintf stat,"%sW_StatsQuantiles",cdf_ev    ; wave stat1=$stat
variable /G  root:Params:RoVibHamiltonian:pd_vmin =stat1[0]
variable /G  root:Params:RoVibHamiltonian:int_min =stat1[0]

variable /G  root:Params:RoVibHamiltonian:pd_vmax =stat1[1]
killwaves /Z stat1
End

//*********************************************************************************************************************
//BATCH NORM

function Batch_norm_disatance(b_nd,pop_bnd,popStr_bnd) : PopupMenuControl	//parameter (eg; distance for wavefunctions)
String b_nd
Variable pop_bnd
String popStr_bnd
string /G root:Params:RoVibHamiltonian:b_norm_distance = popStr_bnd
End

//*********************************************************************************************************************

//BATCH EXP VAL
Function BExpVal_parameter_Selection(bev_sel3,popNum_bev3,popStr_bev3) : PopupMenuControl	//parameter (eg; gamma)
String bev_sel3
Variable popNum_bev3
String popStr_bev3
string /G root:Params:RoVibHamiltonian:b_gwave = popStr_bev3
End

//*********************************************************************************************************************
//BATCH EXP VAL
Function BExpVal_parameter_Selection_d(bev_sel31,popNum_bev31,popStr_bev31) : PopupMenuControl	//parameter (eg; gamma)
String bev_sel31
Variable popNum_bev31
String popStr_bev31
string /G root:Params:RoVibHamiltonian:b_gWave_dis = popStr_bev31
End

//BATCH EXP VAL
Function BExpVal_wfn_distance_Selection(bev_wd,popNum_wd,popStr_wd) : PopupMenuControl	//parameter (eg; gamma)
String bev_wd
Variable popNum_wd
String popStr_wd
string /G root:Params:RoVibHamiltonian:b_wd =popStr_wd
End
//*********************************************************************************************************************
//OtherAtom_PES_R_wave SELECTION
Function pes_rwave__Selection(HPES_R,popNum_P1,popStr_PES_R) : PopupMenuControl	//parameter (eg; gamma)
String HPES_R
Variable popNum_P1
String  popStr_PES_R
string /G  root:Params:RoVibHamiltonian:potential_rwave =  popStr_PES_R
End

Function pes_energy__Selection(HPES_E,popNum_P2,popStr_PES_E) : PopupMenuControl	//parameter (eg; gamma)
String HPES_E
Variable popNum_P2
String  popStr_PES_E
string /G   root:Params:RoVibHamiltonian:potential_energywave =  popStr_PES_E
End

//*********************************************************************************************************************
//Store all constants here; 
//masses of atoms , reduced mass of molecules, etc....

function constants_val()   
newdatafolder /s /O root:constants
variable /G H2_nu = 918.076336945000000000
variable /G D2_nu = 1835.241483925000000000
variable /G HD_nu = 1223.899228923000000000

variable /G mProton = 1836.15267389000000
variable /G mDeuteron = 3670.48296785000000
setdatafolder root:
end

////*******************************************************************************************************************
//GENERATE COEFFICIENTS FOR THE ASYMMETRIC FIRST AND SECOND DERIVATIVES
// ROW 1 AND ROW2  AND MID-

function gen_coef_R1(n)
variable n;
make /O /D /n=(5,5) A=0
make  /O /D /n=5 C=0
C[n]=1 
make /O /D /n=5 D_coef
variable i,j

for(i=0 ; i<5 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<5 ; j=j+1)
			A[j][i]=(i^j) ;
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

// print A

string cmd; sprintf cmd,"matrixop /O D_coefR1= inv(A) x C "
//print cmd
execute cmd
wave D_coefR1 = D_coefR1
// print D_coefR1
end
//=======================================
//ROW-2
function gen_coef_R2(n)
variable n;
make /O /D /n=(5,5) A=0
make  /O /D /n=5 C=0
C[n]=1 
make /O /D /n=5 D_coefR2

variable i,j

for(i=0 ; i<5 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<5 ; j=j+1)
			A[j][i]=(i^j) ;
			if (j>0)
			A[j][i]=((i-1)^j) ;
			endif
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefR2= inv(A) x C "
//print cmd
execute cmd
end
//=======================================
// Mid rows with symmetric derivatives(5 points)
function gen_coef_R3(n)
variable n;
make /O /D /n=(5,5) A=0
make  /O /D /n=5 C=0
C[n]=1 
make /O /D /n=5 D_coefR3

variable i,j

for(i=0 ; i<5 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<5 ; j=j+1)
			A[j][i]=(i^j) ;
			if (j>0)
			A[j][i]=((i-2)^j) ;
			endif
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefR3= inv(A) x C "
//print cmd
execute cmd
end
//=======================================
//Last row: asymmetric

function gen_coef_Rn(n)
variable n;
make /O /D /n=(5,5) A=0
make  /O /D /n=5 C=0
C[n]=1 
make /O /D /n=5 D_coefRn

variable i,j


for(i=0 ; i<5 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<5 ; j=j+1)
			A[j][i]=(i^j) ;
			
			if (j==1 || j==3)
			A[j][i]= -A[j][i] ;
			endif
			
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefRn= inv(A) x C "
//print cmd
execute cmd
end
//=======================================
//Second last row: asymmetric
function gen_coef_Rn2(n)
variable n;
make /O /D /n=(5,5) A=0
make  /O /D /n=5 C=0
C[n]=1 
make /O /D /n=5 D_coefRn2

variable i,j


for(i=0 ; i<5 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<5 ; j=j+1)
			A[j][i]=(i^j) ;
			if (j>0)
			A[j][i]=((i-1)^j) ;
			endif
			A[][0]= abs (A[p][0])
			
			if (j== 1 || j ==3)
			A[j][i]= -A[j][i] ;
			endif
			 A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefRn2= inv(A) x C "
//print cmd
execute cmd
end
//********************************************************************************************************************
//GENERATE COEFFICIENTS FOR THE ASYMMETRIC FIRST AND SECOND DERIVATIVES
// with 3 point numerical differentiation scheme (upto second derivative term in Taylor Series)
// ROW 1
function gen_coef3term_R1(n)
variable n;
make /O /D /n=(3,3) A=0
make  /O /D /n=3 C=0
C[n]=1 
make /O /D /n=3 D_coef
variable i,j

for(i=0 ; i<3 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<3 ; j=j+1)
			A[j][i]=(i^j) ;
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefR1= inv(A) x C "
//print cmd
execute cmd
wave D_coefR1 = D_coefR1
// print D_coefR1
end
//=======================================
//ROW-2 // 4 term approximation, i.e. upto the third derivative.
function gen_coef4term_R2(n)
variable n;
make /O /D /n=(4,4) A=0
make  /O /D /n=4 C=0
C[n]=1 
make /O /D /n=4 D_coefR2

variable i,j

for(i=0 ; i<4 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<4 ; j=j+1)
			A[j][i]=(i^j) ;
			if (j>0)
			A[j][i]=((i-1)^j) ;
			endif
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 
//print A
string cmd; sprintf cmd,"matrixop /O D_coefR2= inv(A) x C "
//print cmd
execute cmd
//print D_coefR2
end
//=======================================
//Second last row: asymmetric
// 4 term approximation for numerical derivative, i.e. upto the 3rd derivative in the Taylor series expansion.
function gen_coef_Rn2_4term(n)
variable n;
make /O /D /n=(4,4) A=0
make  /O /D /n=4 C=0
C[n]=1 
make /O /D /n=4 D_coefRn2

variable i,j


for(i=0 ; i<4 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<4 ; j=j+1)
			A[j][i]=(i^j) ;
			if (j>0)
				A[j][i]=((i-2)^j) ;
			endif
			 A[j][i] = A[j][i] / factorial (j)
		endfor
		
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefRn2= inv(A) x C "
//print cmd
execute cmd
// print A
print D_coefRn2
end
//=======================================
//Second last row: asymmetric
// 4 term approximation for numerical derivative, i.e. upto the 3rd derivative in the Taylor series expansion.
function gen_coef_Rn_3term(n)
variable n;
make /O /D /n=(3,3) A=0
make  /O /D /n=3 C=0
C[n]=1 
make /O /D /n=3 D_coefRn

variable i,j

for(i=0 ; i<3 ; i=i+1)
	A[][i]=i
		for (j=0 ; j<3 ; j=j+1)
			A[j][i]=(i^j) ;
			if (j>0)
				A[j][i]=((i-2)^j) ;
			endif
			A[j][i] = A[j][i] / factorial (j)
		endfor
endfor 

string cmd; sprintf cmd,"matrixop /O D_coefRn= inv(A) x C "
//print cmd
execute cmd
// print A
print D_coefRn
end
//********************************************************************************************************************
//Started :12/13/2015 _ Function to prepare the HamiltonianMatrix for Diatomic molecule.
// Ref. eqn (1) on research noteook, page 7

function MatrixPrep()
variable J,mA,mB
string r_w, pot_wave

//Clear old variables , waves
Killwaves /Z H_matrix,A,C,D_coef,M_inverse

Killwaves /Z A,C,D_coef,M_inverse

prompt J,"Rotational Level"
prompt mA,"Mass atom A"
prompt mB,"Mass atom B"
Prompt r_w,  "Internuclear distance wave", popup WaveList("*", ";", "")
Prompt pot_wave,  "Potential function wave", popup WaveList("*", ";", "")
doprompt "Input data about the diatomic molecule", J,mA,mB,r_w,pot_wave
			if (V_Flag)
			return -1 			// User canceled
			endif
// calculate Nuclear reduced mass
variable nu =1/ (1/mA + 1/mB)
//printf "nuclear reduced mass : %g\r", nu
wave r_wave1 = $r_w ; wave p_wave = $pot_wave ;

printf "Dimensions : %g x %g \r",  dimsize(r_wave1,0), dimsize(r_wave1,0)
make /O /D /n=(dimsize(r_wave1,0), dimsize(r_wave1,0)) H_matrix
wave HM =H_matrix
variable dimn = dimsize(r_wave1,0)

// The J term and potential term
variable i=0, J_val, P_val, C_term
for (i=0  ; i < ((dimsize(r_wave1,0))) ; i =i+1)
J_val=0 ;  P_val=0 ;  C_term=0
J_val = (J*(J+1)) / (2*nu*(r_wave1[i])^2 )
P_val = p_wave [i]
C_term = J_val + P_val
HM[i][i] = C_term
endfor 

//FIRST DERIVATIVE******************************************
variable st,sta,val

// 1. THE FIRST AND LAST POINTS

// FIRST TWO ROWS
variable RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val= -1/(nu*RI)
killwaves /Z  D_coefR1,D_coefR2
gen_coef_R1(1) // Generate first derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<5 ; i=i+1)
HM[0][i]=D_coefR1[i]* (1/st)*val + HM[0][i]

endfor

gen_coef_R2(1) // Generate first derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<5 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st)*val + HM[1][i]
endfor

//LAST TWO ROWS
RI=r_wave1[((dimn)-1)]
val= -1/(nu*RI)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn(1) // Generate first derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st)*val + HM[(dimn-1)][(dimn-1-i)]
endfor

 //Second last row 
 RI=r_wave1[((dimn)-1)]
val= -1/(nu*RI)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn2(1) // Generate first derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st)*val + HM[(dimn-2)][(dimn-1-i)]
endfor
//==============

// mid - terms in the matrix, 5 terms centered around diagonal
variable ista,iend,col,e1,VA
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend

gen_coef_R3(1) // Generate first derv. central difference 5 terms
wave D_coefR3 =D_coefR3

for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)] // STEP SIZE 'h' for the 5 term DERIVATIVES
val= -1/(nu*RI)
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		VA =val*D_coefR3[e1]* (1/st)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor

endfor
//FIRST DERIVATIVE TERMS FINISH HERE ------------------------------------------------------------------------------------------

//SECOND DERIVATIVE**************************************************************************
// 1. THE FIRST AND LAST POINTS
// FIRST TWO ROWS
RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val= -1/(nu*2)
killwaves /Z  D_coefR1,D_coefR2
gen_coef_R1(2) // Generate SECOND  derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<5 ; i=i+1)
HM[0][i]=D_coefR1[i]* (1/st^2)*val + HM[0][i]
endfor

gen_coef_R2(2) // Generate  SECOND  derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<5 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st^2)*val + HM[1][i]
endfor


//LAST TWO ROWS
RI=r_wave1[((dimn)-1)]
val= -1/(nu*2)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn(2) // Generate SECOND  derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st^2)*val + HM[(dimn-1)][(dimn-1-i)]
endfor


 //Second last row 
 RI=r_wave1[((dimn)-1)]
val= -1/(nu*2)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn2(2) // Generate   SECOND  derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st^2)*val + HM[(dimn-2)][(dimn-1-i)]
endfor

//==============

// 2. TERMS IN BETWEEN for SECOND DERIVATIVE - 5 terms 
// mid - terms in the matrix, 5 terms centered around diagonal
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend
killwaves /Z  D_coefR3
gen_coef_R3(2) // Generate SECOND  derv. central difference 5 terms
wave D_coefR3 =D_coefR3


for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)]  ;      // STEP SIZE 'h' for the DERIVATVES 5 term
val= -1/(nu*2)
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		
		VA =val *D_coefR3[e1]* (1/st^2)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor
endfor
//==============



printf "J : %g | mA : %6.11f, mB : %6.11f | nu : %6.11f | H-Matrix prepared\r",J,mA,mB,nu
DFREF cdf1 = GetDataFolderDFR()
newdatafolder /O /S Calc
killwaves /z HMat1,HMat, r_wave
Killwaves /Z M_L_eigenVectors, M_R_eigenVectors
DFREF cdf2 = GetDataFolderDFR()
//setdatafolder cdf1
//duplicate /O  H_matrix , HMat
make /O /D /n=(dimsize(r_wave1,0), dimsize(r_wave1,0))  HMat1 ; wave CHM=HMat1
make /O /D /n=(dimsize(r_wave1,0))   r_wave = r_wave1
CHM = HM ; 
setdatafolder cdf1
end
//********************************************************************************************************************
//********************************************************************************************************************

//12/13/2015 _ Function to prepare the HamiltonianMatrix for a simple harmonic oscillator


function MatrixPrep_SHO()
variable xn
string r_w, pot_wave

//Clear old variables , waves
Killwaves /Z H_matrix, H_matrix_SHO
Killwaves /Z M_L_eigenVectors, M_R_eigenVectors,W_eigenValues,W_MatrixRCONDE,W_MatrixRCONDV
Killwaves /Z A,C,D_coef,M_inverse


prompt xn,"n in (x-n) term "
Prompt r_w,  "Internuclear distance wave", popup WaveList("*", ";", "")
doprompt "Input data about the SHO", xn,r_w
			if (V_Flag)
			return -1 			// User canceled
			endif

wave r_wave1 = $r_w ; 

printf "Dimensions : %g x %g \r",  dimsize(r_wave1,0), dimsize(r_wave1,0)
make /O /D /n=(dimsize(r_wave1,0), dimsize(r_wave1,0)) H_matrix_SHO
wave HM =H_matrix_SHO
variable dimn = dimsize(r_wave1,0)

// The (x-n) term 
variable i=0, J,  C_term
for (i=0  ; i < ((dimsize(r_wave1,0))) ; i =i+1)
J=0 ;   C_term=0
J = r_wave1[i]

C_term =10*((J-xn)^2)
HM[i][i] = C_term
endfor 


//=========================================


//SECOND DERIVATIVE**************************************************************************
// 1. THE FIRST AND LAST POINTS
// FIRST TWO ROWS
variable RI,st,val
RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val=-1/2
killwaves /Z  D_coefR1,D_coefR2
gen_coef_R1(2) // Generate SECOND  derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<5 ; i=i+1)
HM[0][i]=D_coefR1[i]* (1/st^2)*val + HM[0][i]
endfor

gen_coef_R2(2) // Generate  SECOND  derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<5 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st^2)*val + HM[1][i]
endfor


//LAST TWO ROWS
RI=r_wave1[((dimn)-1)]
val=-1/2
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn(2) // Generate SECOND  derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st^2)*val + HM[(dimn-1)][(dimn-1-i)]
endfor


 //Second last row 
 RI=r_wave1[((dimn)-1)]
val=-1/2
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn2(2) // Generate   SECOND  derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st^2)*val + HM[(dimn-2)][(dimn-1-i)]
endfor
//==============
variable ista, iend,e1,VA,col
// 2. TERMS IN BETWEEN for SECOND DERIVATIVE - 5 terms 
// mid - terms in the matrix, 5 terms centered around diagonal
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend

gen_coef_R3(2) // Generate SECOND  derv. central difference 5 terms
wave D_coefR3 =D_coefR3


for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)]     // STEP SIZE 'h' for the DERIVATVES 5 term
val=-1/2
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		
		VA =val *D_coefR3[e1]* (1/st^2)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor
endfor
//==============
//=========================================


//    abort
DFREF cdf1 = GetDataFolderDFR()
newdatafolder /O /S Calc
killwaves /z HMat
DFREF cdf2 = GetDataFolderDFR()
//setdatafolder cdf1
//duplicate /O  H_matrix , HMat
make /O /D /n=(dimsize(r_wave1,0), dimsize(r_wave1,0))  HMat ; wave CHM=HMat
CHM = H_matrix_SHO ; 
setdatafolder cdf1
end
//********************************************************************************************************************

//--------------------------------------------------------------------------------------------------------------------------------
function sorter_eigval_eigvec()
//wave org_eigval, org_eigvec
wave org_eigval = W_eigenValues
wave org_eigvec=M_R_eigenVectors

variable dim= dimsize (org_eigval,0)
wavestats /Q org_eigval
make  /D /O /n=(dim) eval=org_eigval ;  make  /D /O /n=(dim,dim) evec=0;
make  /D /O /n=(dim) evalO=org_eigval
wave eval = eval  ; wave evec = evec 
sort eval,eval
//print dim
variable i,j,O,V
make /D /O /n=(dim) indexw=0;  wave indexw = indexw
	for (i=0 ; i < (dim) ; i=i+1)
		O=eval[i]
		findvalue  /V=(O) /T=0.001  evalO
	//	print V_value,O
		if (V_value != -1)
			V=V_value
			indexw[i] = V
			evec[][i] = org_eigvec[p][V]
		endif
	endfor
end
//--------------------------------------------------------------------------------------------------------------------------------

function auto_solveJ_eigenvalue(x1,x2)
variable x1,x2

variable i; string w="J" ; string df=getdatafolder(1)
// print df
for (i=0 ; i <= (x2) ; i=i+1)
	setdatafolder $df
	sprintf w,"J%g",i
// print w
	setdatafolder $w
	setdatafolder Calc
	printf "Started J%g\r",i
	string cmd = "MatrixEigenV   /R /S=2  HMat1"
	execute cmd
	sorter_eigval_eigvec()
	killwaves /Z HMat1
endfor
setdatafolder $df
print "Done !"
end

//********************************************************************************

function auto_solveJ_eigenval (s7a):ButtonControl
string s7a
variable t0=ticks
nvar x1= root:Params:RoVibHamiltonian:j_min ; nvar x2= root:Params:RoVibHamiltonian:j_max
printf "Eigenvalue-eigenvector solution | J=%g to J=%g\r" x1,x2
variable i=x1; string w="J" ; string df=getdatafolder(1)
// print df
for (i=0 ; i <= (x2) ; i=i+1)
	setdatafolder $df
	
	sprintf w,"J%g",i
// print w
	setdatafolder $w
	wave HM= H_matrix
	setdatafolder Calc
	printf "Started J%g. ",i
	MatrixEigenV   /R  /S=2  HM
	sorter_eigval_eigvec()

endfor
setdatafolder $df
printf "Done ! ( %g sec )\r",(ticks-t0)/60
end
//--------------------------------------------------------------------------------------------------------------------------------

//********************************************************************************
//**********************************************************************
//********************************************************************************
// Defining A and B defines the atom as H or D. 	// RUN FROM COMMAND LINE
//A=1 : proton nuclei  ; =2 : Deuteron nuclei
//B=1 : proton nuclei   ; =2 : Deuteron nuclei
//	x2= highest J level for which matrix is to be generated
function JfoldersPrep_N (x2,A,B)
variable x2,A,B

variable mA,mB ; string mol
if (A==1)
	mA= 1836.15267406000 +1 ; 
	if (B == 1)
		mB= 1836.152673890000 +1 ;  mol="H-H"
	elseif ( B==2 )
		mB = 3670.48296785000000 +1 ;	 mol="H-D"
	endif
endif
	
if (A==2  && B==2 )
		mA = 3670.48296785000000 +1 ; mB = 1836.152673890000 +1 ;	 mol="D-D"
endif

// printf "%s ( %g , %g)\r" mol,mA,mB				// RUN FROM COMMAND LINE
string cdf=getdatafolder(1)
wave r_wave1 = r_wave
wave scaled_y1 = scaled_y
variable dim= dimsize(r_wave1,0)
printf "Dim: %g x %g , %s \r",dim,dim,mol
variable i; string w

for (i=0 ; i <= (x2) ; i=i+1)
sprintf w,"J%g",i
newdatafolder /S $w
make /O /D /n=(dim) r_wave0 =r_wave1 
make /O /D /n=(dim) scaled_y0 =scaled_y1
string r1="r_wave0"   ;  string p1= "scaled_y0"

MatrixPrep_auto(i,  mA,mB ,r1, p1)
setdatafolder $cdf
endfor

// H2 =   1836.152673890000000 ,1836.152673890000000
//D2=   3670.48296785000000   , 3670.48296785000000
//HD =   1836.152673890000000 ,  3670.48296785000000		// RUN FROM COMMAND LINE
end

//--------------------------------------------------------------------------------------------------------------------------------
//********************************************************************************
// NUCLEAR MASS HAVE ERROR ADDED FOR ESTIAMTING ERROR-------------------------- RUN WITH CAUTION
// Defining A and B defines the atom as H or D. 	// RUN FROM COMMAND LINE
//A=1 : proton nuclei  ; =2 : Deuteron nuclei
//B=1 : proton nuclei   ; =2 : Deuteron nuclei
//	x2= highest J level for which matrix is to be generated
function JfoldersPrep_Nmplus (x2,A,B)
variable x2,A,B

variable mA,mB ; string mol
if (A==1)
	mA= 1836.15267406000 +1.0000000000000 ; 
	if (B == 1)
		mB= 1836.152674060000 +1.000000000000 ;  mol="H-H"
	elseif ( B==2 )
		mB = 3670.48296798000000 +1.0000000000000 ;	 mol="H-D"
	endif
endif
	
if (A==2  && B==2 )
		mA = 3670.48296798000000 +1.000000000 ; mB =3670.48296798000000  +1.00000000000000 ;	 mol="D-D"
endif

// printf "%s ( %g , %g)\r" mol,mA,mB				// RUN FROM COMMAND LINE
string cdf=getdatafolder(1)
wave r_wave1 = r_wave
wave scaled_y1 = scaled_y
variable dim= dimsize(r_wave1,0)
printf "Dim: %g x %g , %s \r",dim,dim,mol
variable i; string w

for (i=0 ; i <= (x2) ; i=i+1)
sprintf w,"J%g",i
newdatafolder /S $w
make /O /D /n=(dim) r_wave0 =r_wave1 
make /O /D /n=(dim) scaled_y0 =scaled_y1
string r1="r_wave0"   ;  string p1= "scaled_y0"

MatrixPrep_auto(i,  mA,mB ,r1, p1)
setdatafolder $cdf
endfor

// H2 =   1836.152673890000000 ,1836.152673890000000
//D2=   3670.48296785000000   , 3670.48296785000000
//HD =   1836.152673890000000 ,  3670.48296785000000		// RUN FROM COMMAND LINE
end

//--------------------------------------------------------------------------------------------------------------------------------
//********************************************************************************

// NUCLEAR MASS HAVE ERROR SUBTRACTED FOR ESTIAMTING ERROR-------------------------- RUN WITH CAUTION
// Defining A and B defines the atom as H or D. 	// RUN FROM COMMAND LINE
//A=1 : proton nuclei  ; =2 : Deuteron nuclei
//B=1 : proton nuclei   ; =2 : Deuteron nuclei
//	x2= highest J level for which matrix is to be generated
function JfoldersPrep_Nmsub (x2,A,B)
variable x2,A,B

variable mA,mB ; string mol
if (A==1)
	mA= 1836.15267372000 +1.00000000000 ; 
	if (B == 1)
		mB= 1836.1526737200000 +1.000000000000 ;  mol="H-H"
	elseif ( B==2 )
		mB = 3670.48296772000000 +1.000000000000 ;	 mol="H-D"
	endif
endif
	
if (A==2  && B==2 )
		mA = 3670.48296772000000 +1.000000000 ; mB = 3670.48296772000000 +1.00000000000000 ;	 mol="D-D"
endif

// printf "%s ( %g , %g)\r" mol,mA,mB				// RUN FROM COMMAND LINE
string cdf=getdatafolder(1)
wave r_wave1 = r_wave
wave scaled_y1 = scaled_y
variable dim= dimsize(r_wave1,0)
printf "Dim: %g x %g , %s \r",dim,dim,mol
variable i; string w

for (i=0 ; i <= (x2) ; i=i+1)
sprintf w,"J%g",i
newdatafolder /S $w
make /O /D /n=(dim) r_wave0 =r_wave1 
make /O /D /n=(dim) scaled_y0 =scaled_y1
string r1="r_wave0"   ;  string p1= "scaled_y0"

MatrixPrep_auto(i,  mA,mB ,r1, p1)
setdatafolder $cdf
endfor

// H2 =   1836.152673890000000 ,1836.152673890000000
//D2=   3670.48296785000000   , 3670.48296785000000
//HD =   1836.152673890000000 ,  3670.48296785000000		// RUN FROM COMMAND LINE
end

//--------------------------------------------------------------------------------------------------------------------------------
//********************************************************************************
// Defining A and B defines the atom as H or D.		// RUN FROM PANEL
//A=1 : proton nuclei  ; =2 : Deuteron nuclei
//B=1 : proton nuclei   ; =2 : Deuteron nuclei
// WORKS FROM THE PANEL :::
function JfoldersPrep_N2 (s7b):ButtonControl
string s7b

nvar x2= root:Params:RoVibHamiltonian:pmax; nvar A= root:Params:RoVibHamiltonian:atom1 ; nvar B=root:Params:RoVibHamiltonian:atom2 ; 
variable mA,mB ; string mol
if (A==1)
	mA= ( 1836.15267389000+  1.000000000000) ;  //mass H atom
	if (B == 1)
		mB= (1836.152673890000+  1.000000000000) ;  mol="H-H"
	elseif ( B==2 )
		mB = (3670.48296785000000+  1.000000000000) ;	 mol="H-D"
	endif
endif
	
if (A==2  && B==2 )
		mA = (3670.48296785000000+  1.000000000000000000  ) ; mB = (3670.48296785000000+ 1.000000000000000000) ;	 mol="D-D"
endif

// printf "%s ( %g , %g)\r" mol,mA,mB			// RUN FROM PANEL
string cdf=getdatafolder(1)
wave r_wave1 = r_wave
wave scaled_y1 = scaled_y
variable dim= dimsize(r_wave1,0)
printf "Dim: %g x %g , %s \r",dim,dim,mol
variable i; string w

for (i=0 ; i <= (x2) ; i=i+1)
sprintf w,"J%g",i
newdatafolder /S $w
make /O /D /n=(dim) r_wave0 =r_wave1 
make /O /D /n=(dim) scaled_y0 =scaled_y1
string r1="r_wave0"   ;  string p1= "scaled_y0"

MatrixPrep_auto(i,  mA,mB ,r1, p1)
setdatafolder $cdf
endfor
// H2 =   1836.152673890000000 ,1836.152673890000000			// RUN FROM PANEL
//D2=   3670.48296785000000   , 3670.48296785000000
//HD =   1836.152673890000000 ,  3670.48296785000000			// RUN FROM PANEL
end

//--------------------------------------------------------------------------------------------------------------------------------
//********************************************************************************
// Defining A and B defines the atom as H or D.		// RUN FROM PANEL
//A=1 : proton nuclei  ; =2 : Deuteron nuclei
//B=1 : proton nuclei   ; =2 : Deuteron nuclei
// WORKS FROM THE PANEL :::
function JfoldersPrep_SymmH (s7b):ButtonControl  // symmetric H 
string s7b

nvar x2= root:Params:RoVibHamiltonian:pmax; nvar A= root:Params:RoVibHamiltonian:atom1 ; nvar B=root:Params:RoVibHamiltonian:atom2 ; 
variable mA,mB ; string mol
if (A==1)
	mA= ( 1836.15267389000+  1.000000000000) ;  //mass H atom
	if (B == 1)
		mB= (1836.152673890000+  1.000000000000) ;  mol="H-H"
	elseif ( B==2 )
		mB = (3670.48296785000000+  1.000000000000) ;	 mol="H-D"
	endif
endif
	
if (A==2  && B==2 )
		mA = (3670.48296785000000+  1.000000000000000000  ) ; mB = (3670.48296785000000+ 1.000000000000000000) ;	 mol="D-D"
endif

// printf "%s ( %g , %g)\r" mol,mA,mB			// RUN FROM PANEL
string cdf=getdatafolder(1)
wave r_wave1 = r_wave
wave scaled_y1 = scaled_y
variable dim= dimsize(r_wave1,0)
printf "Dim: %g x %g , %s \r",dim,dim,mol
variable i; string w

for (i=0 ; i <= (x2) ; i=i+1)
sprintf w,"J%g",i
newdatafolder /S $w
make /O /D /n=(dim) r_wave0 =r_wave1 
make /O /D /n=(dim) scaled_y0 =scaled_y1
string r1="r_wave0"   ;  string p1= "scaled_y0"

MatrixPrep_autoSymmetric(i,  mA,mB ,r1, p1)			// for symmetric H 
setdatafolder $cdf
endfor
// H2 =   1836.152673890000000 ,1836.152673890000000			// RUN FROM PANEL
//D2=   3670.48296785000000   , 3670.48296785000000
//HD =   1836.152673890000000 ,  3670.48296785000000			// RUN FROM PANEL
end

//--------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------

//Enter variables to make the hamiltonian matrix for the ro-vibration state

function MatrixPrep_auto(J,mA,mB,r_w, pot_wave)

        variable J,mA,mB
        string r_w, pot_wave

        //Clear old variables , waves
        Killwaves /Z H_matrix,A,C,D_coef,M_inverse
        Killwaves /Z A,C,D_coef,M_inverse

        // calculate Nuclear reduced mass
        variable nu =1/ (1/mA + 1/mB)
        // printf "Reduced mass : %6.8f\r", nu

        wave r_wave1 = $r_w ; wave p_wave = $pot_wave ;

        //	printf "Dimensions : %g x %g \r",  dimsize(r_wave1,0), dimsize(r_wave1,0)
        make /O /D /n=(dimsize(r_wave1,0), dimsize(r_wave1,0)) H_matrix
        wave HM =H_matrix
        variable dimn = dimsize(r_wave1,0)

        // The J term and potential term
        variable i=0, J_val, P_val, C_term
        
        for (i=0  ; i < ((dimsize(r_wave1,0))) ; i =i+1)
        	J_val=0 ;  P_val=0 ;  C_term=0
	        J_val = (J*(J+1)) / (2*nu*(r_wave1[i]^2) )
	        P_val = p_wave [i]
	        C_term = J_val + P_val
	        HM[i][i] = C_term
        endfor 

//FIRST DERIVATIVE******************************************
variable st,sta,val

// 1. THE FIRST AND LAST POINTS

// FIRST TWO ROWS
variable RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val= -1/(nu*RI)
killwaves /Z  D_coefR1,D_coefR2
gen_coef_R1(1) // Generate first derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<5 ; i=i+1)
HM[0][i]=D_coefR1[i]* (1/st)*val + HM[0][i]
//	print "-------------------"
//	print D_coefR1[i]
//	print "-------------------"
endfor

gen_coef_R2(1) // Generate first derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<5 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st)*val + HM[1][i]
endfor

//LAST TWO ROWS
RI=r_wave1[((dimn)-1)]
val= -1/(nu*RI)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn(1) // Generate first derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st)*val + HM[(dimn-1)][(dimn-1-i)]
endfor

 //Second last row 
 RI=r_wave1[((dimn)-1)]
val= -1/(nu*RI)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn2(1) // Generate first derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st)*val + HM[(dimn-2)][(dimn-1-i)]
endfor
//==============

// mid - terms in the matrix, 5 terms centered around diagonal
variable ista,iend,col,e1,VA
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
//	printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_R3(1) // Generate first derv. central difference 5 terms
wave D_coefR3 =D_coefR3

for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)] // STEP SIZE 'h' for the 5 term DERIVATIVES
val= -1/(nu*RI)
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		VA =val*D_coefR3[e1]* (1/st)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor

endfor
//FIRST DERIVATIVE TERMS FINISH HERE ------------------------------------------------------------------------------------------

//SECOND DERIVATIVE**************************************************************************
// 1. THE FIRST AND LAST POINTS
// FIRST TWO ROWS
RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val= -1/(nu*2)
killwaves /Z  D_coefR1,D_coefR2
gen_coef_R1(2) // Generate SECOND  derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<5 ; i=i+1)
HM[0][i]=D_coefR1[i]* (1/st^2)*val + HM[0][i]
endfor

killwaves /z D_coefR2
gen_coef_R2(2) // Generate  SECOND  derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<5 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st^2)*val + HM[1][i]
endfor


//LAST TWO ROWS
RI=r_wave1[((dimn)-1)]
val= -1/(nu*2)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn(2) // Generate SECOND  derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st^2)*val + HM[(dimn-1)][(dimn-1-i)]
endfor


 //Second last row 
 RI=r_wave1[((dimn)-1)]
val= -1/(nu*2)
killwaves /Z  D_coefRn,D_coefRn2,D_coefR3
gen_coef_Rn2(2) // Generate   SECOND  derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<5 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st^2)*val + HM[(dimn-2)][(dimn-1-i)]
endfor

//==============

// 2. TERMS IN BETWEEN for SECOND DERIVATIVE - 5 terms 
// mid - terms in the matrix, 5 terms centered around diagonal
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
//printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend
killwaves /Z  D_coefR3
gen_coef_R3(2) // Generate SECOND  derv. central difference 5 terms
wave D_coefR3 =D_coefR3


for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)]  ;      // STEP SIZE 'h' for the DERIVATVES 5 term
val= -1/(nu*2)
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		
		VA =val *D_coefR3[e1]* (1/st^2)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor
endfor
//==============
killwaves /Z  D_coefR1,D_coefR2,D_coefRn2,D_coefR3,M_inverse,D_coef,A,C


printf "J :  %g | mA :  %6.11f, mB : %6.11f | nu : %6.11f | H-Matrix prepared\r",J,mA,mB,nu
DFREF cdf1 = GetDataFolderDFR()
newdatafolder /O /S Calc
killwaves /z HMat1,HMat, r_wave
Killwaves /Z M_L_eigenVectors, M_R_eigenVectors
DFREF cdf2 = GetDataFolderDFR()
//setdatafolder cdf1
//duplicate /O  H_matrix , HMat
make /O /D /n=(dimsize(r_wave1,0))   r_wave = r_wave1
setdatafolder cdf1
end


//*******************************************************************************************************************
//--------------------------------------------------------------------------------------------------------------------------------

//Enter variables to make the hamiltonian matrix for the ro-vibration state

function MatrixPrep_autoSymmetric(J,mA,mB,r_w, pot_wave)
variable J,mA,mB
string r_w, pot_wave

//Clear old variables , waves
Killwaves /Z H_matrix,A,C,D_coef,M_inverse
Killwaves /Z A,C,D_coef,M_inverse

// calculate Nuclear reduced mass
variable nu =1/ (1/mA + 1/mB)
//printf "nuclear reduced mass : %g\r", nu

wave r_wave1 = $r_w ; wave p_wave = $pot_wave ;

//	printf "Dimensions : %g x %g \r",  dimsize(r_wave1,0), dimsize(r_wave1,0)
make /O /D /n=(dimsize(r_wave1,0), dimsize(r_wave1,0)) H_matrix
wave HM =H_matrix
variable dimn = dimsize(r_wave1,0)

// The J term and potential term
variable i=0, J_val, P_val, C_term
for (i=0  ; i < ((dimsize(r_wave1,0))) ; i =i+1)
J_val=0 ;  P_val=0 ;  C_term=0
J_val = (J*(J+1)) / (2*nu*(r_wave1[i]^2) )
P_val = p_wave [i]
C_term = J_val + P_val
HM[i][i] = C_term
endfor 



//FIRST DERIVATIVE******************************************
variable st,sta,val

// 1. THE FIRST AND LAST POINTS

// FIRST TWO ROWS
variable RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val= -1/(nu*RI)
killwaves /Z  D_coefR1,D_coefR2
gen_coef3term_R1(1) // Generate first derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<3 ; i=i+1)
	HM[0][i]=D_coefR1[i]* (1/st)*val + HM[0][i]
endfor

gen_coef4term_R2(1) // Generate first derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<4 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st)*val + HM[1][i]
endfor

//LAST TWO ROWS

// last row
RI=r_wave1[((dimn)-1)]
val= -1/(nu*RI)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn_3term(1) // Generate first derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<3 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st)*val + HM[(dimn-1)][(dimn-1-i)]
endfor

 //Second last row 
 RI=r_wave1[((dimn)-1)]
val= -1/(nu*RI)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn2_4term(1) // Generate first derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<4 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st)*val + HM[(dimn-2)][(dimn-1-i)]
endfor
//==============

// mid - terms in the matrix, 5 terms centered around diagonal
variable ista,iend,col,e1,VA
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
//	printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_R3(1) // Generate first derv. central difference 5 terms
wave D_coefR3 =D_coefR3

for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)] // STEP SIZE 'h' for the 5 term DERIVATIVES
val= -1/(nu*RI)
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		VA =val*D_coefR3[e1]* (1/st)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor

endfor
//FIRST DERIVATIVE TERMS FINISH HERE ------------------------------------------------------------------------------------------

//SECOND DERIVATIVE**************************************************************************
// 1. THE FIRST AND LAST POINTS
// FIRST TWO ROWS
RI =r_wave1[0]   // first point value 
st = r_wave1[1]- r_wave1[0] // h - value
val= -1/(nu*2)
killwaves /Z  D_coefR1,D_coefR2
gen_coef3term_R1(2) // Generate SECOND  derv. 1st row
wave D_coefR1 =D_coefR1
for(i=0 ; i<3 ; i=i+1)
HM[0][i]=D_coefR1[i]* (1/st^2)*val + HM[0][i]
endfor

killwaves /z D_coefR2
gen_coef_R2(2) // Generate  SECOND  derv. 2nd row
wave D_coefR2 =D_coefR2
for(i=0 ; i<4 ; i=i+1)
HM[1][i]=D_coefR2[i]* (1/st^2)*val + HM[1][i]
endfor


//LAST TWO ROWS

// last row
RI=r_wave1[((dimn)-1)]
val= -1/(nu*2)
killwaves /Z  D_coefRn,D_coefRn2
gen_coef_Rn_3term(2) // Generate SECOND  derv. last  row
wave D_coefRn =D_coefRn
for(i=0 ; i<3 ; i=i+1)
HM[(dimn-1)][(dimn-1-i)]=D_coefRn[i]* (1/st^2)*val + HM[(dimn-1)][(dimn-1-i)]
endfor


 //Second last row 
 RI=r_wave1[((dimn)-1)]
val= -1/(nu*2)
killwaves /Z  D_coefRn,D_coefRn2,D_coefR3
gen_coef_Rn2_4term(2) // Generate   SECOND  derv. second last  row
 wave D_coefRn2 =D_coefRn2
for(i=0 ; i<4 ; i=i+1)
HM[(dimn-2)][(dimn-1-i)]=D_coefRn2[i]* (1/st^2)*val + HM[(dimn-2)][(dimn-1-i)]
endfor

//==============

// 2. TERMS IN BETWEEN for SECOND DERIVATIVE - 5 terms 
// mid - terms in the matrix, 5 terms centered around diagonal
ista= 2 ;  iend = ((dimsize(r_wave1,0))) -3;
//printf "Middle terms (5 term symmetric terms)  range :%g , %g\r", ista, iend
killwaves /Z  D_coefR3
gen_coef_R3(2) // Generate SECOND  derv. central difference 5 terms
wave D_coefR3 =D_coefR3


for (ista=2 ; ista <= iend; ista = ista+1  )
val=0
RI = r_wave1[ista]
st = r_wave1[(ista+1)] - r_wave1[(ista)]  ;      // STEP SIZE 'h' for the DERIVATVES 5 term
val= -1/(nu*2)
col=ista-2

for (e1=0; e1<5 ; e1=e1+1)
		VA=0
		
		VA =val *D_coefR3[e1]* (1/st^2)
	
		HM[(ista)][(col+e1)] =VA + HM[(ista)][(col+e1)]
	//	printf "--  %g , %g, e1 = %g ; val = %g | Coef : %g \r", ista, (col+e1),  e1,VA, D_coefR3[e1]
	
	endfor
endfor
//==============
killwaves /Z  D_coefR1,D_coefR2,D_coefRn2,D_coefR3,M_inverse,D_coef,A,C


printf "J :  %g | mA :  %6.11f, mB : %6.11f | nu : %6.11f | H-Matrix prepared\r",J,mA,mB,nu
DFREF cdf1 = GetDataFolderDFR()
newdatafolder /O /S Calc
killwaves /z HMat1,HMat, r_wave
Killwaves /Z M_L_eigenVectors, M_R_eigenVectors
DFREF cdf2 = GetDataFolderDFR()
//setdatafolder cdf1
//duplicate /O  H_matrix , HMat
make /O /D /n=(dimsize(r_wave1,0))   r_wave = r_wave1
setdatafolder cdf1
end


//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
function norm_wfn  (s7c):Buttoncontrol // (w1,r1,x1,x2,n)
string s7c

svar w1a=root:Params:RoVibHamiltonian:nwave 
svar r1a =root:Params:RoVibHamiltonian:distance
nvar x1= root:Params:RoVibHamiltonian:llim
nvar x2 = root:Params:RoVibHamiltonian:ulim
nvar n =root:Params:RoVibHamiltonian:int_pnts
nvar threshold=root:Params:RoVibHamiltonian:int_thresh
variable th=1*10^(-1*threshold)

string cdf2=GetDataFolder(1)
string s1,s2
sprintf s1,"%s%s",cdf2,w1a
//print s1
wave w1=$s1
sprintf s2,"%s%s",cdf2,r1a
print s1, s2
wave r1=$s2

variable t0=ticks
string wname=NameOfWave(w1)

wavestats /Q w1
if (V_min < -1e-5)
print "inv"
w1=-1*(w1[p])
endif

make /O /D /n=(dimsize(w1,0)) amp_wave=w1[p]*w1[p]   * ((r1[p])^2 )  ;
wave amp_wave =amp_wave
customsplineU_auto(amp_wave,r1)
wave amp_wave_a=amp_wave_a
glnumerical_integrate(x1,x2,n,amp_wave_a)
nvar res = intResult
variable N1=(1/(sqrt(res))) ; //print N1 , res
printf "normalization factor : %g\r",N1

make /O /D /n=(dimsize(w1,0)) norm_wave1=w1[p]*N1
make /O /D /n=(dimsize(w1,0)) test_wave=norm_wave1[p]*norm_wave1[p] * ((r1[p])^2 )


// testing for correct normalization
customsplineU_auto(test_wave,r1)
wave test_wave_a=test_wave_a
glnumerical_integrate(x1,x2,n,test_wave_a)
nvar res=intResult; 
variable delta =abs((1.0000000000000000- res))
//printf "Delta : %6.14f | time : %g sec \r",(delta),(ticks-t0)/60
if (delta < (th))
	killwaves /Z test_wave_a,test_wave,GLpoints1,amp_wave_a,amp_wave
	duplicate /O /D norm_wave1, $wname+"_norm" ; killwaves /Z norm_wave1
	printf "delta < threshold (%e)| Finished:%g sec\r",(th),(ticks-t0)/60
	else
	printf "delta >(threshold  %e) :delta=%6.14f\r",th,delta
	killwaves /Z test_wave_a,test_wave,GLpoints1,amp_wave_a,amp_wave
endif
end



//*******************************************************************************************************************
//*******************************************************************************************************************


//*******************************************************************************************************************
//*******************************************************************************************************************

//**************************************************************************************************************
function norm_wfn_batch (s1g):buttoncontrol
string s1g // (w1,r1,x1,x2,n)

svar r1a =root:Params:RoVibHamiltonian:b_norm_distance
nvar x1= root:Params:RoVibHamiltonian:llim
nvar x2 = root:Params:RoVibHamiltonian:ulim
nvar n =root:Params:RoVibHamiltonian:int_pnts
nvar threshold=root:Params:RoVibHamiltonian:int_thresh
variable th=1*10^(-1*threshold)
//--------------------------------
killvariables /z intResult
//--------------------------------
variable we1; string nam5
for (we1=0; we1<40; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

//print nam5
string w1a=nam5

string cdf2=GetDataFolder(1)
string s1,s2

wave w1=$nam5			//wavefunction wave
sprintf s2,"%s%s",cdf2,r1a 	// distance wave
printf "------ input: %s,  %s ------\r" nam5, s2
wave r1=$s2

variable t0=ticks
string wname=NameOfWave(w1)

wavestats /Q w1
if (V_min < -1e-5)
print "inv"
w1=-1*(w1[p])
endif

make /O /D /n=(dimsize(w1,0)) amp_wave=w1[p]*w1[p] * (r1[p]^2)  ;
wave amp_wave =amp_wave
customsplineU_auto(amp_wave,r1)
wave amp_wave_a=amp_wave_a
glnumerical_integrate(x1,x2,n,amp_wave_a)
nvar res = intResult
variable N1=1/sqrt(res) ; //print N1 , res
make /O /D /n=(dimsize(w1,0)) norm_wave1=w1[p]*N1
make /O /D /n=(dimsize(w1,0)) test_wave=norm_wave1[p]*norm_wave1[p] * (r1[p]^2)

// testing for correct normalization
customsplineU_auto(test_wave,r1)
wave test_wave_a=test_wave_a
glnumerical_integrate(x1,x2,n,test_wave_a)
nvar res=intResult; 
variable delta =abs((1.0000000000000000 - res))
//printf "Delta : %6.14f | time : %g sec \r",(delta),(ticks-t0)/60
if (delta < (th))
	killwaves /Z test_wave_a,test_wave,GLpoints1,amp_wave_a,amp_wave
	duplicate /O /D norm_wave1, $wname+"_norm" ; killwaves /Z norm_wave
	printf "delta < threshold (%e)| Finished:%g sec\r",(th),(ticks-t0)/60
	else
	printf "delta >(threshold  %e) :delta=%6.14f\r",th,delta
	killwaves /Z test_wave_a,test_wave,GLpoints1,amp_wave_a,amp_wave
endif
printf "N1 factor : %f  ----------------------------\r",N1
endfor

end
//**************************************************************************************************************
// function to extract the wavefunctions for the different J levels within the ZEROTH vibrational state.
function get_wfns_v0(s24):buttoncontrol
string s24
string saveDFR = GetDataFolder(1)

variable i=0, x1,x2

// v=0 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
string s2
sprintf s2,"%sJ%g",saveDFR,i

if (datafolderexists (s2))
string s1,s11,s12,s13
sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
sprintf s11,"%sJ%g:Calc:evec[p][0]"saveDFR,i
sprintf s12,"v0J%g",i
sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
print s11, s12,s13

wave w1=$s1 ; wave w2= $s13
//print nameofwave(w1)

x1=dimsize (w1,0) ; x2= dimsize(w1,1)

//print x1,x2
newdatafolder /O /s wfns
make /d /O /n=(x1) r_wave =w2[p][0]
make /d /O /n=(x1) $s12 =w1[p][0]	// 0th vibrational state
setdatafolder saveDFR
else 
abort
endif
endfor

end
//**************************************************************************************************************
//**************************************************************************************************************
// function to extract the wavefunctions for the different J levels within the first vibrational state.
function get_wfns_v1(s24):buttoncontrol
string s24
string saveDFR = GetDataFolder(1) // present data folder

variable i=0, x1,x2
//------------------------------------------------------------------------
// v=1 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
	string s2
	sprintf s2,"%sJ%g",saveDFR,i

	if (datafolderexists (s2))
		string s1,s11,s12,s13
		sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
		sprintf s11,"%sJ%g:Calc:evec[p][1]"saveDFR,i
		sprintf s12,"v1J%g",i
		sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
		print s11, s12,s13

		wave w1=$s1 ; wave w2= $s13
		//print nameofwave(w1)

		x1=dimsize (w1,0) ; x2= dimsize(w1,1)

		//print x1,x2
		newdatafolder /O /s wfns
		make /d /O /n=(x1) r_wave =w2[p][0]
		make /d /O /n=(x1) $s12 =w1[p][1]	// 1st vibrational state
		setdatafolder saveDFR
	else 
		abort
	endif
endfor
//------------------------------------------------------------------------
end
//**************************************************************************************************************
//**************************************************************************************************************
// function to extract the wavefunctions for the different J levels within the SECOND vibrational state.
function get_wfns_v2(s24):buttoncontrol
string s24
string saveDFR = GetDataFolder(1) // present data folder

variable i=0, x1,x2
//------------------------------------------------------------------------
// v=1 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
	string s2
	sprintf s2,"%sJ%g",saveDFR,i

	if (datafolderexists (s2))
		string s1,s11,s12,s13
		sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
		sprintf s11,"%sJ%g:Calc:evec[p][2]"saveDFR,i
		sprintf s12,"v2J%g",i
		sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
		print s11, s12,s13

		wave w1=$s1 ; wave w2= $s13
		//print nameofwave(w1)

		x1=dimsize (w1,0) ; x2= dimsize(w1,1)

		//print x1,x2
		newdatafolder /O /s wfns
		make /d /O /n=(x1) r_wave =w2[p][0]
		make /d /O /n=(x1) $s12 =w1[p][2]	// 2nd vibrational state
		setdatafolder saveDFR
	else 
		abort
	endif
endfor
//------------------------------------------------------------------------
end
//**************************************************************************************************************
//**************************************************************************************************************
// function to extract the wavefunctions for the different J levels within the THIRD vibrational state.
function get_wfns_v3(s24):buttoncontrol
string s24
string saveDFR = GetDataFolder(1) // present data folder

variable i=0, x1,x2
//------------------------------------------------------------------------
// v=1 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
	string s2
	sprintf s2,"%sJ%g",saveDFR,i

	if (datafolderexists (s2))
		string s1,s11,s12,s13
		sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
		sprintf s11,"%sJ%g:Calc:evec[p][2]"saveDFR,i
		sprintf s12,"v3J%g",i
		sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
		print s11, s12,s13

		wave w1=$s1 ; wave w2= $s13
		//print nameofwave(w1)

		x1=dimsize (w1,0) ; x2= dimsize(w1,1)

		//print x1,x2
		newdatafolder /O /s wfns
		make /d /O /n=(x1) r_wave =w2[p][0]
		make /d /O /n=(x1) $s12 =w1[p][3]	// 3rd vibrational state
		setdatafolder saveDFR
	else 
		abort
	endif
endfor
//------------------------------------------------------------------------
end
//**************************************************************************************************************
//**************************************************************************************************************
// function to extract the wavefunctions for the different J levels within the FOURTH vibrational state.
function get_wfns_v4(s24):buttoncontrol
string s24
string saveDFR = GetDataFolder(1) // present data folder

variable i=0, x1,x2
//------------------------------------------------------------------------
// v=1 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
	string s2
	sprintf s2,"%sJ%g",saveDFR,i

	if (datafolderexists (s2))
		string s1,s11,s12,s13
		sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
		sprintf s11,"%sJ%g:Calc:evec[p][2]"saveDFR,i
		sprintf s12,"v4J%g",i
		sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
		print s11, s12,s13

		wave w1=$s1 ; wave w2= $s13
		//print nameofwave(w1)

		x1=dimsize (w1,0) ; x2= dimsize(w1,1)

		//print x1,x2
		newdatafolder /O /s wfns
		make /d /O /n=(x1) r_wave =w2[p][0]
		make /d /O /n=(x1) $s12 =w1[p][4]	// 4th vibrational state
		setdatafolder saveDFR
	else 
		abort
	endif
endfor
//------------------------------------------------------------------------
end
//**************************************************************************************************************
//**************************************************************************************************************
// function to extract the wavefunctions for the different J levels within the FIFTH vibrational state.
function get_wfns_v5(s24):buttoncontrol
string s24
string saveDFR = GetDataFolder(1) // present data folder

variable i=0, x1,x2
//------------------------------------------------------------------------
// v=1 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
	string s2
	sprintf s2,"%sJ%g",saveDFR,i

	if (datafolderexists (s2))
		string s1,s11,s12,s13
		sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
		sprintf s11,"%sJ%g:Calc:evec[p][2]"saveDFR,i
		sprintf s12,"v5J%g",i
		sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
		print s11, s12,s13

		wave w1=$s1 ; wave w2= $s13
		//print nameofwave(w1)

		x1=dimsize (w1,0) ; x2= dimsize(w1,1)

		//print x1,x2
		newdatafolder /O /s wfns
		make /d /O /n=(x1) r_wave =w2[p][0]
		make /d /O /n=(x1) $s12 =w1[p][5]	// 5th vibrational state
		setdatafolder saveDFR
	else 
		abort
	endif
endfor
//------------------------------------------------------------------------
end
//**************************************************************************************************************

function get_wfns_v_level(vib_level)
variable vib_level
string s24
string saveDFR = GetDataFolder(1)

variable i=0, x1,x2
//------------------------------------------------------------------------
// v=1 level, and iterate over all available J levels

for (i=0 ; i< 20 ; i=i+1)
string s2
sprintf s2,"%sJ%g",saveDFR,i

sprintf s2,"%sJ%g",saveDFR,i

if (datafolderexists (s2))
string s1,s11,s12,s13
sprintf s1,"%sJ%g:Calc:evec"saveDFR,i
sprintf s11,"%sJ%g:Calc:evec[p][%g]"saveDFR,i, vib_level
sprintf s12,"v%gJ%g",vib_level,i
sprintf s13,"%sJ%g:Calc:r_wave"saveDFR,i
print s11, s12,s13

wave w1=$s1 ; wave w2= $s13
//print nameofwave(w1)

x1=dimsize (w1,0) ; x2= dimsize(w1,1)

//print x1,x2
newdatafolder /O /s wfns
make /d /O /n=(x1) r_wave =w2[p][0]
make /d /O /n=(x1) $s12 =w1[p][vib_level]
setdatafolder saveDFR
else 
abort
endif
endfor
//------------------------------------------------------------------------


end
//**************************************************************************************************************
function batch_expVal_gamma(s25):buttoncontrol
string s25
svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet

variable i1=0, i2=0

i2=((jf-ji)+1)
	//print i2
string s6
sprintf s6,"%s_ResultSet",g_wave
print s6
make /d /O /n=(i2,5) $s6=0
wave gjs=$s6

	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",i1
	sprintf n4,"v0J%g_norm",(i1+2)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	printf  "*****  %s <-> %s (%s, %s)  ***** \r",wav1,wav2,wav4,wav3
	//------------------------------------------------------------------------
	
	
//**CLEAR OLD WAVES *****************
killwaves /Z int_x,scaled_y

//*********************************************
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param : waves input

string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

s4=nameofwave (param)			//parameter wave


variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq


//print d_wf,d_p
//RANGES------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:Bint_min
nvar b= root:Params:RoVibHamiltonian:Bint_max
nvar h= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a,b,h
variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------

//1. SPLINES
customsplineu_auto(w1a,wf_x)
string array1=nameofwave(w1a)+"_a";
//print array1 ; 
s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

customsplineu_auto(w2a,wf_x)
array1=nameofwave(w2a)+"_a"; 	s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2


customsplineu_auto(param,param_x)
array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3


customsplineu_auto(rsq,wf_x)		//r_square term
array1=nameofwave(rsq)+"_a";	sq=nameofwave(rsq)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $sq
wave r_sq=$sq


//remove USED waves  for the three splines
killwaves /Z $array1,scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
//print s5
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] * (r_sq[p])
wave integral=$s5

killwaves /Z $s2,$s3

//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
//print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1

//------------------------------------------------------------------------
nvar res= intResult

//RESULT ASSIGNMENT
gjs[i1][0]=i1
gjs[i1][1]=(i1+2)
gjs[i1][2]=(res)
gjs[i1][3]=(res * (1.4818e-25))
gjs[i1][4]=((res * (1.4818e-25))^2)

//print /D res*1.4818e-25
//variable /G gamma_cm_sq=(res*1.4818e-25)^2
//variable /G gamma_au=(res)
//print /D   gamma_cm_sq
printf "---- SOLVED FOR J=%g <-> J=%g ----\r",i1,(i1+2)
killwaves /Z GLpoints1	,int_x,scaled_y
	//------------------------------------------------------------------------
	nvar res= intResult
//	printf "res(a.u.) = %1.5g | res(cm^-3) = %1.5g \r"res, (res*1.4818e-25)
//	printf "res(cm^-6) = %1.5g\r",(res*1.4818e-25)^2
//----------------------------------------------------------------------------
killwaves /Z $array1, integral, $s3
string kw=nameofwave(w1a)+"_sc_intx"
string kw2=nameofwave(w2a)+"_sc_intx"
string kw3=s3
printf "%s | %s | %s\r", kw,kw2,kw3
killwaves /z $kw,$kw2,$kw3
//----------------------------------------------------------------------------
endfor

// clean the used waves ----------------------------
killwaves /z int_x, $array1 ;  killwaves /z GLpoints1
//-----------------------------------------------------------
print "Done !"
end


//*******************************************************************************************************************
//*******************************************************************************************************************
//**************************************************************************************************************
function batch_expVal_gamma_opt(s25):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// spline wave stored for run // for each wave //
string s25
svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet

variable i1=0, i2=0

i2=((jf-ji)+1)
	//print i2
string s6
sprintf s6,"%s_ResultSet",g_wave
print s6
make /d /O /n=(i2,5) $s6=0
wave gjs=$s6

variable g0,g1,g2=0 // special variable for running for matrix element for <00|gamma|00>
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
	g0=i1
	if(g0==0) 
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",(g0)
	sprintf n4,"v0J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	printf  "*****  %s <-> %s (%s, %s)  ***** \r",wav1,wav2,wav4,wav3
	//------------------------------------------------------------------------
	
	else
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",(i1)
	sprintf n4,"v0J%g_norm",(i1+2)
//	print i1, n3,n4
	 wav1=n3		;  wav2=n4
	 wav3=g_wave	;  wav4 = g_dis_wave
	 wav5=wfn_dis
	printf  "*****  %s <-> %s (%s, %s)  ***** \r",wav1,wav2,wav4,wav3
	//------------------------------------------------------------------------
	
	endif
	
//**CLEAR OLD WAVES *****************
killwaves /Z int_x,scaled_y

//*********************************************
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param : waves input

string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

s4=nameofwave (param)			//parameter wave


variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq


//print d_wf,d_p
//RANGES------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:Bint_min
nvar b= root:Params:RoVibHamiltonian:Bint_max
nvar h= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a,b,h
variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------

//1. SPLINES
string array1=nameofwave(w1a)+"_a";
variable w1
w1= exists (array1)
//print w1
if (w1==1)
print "array wave found"
else
 print "Array wave not found."
 customsplineu_auto(w1a,wf_x)
endif 

s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

//---------------------------------------------------
array1=nameofwave(w2a)+"_a";
variable w2
w2= exists (array1)
if (w2==1)
print "array wave found"
else
 print "Array wave not found."
customsplineu_auto(w2a,wf_x)
endif 
s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2
//---------------------------------------------------
array1=nameofwave(rsq)+"_a";
variable w3
w3= exists (array1)
if (w3==1)
print "array wave found"
else
 print "Array wave not found."
customsplineu_auto(rsq,wf_x)
endif 
sq=nameofwave(rsq)+"_sc_intx"
ScaledY_gen(int_x,$array1)
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $sq
wave r_sq=$sq
//--------------------------------------------------
	customsplineu_auto(rsq,wf_x)		//r_square term
	array1=nameofwave(rsq)+"_a";	sq=nameofwave(rsq)+"_sc_intx"
	ScaledY_gen(int_x,$array1)
	killwaves /Z $array1
	duplicate /O /D scaled_y, $sq
	wave r_sq=$sq
//---------------------------------------------------

customsplineu_auto(param,param_x)
array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3





//remove USED waves  for the three splines
killwaves /Z $array1,scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
//print s5
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] * (r_sq[p])
wave integral=$s5

killwaves /Z $s2,$s3

//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
//print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1

//------------------------------------------------------------------------
nvar res= intResult

//RESULT ASSIGNMENT
gjs[i1][0]=i1
gjs[i1][1]=(i1+2)
gjs[i1][2]=(res)
gjs[i1][3]=(res * (1.4818e-25))
gjs[i1][4]=((res * (1.4818e-25))^2)

//print /D res*1.4818e-25
//variable /G gamma_cm_sq=(res*1.4818e-25)^2
//variable /G gamma_au=(res)
//print /D   gamma_cm_sq
printf "---- SOLVED FOR J=%g <-> J=%g ----\r",i1,(i1+2)
killwaves /Z GLpoints1	,int_x,scaled_y
	//------------------------------------------------------------------------
	nvar res= intResult
//	printf "res(a.u.) = %1.5g | res(cm^-3) = %1.5g \r"res, (res*1.4818e-25)
//	printf "res(cm^-6) = %1.5g\r",(res*1.4818e-25)^2
//----------------------------------------------------------------------------
killwaves /Z $array1, integral, $s3
string kw=nameofwave(w1a)+"_sc_intx"
string kw2=nameofwave(w2a)+"_sc_intx"
string kw3=s3
printf "%s | %s | %s\r", kw,kw2,kw3
killwaves /z $kw,$kw2,$kw3
//----------------------------------------------------------------------------
endfor

// clean the used waves ----------------------------
killwaves /z int_x, $array1 ;  killwaves /z GLpoints1
//-----------------------------------------------------------
print "Done !"
end

//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
function cleanCustomSpline(s2h):buttoncontrol
string s2h
string cdf=getdatafolder(1)
string folder="root:Params:Customspline:"
setdatafolder folder
killwaves /a /z
setdatafolder cdf

end

//*******************************************************************************************************************
//*******************************************************************************************************************

function clean_H_Matrix(s2i):buttoncontrol
string s2i
string cdf=getdatafolder(1)
string fname,j
nvar jm=root:Params:RoVibHamiltonian:pmax
//print jm
variable i=0
for (i=0 ; i<(jm+1) ; i=i+1)
sprintf fname,"J%g",i
//print fname
//print cdf
string ndf
sprintf ndf,"%s%s",cdf,fname
//print ndf

setdatafolder ndf
killwaves /Z  H_matrix

sprintf ndf,"%s:Calc:",ndf
setdatafolder ndf
killwaves /Z M_R_eigenVectors
//print ndf

endfor


setdatafolder cdf
print "H-matrices removed !"
end


//*******************************************************************************************************************
//*******************************************************************************************************************
function print_eval_J0(s2k):buttoncontrol
string s2k
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Calculated Eigenvalues ---"
string cdf=getdatafolder(1)
variable i=0
for (i=0 ; i<(15) ; i=i+1)

sprintf fname,"J0",i
string ndf
sprintf ndf,"%s%s:Calc",cdf,fname
//print ndf
setdatafolder ndf
wave eval=eval

nvar asymp_corr = root:Params:RoVibHamiltonian:AsympCorr ; nvar valp=root:Params:RoVibHamiltonian:hartree_asym
variable corr=asymp_corr ; variable p=valp
 //119.5320000  // CORRECTION OFFSET FOR INFINITY (119.532+0.542-2.9218)

printf "%6.7f\t%6.8f\t %6.5f\r",  (( eval [(i)] * 219474.631370500000 )-(((p)*219474.63137050000)+ (corr))),  eval [(i)], (eval [(i)] * 219474.631370500000)

endfor
setdatafolder cdf
printf "Asymptote of potential: %7.7f | Correction : %7.7f ",p,corr
print "----------------------------"
end
//*******************************************************************************************************************
//*******************************************************************************************************************

//********************************************************************************

// WORKS FROM THE PANEL :::
function JfoldersPrep_OtherAtoms (s7b):ButtonControl
        string s7b

        nvar x2= root:Params:RoVibHamiltonian:pmax; 
        nvar mA= root:Params:RoVibHamiltonian:massA ; 
        nvar mB= root:Params:RoVibHamiltonian:massB ;
        nvar neA = root:Params:RoVibHamiltonian:n_electronA
        nvar neB=root:Params:RoVibHamiltonian:n_electronB ; 

        svar  PES_Ewave =root:Params:RoVibHamiltonian:potential_energywave
        svar PES_rwave = root:Params:RoVibHamiltonian:potential_rwave

        variable uA,uB
        uA=mA+neA ; 
        uB=mB+neB ; 

        // printf "%s ( %g , %g)\r" mol,mA,mB
        string cdf=getdatafolder(1)

        wave r_wave1 =$PES_rwave 
        wave scaled_y1 = $PES_Ewave
        variable dim= dimsize(r_wave1,0)
        printf "Dim: %g x %g , %g , %g \r",dim,dim, uA,uB
        variable i; string w

        for (i=0 ; i <= (x2) ; i=i+1)
	        sprintf w,"J%g",i
	        newdatafolder /S $w
	        make /O /D /n=(dim) r_wave0 =r_wave1 
	        make /O /D /n=(dim) scaled_y0 =scaled_y1
	        string r1="r_wave0"   ;  string p1= "scaled_y0"

	        MatrixPrep_auto(i,  uA,uB ,r1, p1)
	        setdatafolder $cdf
        endfor
        
        // H2 =   1836.152673890000000 ,1836.152673890000000
        //D2=   3670.48296785000000   , 3670.48296785000000
        //HD =   1836.152673890000000 ,  3670.48296785000000
end

//********************************************************************************
Window Prep_HMatrixOther() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(539,54,802,500)
	SetDrawLayer UserBack
	DrawLine 2,244,257,244
	DrawLine 4,293,259,293
	DrawLine 4,343,259,343
	DrawLine 2,414,257,414
	TitleBox Prepare_H_matrix,pos={1.00,4.00},size={108.00,23.00},title="Prepare_H_matrix"
	TitleBox Prepare_H_matrix,fStyle=1
	TitleBox specificMasses,pos={9.00,33.00},size={143.00,23.00},title="Entering_specific_ masses"
	SetVariable NuclearMass_A,pos={1.00,64.00},size={259.00,18.00},title="NuclearMass_A"
	SetVariable NuclearMass_A,format="%7.12f",valueColor=(65535,0,0)
	SetVariable NuclearMass_A,limits={0,inf,0},value= root:Params:RoVibHamiltonian:massA
	SetVariable NuclearMass_B,pos={1.00,85.00},size={259.00,18.00},title="NuclearMass_B"
	SetVariable NuclearMass_B,format="%7.12f"
	SetVariable NuclearMass_B,limits={-inf,inf,0},value= root:Params:RoVibHamiltonian:massB
	SetVariable n_electron,pos={9.00,108.00},size={222.00,18.00},title="NumberOf_electron(Atom-A)"
	SetVariable n_electron,limits={1,30,1},value= root:Params:RoVibHamiltonian:n_electronA
	SetVariable j_max,pos={135.00,150.00},size={105.00,18.00},title="J_max"
	SetVariable j_max,value= root:Params:RoVibHamiltonian:pmax
	Button j_folders_prep_heavy,pos={39.00,219.00},size={165.00,19.00},proc=JfoldersPrep_OtherAtoms,title="Make_H_(J_folders)"
	Button j_folders_prep_heavy,fColor=(65535,32764,16385)
	SetVariable n_electronB,pos={9.00,129.00},size={222.00,18.00},title="NumberOf_electron(Atom-B)"
	SetVariable n_electronB,limits={1,30,1},value= root:Params:RoVibHamiltonian:n_electronB
	PopupMenu r_wave_select,pos={40.00,174.00},size={167.00,19.00},proc=pes_rwave__Selection,title="Potential_r_wave"
	PopupMenu r_wave_select,mode=3,popvalue="r_wave_au",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu PES_wave_select1,pos={9.00,195.00},size={228.00,19.00},proc=pes_energy__Selection,title="Potential_Energy_wave"
	PopupMenu PES_wave_select1,mode=2,popvalue="O2_PES_Hartree",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	Button Rot_Transition_energies,pos={7.00,270.00},size={105.00,19.00},proc=print_transition_rotn_even,title="RotnTransm even"
	Button Rot_Transition_energies,fColor=(32768,40777,65535)
	Button VibrationalTrans,pos={129.00,270.00},size={129.00,19.00},proc=print_Vib_transitionEnergy,title="VibrationalTransitions"
	Button VibrationalTrans,fColor=(1,52428,26586)
	TitleBox title01,pos={115.00,246.00},size={144.00,23.00},title="Obtain transition energies"
	TitleBox Normalization,pos={3.00,295.00},size={82.00,23.00},title="Norm_Int(x)dr"
	Button button0,pos={90.00,297.00},size={105.00,19.00},proc=norm_wfn_wtr2,title="Wfn_normalize"
	Button button1,pos={93.00,319.00},size={150.00,19.00},proc=expectation_val_wtr2,title="Calc_expectationVal"
	Button O1,pos={0.00,345.00},size={82.00,19.00},proc=print_TE_O1,title="O1_transition"
	Button Q1,pos={87.00,343.00},size={82.00,19.00},proc=print_TE_Q1,title="Q1_transition"
	Button S1,pos={175.00,342.00},size={82.00,19.00},proc=print_TE_S1,title="S1_transition"
	Button O1_q,pos={0.00,367.00},size={39.00,19.00},proc=print_TE_O1_q,title="O1_q"
	Button Q1_q,pos={43.00,367.00},size={39.00,19.00},proc=print_TE_Q1_q,title="Q1_q"
	Button S1_q,pos={87.00,367.00},size={39.00,19.00},proc=print_TE_S1_q,title="S1_q"
	Button del_G,pos={132.00,367.00},size={39.00,19.00},proc=print_TE_delG_q,title="del_G"
	Button R1_,pos={174.00,366.00},size={28.00,19.00},proc=print_TE_R1_q,title="R1"
	Button P1,pos={211.00,366.00},size={27.00,19.00},proc=print_TE_P1_q,title="P1"
	Button Rotn_Energy,pos={3.00,394.00},size={120.00,19.00},proc=print_TE_EJ,title="Rotn_Energy(cm-1)"
	Button Rot_Transition_energies1,pos={4.00,248.00},size={98.00,19.00},proc=print_transitionEnergy,title="RotnTransitions"
	Button Rot_Transition_energies1,fColor=(32768,40777,65535)
	Button O1_q1,pos={5.00,416.00},size={67.00,19.00},proc=print_TE_O1_q_alt,title="O1_q_alt"
	Button O1_q1,fColor=(16385,28398,65535)
	Button S1_q1,pos={88.00,417.00},size={63.00,19.00},proc=print_TE_S1_q_alt,title="S1_q_alt"
	Button S1_q1,fColor=(3,52428,1)
EndMacro


//********************************************************************************
//*******************************************************************************************************************
function print_transitionEnergy(s2m):buttoncontrol
string s2m
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Rotational transition energies ---"
string cdf=getdatafolder(1)
print cdf
variable i=0
for (i=0 ; i<(jm-1) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,i
sprintf addr2,"%sJ%g:Calc",cdf,(i+2)
//print addr1, addr2 
//a=:Calc:eval[0]
setdatafolder addr1
wave eval=eval
a=eval[0]

setdatafolder addr2
wave eval=eval
b=eval[0]


printf "%6.8f\r",( (b-a )*219474.6313702 )

endfor
setdatafolder cdf

print "----------------------------------------------"
end

//********************************************************************************
//*******************************************************************************************************************
function print_transition_rotn_even (s2m):buttoncontrol
string s2m
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Rotational transition energies ---"
string cdf=getdatafolder(1)
print cdf
variable i=0
for (i=1 ; i<(jm-1) ; i=i+2)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,i
sprintf addr2,"%sJ%g:Calc",cdf,(i+2)
//print addr1, addr2 
//a=:Calc:eval[0]
setdatafolder addr1
wave eval=eval
a=eval[0]

setdatafolder addr2
wave eval=eval
b=eval[0]


printf "%g -> %g, %6.8f\r",  i  , (i+2) ,( (b-a )*219474.6313702 )

endfor
setdatafolder cdf

print "----------------------------------------------"
end

//*******************************************************************************************************************
//*******************************************************************************************************************
function print_Vib_transitionEnergy(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Vibrational transition energies ---"
string cdf=getdatafolder(1)
variable i=0
for (i=1 ; i<(6) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ0:Calc",cdf,i

//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(0)]


wave eval=eval
b=eval[(i)]


printf "v=%g -> v=%g = %6.8f\r",0,(i), ( (b-a )*219474.6310000)

endfor
setdatafolder cdf

print "----------------------------------------------"
end
//*******************************************************************************************************************
//*******************************************************************************************************************
function expectation_val_wtr2 (s72d):buttoncontrol
string s72d
Print "Expectation value (wt r^2) for selected parameter:----------"
svar wav1 = root:Params:RoVibHamiltonian:iwv1
svar wav2 =root:Params:RoVibHamiltonian:iwv2
svar wav3 =root:Params:RoVibHamiltonian:iwv3
svar wav5 = root:Params:RoVibHamiltonian:int_distance //wavefunction distance
svar wav4 =root:Params:RoVibHamiltonian:iwv4   // parameter_wave_distance
nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
//**CLEAR OLD WAVES *****************
killwaves /Z int_x,scaled_y

//*********************************************
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param : waves input



string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /o /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq
s4=nameofwave (param)			//parameter wave

print d_wf,d_p
//RANGES-------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:int_min
nvar b= root:Params:RoVibHamiltonian:int_max
nvar h= root:Params:RoVibHamiltonian:step_h

variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------

//1. SPLINES
customsplineu_auto(w1a,wf_x)
string array1=nameofwave(w1a)+"_a";
//print array1 ; 
s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

customsplineu_auto(w2a,wf_x)
array1=nameofwave(w2a)+"_a"; 	s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2


//customsplineu_auto(param,param_x)
//array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
//ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
//duplicate /O /D scaled_y, $s3
//wave param_sc=$s3







customsplineu_auto(param,param_x)
array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3



//remove USED waves  for the three splines
killwaves /Z $array1,scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
print s5
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] 
wave integral=$s5
//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1
nvar res= intResult
printf "Gamma (au): %g | (cm^3): %g \r ",res, res*1.4818e-25
variable /G gamma_cm_sq=(res*1.48185e-25)^2
variable /G gamma_au=(res)
printf "Gamma_sq (cm^6): %g\r",gamma_cm_sq

// clean the used waves ----------------------------
killwaves /z int_x, $array1, rsq, $s4+"_expVal" , $s4+"_sc_intx", $sq+"_sc_intx",$s2+"_sc_intx",$s3+"_sc_intx" , $s1+"_sc_intx"  ;
 killwaves /z GLpoints1
//-----------------------------------------------------------
end


//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
function norm_wfn_wtr2  (s72c):Buttoncontrol // (w1,r1,x1,x2,n)
string s72c
Print "Normalize wavefuntion (wt r^2): ----------------------------"
svar w1a=root:Params:RoVibHamiltonian:nwave 
svar r1a =root:Params:RoVibHamiltonian:distance
nvar x1= root:Params:RoVibHamiltonian:llim
nvar x2 = root:Params:RoVibHamiltonian:ulim
nvar n =root:Params:RoVibHamiltonian:int_pnts
nvar threshold=root:Params:RoVibHamiltonian:int_thresh
variable th=1*10^(-1*threshold)

string cdf2=GetDataFolder(1)
string s1,s2
sprintf s1,"%s%s",cdf2,w1a
//print s1
wave w1=$s1
sprintf s2,"%s%s",cdf2,r1a
print s1, s2
wave r1=$s2

variable t0=ticks
string wname=NameOfWave(w1)

wavestats /Q w1
if (V_min < -1e-5)
print "inv"
w1=-1*(w1[p])
endif

make /O /D /n=(dimsize(w1,0)) amp_wave=w1[p]*w1[p]   ;
wave amp_wave =amp_wave
customsplineU_auto(amp_wave,r1)
wave amp_wave_a=amp_wave_a
glnumerical_integrate(x1,x2,n,amp_wave_a)
nvar res = intResult
variable N1=(1/(sqrt(res))) ; //print N1 , res
printf "normalization factor : %g\r",N1

make /O /D /n=(dimsize(w1,0)) norm_wave1=w1[p]*N1
make /O /D /n=(dimsize(w1,0)) test_wave=norm_wave1[p]*norm_wave1[p] ;


// testing for correct normalization
customsplineU_auto(test_wave,r1)
wave test_wave_a=test_wave_a
glnumerical_integrate(x1,x2,n,test_wave_a)
nvar res=intResult; 
variable delta =abs((1.0000000000000000- res))
//printf "Delta : %6.14f | time : %g sec \r",(delta),(ticks-t0)/60
if (delta < (th))
	killwaves /Z test_wave_a,test_wave,GLpoints1,amp_wave_a,amp_wave
	duplicate /O /D norm_wave1, $wname+"_norm" ; killwaves /Z norm_wave1
	printf "delta < threshold (%e)| Finished:%g sec\r",(th),(ticks-t0)/60
	else
	printf "delta >(threshold  %e) :delta=%6.14f\r",th,delta
	killwaves /Z test_wave_a,test_wave,GLpoints1,amp_wave_a,amp_wave
endif
end



//*******************************************************************************************************************

function MaxValNormalize (s11a):buttoncontrol                     //
string s11a
string cdf=getdatafolder(1)
print cdf
wave wl
string option,trail_text
trail_text="_MxVal_norm"
prompt option,"Selec the wave to be normalized (to max y-value of 1) ", popup WaveList("*", ";", "")
prompt trail_text,"Trailing text for the normalized spectra"
DoPrompt "Normalizing white light spectra", option,trail_text
 				if( V_Flag )
     				return 0         		// user canceled
   				endif
wave wl=$option   				

variable d1
d1=dimsize(wl,0)
print d1
wavestats /Q wl
print V_max

string currentdf 
string part1=nameofwave(wl)
string wname ; sprintf wname,"%s%s",part1,trail_text
currentdf= GetDataFolder(1)		//Save the current data folder 
//printf "Current Data folder is : %s\r",currentdf

//NewDataFolder /O root:D_Sub			//Data folder where the normalized wl wave will be made
variable maxvalue=V_max
print maxvalue
	//	SetDataFolder root:D_Sub						//Current data folder is changed to root:D_Sub
		make /o /d  /n=(d1) $wname=wl / maxvalue			//New wave is made there 
		
		display $part1+trail_text
		TextBox/C/N=text0/A=LT "Normalized wave  \r  " + NameOfWave(wl);
		SetDataFolder currentdf						// Restore current data folder
		make /o /n=(d1) $part1+trail_text=wl / maxvalue			//Also makes wl norm in current folder
KillVariables /Z maxvalue
end


//***********************************************************************************************************
//O1 branch
function print_TE_O1(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- O1 branch transition energies ---"
string cdf=getdatafolder(1) ; // print cdf ;
variable i=0
for (i=0 ; (i+2)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i)
sprintf addr2,"%sJ%g:Calc",cdf,(i+2)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "O1-branch"
printf "v=1,J=%g <-> v=0,J=%g  --->   %6.4f\r",(i),(i+2), ( (a-b)  *219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   O1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//Q1 branch
function print_TE_Q1(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Q1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+2)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "Q1-branch"
printf "v=1,J=%g <-> v=0,J=%g  --->  %6.8f\r",(i),(i), ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   Q1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//S1 branch
function print_TE_S1(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- S1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+2)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i+2)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "S1-branch"
printf "v=1,J=%g <-> v=0,J=%g  --->   %6.8f\r",(i+2),(i), ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   S1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//O1 branch
function print_TE_O1_q(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- O1 branch transition energies ---"
string cdf=getdatafolder(1) ; // print cdf ;
variable i=0
for (i=0 ; (i+2)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i)
sprintf addr2,"%sJ%g:Calc",cdf,(i+2)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "O1-branch"
printf "%g, %6.4f\r", i, ( (a-b)  *219474.6313702000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   O1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//O1 branch - even


//***********************************************************************************************************
//***********************************************************************************************************
//Q1 branch
function print_TE_Q1_q(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Q1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "Q1-branch"
printf "%6.4f\r", ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   Q1 branch
end



//***********************************************************************************************************
//***********************************************************************************************************
//S1 branch
function print_TE_S1_q(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- S1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+2)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i+2)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "S1-branch"
printf "%6.4f\r", ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   S1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//S1 branch - even
function print_TE_S1_q_alt (s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- S1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=1 ; (i+2)<(jm) ; i=i+2)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i+2)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "S1-branch"
printf "%6.4f\r", ( (a-b)*219474.631370200000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   S1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//S(1) branch relative to the v1 i.e. q1 
function print_TE_S1_relv1(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- S1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+2)<(jm) ; i=i+1)

variable a,b,c, ge
string addr1,addr2, addr3
sprintf addr1,"%sJ%g:Calc",cdf,(i+2)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
sprintf addr3,"%sJ0:Calc",cdf
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

setdatafolder addr3
wave eval=eval
c=eval[(1)]
ge=eval[(0)]

//print "S1-branch"
printf "%6.5f, %6.5f\r", ( (a-b)*219474.6310000000000000000000) , ( (a-b)*219474.6310000000000000000000) - ((c-ge)*219474.631)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   S1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************

//Rotational energies with respect to the ground state

function print_TE_EJ(s2o):buttoncontrol
string s2o
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- E(j)  rotational energies (from j=1) with respect to the ground state (j=0) ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+1)<(jm) ; i=i+1)
variable k=0
variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i+1)
sprintf addr2,"%sJ%g:Calc",cdf,(k)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(0)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "S1-branch"
printf "%6.10f\r", ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   S1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//del_G_difference between vibrational_levels_

//***********************************************************************************************************
//***********************************************************************************************************
//R1 branch - from the first vibrational level
function print_TE_R1_q(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- R1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+1)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i+1)
sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(2)]

setdatafolder addr2
wave eval=eval
b=eval[(1)]

//print "R1-branch"
printf "%6.8f\r", ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   R1 branch
end

//***********************************************************************************************************
//***********************************************************************************************************
//P1 branch - from the first vibrational level
function print_TE_P1_q(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- P1 branch transition energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+1)<(jm) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i)
sprintf addr2,"%sJ%g:Calc",cdf,(i+1)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "P1-branch"
printf "%6.8f\r", ( (a-b)*219474.6310000000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   P1 branch
end

//***********************************************************************************************************
//*******************************************************************************************************************
//**************************************************************************************************************
function batch_expVal_gamma_2(wave1,wave2_xAxis)  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// spline wave stored for run // for each wave //
wave wave1, wave2_xAxis

string g_wave=nameofwave(wave1)
string g_dis_wave=nameofwave(wave2_xAxis)


//svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
//svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet

variable i1=0, i2=0

i2=((jf-ji)+1)
	//print i2
string s6
sprintf s6,"%s_ResultSet",g_wave
print s6
make /d /O /n=(i2,5) $s6=0
wave gjs=$s6

	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",i1
	sprintf n4,"v0J%g_norm",(i1+2)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	printf  "*****  %s <-> %s (%s, %s)  ***** \r",wav1,wav2,wav4,wav3
	//------------------------------------------------------------------------
	
	
//**CLEAR OLD WAVES *****************
killwaves /Z int_x,scaled_y

//*********************************************
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param : waves input

string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

s4=nameofwave (param)			//parameter wave


variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq


//print d_wf,d_p
//RANGES------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:Bint_min
nvar b= root:Params:RoVibHamiltonian:Bint_max
nvar h= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a,b,h
variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------

//1. SPLINES
string array1=nameofwave(w1a)+"_a";
variable w1
w1= exists (array1)
//print w1
if (w1==1)
print "array wave found"
else
 print "Array wave not found."
 customsplineu_auto(w1a,wf_x)
endif 

s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

//---------------------------------------------------
array1=nameofwave(w2a)+"_a";
variable w2
w2= exists (array1)
if (w2==1)
print "array wave found"
else
 print "Array wave not found."
customsplineu_auto(w2a,wf_x)
endif 
s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2
//---------------------------------------------------
array1=nameofwave(rsq)+"_a";
variable w3
w3= exists (array1)
if (w3==1)
print "array wave found"
else
 print "Array wave not found."
customsplineu_auto(rsq,wf_x)
endif 
sq=nameofwave(rsq)+"_sc_intx"
ScaledY_gen(int_x,$array1)
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $sq
wave r_sq=$sq
//--------------------------------------------------
// 	customsplineu_auto(rsq,wf_x)		//r_square term
//	array1=nameofwave(rsq)+"_a";	sq=nameofwave(rsq)+"_sc_intx"
//	ScaledY_gen(int_x,$array1)
//	killwaves /Z $array1
//	duplicate /O /D scaled_y, $sq
//	wave r_sq=$sq
//---------------------------------------------------

customsplineu_auto(param,param_x)
array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3





//remove USED waves  for the three splines
killwaves /Z $array1,scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
//print s5
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] * (r_sq[p])
wave integral=$s5

killwaves /Z $s2,$s3

//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
//print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1

//------------------------------------------------------------------------
nvar res= intResult

//RESULT ASSIGNMENT
gjs[i1][0]=i1
gjs[i1][1]=(i1+2)
gjs[i1][2]=(res)
gjs[i1][3]=(res * (1.4818e-25))
gjs[i1][4]=((res * (1.4818e-25))^2)

//print /D res*1.4818e-25
//variable /G gamma_cm_sq=(res*1.4818e-25)^2
//variable /G gamma_au=(res)
//print /D   gamma_cm_sq
printf "---- SOLVED FOR J=%g <-> J=%g ----\r",i1,(i1+2)
killwaves /Z GLpoints1	,int_x,scaled_y
	//------------------------------------------------------------------------
	nvar res= intResult
//	printf "res(a.u.) = %1.5g | res(cm^-3) = %1.5g \r"res, (res*1.4818e-25)
//	printf "res(cm^-6) = %1.5g\r",(res*1.4818e-25)^2
//----------------------------------------------------------------------------
killwaves /Z $array1, integral, $s3
string kw=nameofwave(w1a)+"_sc_intx"
string kw2=nameofwave(w2a)+"_sc_intx"
string kw3=s3
printf "%s | %s | %s\r", kw,kw2,kw3
killwaves /z $kw,$kw2,$kw3
//----------------------------------------------------------------------------
endfor

// clean the used waves ----------------------------
killwaves /z int_x, $array1 ;  killwaves /z GLpoints1
//-----------------------------------------------------------
print "Done !"
end

//************************************************************************************************************************
//************************************************************************************************************************
//Select the parameter wave , gamma or mean polarizability and run it.

function run_batch_gamma_fullSet()
variable we1; string nam5
for (we1=0; we1<30; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

//print nam5
string w1a=nam5
print w1a
string distance_wave="root:HD_PES:wfns:g_distance" 
batch_expVal_gamma_2($w1a,  $distance_wave)

endfor
end
//************************************************************************************************************************
//************************************************************************************************************************
// processing the 'ResultSet' for wavelength dependent polarizability (anisotropy, gamma).
function table_process_anisotropy(st01):buttoncontrol
string st01
string cdf, ocdf; 
ocdf=getdatafolder(1)
newdatafolder /o /s  processed
cdf=getdatafolder(1)
// print cdf
variable we1; string nam5
for (we1=0; we1<75; we1=we1+1)

	nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
		break
	endif

	//print nam5
	string w1a=nam5
	// print w1a
	string n_wave
	sprintf  n_wave,"%s_t",w1a 
	//print n_wave
	wave data_matrix = $w1a
	matrixop $n_wave = (abs(data_matrix))^t

	string name; name=nameofwave ($n_wave) ; print name ; 
	string fwave; sprintf fwave,"%s%s",cdf,name
	print fwave
	movewave $n_wave, $fwave

	wave data_matrix = $fwave
	DeletePoints 0,2, data_matrix
	DeletePoints 1,2, data_matrix
endfor

concat_gamma_matrix()		//	concatenate the waves to a 2D matrix

// wave definitions ----------
string outputwName="Gamma_all_DataSet"
wave output = $outputwName
wave omega_originalWave=root:data:omega_nm
wave omega_finalWave = root:data:omega_final
//----------------------------------

// copy the output to a new folder 'scaling'
newdatafolder /O /S scaling

duplicate /d /o output, Gamma_all_DataSet
duplicate /d /o omega_originalWave , Omega_Original
duplicate /d /o omega_finalWave , Omega_Final
//--------------------------------------------------

setdatafolder  ocdf
end

//************************************************************************************************************************
//************************************************************************************************************************
// processing the 'ResultSet' for wavelength dependent polarizability (Isotropy, mean polarizability).
function table_process_isotropy(st02):buttoncontrol
string st02
string cdf, ocdf; 
ocdf=getdatafolder(1)
newdatafolder /o /s  processed
 cdf=getdatafolder(1)
// print cdf
variable we1; string nam5
for (we1=0; we1<52; we1=we1+1)

	nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
		break
	endif

	//print nam5
	string w1a=nam5
	// print w1a
	string n_wave
	sprintf  n_wave,"%s_t",w1a 
	//print n_wave
	wave data_matrix = $w1a
	matrixop $n_wave = (abs(data_matrix))^t

	string name; name=nameofwave ($n_wave) ; print name ; 
	string fwave; sprintf fwave,"%s%s",cdf,name
	print fwave
	movewave $n_wave, $fwave

	wave data_matrix = $fwave
	DeletePoints 0,4, data_matrix
	DeletePoints 1,2, data_matrix
endfor

concat_isotropy_matrix()		//	concatenate the waves to a 2D matrix
// wave definitions ----------
string outputwName="Isotropy_all_DataSet"
wave output = $outputwName
wave omega_originalWave=root:data:omega_nm
wave omega_finalWave = root:data:omega_final
//----------------------------------

// copy the output to a new folder 'scaling'
newdatafolder /O /S scaling
duplicate /d /o output, Isotropy_all_DataSet
duplicate /d /o omega_originalWave , Omega_Original
duplicate /d /o omega_finalWave , Omega_Final
//--------------------------------------------------
setdatafolder  ocdf
end

//************************************************************************************************************************
//************************************************************************************************************************
// convert all the elements in the selected waves to positive.
function wave_abs_process()
string cdf, ocdf; 
variable we1; string nam5
for (we1=0; we1<30; we1=we1+1)
	 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
		break
	endif
	wave input=$nam5
	input=abs(input)
endfor
end

//************************************************************************************************************************
//************************************************************************************************************************

function remove_rows_data_matrix ()
variable we1; string nam5
for (we1=0; we1<30; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

//print nam5
string w1a=nam5
wave data_matrix = $w1a
DeletePoints 0,2, data_matrix
DeletePoints 1,2, data_matrix
endfor
end
//************************************************************************************************************************
//************************************************************************************************************************

// concatenate the isotropy data waves.
function concat_isotropy_matrix()
string sel
sel=wavelist("isotropy*_t",";","")
print sel
concatenate /NP =0 sel, Isotropy_all_DataSet
end
//************************************************************************************************************************
//************************************************************************************************************************

// concatenate the anisotropy, gamma waves.
function concat_gamma_matrix()
string sel
sel=wavelist("gamma*_t",";","")
print sel
concatenate /NP =0 sel, Gamma_all_DataSet
end

//************************************************************************************************************************
//************************************************************************************************************************
//Select the parameter wave , gamma or mean polarizability and run it.

function run_batch_meanpol_fullSet()
variable we1; string nam5
for (we1=0; we1<30; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

//print nam5
string w1a=nam5
print w1a
string distance_wave="g_distance" 
batch_expVal_SPmeanPol_opt ($w1a,  $distance_wave)

endfor
end
//************************************************************************************************************************

function batch_expVal_SPmeanPol_opt(wave_a, x_axis_wave):buttoncontrol  
//Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// spline wave stored for run // for each wave //
wave wave_a, x_axis_wave

string g_wave=nameofwave(wave_a)
string g_dis_wave=nameofwave(x_axis_wave)

string s25
//svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
//svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet

variable i1=0, i2=0

i2=((jf-ji)+1)
	//print i2
string s6
sprintf s6,"%s_ResultSet",g_wave
print s6
make /d /O /n=(i2,5) $s6=0
wave gjs=$s6

	
	for (i1=ji ; i1<(9) ; i1=i1+1)
	
	if (i1==0)
		n3="v0J0_norm"
		n4="v0J0_norm"	
	elseif(i1==1)
		n3="v1J0_norm"
		n4="v1J0_norm"
	elseif (i1==2)
		n3="v1J0_norm"
		n4="v0J0_norm"
	elseif (i1==3)
		n3="v1J1_norm"
		n4="v0J1_norm"
	elseif(i1==4)
		n3="v1J2_norm"
		n4="v0J2_norm"
	elseif(i1==5)	
		n3="v1J3_norm"
		n4="v0J3_norm"
	elseif (i1==6)	
		n3="v1J4_norm"
		n4="v0J4_norm"
	elseif (i1==7)
		n3="v1J5_norm"
		n4="v0J5_norm"
	elseif (i1==8)	
		n3="v1J6_norm"
		n4="v0J6_norm"
	endif
	
	
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
//	sprintf n3,"v0J%g_norm",i1
//	sprintf n4,"v0J%g_norm",(i1+2)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	printf  "*****  %s <-> %s (%s, %s)  ***** \r",wav1,wav2,wav4,wav3
	//------------------------------------------------------------------------
	
	
//**CLEAR OLD WAVES *****************
killwaves /Z int_x,scaled_y

//*********************************************
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param : waves input

string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

s4=nameofwave (param)			//parameter wave


variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq


//print d_wf,d_p
//RANGES------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:Bint_min
nvar b= root:Params:RoVibHamiltonian:Bint_max
nvar h= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a,b,h
variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------

//1. SPLINES
string array1=nameofwave(w1a)+"_a";
variable w1
w1= exists (array1)
//print w1
if (w1==1)
printf "array wave found: %s\r", array1
else
printf "array wave NOT found: %s\r", array1
 customsplineu_auto(w1a,wf_x)
endif 

s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

//---------------------------------------------------
array1=nameofwave(w2a)+"_a";
variable w2
w2= exists (array1)
if (w2==1)
printf "array wave found: %s\r", array1
else
printf "array wave NOT found: %s\r", array1
customsplineu_auto(w2a,wf_x)
endif 
s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2
//---------------------------------------------------
array1=nameofwave(rsq)+"_a";
variable w3
w3= exists (array1)
if (w3==1)
printf "array wave found: %s\r", array1
else
printf "array wave NOT found: %s\r", array1
customsplineu_auto(rsq,wf_x)
endif 
sq=nameofwave(rsq)+"_sc_intx"
ScaledY_gen(int_x,$array1)
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $sq
wave r_sq=$sq
//--------------------------------------------------
//---------------------------------------------------

customsplineu_auto(param,param_x)
array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3


//remove USED waves  for the three splines
killwaves /Z $array1,scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
//print s5
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] * (r_sq[p])
wave integral=$s5

killwaves /Z $s2,$s3

//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
//print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1

//------------------------------------------------------------------------
nvar res= intResult

//RESULT ASSIGNMENT
gjs[i1][0]=i1
gjs[i1][1]=(i1+2)
gjs[i1][2]=(res)
gjs[i1][3]=(res * (1.4818e-25))
gjs[i1][4]=((res * (1.4818e-25))^2)

//print /D res*1.4818e-25
//variable /G gamma_cm_sq=(res*1.4818e-25)^2
//variable /G gamma_au=(res)
//print /D   gamma_cm_sq
printf "---- SOLVED FOR J=%g <-> J=%g ----\r",i1,(i1+2)
killwaves /Z GLpoints1	,int_x,scaled_y
	//------------------------------------------------------------------------
	nvar res= intResult
//	printf "res(a.u.) = %1.5g | res(cm^-3) = %1.5g \r"res, (res*1.4818e-25)
//	printf "res(cm^-6) = %1.5g\r",(res*1.4818e-25)^2
//----------------------------------------------------------------------------
killwaves /Z $array1, integral, $s3
string kw=nameofwave(w1a)+"_sc_intx"
string kw2=nameofwave(w2a)+"_sc_intx"
string kw3=s3
printf "%s | %s | %s\r", kw,kw2,kw3
killwaves /z $kw,$kw2,$kw3
//----------------------------------------------------------------------------
endfor

// clean the used waves ----------------------------
killwaves /z int_x, $array1 ;  killwaves /z GLpoints1
//-----------------------------------------------------------
print "Done !"
end

//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************

// A modular function which accepts 5 waves and calculates the integral using custom spline and Gauss_Legendre quadrature.

// wf1 = wavefunction 1
// wf2 = wavefunction 2
// parameter = parameter like alpha or gamma or dipole moment or etc....
// wfx = x-axis of the wavefunction
// paramx= xaxis of the parameter
// init = initial J number

function matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
wave wf1,wf2, parameter,paramx,wfx
variable init

string wav1=NameOfWave(wf1 )
string wav2=NameOfWave( wf2)
string wav3=NameOfWave(parameter )
string wav4=NameOfWave( paramx)
string wav5=NameOfWave( wfx)


nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts

//---------------------------------------------
variable i1=init
//**CLEAR OLD WAVES *****************
killwaves /Z int_x,scaled_y

//*********************************************
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param : waves input

string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

s4=nameofwave (param)			//parameter wave


variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /o /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq
  


//print d_wf,d_p
//RANGES------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:Bint_min
nvar b= root:Params:RoVibHamiltonian:Bint_max
nvar h= root:Params:RoVibHamiltonian:step_h
// printf "int min=%g , int_max=%g , int_step=%g \r",a,b,h
variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------

//1. SPLINES
string array1=nameofwave(w1a)+"_a";
variable w1
w1= exists (array1)
//print w1
if (w1==1)
	// print "array wave found"
else
	// print "Array wave not found."
	 customsplineu_auto(w1a,wf_x)
endif 

s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

//---------------------------------------------------
array1=nameofwave(w2a)+"_a";
variable w2
w2= exists (array1)
if (w2==1)
	// print "array wave found"
else
	// print "Array wave not found."
	customsplineu_auto(w2a,wf_x)
endif 
s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2
//---------------------------------------------------
array1=nameofwave(rsq)+"_a";
variable w3
w3= exists (array1)
if (w3==1)
	// print "array wave found"
else
	// print "Array wave not found."
	customsplineu_auto(rsq,wf_x)
endif 
sq=nameofwave(rsq)+"_sc_intx"

ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $sq
wave r_sq=$sq
//--------------------------------------------------

//---------------------------------------------------
// PARAMETER LIKE ANISOTROPY OR MEAN POLARIZABILITY FOR INTEGRATION

array1=nameofwave(param)+"_a";
variable w4
w4= exists (array1)
if (w4==1)
// print "array wave found"
else

// print "Array wave not found."
customsplineu_auto(param,param_x)
endif



array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
// killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3



//remove USED waves  for the three splines
killwaves /Z scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
//print s5
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] * (r_sq[p])
wave integral=$s5

killwaves /Z $s2

//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
//print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
//printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1

//------------------------------------------------------------------------
nvar res= intResult
variable me_au
//RESULT ASSIGNMENT
me_au=(res)
//variable /G gamma_au=(res)
killwaves /Z GLpoints1	,int_x,scaled_y

//----------------------------------------------------------------------------
killwaves /Z $array1, integral //, $s3  // keep the array for later use in subsequent iterations
string kw=nameofwave(w1a)+"_sc_intx"
string kw2=nameofwave(w2a)+"_sc_intx"
string kw3=s3
//printf "<%s | %s | %s>\r", kw,kw2,kw3
killwaves /z $kw,$kw2,$kw3
//----------------------------------------------------------------------------
return (me_au)
print me_au
end

//*******************************************************************************************************************
//*******************************************************************************************************************
function batch_integral_me_gamma(s55):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s55
svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)

variable t0=ticks

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet

variable i1=0, i2=0

i2=((jf-ji)+1)
	//print i2

//---making wave to keep result set --//
string s6
sprintf s6,"%s_ResultSet",g_wave
print s6
make /d /O /n=(i2+1,5) $s6=0
wave gjs=$s6
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",(g0)
	sprintf n4,"v0J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1][0]=0
	gjs[i1][1]=0
	gjs[i1][2]=(result)
	gjs[i1][3]=(result * (1.4818e-25))
	gjs[i1][4]=((result * (1.4818e-25))^2)

	//------------------------------------------------------------------------
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v0J%g_norm",(i1)
	sprintf n4,"v0J%g_norm",(i1+2)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=(i1)
	gjs[i1+1][1]=(i1+2)
	gjs[i1+1][2]=(result)
	gjs[i1+1][3]=(result * (1.4818e-25))
	gjs[i1+1][4]=((result * (1.4818e-25))^2)
//------------------------------------------------------------------------------	
	

endfor
nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array
printf "Done ! ( %g sec )\r",(ticks-t0)/60

end
//*******************************************************************************************************************
//*******************************************************************************************************************

function batch_integral_me_meanpol(s55):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s55
svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)

variable t0=ticks

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet

variable i1=0, i2=0

i2=((jf-ji)+1)
	//print i2

//---making wave to keep result set --//
string s6
sprintf s6,"%s_ResultSet",g_wave
print s6
make /d /O /n=(i2+2,7) $s6=0
wave gjs=$s6

printf "Started : matrix elements bar(alpha) %s -----------\r",g_wave
//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",(g0)
	sprintf n4,"v0J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	
	gjs[i1][0]=0
	gjs[i1][1]=0
	gjs[i1][2]=0
	gjs[i1][3]=0
	gjs[i1][4]=(result)
	gjs[i1][5]=(result * (1.4818e-25))
	gjs[i1][6]=((result * (1.4818e-25))^2)

	//------------------------------------------------------------------------
	//-----------v=1,J=0;v=1,J=0----------------------------------------------------	
	sprintf n3,"v1J%g_norm",(g0)
	sprintf n4,"v1J%g_norm",(g0)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=1
	gjs[i1+1][1]=0
	gjs[i1+1][2]=1
	gjs[i1+1][3]=0
	gjs[i1+1][4]=(result)
	gjs[i1+1][5]=(result * (1.4818e-25))
	gjs[i1+1][6]=((result * (1.4818e-25))^2)
//------------------------------------------------------------------------------	
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v1J%g_norm",(i1)
	sprintf n4,"v0J%g_norm",(i1)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+2][0]=1
	gjs[i1+2][1]=i1
	gjs[i1+2][2]=0
	gjs[i1+2][3]=i1
	gjs[i1+2][4]=(result)
	gjs[i1+2][5]=(result * (1.4818e-25))
	gjs[i1+2][6]=((result * (1.4818e-25))^2)
//------------------------------------------------------------------------------	
	


endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array


nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0

printf "Done ! ( %g sec )\r",(ticks-t0)/60
end
//*******************************************************************************************************************
//*******************************************************************************************************************
// Following function is an extension of above function (for gamma) for calculating att the matrix elements in one go, i.e. for all wavelengths...
// This is only for the GROUND VIBRATIONAL STATE (v=0).
// all gamma_waves should have same x-axis wave which is selected from the `g_dis_wave` parameter.

function Set_integral_me_gamma_v0(s55):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s55
// svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)
newdatafolder/o /s result_sets
string cdf2=getdatafolder(1)
//print cdf1,cdf2
setdatafolder cdf1
variable t0=ticks

variable we1; string nam5 // main LOOP runnning over the gamma_wavelength waves ------------------------
for (we1=0; we1<75 ; we1=we1+1)
	 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
		break
	endif
string g_wave=nam5
//---------------------------------------------------------
// Generating result set wave to keep calculated matrix elements:
variable i1=0, i2=0

i2=((jf-ji)+1)

//---making wave to keep result set --//
string s6
wave parameter_wave =$g_wave
string sname2=nameofwave(parameter_wave)
sprintf s6,"%s%s_ResultSet",cdf2,sname2
print s6
make /d /O /n=(i2+1,5) $s6=0
wave gjs=$s6
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------------------

string n1,n2
string n3,n4

//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",(g0)
	sprintf n4,"v0J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1][0]=0
	gjs[i1][1]=0
	gjs[i1][2]=(result)
	gjs[i1][3]=(result * (1.4818e-25))
	gjs[i1][4]=((result * (1.4818e-25))^2)

	//------------------------------------------------------------------------
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v0J%g_norm",(i1)
	sprintf n4,"v0J%g_norm",(i1+2)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=(i1)
	gjs[i1+1][1]=(i1+2)
	gjs[i1+1][2]=(result)
	gjs[i1+1][3]=(result * (1.4818e-25))
	gjs[i1+1][4]=((result * (1.4818e-25))^2)
//------------------------------------------------------------------------------	
	

endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array



endfor

nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0
printf "Done ! ( %g sec )\r",(ticks-t0)/60

setdatafolder cdf1
end
//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
// Following function is an extension of above function (for gamma) for calculating att the matrix elements in one go, i.e. for all wavelengths...
// This is only for the FIRST VIBRATIONAL STATE (v=1).
// all gamma_waves should have same x-axis wave which is selected from the `g_dis_wave` parameter.

function Set_integral_me_gamma_v1(s54):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s54
// svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)
newdatafolder/o /s result_sets
string cdf2=getdatafolder(1)
//print cdf1,cdf2
setdatafolder cdf1
variable t0=ticks

variable we1; string nam5 // main LOOP runnning over the gamma_wavelength waves ------------------------
for (we1=0; we1<75; we1=we1+1)
	 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
		break
	endif
string g_wave=nam5
//---------------------------------------------------------
// Generating result set wave to keep calculated matrix elements:
variable i1=0, i2=0

i2=((jf-ji)+1)

//---making wave to keep result set --//
string s6
wave parameter_wave =$g_wave
string sname2=nameofwave(parameter_wave)
sprintf s6,"%s%s_ResultSet",cdf2,sname2
print s6
make /d /O /n=(i2+1,5) $s6=0
wave gjs=$s6
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------------------

string n1,n2
string n3,n4

//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v1J%g_norm",(g0)
	sprintf n4,"v1J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1][0]=0
	gjs[i1][1]=0
	gjs[i1][2]=(result)
	gjs[i1][3]=(result * (1.481847096e-25))
	gjs[i1][4]=((result * (1.481847096e-25))^2)

	//------------------------------------------------------------------------
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v1J%g_norm",(i1)
	sprintf n4,"v1J%g_norm",(i1+2)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=(i1)
	gjs[i1+1][1]=(i1+2)
	gjs[i1+1][2]=(result)
	gjs[i1+1][3]=(result * (1.481847096e-25))
	gjs[i1+1][4]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	
	

endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array



endfor

nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0
printf "Done ! ( %g sec )\r",(ticks-t0)/60

setdatafolder cdf1
end
//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
// Following function is an extension of above function (for gamma) for calculating att the matrix elements in one go, i.e. for all wavelengths...
// This is only for the SECOND VIBRATIONAL STATE (v=2).
// all gamma_waves should have same x-axis wave which is selected from the `g_dis_wave` parameter.

function Set_integral_me_gamma_v2(s54):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run for each wave //
string s54
// svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final
string cdf1=getdatafolder(1)
newdatafolder/o /s result_sets
string cdf2=getdatafolder(1)
//print cdf1,cdf2
setdatafolder cdf1
variable t0=ticks

variable we1; string nam5 // main LOOP runnning over the gamma_wavelength waves ------------------------
for (we1=0; we1<75; we1=we1+1)
	 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
		break
	endif
string g_wave=nam5
//---------------------------------------------------------
// Generating result set wave to keep calculated matrix elements:
variable i1=0, i2=0

i2=((jf-ji)+1)

//---making wave to keep result set --//
string s6
wave parameter_wave =$g_wave
string sname2=nameofwave(parameter_wave)
sprintf s6,"%s%s_ResultSet",cdf2,sname2
print s6
make /d /O /n=(i2+1,5) $s6=0
wave gjs=$s6
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------------------

string n1,n2
string n3,n4

//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v2J%g_norm",(g0)
	sprintf n4,"v2J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1][0]=0
	gjs[i1][1]=0
	gjs[i1][2]=abs(result)
	gjs[i1][3]=abs(result * (1.481847096e-25))
	gjs[i1][4]=((result * (1.481847096e-25))^2)

	//------------------------------------------------------------------------
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v2J%g_norm",(i1)
	sprintf n4,"v2J%g_norm",(i1+2)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=(i1)
	gjs[i1+1][1]=(i1+2)
	gjs[i1+1][2]=abs(result)
	gjs[i1+1][3]=abs(result * (1.481847096e-25))
	gjs[i1+1][4]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	
	

endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array



endfor

nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0
printf "Done ! ( %g sec )\r",(ticks-t0)/60

setdatafolder cdf1
end
//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
// Following function is an extension of above function (for mean polarizability) for calculating att the matrix elements in one go, i.e. for all wavelengths...

// all isotropy_waves should have same x-axis wave which is selected from the `g_dis_wave` parameter.

function Set_integral_me_meanPolv0v1 (s55):buttoncontrol  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s55
// svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final

string cdf1=getdatafolder(1)
newdatafolder/o /s result_sets
string cdf2=getdatafolder(1)
//print cdf1,cdf2
setdatafolder cdf1

variable t0=ticks



variable we1; string nam5 // main LOOP runnning over the gamma_wavelength waves ------------------------
for (we1=0; we1<75; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

print nam5
string g_wave=nam5
//-------------------------------------------------------
// Generating wave to keep the result set: 
variable i1=0, i2=0

i2=((jf-ji)+1)

//---making wave to keep result set --//
string s6
wave parameter_wave =$g_wave
string sname2=nameofwave(parameter_wave)
sprintf s6,"%s%s_ResultSet",cdf2,sname2
print s6
make /d /o /n=(i2+2,7) $s6=0
wave gjs=$s6
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------------------

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet
//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v0J%g_norm",(g0)
	sprintf n4,"v0J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	
	gjs[i1][0]=0
	gjs[i1][1]=0
	gjs[i1][2]=0
	gjs[i1][3]=0
	gjs[i1][4]=(result)
	gjs[i1][5]=(result * (1.481847096e-25))
	gjs[i1][6]=((result * (1.481847096e-25))^2)

	//------------------------------------------------------------------------
	//-----------v=1,J=0;v=1,J=0----------------------------------------------------	
	sprintf n3,"v1J%g_norm",(g0)
	sprintf n4,"v1J%g_norm",(g0)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=1
	gjs[i1+1][1]=0
	gjs[i1+1][2]=1
	gjs[i1+1][3]=0
	gjs[i1+1][4]=(result)
	gjs[i1+1][5]=(result * (1.481847096e-25))
	gjs[i1+1][6]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v1J%g_norm",(i1)
	sprintf n4,"v0J%g_norm",(i1)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+2][0]=1
	gjs[i1+2][1]=i1
	gjs[i1+2][2]=0
	gjs[i1+2][3]=i1
	gjs[i1+2][4]=(result)
	gjs[i1+2][5]=(result * (1.481847096e-25))
	gjs[i1+2][6]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	
	


endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array


endfor // main loop
nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0
end
//*******************************************************************************************************************
//*******************************************************************************************************************
//*******************************************************************************************************************
// Following function is an extension of above function (for mean polarizability) for calculating all the matrix elements in one go, i.e. for all wavelengths...
// For TRANSITIONS CONNECTING v=1 and v=2, mean polarizability
// all isotropy_waves should have same x-axis wave which is selected from the `g_dis_wave` parameter.

function Set_integral_me_meanPolv1v2(s55):buttoncontrol  //Batch function for expVal_anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s55
// svar g_wave =root:Params:RoVibHamiltonian:b_gwave 
svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final

string cdf1=getdatafolder(1)
newdatafolder/o /s result_sets
string cdf2=getdatafolder(1)
//print cdf1,cdf2
setdatafolder cdf1

variable t0=ticks		// for timer

variable we1; string nam5 // main LOOP runnning over the gamma_wavelength waves ------------------------

for (we1=0; we1<75; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

print nam5
string g_wave=nam5
//-------------------------------------------------------
// Generating wave to keep the result set: 
variable i1=0, i2=0

i2=((jf-ji)+1)

//---making wave to keep result set --//
string s6
wave parameter_wave =$g_wave
string sname2=nameofwave(parameter_wave)
sprintf s6,"%s%s_ResultSet",cdf2,sname2
print s6
make /d /o /n=(i2+1,7) $s6=0
wave gjs=$s6
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------------------

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet
//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result


 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
	//------------------------------------
	//	printf "Solving for all pairs for %s",g_wave
	sprintf n3,"v2J%g_norm",(g0)
	sprintf n4,"v2J%g_norm",(g0)
//	print i1, n3,n4
	string wav1=n3		; string wav2=n4
	string wav3=g_wave	; string wav4 = g_dis_wave
	string wav5=wfn_dis
	
	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
		
	gjs[i1][0]=2
	gjs[i1][1]=0
	gjs[i1][2]=2
	gjs[i1][3]=0
	gjs[i1][4]=abs(result)
	gjs[i1][5]=abs(result * (1.481847096e-25))
	gjs[i1][6]=((result * 1.481847096e-25)^2)

	//------------------------------------------------------------------------
	
for (i1=ji ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v2J%g_norm",(i1)
	sprintf n4,"v1J%g_norm",(i1)
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  \r",wav1,wav3,wav4,wav2
	
	gjs[i1+1][0]=2
	gjs[i1+1][1]=i1
	gjs[i1+1][2]=1
	gjs[i1+1][3]=i1
	gjs[i1+1][4]=abs(result)
	gjs[i1+1][5]=abs(result * (1.481847096e-25))
	gjs[i1+1][6]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	
	
endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array


endfor // main loop
nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0
end
//*******************************************************************************************************************
//*******************************************************************************************************************

function print_TE_delG_q(s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- Delta G energies ---"
string cdf=getdatafolder(1) ;// print cdf ;
variable i=0
for (i=0 ; (i+1)<(10) ; i=i+1)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,0
//sprintf addr2,"%sJ%g:Calc",cdf,(i)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(i)]

wave eval=eval
b=eval[(i+1)]

//print "S1-branch"
printf "%6.8f\r", ( (b-a)*219474.6313702000000000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   S1 branch
end
//*******************************************************************************************************************


//*******************************************************************************************************************
//*******************************************************************************************************************

// Following function is an extension of above function (for gamma) for calculating att the matrix elements in one go, i.e. for all wavelengths...
// all anisotropy_waves should have same x-axis wave which is selected from the `g_dis_wave` parameter.

function Set_integral_me_v0v1_batch ( )  //Batch function for expVal_gamma or anisotropy_OPTIMIZED for speed.
// with modular implementation of the integration.
// spline wave stored for run // for each wave //
string s55


svar g_dis_wave =root:Params:RoVibHamiltonian:b_gWave_dis
svar wfn_dis=root:Params:RoVibHamiltonian:b_wd

nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts
nvar ji =root:Params:RoVibHamiltonian:J_ini
nvar  jf= root:Params:RoVibHamiltonian:J_final

string cdf1=getdatafolder(1)
newdatafolder/o /s result_sets
string cdf2=getdatafolder(1)
//print cdf1,cdf2
setdatafolder cdf1

variable t0=ticks

variable we1; string nam5 // main LOOP runnning over the gamma_wavelength waves ------------------------
for (we1=0; we1<75; we1=we1+1)

 nam5 = getbrowserselection(we1)
  	if (strlen(nam5) == 0)
	break
	endif

print nam5
string g_wave=nam5
//-------------------------------------------------------
// Generating wave to keep the result set: 
variable i1=0, i2=0

i2=((jf-ji)+1)

//---making wave to keep result set --//
string s6, s7
wave parameter_wave =$g_wave
string sname2=nameofwave(parameter_wave)
sprintf s6,"%s%s_ResultSet",cdf2,sname2
sprintf s7,"%s%s_ResultSet_2",cdf2,sname2
print s6, s7
make /d /o /n=(i2+2,7) $s6=0
make /d /o /n=(i2+2,7) $s7=0
wave gjs=$s6
wave gjs2=$s7
printf "Started : matrix elements gamma %s -----------\r",g_wave
//---------------------------------------------------------

string n1,n2
string n3,n4

killwaves /z GammaJ_ResultSet
//---------------------------------------------//
variable g0=0 // special variable for running for matrix element for <00|gamma|00>
variable result

//-------------------------------------------------------------------
string wav1=n3		; string wav2=n4
string wav3=g_wave	; string wav4 = g_dis_wave
string wav5=wfn_dis
//-------------------------------------------------------------------

 //    matrix_element_integral_mod(wf1,parameter, wf2, wfx, paramx,init)
jf=13
for (i1=0 ; i1<=(jf) ; i1=i1+1)
	
//------------------------------------------------------------------------------	
	sprintf n3,"v0J%g_norm",(i1)
	sprintf n4,"v1J%g_norm",(i1+2)
	
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  %g, %g \r",wav1,wav3,wav4,wav2, i1, i1+2
	
	gjs[i1][0]=0
	gjs[i1][1]=i1
	gjs[i1][2]=1
	gjs[i1][3]=i1+2
	gjs[i1][4]=(result)
	gjs[i1][5]=(result * (1.481847096e-25))
	gjs[i1][6]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	

//------------------------------------------------------------------------------	
	sprintf n3,"v0J%g_norm",(i1+2)
	sprintf n4,"v1J%g_norm",(i1)
	
	 wav1=n3		; wav2=n4
	 wav3=g_wave	; wav4 = g_dis_wave
	 wav5=wfn_dis

	wave wf1=$n3
	wave param=$wav3
	wave wf2=$n4
	wave paramx=$wav5
	wave wfx=$wav4
	
	result=matrix_element_integral_mod(wf1,param,wf2,paramx,wfx,0)
	printf  " < %s | (%s, %s)  | %s >  %g, %g \r",wav1,wav3,wav4,wav2, i1, i1+2
	
	gjs2[i1][0]=0
	gjs2[i1][1]=i1+2
	gjs2[i1][2]=1
	gjs2[i1][3]=i1
	gjs2[i1][4]=(result)
	gjs2[i1][5]=(result * (1.481847096e-25))
	gjs2[i1][6]=((result * (1.481847096e-25))^2)
//------------------------------------------------------------------------------	
	


endfor

string param_array=nameofwave(param)+"_a";
killwaves /Z $param_array


endfor // main loop
nvar a0 =root:Params:RoVibHamiltonian:Bint_min
nvar b0= root:Params:RoVibHamiltonian:Bint_max
nvar h0= root:Params:RoVibHamiltonian:step_h
printf "int min=%g , int_max=%g , int_step=%g \r",a0,b0,h0
end
//*******************************************************************************************************************
//*******************************************************************************************************************


//****************************************************************************************************

// Extract the energy in (cm-1) for the rovibrational states of the molecule for which the 
// computation has been performed.
// Function will generate two waves in the folder 'root:Params:RoVibHamiltonian:'
// with the name  'prefix_string'eV0, 'prefix_string'eV1 for the energy of the rotational states within the ground and the first vibrational state.

function extract_rovibstates( prefix_string)
//	Params
//	prefix_string		=	string to be appended to the name of the energy wave
string prefix_string

nvar jm=root:Params:RoVibHamiltonian:j_max
variable v0Jmax
variable v0J0
variable v1J0
string addr1,addr2
variable ej
variable i,j
variable a,b
print "--- Energies of the states (relative to the ground state ) ---"

string cdf=getdatafolder(1)

sprintf addr1,"%sJ0:Calc",cdf
setdatafolder addr1
print addr1
wave eval=eval
v0J0=eval[(0)]
v1J0=eval[(1)]
v0Jmax = (v1J0-v0J0)
//	print /d v1J0,v0J0,(v1J0-v0J0)*219474.63137020


// Find out the max J state possible for the ground vib state.
// For the ground state ==> Jmax is 15
// For the first vibrational state ==> Jmax is 11

//----------------------------------------------
// make the waves to keep data
setdatafolder root:Params:RoVibHamiltonian:
string wname
sprintf wname, "%seV0",prefix_string
//	print wname

make /o /d /n=(16) $wname=0
wave eV0=$wname
sprintf wname, "%seV1",prefix_string
make /o /d /n=(12) $wname=0
wave eV1=$wname
//----------------------------------------------
	
// Extract the energy of the rotational state for the ground vib state	
for (i=0 ; i < 16 ; i=i+1)

	sprintf addr1,"%sJ%g:Calc",cdf, 0
	sprintf addr2,"%sJ%g:Calc",cdf,i
	
	setdatafolder addr1
	wave eval=eval
	a=eval[(0)]

	setdatafolder addr2
	wave eval=eval
	b=eval[(0)]
	printf "%6.8f\r",((b-a )* 219474.63137020)
	eV0[i]=((b-a )* 219474.63137020)
endfor

// Extract the energy of the rotational states within the first vib state
for (i=0 ; i < 12 ; i=i+1)

	sprintf addr1,"%sJ%g:Calc",cdf, 0
	sprintf addr2,"%sJ%g:Calc",cdf,i
	
	setdatafolder addr1
	wave eval=eval
	a=eval[(0)]

	setdatafolder addr2
	wave eval=eval
	b=eval[(1)]
	printf "%6.8f\r",((b-a )* 219474.63137020)
	eV1[i]=((b-a )* 219474.63137020)
endfor

setdatafolder cdf

print "----------------------------------------------"
end

//****************************************************************************************************

function expectation_val (s7d):buttoncontrol
string s7d
svar wav1 = root:Params:RoVibHamiltonian:iwv1			// wavefunction
svar wav2 =root:Params:RoVibHamiltonian:iwv2			// wavefunction
svar wav3 =root:Params:RoVibHamiltonian:iwv3			// parameter selected
svar wav5 = root:Params:RoVibHamiltonian:int_distance	// wavefunction distance
svar wav4 =root:Params:RoVibHamiltonian:iwv4			// parameter_wave_distance
nvar int_pnts1=root:Params:RoVibHamiltonian:int_pnts	// integration points ( Quadrature)
//**CLEAR OLD WAVES -----------------------
killwaves /Z int_x,scaled_y
// --------------------------------------------------------
string cdf2=GetDataFolder(1)
string s1,s2,s3,s4,sq
sprintf s1,"%s%s",cdf2,wav1
wave w1a=$s1

sprintf s2,"%s%s",cdf2,wav2
wave w2a=$s2
sprintf s3,"%s%s",cdf2,wav3
wave param=$s3				// w1a , w2a, param , r_sq: waves input


string d_wf,d_p
sprintf d_p,"%s%s",cdf2,wav4
wave param_x=$d_p

sprintf d_wf,"%s%s",cdf2,wav5
wave wf_x=$d_wf

variable dimx = dimsize(wf_x,0)		// r_square term in the integrand.
make /o /d /n=(dimx) rsq= (wf_x[p])^2 ; wave rsq=rsq
s4=nameofwave (param)			//parameter wave

	print "------------------------------------------------------------------------------------------------------" 
//	print  d_wf,d_p
//RANGES-------------------------------------------------------------------------------------------
nvar a =root:Params:RoVibHamiltonian:int_min		//	int_min
nvar b= root:Params:RoVibHamiltonian:int_max		// int_max
nvar h= root:Params:RoVibHamiltonian:step_h		// step ( to be used for the final waves which go for multiplication)

variable f1=a-h, num=(((b-a)/h)+2)

make /O /D /n=(num) int_x=a							// int_x = integration x axis, of defined range and step size
wave int_x=int_x; variable x2=dimsize(int_x,0)
variable i=0
for (i=0; i < (x2) ; i=i+1)
	int_x[i]=int_x[i]+(h*i)
endfor 

//-------------------------------------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------------------

//1. SPLINES
string array1=nameofwave(w1a)+"_a";			// wfn 1
variable w1
w1= exists (array1)
//print w1
if (w1==1)
	printf "array wave found = %s\r", array1
else
	printf "array wave NOT found = %s\r", array1
 customsplineu_auto(w1a,wf_x)
endif 

s1=nameofwave(w1a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
wave scaled_y = scaled_y ; 
duplicate /O /D scaled_y, $s1 
wave wf1=$s1

//---------------------------------------------------
array1=nameofwave(w2a)+"_a";					// wfn 2
variable w2
w2= exists (array1)
if (w2==1)
	printf "array wave found = %s\r", array1
else
	printf "array wave NOT found = %s\r", array1
customsplineu_auto(w2a,wf_x)
endif 
s2=nameofwave(w2a)+"_sc_intx"
ScaledY_gen(int_x,$array1)
//killwaves /Z $array1
duplicate /O /D scaled_y, $s2
wave wf2=$s2
//---------------------------------------------------
array1=nameofwave(rsq)+"_a";					// r square
variable w3
w3= exists (array1)
if (w3==1)
	printf "array wave found = %s\r", array1
else
	printf "array wave NOT found = %s\r", array1
customsplineu_auto(rsq,wf_x)
endif 
sq=nameofwave(rsq)+"_sc_intx"
ScaledY_gen(int_x,$array1)
duplicate /O /D scaled_y, $sq
wave r_sq=$sq
//--------------------------------------------------
// parameter
customsplineu_auto(param,param_x)				// parameter
array1=nameofwave(param)+"_a";	s3=nameofwave(param)+"_sc_intx"
ScaledY_gen(int_x,$array1)
killwaves /Z $array1
duplicate /O /D scaled_y, $s3
wave param_sc=$s3

//-------------------------------------------------------------------
//remove USED waves  for the three splines
//killwaves /Z $array1,scaled_y
//-----INTEGRAL  -------------------------------------------------------------------------------------------
string s5=(s4)+"_expVal"
//	print s5
// print x2
make /O /D /n=(x2) $s4+"_expVal" =wf1[p]*wf2[p]*param_sc[p] * (r_sq[p])
wave integral=$s5
//--INTEGRAL SPLINE------------------------------------------------------------------------------------
customsplineu_auto(integral,int_x)
array1=nameofwave(integral)+"_a" //;s3=nameofwave(integral)+"_sc_intx"
//wave nwave=$s4+"_expVal"
//	print array1
glnumerical_integrate(a,b,(int_pnts1),$array1)
printf "< %s |%s |%s >\r", nameofwave(w1a),  s5, nameofwave(w2a)
printf "limits : %g - %g | points : %g | %s\r",a,b,int_pnts1,array1
nvar res= intResult
printf "ME (au): %5.8f | (cm^3): %g \r ",res, res*1.481847096e-25
variable /G gamma_cm_sq=(res*1.481847096e-25)^2
variable /G gamma_au=(res)
printf "ME_sq (cm^6): %g -----------------------------------------------------------    \r",gamma_cm_sq

// clean the used waves ----------------------------
killwaves /z int_x, $array1, rsq, $s4+"_expVal" , $s4+"_sc_intx", $sq+"_sc_intx",$s2+"_sc_intx",$s3+"_sc_intx" , $s1+"_sc_intx"  ;
 killwaves /z GLpoints1
//-----------------------------------------------------------
end

//************************************************************************************************************************

function print_TE_O1_q_alt (s2n):buttoncontrol
string s2n
nvar jm=root:Params:RoVibHamiltonian:j_max
string fname,j
print "--- O1 branch transition energies ---"
string cdf=getdatafolder(1) ; // print cdf ;
variable i=0
for (i=1 ; (i+2)<(jm) ; i=i+2)

variable a,b
string addr1,addr2
sprintf addr1,"%sJ%g:Calc",cdf,(i)
sprintf addr2,"%sJ%g:Calc",cdf,(i+2)
//print addr1, addr2 

setdatafolder addr1
wave eval=eval
a=eval[(1)]

setdatafolder addr2
wave eval=eval
b=eval[(0)]

//print "O1-branch"
printf "%g->%g, %6.4f\r", (i+2) , (i), ( (a-b)  *219474.6313702000000000)

endfor
setdatafolder cdf

print "----------------------------------------------"		////   O1 branch
end
