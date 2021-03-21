#pragma rtGlobals=3		// Use modern global access method and strict wave access.

//*********************************************************************************************

function matrix_xx(s2d): buttoncontrol
string s2d
string sel
// SELECT FOLDER THEN RUN-------------------------
string cdf=getbrowserselection(0)
setdatafolder $cdf
//------------------------------------

sel=wavelist("wave*",";","")
killwaves /z  matrix_xx1
concatenate sel,matrix_xx1
setdatafolder root:
newdatafolder /O data
duplicate /O matrix_xx1,  root:data:matrix_xx1
matrixtranspose root:data:matrix_xx1
end
//*********************************************************************************************
function matrix_zz(s2e):buttoncontrol
string s2e
string sel
// SELECT FOLDER THEN RUN-------------------------
string cdf=getbrowserselection(0)
setdatafolder $cdf
//------------------------------------

sel=wavelist("wave*",";","")
killwaves /z matrix_zz1
concatenate sel,matrix_zz1
setdatafolder root:
newdatafolder /O data
duplicate /O matrix_zz1,  root:data:matrix_zz1
matrixtranspose root:data:matrix_zz1
end

//*********************************************************************************************

function prep_gamma(s2a) : ButtonControl
string s2a
setdatafolder root:data
wave m1=matrix_zz1 ; wave m2=matrix_xx1
variable x1,x2
x1=dimsize(m1,0)
x2=dimsize(m1,1)
make /O /D /n=(x1,x2) gamma_data=m1-m2;
make /O /D /n=(x1,x2)  isotropy= (m1+2*m2)/3;

//print x1,x2
end

//*********************************************************************************************

function ini_alpha()
DFREF cdf1 = GetDataFolderDFR()
 variable checkf=datafolderexists ("root:Params:alpha")
 if (checkf == 0)
 	setdatafolder  root:
 	newdatafolder /O /s Params
 	newdatafolder /O Customspline
 	newdatafolder /O /s alpha
 	
 	string /G xx_folder, zz_folder
 	string /G root:Params:Alpha:g_distance
 	string /G root:Params:Alpha:s_Omega_nm 
	variable /G root:Params:Alpha:v_Omega_nm
	string /G root:Params:Alpha:s_save_name
	string /G root:Params:Alpha:g_distance
	string /G root:Params:Alpha:scalingOmegaFinal
	string /G root:Params:Alpha:scalingOmegaOriginal
	string /G root:Params:Alpha:scalingME2D	// for scaling the matrix of MEs
	string /G root:Params:Alpha:scalingME2D_g	// for scaling to one wavelength, gamma
	string /G root:Params:Alpha:scalingME2D_i	// for scalng to one wavelength, isotropy
	variable /G root:Params:Alpha:scalingOmega_wavelength 	// wavelength	
endif
 
 	setdatafolder root: 
 	string cmd ="Alpha_analysis()"
 	execute cmd
 	setdatafolder cdf1
end
//*********************************************************************************************


Window Alpha_analysis() : Panel
	PauseUpdate; Silent 1		// building window...
	NewPanel /W=(1437,54,1700,542)
	SetDrawLayer UserBack
	DrawLine 5,145,255,145
	DrawLine -2,234,255,234
	DrawLine 3,97,251,97
	DrawLine -3,263,254,263
	DrawLine -1,309,256,309
	DrawLine -1,374,256,374
	Button Gen_Gamma,pos={31.00,76.00},size={95.00,20.00},proc=prep_gamma,title="Gen_Gamma"
	Button Gen_Gamma,fColor=(32768,32768,65280)
	TitleBox Alpha_title,pos={3.00,2.00},size={159.00,23.00},title="Alpha: Polarizability analysis"
	TitleBox Alpha_title,fSize=12,fStyle=1
	PopupMenu Gamma_,pos={2.00,125.00},size={148.00,19.00},proc=gamma_wave_Selection,title="Gamma_wave"
	PopupMenu Gamma_,mode=6,popvalue="omega_h",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	Button Gen_omega_nm,pos={152.00,106.00},size={95.00,20.00},proc=gen_gamma_nm,title="Gen_omega_nm"
	Button Gen_omega_nm,fColor=(32768,54528,65280)
	PopupMenu Omega_nm_val,pos={2.00,172.00},size={126.00,19.00},proc=omega_nm_Selection,title="Omega_nm"
	PopupMenu Omega_nm_val,mode=2,popvalue="182.253",value= #"root:Params:gen_popupString"
	Button Gen_waves_Gamma,pos={5.00,213.00},size={185.00,20.00},proc=generate_gamma_nm_waves,title="Gen_waves_anisotropy isotropy"
	Button Gen_waves_Gamma,fColor=(0,52224,26368)
	TitleBox Specific_freq_dependent_waves,pos={0.00,149.00},size={187.00,23.00},title="Specific_freq_dependent_waves"
	TitleBox Specific_freq_dependent_waves,fStyle=1
	TitleBox Convert_unit,pos={3.00,100.00},size={130.00,23.00},title="Convert_unit(H -> nm)"
	PopupMenu popup_distance_wave,pos={5.00,192.00},size={190.00,19.00},proc=gamma_rwave_sel,title="R_wave (for gamma)"
	PopupMenu popup_distance_wave,mode=7,popvalue="g_distance",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	Button make_alpha_xx,pos={34.00,28.00},size={95.00,20.00},proc=matrix_xx,title="Make_alpha_xx"
	Button make_alpha_xx,fColor=(65280,48896,48896)
	Button make_alpha_zz,pos={32.00,53.00},size={95.00,20.00},proc=matrix_zz,title="Make_alpha_zz"
	Button make_alpha_zz,fColor=(65280,21760,0)
	TitleBox title0,pos={132.00,30.00},size={123.00,23.00},title="Select folder then USE"
	TitleBox title1,pos={132.00,52.00},size={123.00,23.00},title="Select folder then USE"
	TitleBox title2,pos={13.00,28.00},size={14.00,23.00},title="1"
	TitleBox title3,pos={12.00,52.00},size={14.00,23.00},title="2"
	TitleBox title4,pos={10.00,76.00},size={14.00,23.00},title="3"
	TitleBox Specific_freq_dependent_Gamma1,pos={3.00,236.00},size={118.00,23.00},title="Save_output_waves"
	TitleBox Specific_freq_dependent_Gamma1,fStyle=1
	Button button_save,pos={139.00,236.00},size={100.00,20.00},proc=Gamma_save_itx,title="Save_itx_waves"
	Button button_save,fColor=(65280,43520,0)
	Button Arrange_Gamma_sq_val,pos={6.00,267.00},size={190.00,20.00},proc=Arrange_gamma_square,title="Arrange_Gamma_square_Values"
	Button Arrange_Gamma_sq_val,fColor=(16385,49025,65535)
	Button Process_gamma,pos={6.00,287.00},size={115.00,20.00},proc=table_process_anisotropy,title="Process_anisotropy"
	Button Process_gamma,fColor=(65535,54607,32768)
	Button Process_gamma1,pos={124.00,287.00},size={115.00,20.00},proc=table_process_isotropy,title="Process_Isotropy"
	Button Process_gamma1,fColor=(26205,52428,1)
	Button Gen_waves_Gamma1,pos={134.00,171.00},size={115.00,20.00},proc=generate_gamma_nm_waves_ALL,title="Gen_nm_waves_all"
	Button Gen_waves_Gamma1,fColor=(65535,49157,16385)
	TitleBox Scaling_ME_matrix,pos={1.00,311.00},size={72.00,23.00},title="Scaling MEs"
	TitleBox Scaling_ME_matrix,fStyle=1
	PopupMenu ME_matrix,pos={75.00,313.00},size={176.00,19.00},proc=ME2D_wave_Selection,title="MEs:2D"
	PopupMenu ME_matrix,help={"select 2D wave containing MEs along wavelength ( increasing order of wavelength)"}
	PopupMenu ME_matrix,mode=2,popvalue="Gamma_full_data_set",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu Omega_nm_original,pos={75.00,334.00},size={169.00,19.00},proc=OmegaOriginal_wave_Selection,title="Omega:Original"
	PopupMenu Omega_nm_original,help={"select 2D wave containing MEs along wavelength ( increasing order of wavelength)"}
	PopupMenu Omega_nm_original,mode=4,popvalue="omega_nm",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu Omega_nm_final,pos={63.00,354.00},size={191.00,19.00},proc=OmegaFinal_wave_Selection,title="Omega:Final(array)"
	PopupMenu Omega_nm_final,help={"select 2D wave containing MEs along wavelength ( increasing order of wavelength)"}
	PopupMenu Omega_nm_final,mode=5,popvalue="omega_final",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	Button Scale_MEmatrix,pos={0.00,335.00},size={70.00,22.00},proc=scale_MEmatrix_to_wavelength,title="Scale_MEs"
	Button Scale_MEmatrix,fColor=(16385,49025,65535)
	TitleBox Scaling_ME_matrix2,pos={95.00,376.00},size={166.00,23.00},title="Scaling MEs to 1 wavelength"
	TitleBox Scaling_ME_matrix2,fStyle=1
	PopupMenu ME_matrix1,pos={79.00,398.00},size={183.00,19.00},proc=ME2D_Gamma_wave_Selection,title="MEs:2Dg"
	PopupMenu ME_matrix1,help={"select 2D wave containing MEs along wavelength ( increasing order of wavelength)"}
	PopupMenu ME_matrix1,mode=2,popvalue="Gamma_full_data_set",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu ME_matrix2,pos={79.00,418.00},size={183.00,19.00},proc=ME2D_Iso_wave_Selection,title="MEs:2D i"
	PopupMenu ME_matrix2,help={"select 2D wave containing MEs along wavelength ( increasing order of wavelength)"}
	PopupMenu ME_matrix2,mode=3,popvalue="Isotropy_full_data_set",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	PopupMenu Omega_nm_original1,pos={71.00,440.00},size={169.00,19.00},proc=OmegaOriginal_wave_Selection,title="Omega:Original"
	PopupMenu Omega_nm_original1,help={"select 2D wave containing MEs along wavelength ( increasing order of wavelength)"}
	PopupMenu Omega_nm_original1,mode=4,popvalue="omega_nm",value= #"\"_none_;\" + WaveList(\"*\",\";\",\"\")"
	SetVariable Set_omega_nm,pos={35.00,462.00},size={175.00,18.00},bodyWidth=80,title="Omega_final: nm"
	SetVariable Set_omega_nm,format="%7.6f",valueBackColor=(49151,60031,65535)
	SetVariable Set_omega_nm,limits={5,3000,1e-007},value= root:Params:alpha:scalingOmega_wavelength,live= 1
	Button Scale_MEmatrix1,pos={4.00,379.00},size={64.00,45.00},proc=scale_MEmatrix_to_OneWavelength,title="Scale_ME_\rto final"
	Button Scale_MEmatrix1,fColor=(32768,40777,65535)
EndMacro


//*********************************************************************************************
Function gamma_wave_Selection(gamma_sel1,pop_gw1,popStr_gw1) : PopupMenuControl	
String gamma_sel1
Variable pop_gw1
String popStr_gw1
string /G root:Params:Alpha:gwave = popStr_gw1
End

//*********************************************************************************************
function gen_popupwave_string(wv1)
wave wv1
variable x2= dimsize (wv1,0)
print x2
variable i1; string s1="0"
for (i1=0; i1 < x2 ; i1=i1+1)

sprintf s1,"%s;%g",s1,wv1[i1]
//print s1
string /G root:Params:gen_popupString=s1
endfor
end
//*********************************************************************************************

function gen_gamma_nm(s2b) : buttoncontrol
string s2b
svar gamma_wave=root:Params:alpha:gwave
wave gamma_wave1=$gamma_wave
variable x1=dimsize(gamma_wave1,0)
make /d /n=(x1) omega_nm=(1e7/(  gamma_wave1[p]*219474.6313702 )) 
gen_popupwave_string(omega_nm)
end
//*********************************************************************************************
Function omega_nm_Selection(omega_nm_sel1,pop_on1,popStr_on1) : PopupMenuControl	
String  omega_nm_sel1
Variable pop_on1
String popStr_on1

string /G root:Params:Alpha:s_Omega_nm = popStr_on1
variable /G root:Params:Alpha:v_Omega_nm = str2num(popStr_on1)
//variable /G  = popStr_on1
//string cmd="variable  root:Params:Alpha:v_Omega_nm =popStr_on1"
//execute cmd
End



//*********************************************************************************************

//*********************************************************************************************

Function gamma_rwave_sel(rw_sel1,pop_rw1,popStr_rw1) : PopupMenuControl	
String  rw_sel1
Variable pop_rw1
String  popStr_rw1
string /G root:Params:Alpha:g_distance =popStr_rw1
End

//*********************************************************************************************
function generate_gamma_nm_waves(s2c):buttoncontrol
string s2c
string cdf=getdatafolder(1)
setdatafolder root:data
nvar num=root:Params:alpha:v_Omega_nm
svar sname =root:Params:alpha:s_Omega_nm 
wave omw=root:data:omega_nm
svar rwave=root:Params:alpha:g_distance
wave r_wave= $rwave
//print num
string cdf1=getdatafolder(1)
findvalue /V=(num) /T=0.1 omw 
variable n1= V_value
wave gd1=root:data:gamma_data
variable d1=dimsize(gd1,0)
variable d2=dimsize(gd1,1)
string wname ; sprintf wname,"gamma_%s",sname
string nname= ReplaceString(".",wname, "_")
//print nname;
setdatafolder root:
newdatafolder /O /S Output

make /D /n=(d1) /O  $nname=gd1[p][n1]
make /D /n=(d1) /O  g_distance=r_wave[p]

string /G root:Params:Alpha:s_save_name =nname
//print d1, n1
SetDataFolder cdf1
//printf "Frequency dep. gamma wave generated. for %g nm\r",num
generate_isotropy_nm_waves()
end
//*********************************************************************************************
//*********************************************************************************************
function generate_isotropy_nm_waves()
string cdf=getdatafolder(1)
setdatafolder root:data
nvar num=root:Params:alpha:v_Omega_nm
svar sname =root:Params:alpha:s_Omega_nm 
wave omw=root:data:omega_nm
svar rwave=root:Params:alpha:g_distance
wave r_wave= $rwave
//print num
string cdf1=getdatafolder(1)
findvalue /V=(num) /T=0.1 omw 
variable n1= V_value
wave gd1=root:data:isotropy
variable d1=dimsize(gd1,0)
variable d2=dimsize(gd1,1)
string wname ; sprintf wname,"isotropy_%s",sname
string nname= ReplaceString(".",wname, "_")
//print nname;
setdatafolder root:
newdatafolder /O /S Output

make /D /n=(d1) /O  $nname=gd1[p][n1]
make /D /n=(d1) /O  g_distance=r_wave[p]

string /G root:Params:Alpha:s_save_name =nname
//print d1, n1
SetDataFolder cdf1
printf "Anisotropy & isotropy wave generated for %g nm.\r",num
end
//*********************************************************************************************
//*********************************************************************************************
function Gamma_save_itx(s2f):buttoncontrol
string s2f
svar sname = root:Params:Alpha:s_save_name
string dname = "g_distance"
string nname= ReplaceString(".",sname, "_")
string a,b
sprintf a,"root:Output:'%s'",sname
sprintf b,"root:Output:%s",dname
print a,b
wave a1=$a ; wave b1=$b ;
//print removeending (sname,".")

Save/T a1,b1 as nname+".itx"
end
//*********************************************************************************************
//*********************************************************************************************
function Arrange_gamma_square(s2g):buttoncontrol
string s2g
string s1,s2
string cdf1=getdatafolder(1)

killwaves /Z dest
string sel
sel=wavelist("gamma_*_ResultSet",";","")
print sel
concatenate /DL  /np sel,gammaSq_AllData
duplicate gammaSq_AllData,gammaSq_AD_T
wave gw=gammaSq_AD_T
matrixtranspose gw
end
//*********************************************************************************************
//*********************************************************************************************

function generate_gamma_nm_waves_ALL(s3c):buttoncontrol
string s3c
string cdf=getdatafolder(1)
setdatafolder root:data

svar sname =root:Params:alpha:s_Omega_nm 
wave omw=root:data:omega_nm
svar rwave=root:Params:alpha:g_distance
wave r_wave= $rwave
//---------------------------------------------------------------------

variable p1,p2
p1=dimsize(omw,0) ; p2=dimsize(omw,1) 
variable i1
variable num
for (i1=0 ; i1<p1 ;i1=i1+1)
num=omw[i1]
print i1,num
findvalue /V=(num) /T=0.1 omw 
variable n1= V_value
wave gd1=root:data:gamma_data
variable d1=dimsize(gd1,0)
variable d2=dimsize(gd1,1)

variable lambda=num
string wname ; sprintf wname,"gamma_%g",lambda
string nname= ReplaceString(".",wname, "_")

setdatafolder root:
newdatafolder /O /S Output

make /D /n=(d1) /O  $nname=gd1[p][n1]
make /D /n=(d1) /O  g_distance=r_wave[p]

string /G root:Params:Alpha:s_save_name =nname

SetDataFolder cdf1
endfor


//--------------------------------------------------------
for (i1=0 ; i1<p1 ;i1=i1+1)
num=omw[i1]
print i1,num
findvalue /V=(num) /T=0.1 omw 
 n1= V_value
wave gd1=root:data:isotropy
d1=dimsize(gd1,0)
d2=dimsize(gd1,1)

lambda=num
sprintf wname,"isotropy_%g",lambda
nname= ReplaceString(".",wname, "_")

setdatafolder root:
newdatafolder /O /S Output

make /D /n=(d1) /O  $nname=gd1[p][n1]
make /D /n=(d1) /O  g_distance=r_wave[p]

string /G root:Params:Alpha:s_save_name =nname

SetDataFolder cdf1

endfor
//---------------------------------

end

//*********************************************************************************************
//	function to scale the matrix elements to other wavelengths.
//	function generates a new matrix of elements which are scaled to new wavelengths (defined by new_wavelength 1D wave).
//	scaling is done by using spline interpolation

function scale_MEmatrix_to_wavelength(st01):buttoncontrol
string st01
// parameters :
//	ME_matrix	=	2D wave of matrix elements ( transition polarizability)
//	old_wavelength	=	1D wave of the wavelength at which the matrix elements are available
//	new_wavelength	=	1D wave of wavelength to   which the matrix elements will be scaled.

svar MEwName	=	root:Params:Alpha:scalingME2D
svar OmegawName	=	root:Params:Alpha:scalingOmegaOriginal
svar OmegaFinalwName	=	root:Params:Alpha:scalingOmegaFinal

wave ME_matrix = $MEwName

wave oldnm= $OmegawName	// should be in increasing order
wave newnm= $OmegaFinalwName	// should be in increasing order

variable x0,x1,x2
x0=dimsize(ME_matrix,0)
x1=dimsize(newnm,0)
x2=dimsize(ME_matrix,1)

string newname
sprintf newname,"%s_scaled",nameofwave(ME_matrix)
print newname, x0, x1,x2

make /o /d /n=( (x1+1) , x2) $newname=0
wave scaledME=$newname

variable i,j, nm, value
// varibale i goes over cols
// varibale j goes over individual nm ( new_wavelength)

// loop over cols, subloop over new nm
for (i=0 ; i<x2 ; i=i+1)

	make /o /d /n=(x0) input_wave=ME_matrix[p][i]
	wave input_wave=input_wave
	deletepoints /M=0 (x0-1),1,input_wave	//	one point is removed, which is static polarizability

	//	CustomSplineU_auto(source,xscale)
	CustomSplineU_auto(input_wave , oldnm)

	string array1=nameofwave(input_wave)+"_a";
	wave array=$array1
	//	print nameofwave(array)

	for (j=0 ; j<x1 ; j=j+1)
		nm=newnm[j]
		value=findY2(nm, array)
		scaledME[j][i]=value
		//	print i,j
	endfor
	
endfor

// assign the last row as the static polarizability from the original file, where it is the last row.
scaledME[(x1)][]=ME_matrix[x0-1][q]
//	print x0,x1,x2

killwaves /Z input_wave, input_wave_a
print "Last row is for static. Transpose scaled matrix to arrange MEs vertically.\r"
printf "Done.\r"
end
//--------------------------------------------------------------------------------------------

//*********************************************************************************************
//	function to scale the matrix elements to other wavelengths.
//	function generates a new matrix of elements which are scaled to new wavelengths (defined by new_wavelength 1D wave).
//	scaling is done by using spline interpolation

function scale_MEmatrix_to_OneWavelength(st01):buttoncontrol
string st01
// parameters :
//	ME_matrix	=	2D wave of matrix elements ( transition polarizability)
//	old_wavelength	=	1D wave of the wavelength at which the matrix elements are available
//	new_wavelength	=	1D wave of wavelength to   which the matrix elements will be scaled.

svar MEwName_g	=	root:Params:Alpha:scalingME2D_g	// for gamma
svar MEwName_i	=	root:Params:Alpha:scalingME2D_i	// for isotropy
svar OmegawName	=	root:Params:Alpha:scalingOmegaOriginal
nvar OmegaFinal_nm	=	root:Params:Alpha:scalingOmega_wavelength 	// wavelength

wave ME_matrix = $MEwName_g
wave oldnm = $OmegawName	// should be in increasing order

string cdf
string newname
string path
string list
string a,b
string number
string regexp
string array1
variable x0,x1
variable i,j, nm, value

cdf=getdatafolder(1)
x0=dimsize(ME_matrix,0)
x1=dimsize(ME_matrix,1)

sprintf newname,"ME_gamma_%7.6f", OmegaFinal_nm
print cdf
print newname, x0, x1
make /o /d /n=( 1 , x1) $newname=0
wave scaledME=$newname

printf "Wavelength asked :%7.6f\r",OmegaFinal_nm

// -------------------------------------------------------------
//     Gamma matrix elements

// loop over cols 
for (i=0 ; i < x1 ; i=i+1)

	make /o /d /n=(x0) input_wave=ME_matrix[p][i]
	wave input_wave=input_wave
	deletepoints /M=0 (x0-1),1,input_wave	//	one point is removed, which is static polarizability

	//	CustomSplineU_auto(source,xscale)
	CustomSplineU_auto(input_wave , oldnm)

	array1=nameofwave(input_wave)+"_a";
	wave array=$array1
	//	print nameofwave(array)

	value=findY2( OmegaFinal_nm , array)
	scaledME[0][i]=value
endfor

// Code to assign the first two cols as J indices ---
// transpose
matrixop /o scaledME = scaledME^t

// move to directory above the current directory
path="::"
setdatafolder $path
// print getdatafolder(1)
 
number=num2str(oldnm[0])
// print number
regexp = "([[:digit:]]+).([[:digit:]]+)"
splitstring /E=(regexp) number,a,b
print a,b
sprintf regexp,"gamma_%s*_ResultSet",a

// generate the new match string for wavelist
list=wavelist (regexp,"","")

wave ResultSetWave=$list

InsertPoints/M=1 0,2, scaledME

// assign the first two cols for J indices
scaledME[][0]=ResultSetWave[p][0]
scaledME[][1]=ResultSetWave[p][1]
//------------------------------------------------------------
setdatafolder cdf
//---------------------------------------------------------------------------
//     Isotropy matrix elements

wave ME_matrix = $MEwName_i

x0=dimsize(ME_matrix,0)
x1=dimsize(ME_matrix,1)

sprintf newname,"ME_iso_%7.6f", OmegaFinal_nm
print newname, x0, x1
make /o /d /n=( 1 , x1) $newname=0
wave scaledME=$newname


// loop over cols 
for (i=0 ; i < x1 ; i=i+1)

	make /o /d /n=(x0) input_wave=ME_matrix[p][i]
	wave input_wave=input_wave
	deletepoints /M=0 (x0-1),1,input_wave	//	one point is removed, which is static polarizability

	//	CustomSplineU_auto(source,xscale)
	CustomSplineU_auto(input_wave , oldnm)

	array1=nameofwave(input_wave)+"_a";
	wave array=$array1
	//	print nameofwave(array)

	value=findY2( OmegaFinal_nm , array)
	scaledME[0][i]=value
endfor


// Code to assign the first two cols as J indices ---
// transpose
matrixop /o scaledME = scaledME^t

// move to directory above the current directory
path="::"
setdatafolder $path
// print getdatafolder(1)
 
number=num2str(oldnm[0])
// print number
regexp = "([[:digit:]]+).([[:digit:]]+)"
splitstring /E=(regexp) number,a,b

sprintf regexp,"isotropy_%s*_ResultSet",a

// generate the new match string for wavelist
list=wavelist (regexp,"","")

wave ResultSetWave=$list

InsertPoints/M=1 0,2, scaledME

// assign the first two cols for J indices
scaledME[][0]=ResultSetWave[p][0]
scaledME[][1]=ResultSetWave[p][1]
//------------------------------------------------------------
setdatafolder cdf
//---------------------------------------------------------------------------
killwaves /Z input_wave, input_wave_a
printf "Done.\r"
end
//--------------------------------------------------------------------------------------------

//*********************************************************************************************
//*********************************************************************************************

Function ME2D_wave_Selection(ME2D_sel1,pop_gw1,popStr_gw1) : PopupMenuControl	
String ME2D_sel1
Variable pop_gw1
String popStr_gw1
string /G root:Params:Alpha:scalingME2D = popStr_gw1
End
//*********************************************************************************************
Function OmegaOriginal_wave_Selection(Omega_sel1,pop_gw1,popStr_gw1) : PopupMenuControl	
String Omega_sel1
Variable pop_gw1
String popStr_gw1
string /G root:Params:Alpha:scalingOmegaOriginal = popStr_gw1
End
//*********************************************************************************************
Function OmegaFinal_wave_Selection(OmegaF_sel1,pop_gw1,popStr_gw1) : PopupMenuControl	
String OmegaF_sel1
Variable pop_gw1
String popStr_gw1
string /G root:Params:Alpha:scalingOmegaFinal = popStr_gw1
End
//*********************************************************************************************

//*********************************************************************************************

Function ME2D_Gamma_wave_Selection(ME2D_sel1,pop_gw1,popStr_gw1) : PopupMenuControl	
String ME2D_sel1
Variable pop_gw1
String popStr_gw1
string /G root:Params:Alpha:scalingME2D_g = popStr_gw1
End
//*********************************************************************************************
Function ME2D_Iso_wave_Selection(ME2D_sel1,pop_gw1,popStr_gw1) : PopupMenuControl	
String ME2D_sel1
Variable pop_gw1
String popStr_gw1
string /G root:Params:Alpha:scalingME2D_i = popStr_gw1
End
//*********************************************************************************************

