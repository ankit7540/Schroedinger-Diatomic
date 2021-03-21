#pragma TextEncoding = "Windows-1252"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#include <All Gizmo Procedures>
// FUNCTIONS INSIDE THIS FILE..

// 	ScaledY_gen(source_x-wave, result_array)
// 	CustomsplineU()
//	GLNumerical_integrate(a,b,n,array)
//***************************************************************************************************************
function waveConcatenate()
//setdatafolder root:
string sel
print "concatenating waves with name wave* "
sel=wavelist("wave*",";","")
concatenate sel,matrix1
end
//***************************************************************************************************************
// concatenates wave with name  "wave" to make matrix1

function MatReady()
waveConcatenate()
wave matrix1=matrix1
duplicate matrix1, matrix2
duplicate matrix2, matrix3
DeletePoints 0,1, matrix3
matrixtranspose matrix3
end
//***************************************************************************************************************
function MatReady_1()
waveConcatenate()
wave matrix1=matrix1
duplicate matrix1, matrix2
duplicate matrix2, matrix3
matrixtranspose matrix3
end

//***************************************************************************************************************

// WaveSum(w)
// Returns the sum of the entire wave, just like Igor’s sum function.
Function WaveSum(w)
Wave w
Variable i, n=numpnts(w), total=0
for(i=0;i<n;i+=1)
total += w[i]
endfor
return total
End
//***************************************************************************************************************
//5 point asymmetric derivative on the first 4 points of the wave
// Used by customSplineU procedure
function gen_coef_wr1(inp_x_wave, inp_y_wave)
wave inp_x_wave,inp_y_wave
variable n;
Killwaves /Z A,C, D_coefR1
make /O /D /n=(4,4) A=0
make  /O /D /n=4 C=0

for (n=0 ;n<4 ; n=n+1)
	C[n] = inp_y_wave[n]
endfor 

make /O /D /n=4 D_coefR1
variable i,j

for(i=0 ; i<4 ; i=i+1)
	A[i][]=(inp_x_wave[i])
		for (j=0 ; j<4 ; j=j+1)
			A[i][j]=(A[i][j]^j)/ factorial (j);
		endfor
endfor 

matrixop /O  D_coefR1=inv(A) x C 

wave D_coefR1 = D_coefR1
end


//***************************************************************************************************************


function auto_curveFit(n,s, inp_xwave, inp_ywave)
variable n,s
wave inp_xwave, inp_ywave

//-----------------------------------------
killwaves /Z pcoefs, V_chisq,W_sigma
//-----------------------------------------
make  /O /D /n=(n) yd=0
make  /O /D /n=(n) xd=0

variable j,k
j=dimsize( inp_ywave,0)
//print j
variable i1
if (s==0)
for (i1=0 ;i1<(n) ; i1=i1+1)
	yd[i1] = inp_ywave[i1]
	xd[i1] = inp_xwave[i1]
endfor 
//print xd
//print yd
make /o /D /n=4 pcoefs=0
wave coefs=pcoefs
CurveFit /B=2 /NTHR=0 /Q poly 4, kwCWave=coefs, yd /X=xd /D /R /A=0 

endif

if (s==1)

for (i1=0 ;i1<n ; i1=i1+1)
//print j, i1
	yd[i1] = inp_ywave[ j-1 -i1]
	xd[i1] = inp_xwave[j-1 -i1]
endfor 
//print xd
//print yd
make /o /D /n=4 pcoefs=0
wave coefs=pcoefs

CurveFit  /B=2  /NTHR=0 /Q poly 4, kwCWave=coefs, yd /X=xd /D /R /A=0 


endif

killwaves /Z V_chisq,W_sigma,W_ParamConfidenceInterval
end
//***************************************************************************************************************

function CustomSplineU()
	string source,xscale
	variable J0,x1,y1,index
	index=0;	J0=4; 
		Prompt source,  "source wave  ", popup WaveList("*", ";", "")
		Prompt xscale,  "x scale wave  ", popup WaveList("*", ";", "")
		//Prompt  J0,  "Number of splines	:  "
		Prompt index,"wave in the matrix to be selected"
				DoPrompt "input waves  ", source,xscale,index
 				if( V_Flag )
     				return 0         		// user canceled
				endif
	wave  inputwave = $source ; 
	wave x_scale = $xscale	;
	variable t0= ticks
	x1=dimsize(inputwave,0); y1=dimsize(inputwave,1) ; printf "Dimension : %g x %g  | Splines asked : %g \r" x1, y1,x1-1 ; 
	J0=x1-1
	if (index >x1)
	index=0
	endif
	
	//clear previous versions of matrices etc...
	KillWaves /Z P1,Q1,R1,P2,P3,mat1,xw1,soln,product,inp,inp_x,value_x,value_y,Yval
	//
	make /o /D /n=(4*J0,4*J0) P1=0
	make /o /D /n=(4*J0) R1=0
	
	wave P1=P1
	wave Q1=Q1
	wave R1=R1
	
	//PRE-CALCULATION=====================================================================
	variable div=x1/J0; div=round(div)

	// Setting the intervals in the x wave 
	variable n=0; variable k=1;
	for (n=0;n<J0;n=n+1)
		k=div*n
	endfor
	
	make /o /D /n=(n+1) value_y ; wave value_y = value_y
	make /o /D  /n=(n+1) value_x ; wave value_x = value_x
	make /o /D /n=(n+1) index_x1 ; wave index_x1 = index_x1
	
	variable k2
	for  (  k2=0;k2 < (n+1) ; k2=k2+1)                                     
		k=div*k2
		 value_y [(k2)]=inputwave[(k)]
		 index_x1[(k2)]=k
		 value_x[k2]=x_scale[(k)]
	endfor
	//********************************************************
	variable i2=0, i=0
	//********************************************************
	//MATRIX VALUES
	// 1. PART A
	variable e2=-1
	variable j=1,counter=-1,dx=0
	for (j=1; j <= (J0) ; j=j+1 )
	
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		variable e1=0, val=0
		for (e1=0 ; e1 < 4 ; e1=e1+1)
	
			val = dx^e1
			e2=e2+1
			P1[(counter)][(e2)] = val
		//	printf "e2 : %g ; val : %g ; e1: %g\r",e2,val,e1
		endfor
	endfor
	variable fd=e2
	//**************************************************************************************
	e2=-1;
	//2. PART B // CONDITION 1
	 j=1;dx=0
	for (j=1; j <(J0) ; j=j+1 )

		//printf "counter B  :%g\r",counter
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		e1=0; val=0
		for (e1=0 ; e1 < 5 ; e1=e1+1)
	
			val = dx^e1
				if (e1==4)
				val = -1
				endif
				
			e2=e2+1
			//print counter
				if (e1==0 && j>1)
					e2=e2-1
				endif
			P1[(counter)][(e2)] = val
			//e2=e2-1
			//printf "val : %g : %g\r",e1,val
		endfor
	endfor
	//printf "B rows :%g\r",j
	//printf "counter :%g\r",counter
	//**************************************************************************************
	
	e2=-1;
	//2. PART C // CONDITION 2 : FIRST DERIVATIVE
	 j=1;dx=0
	for (j=1; j <(J0) ; j=j+1 )

		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		e1=0; val=0
		for (e1=0 ; e1 < 6 ; e1=e1+1)
			val = e1*dx^(e1-1)
			if (e1==4)
			val=0
			elseif (e1==5)
			val=-1
			endif
			e2=e2+1
			if (e1==0 && j>1)
			e2=e2-2
			endif
			//print counter
			P1[(counter)][(e2)] = val
		//	printf "val : %g : %g\r",e1,val
		endfor
	endfor

	//**************************************************************************************
	
	e2=-1;
	//2. PART D // CONDITION 3 : SECOND DERIVATIVE
	 j=1;dx=0
	for (j=1; j <(J0) ; j=j+1 )
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		e1=0; val=0
		for (e1=0 ; e1 < 7 ; e1=e1+1)
			
			val = e1*dx^(e1-2)
			if (e1==1) // second term
				val=0
			endif
			
			if (e1==3) // fourth term
 				val=val*2
			endif
			if (e1==4 || e1 ==5 )
				val=0
			endif
			if(e1==6)
				val=-2
			endif			
			e2=e2+1    // moves the column number to next 
			if (e1==0 && j>1)
				e2=e2-3
			endif

			P1[(counter)][(e2)] = val
		//	printf "counter : %g | val : %g : %g\r",counter, e1,val
		endfor
	endfor
	//printf "D rows :%g\r",j
	//printf "counter :%g\r",counter
	//**************************************************************************************

	e2=-1;
	//2. PART E // ADDITIONAL THREE ROWS TO BE FILLED UP HERE
	
	j=0;dx=0
	for (j=0; j < 2 ; j=j+1 )
		counter =counter+1
		e1=0; val=0
		for (e1=0 ; e1 < 4 ; e1=e1+1)
		
		val=e1
	//	print e1, val
		
			if (j==0)
				if (e1==1)
					val=1
				else
					val=0
				endif
			endif
			
			e2=e2+1
			
			if (e1==0 && j==1)
				e2=(4*J0)  -  4    // second last row , first derivative assignment in matrix
			endif
			
			P1[(counter)][(e2)] = val
		//	printf "e1: %g | val: %g\r",e1,val
		endfor
	endfor

	// last row
	e2=-1
	 j=1;dx=0
	for (j=1; j <2; j=j+1 )
	counter =counter+1
		dx=0;
		e1=0; val=0
		for (e1=0 ; e1 < 4 ; e1=e1+1)
			
			if (j==1 && e1==0)
				val=1.0000000
			else
				val=00
			endif
			e2=e2+1
			P1[(counter)][(e2)] = val
		//	print counter
		//	printf "val : %g : %g\r",e1,val
		endfor
	endfor
	
//	abort
//===============================================================================
// ASSIGN Y VALUES:::

for (j=1; j <=J0; j=j+1 )
	val=value_y[(j)]
	R1[j-1]=val
endfor

R1[(4*J0)-1] = value_y[(0)] 		// last point in R1
//	R1[(4*J0)-2] = -5;  //  for testing purpose
variable s1=0
// print /D s1, (4*J0)-3
R1[(4*J0)-3] = s1

//VALUE OF  DERIVATIVE

auto_curveFit(4,0, value_x,value_y) //first point
wave D_coefR1 =pcoefs
s1=D_coefR1[1]+2*D_coefR1[2]* (value_x[0])   + 3*D_coefR1[3]*(value_x[0]^2)
R1[(4*J0)-3] = s1

auto_curveFit(4,1, value_x,value_y)  //last point
wave D_coefR1 =pcoefs
s1=D_coefR1[1]+2*D_coefR1[2]* (value_x[x1-1])   + 3*D_coefR1[3]*(value_x[x1-1]^2)
R1[(4*J0)-2] = s1
// print /D s1, (4*J0)-3,  (value_x[x1-1])
// 	abort
//------------------------------------------------
//MAKE X WAVE FOR TESTING 
// print "making Test wave "
make /o /D /n=(J0+1,4*J0) xw1=0
	
	wave xw1=xw1
	e2=-1
	j=1;counter=0;dx=0
	for (j=0; j <= (J0) ; j=j+1 )
		counter = j
		
		if (j==0)
		dx=0
		elseif (j>0)
		dx=value_x[(j)] - value_x[(j-1)]  
		endif
		//printf "dx : %g\r",dx
		e1=0; val=0
		if (counter==1)
			e2=0
			endif
		for (e1=0 ; e1 < 4 ; e1=e1+1)
			
			val = dx^e1
			e2=e2+1
			if (e1==0 && counter==1)
				e2=0
			endif
			xw1[(counter)][(e2)] = val
			endfor
	endfor
//------------------------------------------------

//MATRIX SOLVE

//Check IDENTITY MATRIX OBTAINED OR NOT
// matrixop /o mat1=P2 x P1

// ACTUAL SOLUTION
MatrixOp /O product = Inv(P1) x R1
duplicate /o product, soln
killwaves product //remove product wave !
//TESTING
//print "Running test"
matrixop /o NewY=xw1  x soln
//duplicate /o NewY, testy

//**************************************************************************************
//ERROR CHECK
//print J0
make /o /D /n=(J0+1) diff = NewY[p] - value_y[p]
Display /K=1 value_y,NewY vs value_x
ModifyGraph mode=3,marker(value_y)=17
ModifyGraph marker(value_y)=23,rgb(value_y)=(16384,16384,65280)

AppendToGraph diff vs value_x
RemoveFromGraph diff
AppendToGraph/L=to11 diff vs value_x
ModifyGraph axisEnab(left)={0,0.65},axisEnab(to11)={0.7,1},freePos(to11)=0
ModifyGraph mode(diff)=1
ModifyGraph mode=3,marker(diff)=26,rgb(diff)=(21760,21760,21760)
ModifyGraph marker(diff)=18,msize(diff)=1.5 ; doupdate
ModifyGraph zero(to11)=3; doupdate
variable error = WaveSum (diff)
 //**************************************************************************************
 //MAKE RESULTING ARRAY
//printf "J0 : %g\r", J0
make /O /D /n=(J0,6) $source+"_a"
wave res = $source+"_a"

 i=0
for (i=0; i<J0 ;  i=i+1)
res[(i)][0]= value_x[(i)]	 	

res [(i)][1] = value_x[(i+1)]
res[i][2]=soln[4*i]
res[i][3]=soln[(4*i)+1]
res[i][4]=soln[(4*i)+2]
res[i][5]=soln[(4*i)+3]

endfor	
//**************************************************************************************

printf "Finished in :  %g sec  |  Sum of errors at the points : %g\r",(ticks-t0)/60,  abs (error)

 //**************************************************************************************
 variable next	,points
		Prompt next,  "enter 0 or 1   "
		//Prompt xscale,  "x scale wave  ", popup WaveList("*", ";", "")
		DoPrompt "Next step  (0 = stop here ; 1 = insert x-points ) ",next
		if( V_Flag )
     				return 0         		// user canceled
				endif
				
				
			if (next==1)
			Prompt points,  "Enter number of points per spline"
			DoPrompt " Points per spline  ", points
			if( V_Flag )
     				return 0         		// user canceled
				endif
				KillWaves /Z inp, inp_x
				make /O /D /n=(J0*(points+1), J0*4) inp
				make  /O /D  /n=(J0*(points+1)) inp_x
				
			 i=0 ; counter=-1
				for (i=1 ; i<=J0 ; i=i+1 )
				variable range=value_x[i]-value_x[i-1]
				//print i , range
				
				j=range / points
			//	printf "i = %g | range = %g | minimum step (dx) = %g \r", i,range,j
				k=0 ; div=0 ; e2=-1
					for (k=0; k<=points ; k=k+1)
						div=k*j
						//print k, div
						counter=counter+1
						//print counter
				     	 	for (e1=0 ; e1 < 4 ; e1=e1+1)
							val = div^e1
							
							e2=e2+1
							//print counter
							if (e1==0 && k<=points)
								e2=4*(i-1)
							endif
							inp[(counter)][(e2)] = val
							inp_x [(counter)] = value_x[(i-1)] + (div)
							
							//printf "val : %g : %g\r",e1,val		
						endfor
				 
					endfor
				
				endfor
			matrixop /o Yval=inp  x soln
			
			Display /K=1 Yval vs inp_x
			AppendToGraph inputwave[][index] vs x_scale
			ModifyGraph mode=3,marker(Yval)=19,rgb(Yval)=(16384,16384,65280);DelayUpdate
			ModifyGraph marker(matrix3)=17 ; ModifyGraph mode=4 ; ModifyGraph msize=1 ; 
			Legend/C/N=text0/A=MCdoupdate; doupdate
			
			print "second step finished !"	
			endif
			
			killwaves /Z Q1,value_x, value_y, value_x2,D_coef
			killwaves /Z value_y2, P2,soln,inp,inp_x,inp_x1,xw1,index_x1,index_x12
			killwaves /Z M_inverse,A,C,D_coefR1
			Killwaves /Z P1,R1
			
			//kill waves used for the fit
			killwaves /z yd,xd,fit_yd,Res_yd
			
end


//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************
//***************************************************************************************************************

 
function findY(xn)
 variable xn
 wave Wname ; 
 wave Wname =root:fitting:Rych:mp0000_Rych_a  ; // matrix array which has coeffiecients of the polynomial
 
 variable i=0,x1=0,y1=0,val1=0,val2=0, row

// x1=dimsize(Wname,0); y1=dimsize(Wname,1) ; print x1, y1, xn

if (xn == Wname[0][0] )
print "found in 0" ; row=0
endif

  for (i=0; i<x1 ; i=i+1)
  val1=0 ; val2=0
 // FIND THE SEGMENT 
val1= Wname[i][0]		; 	val2= Wname[i][1] ; //print val1,val2

 if   (((xn < val2)))
	if ((xn > val1 ))
		printf "found in A : %g\r",i ; 
		row=i;
	endif
endif

endfor 

// print row
variable dx = xn - Wname[(row)][0] ; //print  Wname[(row)][0]
variable yval = Wname[(row)][2] + Wname[(row)][3]*dx + Wname[(row)][4]*dx^2 +  Wname[(row)][5]*dx^3 
 
if ( numtype(row) != 2 )
// print dx, yval
return yval // GIVING VALUE TO  THE INTEGRAL FUNCTION
else
print "not found !"
endif
 
 end 
 
//*************************************************************************************
 //**************************************************************************************
 // DEFINE THE MATRIX RESULT ARRAY SIMULTANEOUSLY
function findY2(xn,array)
variable xn ; wave array
//print array
wave Wname = array ; // matrix array which has coefficients of the polynomial

 
 variable i=0,x1=0,val1=0,val2=0, row

x1=dimsize(Wname,0);  // print x1, y1, xn

if (xn == Wname[0][0] )
// print "found in 0" ; row=0
endif

  for (i=0; i<x1 ; i=i+1)
  val1=0 ; val2=0
 // FIND THE SEGMENT 
val1= Wname[i][0]		; 	val2= Wname[i][1] ; //print val1,val2
//print "h"
 if   (((xn <= val2)))
	if ((xn >= val1 ))
	//print "hi"
	//	printf "found in A : %g\r",i ; 
		row=i;
	endif
endif

endfor 

// print row
variable dx = xn - Wname[(row)][0] ; //print  Wname[(row)][0]
variable yval = Wname[(row)][2] + Wname[(row)][3]*(dx) + Wname[(row)][4]*(dx^2) +  Wname[(row)][5]*(dx^3 )
 
if ( numtype(row) != 2 )
// print dx, yval
return yval // GIVING VALUE TO  THE INTEGRAL FUNCTION
else
print "not found !"
endif
 
 end 
 
//*************************************************************************************

function t1(x1)
variable x1
variable y1= sin(5*x1)
return y1
end
//*************************************************************************************	
function GCNumerical_integrate(a,b,n)
variable a,b,n


wave coef= abscissa ;//  wave weight= weight  ;


variable j ;  variable dx,xm,xr,ss ;
  ss=0
 xm=0.5*(b+a); xr = 0.5*(b-a);
 for (j=0 ; j<n ; j=j+1)
 //print coef[j][(n-2)]
 dx=xr*coef[j][(n-2)] ; // print dx ; 
 ss=ss+( 1 *( t1(xm+dx)  + t1(xm-dx)))
 //print t1(xm+dx) ;
 endfor
//print ss 
ss=xr*ss*(1/n)

print ss;
 
 end
  
//*************************************************************************************

Threadsafe Function/wave gauss_legendre(N)
	///routine from Numerical Recipes to calculate weights and abscissae for 
	// Gauss-Legendre quadrature.

 
	variable N
 
	Variable x1, x2
	variable m, jj, ii
	variable z1, z, xm, xl, pp, p3, p2, p1
	Variable eps = 3e-11
 
	make /D  /O /free/n=(N, 2)/d GLpoints
 
	x1 = -1
	x2 = 1
 
	m = (N+1)/2	
	xm = 0.5 * (x2 + x1)
	xl = 0.5 * (x2 - x1)
 
	for (ii = 1; ii <= m; ii += 1) 
		z = cos(pi * (ii - 0.25)/(N + 0.5))    					   //approximation of root stored in z
		do 
			p1 = 1.0
			p2 = 0.0
			for (jj = 1;jj <= N;jj += 1) 
				p3 = p2
				p2 = p1
				p1 = ((2.0 * jj - 1.0) * z * p2 - (jj - 1.0) * p3) / jj    //polynomial eval at z
			endfor
			pp = N * (z * p1 - p2) / (z * z -1.0)   // derivative
			z1 = z
			z = z1 - p1 / pp
		while (abs(z - z1) > EPS)
		GLpoints[ii - 1][0] = xm - xl * z					   //abscissa
		GLpoints[N  - ii][0] = xm + xl * z					   //abscissa symmetric counterpart
		GLpoints[ii - 1][1] = 2.0 * xl / ((1.0 - z * z) * pp * pp) //weight
		GLpoints[N  - ii][1] = GLpoints[ii - 1][1] 		          // symmetric counterpart of weight
	Endfor
 	make /D /O  /n=(N, 2)/d GLpoints1=GLpoints                      //Copy stored in current folder for use.
	return GLpoints
	
End	
//*************************************************************************************	
function GLNumerical_integrate(a,b,n,array)
variable a,b,n
wave array
gauss_legendre(n) // generating Gauss-Legendre Abscissa and weights
wave coef= GLpoints1;  wave weight= GLpoints1 ;

variable j ;  variable dx,xm,xr,ss ;
ss=0    //sum of area in segment

xm=0.5*(b+a); xr = 0.5*(b-a);

for (j=0 ; j<n ; j=j+1)
	
	 dx=xr*coef[j][0] ;
	 ss=ss+  weight[j][1] *(( findY2 (xm+dx,array) )/2 + (( findY2(xm-dx,array))/2))
	
endfor

ss=xr*ss
// print /D ss;

variable /D /G intResult=ss
return ss
end
  
//*************************************************************************************	

// NUMERICAL INTEGRATION TESTING FUNCTIONS :::::::::::::::::
// POLYNOMIALS
function pol1(x1)
variable x1
variable y=6*x1^3 + 5*x1^2 + 4*x1 + 3
return y
end 

function pol2(x1)
variable x1
variable y=x1^2 * (x1-45)*(x1-90)
return y
end 

function pol3(x1)
variable x1
variable y=(x1-3)*(x1+3)
return y
end 



//*************************************************************************************	

// Generate y wave from a given x-wave (source wave) and the spline fit  array made previously.
// source = xaxis wave
// ra= result array from spline 
//using the customspline

function  ScaledY_gen( newX , rarray )
	wave newX			//	new xaxis wave, 1D
	wave rarray			//	coefs array wave, 2D

	wave res_array = rarray
	variable x2
	x2=dimsize( newX, 0)
	print x2
	make /D  /O /n=(x2) scaled_y
	wave scy= scaled_y

	variable i=0, y2
	
	for (i=0 ;i<x2;  i=i+1 )
		y2 = newX [i]
		scy[i]=findY2(y2, res_array)
		//print  y2, findY2(y2, res_array) 
	endfor

end

//*************************************************************************************	
function CustomSplineU_auto(source,xscale)
	wave source,xscale
	DFREF cdf1 = GetDataFolderDFR()
	variable J0,x1,y1,adp,index
	setdatafolder root:Params:Customspline:
	
	string name = NameOfWave(source) , name_rwave =  NameOfWave(xscale)
	wave  inputwave = source ; 
	wave x_scale = xscale	;
	
	x1=dimsize(inputwave,0); y1=dimsize(inputwave,1) ; 
	J0=(x1-1);  // set number of spline as one less than number of points
	//clear previous versions of matrices etc...
	KillWaves /Z P1,Q1,R1,P2,P3,mat1,xw1,soln,product,inp,inp_x,value_x,value_y,Yval
	//
	make /o /D /n=(4*J0,4*J0) P1=0
	make /o /D /n=(4*J0) R1=0
	
	wave P1=P1
	wave R1=R1
	
	//PRE-CALCULATION=====================================================================
	variable div=x1/J0; div=round(div)
	// Setting the intervals in the x wave 
	variable n=0; variable k=1;
	for (n=0;n<J0;n=n+1)
		k=div*n
	endfor
	make /o /D /n=(n+1) value_y ; wave value_y = value_y
	make /o /D  /n=(n+1) value_x ; wave value_x = value_x
	make /o /D /n=(n+1) index_x1 ; wave index_x1 = index_x1
	
	variable k2
	for  (  k2=0;k2 < (n+1) ; k2=k2+1)                                     
		k=div*k2
		 value_y [(k2)]=inputwave[(k)]
		 index_x1[(k2)]=k
		 value_x[k2]=x_scale[(k)]
	endfor
	
	//******************************************************
	variable i2=0, i=0
	//******************************************************
	//MATRIX VALUES
	// 1. PART A
	variable e2=-1
	variable j=1,counter=-1,dx=0
	for (j=1; j <= (J0) ; j=j+1 )
	
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		variable e1=0, val=0
		for (e1=0 ; e1 < 4 ; e1=e1+1)
	
			val = dx^e1
			e2=e2+1
			P1[(counter)][(e2)] = val
		//	printf "e2 : %g ; val : %g ; e1: %g\r",e2,val,e1
			
		endfor
	endfor
	variable fd=e2
	//**************************************************************************************
	e2=-1;
	//2. PART B // CONDITION 1
	 j=1;dx=0
	for (j=1; j <(J0) ; j=j+1 )
	
	//print j
		//printf "counter B  :%g\r",counter
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		e1=0; val=0
		for (e1=0 ; e1 < 5 ; e1=e1+1)
			val = dx^e1
				if (e1==4)
				val = -1
				endif
				
			e2=e2+1

				if (e1==0 && j>1)
					e2=e2-1
				endif
			P1[(counter)][(e2)] = val
			//e2=e2-1
			//printf "val : %g : %g\r",e1,val
		endfor
	endfor

	//**************************************************************************************
	
	e2=-1;
	//2. PART C // CONDITION 2 : FIRST DERIVATIVE
	 j=1;dx=0
	for (j=1; j <(J0) ; j=j+1 )
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		e1=0; val=0
		for (e1=0 ; e1 < 6 ; e1=e1+1)	
			val = e1*dx^(e1-1)
			if (e1==4)
				val=0
			elseif (e1==5)
				val=-1
			endif
			e2=e2+1
			if (e1==0 && j>1)
				e2=e2-2
			endif
			//print counter
			P1[(counter)][(e2)] = val
		//	printf "val : %g : %g\r",e1,val
		endfor
	endfor

	//**************************************************************************************
	e2=-1;
	//2. PART D // CONDITION 3 : SECOND DERIVATIVE
	 j=1;dx=0
	for (j=1; j <(J0) ; j=j+1 )
		counter =counter+1
		dx=value_x[(j)] - value_x[(j-1)]
		//printf "dx : %g\r",dx
		e1=0; val=0
		for (e1=0 ; e1 < 7 ; e1=e1+1)
			
			val = e1*dx^(e1-2)
			if (e1==1) // second term
				val=0
			endif
			if (e1==3) // fourth term
 				val=val*2
			endif
			if (e1==4 || e1 ==5 )
				val=0
			endif
			if(e1==6)
				val=-2
			endif			
			e2=e2+1    // moves the column number to next 
			if (e1==0 && j>1)
				e2=e2-3
			endif

			P1[(counter)][(e2)] = val
		//	printf "counter : %g | val : %g : %g\r",counter, e1,val
		endfor
	endfor
	//printf "D rows :%g\r",j
	//printf "counter :%g\r",counter
	//**************************************************************************************
	e2=-1;
	//2. PART E // ADDITIONAL THREE ROWS TO BE FILLED UP HERE
	
	j=0;dx=0
	for (j=0; j < 2 ; j=j+1 )
		counter =counter+1
		e1=0; val=0
		for (e1=0 ; e1 < 4 ; e1=e1+1)
		
		val=e1
	//	print e1, val
		
			if (j==0)
				if (e1==1)
					val=1
				else
					val=0
				endif
			endif
			
			e2=e2+1
			
			if (e1==0 && j==1)
				e2=(4*J0)  -  4    // second last row , first derivative assignment in matrix
			endif
			P1[(counter)][(e2)] = val
		//	printf "e1: %g | val: %g\r",e1,val
		endfor
	endfor

	// last row
	e2=-1
	 j=1;dx=0
	for (j=1; j <2; j=j+1 )
	counter =counter+1
		dx=0;
		e1=0; val=0
		for (e1=0 ; e1 < 4 ; e1=e1+1)
			if (j==1 && e1==0)
				val=1.0000000000000000000000000
			else
				val=00
			endif
			e2=e2+1
			P1[(counter)][(e2)] = val
		//	print counter
		//	printf "val : %g : %g\r",e1,val
		endfor
	endfor
	
//	abort
//===============================================================================
// ASSIGN Y VALUES:::

for (j=1; j <=J0; j=j+1 )
	val=value_y[(j)]
	R1[j-1]=val
endfor

R1[(4*J0)-1] = value_y[(0)] 		// last point in R1

//	R1[(4*J0)-2] = -5;  //  for testing purpose
variable s1=0

// print /D s1, (4*J0)-3
R1[(4*J0)-3] = s1

//VALUE OF  DERIVATIVEs
auto_curveFit(4,0, value_x,value_y)	//FIRST POINT
wave D_coefR1 =pcoefs
s1=D_coefR1[1]+2*D_coefR1[2]* (value_x[0])   + 3*D_coefR1[3]*(value_x[0]^2)
R1[(4*J0)-3] = s1

auto_curveFit(4,1, value_x,value_y)	//LAST POINT
wave D_coefR1 =pcoefs
s1=D_coefR1[1]+2*D_coefR1[2]* (value_x[x1-1])   + 3*D_coefR1[3]*(value_x[x1-1]^2)
R1[(4*J0)-2] = s1
// print /D s1, (4*J0)-3,  (value_x[x1-1])
// 	abort
//------------------------------------------------
//------------------------------------------------
//MAKE X WAVE FOR TESTING 
// print "making Test wave "
make /o /D /n=(J0+1,4*J0) xw1=0
	
	wave xw1=xw1
	e2=-1
	j=1;counter=0;dx=0
	for (j=0; j <= (J0) ; j=j+1 )
		counter = j
		
		if (j==0)
		dx=0
		elseif (j>0)
		dx=value_x[(j)] - value_x[(j-1)]  
		endif
		//printf "dx : %g\r",dx
		e1=0; val=0
		if (counter==1)
			e2=0
			endif
		for (e1=0 ; e1 < 4 ; e1=e1+1)
			
			val = dx^e1
			e2=e2+1
			if (e1==0 && counter==1)
				e2=0
			endif
			xw1[(counter)][(e2)] = val
			endfor
	endfor
//------------------------------------------------
//MATRIX SOLVE
MatrixOp /O product = Inv(P1) x R1;
wave product=product
matrixop /o NewY=xw1  x product
//**************************************************************************************
//ERROR CHECK
make /o /D /n=(J0+1) diff = NewY[p] - value_y[p]
variable esum=0
for (i=0; i<J0+1 ;  i=i+1)
	esum=esum+(diff[i]*diff[i])
endfor	

variable error = sqrt(esum)
 //**************************************************************************************
 //MAKE RESULTING ARRAY
make /O /D /n=(J0,6) $name +"_a"
wave res = $name +"_a"
 i=0
for (i=0; i<J0 ;  i=i+1)
res[(i)][0]= value_x[(i)]	 	
res [(i)][1] = value_x[(i+1)]
res[i][2]=product[4*i]
res[i][3]=product[(4*i)+1]
res[i][4]=product[(4*i)+2]
res[i][5]=product[(4*i)+3]
endfor	
//**************************************************************************************
//Clear off the generated waves for calculation 
killwaves /Z 	value_x, value_y, value_x2,D_coef
			killwaves /Z value_y2, P2,soln,inp,inp_x,inp_x1,xw1,index_x1,index_x12
			killwaves /Z M_inverse,A,C, D_coefR1
			Killwaves /Z P1,R1,product
printf "Dimension : %g x %g | Splines : %g | Sum of errors at the points : %g\r", x1, y1,J0,error ; 		
setdatafolder cdf1	
make /O /D /n=(J0,6) $name +"_a" =res
end

//**************************************************************************************
//**************************************************************************************
//**************************************************************************************
function findY2_a(xn,array)
variable xn ; wave array
//print array
wave Wname = array ; // matrix array which has coeffiecients of the polynomial

 
 variable i=0,x1=0,val1=0,val2=0, row

x1=dimsize(Wname,0);  // print x1, y1, xn

if (xn == Wname[0][0] )
// print "found in 0" ; row=0
endif

  for (i=0; i<x1 ; i=i+1)
  val1=0 ; val2=0
 // FIND THE SEGMENT 
val1= Wname[i][0]		; 	val2= Wname[i][1] ; //print val1,val2

 if   (((xn <= val2)))
	if ((xn >= val1 ))
	//	printf "found in A : %g\r",i ; 
		row=i;
	endif
endif

endfor 

// print row
variable dx = xn - Wname[(row)][0] ; //print  Wname[(row)][0]
variable yval = Wname[(row)][2] + Wname[(row)][3]*dx + Wname[(row)][4]*dx^2 +  Wname[(row)][5]*dx^3 
// printf "%7.7f\r", yval
if ( numtype(row) != 2 )
// print dx, yval
return yval // giving value to the called function 
else
print "not found !"
endif
 
 printf "%7.7f\r", yval
 end 
 
 //**************************************************************************************
//**************************************************************************************
//--------------------------------------------------------------------------------------------
function wavemaker_auto(x1,x2,step)
variable x1,x2,step
	variable range=(x2-x1)/step
	variable i1
	print x1,x2,step,range
	make  /O /d /n=(range+1) x_wave_n
	wave xwn=x_wave_n
	for (i1=0; i1<(range+1); i1=i1+1) 
		xwn[i1]=x1+i1*step
	endfor
	
	
end
//--------------------------------------------------------------------------------------------