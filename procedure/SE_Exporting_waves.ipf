
//***************************************************************************************************************

// function to save the wavefunctions in the folder as txt file to a directory in the OS file system
// Define the OS file system manually.

function export_wavefunctions(prefixString, firstNum, lastNum)
	String prefixString		// The part of the name that is common to all waves.
	Variable firstNum		// Number of the first wave in the series.
	Variable lastNum		// Number of the last wave in the series.
	
	// Example :
	//		export_wavefunctions ("v1J", 0, 10)
	
	string wave_ref
	string FileName
	string file
	Variable n
	String currentWaveName
 
	//------------------------------------------------------------------------

	
	For (n=firstNum; n<=lastNum; n+=1)
		sprintf currentWaveName, "%s%d_norm", prefixString, n	// "_norm is appended"
 		
 		//-----------------------------------------------------------
		WAVE currentWave = $(currentWaveName)
		printf "%s, %d\r" nameofwave(currentWave), dimsize(currentWave,0)   
		wave_ref=GetWavesDataFolder(currentWave,1)
		file=nameofwave(currentWave)
		sprintf FileName,"%s.txt",file
		Save /O /J  currentWave as FileName
		
		//-----------------------------------------------------------
		
	EndFor
	
	Wave rwave=r_wave 
	sprintf FileName,"r_wave.txt",file
	Save /O /J  rwave as FileName
	print "processed"
End
//***************************************************************************************************************
