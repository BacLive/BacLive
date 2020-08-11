// Set Environment and collect reconstruction info
	run("Set Measurements...", "mean min redirect=None decimal=3");
	SourceID = getImageID();
	SourceName = getTitle(); //needed for a specific step when source ID doesnt work
	setBatchMode("false");  //aim to take advantage of selective hiding of images to improve performance without weird behaviour
	zStart = 0;  //default values for stack pre-processing
	zEnd = 0;
	
// Get basic image info
	Stack.getDimensions(PXWidth,PXHeight,TotalChannels,TotalSlices,TotalTimepoints);

//Get Lut information from all channels and store in CHluts array	
	CHluts = lutGetter(TotalChannels);

// New position for pre-analysis option - with idea of performing analysis before making any other choices

runGraphquestion = getBoolean("Plot highest intensity Z-slices over time? (recommended) \n \n This can help choosing the best channel, reference timepoints and stack preprocessing options", "Yes", "Skip");
	if (runGraphquestion == 1){
		selectImage(SourceID);
		setBatchMode("hide");
		containR = scanR(TotalChannels,SourceID,TotalSlices,TotalTimepoints,PXWidth,PXHeight);
		graphR(containR,TotalChannels,TotalTimepoints);
		probeR = 0;
		for (c = 1; c < TotalChannels+1; c++){  // aim to get calculate average max. intensity Z-slices during the first 10 timepoints for all the channels
			for (i = 0; i < 10; i++) {
				probeR = probeR + containR [i*c];
			}	}
		probeR = probeR / (TotalChannels*10);
		if (probeR < (TotalSlices/4)){  //if average Z-position for first 10 TPs is in the first 25% of z-slices we add a default 25% slices to start
			zStart = round(TotalSlices/4);
			}
		if (probeR > (TotalSlices-(TotalSlices/4))){  //if average Z-position for first 10 TPs is in the last 25% of z-slices we add a default 25% slices to end
			zEnd = round(TotalSlices/4);
			}
		setBatchMode("show");		
		waitForUser("Analysis complete \n \n  Save, close or repostion graph output then click to continue");
		} 

// Option for manually-capturing Z-positions, posibility using a button menu
	selectImage(SourceID);
	Dialog.create("Bac-Live Analysis");
	Dialog.addChoice("Calculation Type:", newArray("Manual Z-Selection", "Semi-Automatic Z-detection"));
	Dialog.addMessage("\n If structures of interest are close to the top or bottom z-slice positions, \n it is recommended to add additional slices to prevent crop volumes exceeding stack limits");   
	Dialog.addSlider("Add extra slices to Z-stack starting position", 0, 50, zStart); 
	Dialog.addSlider("Add extra slices to Z-stack ending position", 0, 50, zEnd);
	Dialog.show();
	analysisType = Dialog.getChoice();
	runGraphquestion = Dialog.getCheckbox();
	AddZSliceStart = Dialog.getNumber();
    AddZSliceEnd = Dialog.getNumber();


// Z-stack preprocessing 

if (AddZSliceStart > 0) {
	for (i = 0; i < AddZSliceStart; i++) {
		Stack.setPosition(1, 1, 1);
		run("Add Slice", "add=slice prepend");
		TotalSlices = TotalSlices + 1;
		  }	}
if (AddZSliceEnd > 0) {
	for (i = 0; i < AddZSliceEnd; i++) {
		Stack.setPosition(1, TotalSlices, 1);
		run("Add Slice", "add=slice");
		TotalSlices = TotalSlices + 1;
		  } }

// Perform manual slice selction
	if(analysisType=="Manual Z-Selection"){
		setBatchMode(false); 
		selectImage(SourceID);	
		Stack.setPosition(1, 1, 1);
		waitForUser("Navigate to reference start Z-position");
		Stack.getPosition(MRefCh, MStartslice, MStartTP); 
		waitForUser("Navigate to reference end Z-position");
		selectImage(SourceID);	
		Stack.getPosition(MRefCh, MEndslice, MEndTP); 
		ManualDV=(MEndslice-MStartslice)/(MEndTP-MStartTP);
		showMessage("Calculated slice displacement is " + ManualDV + " slices per timepoint");
		ADRefCh = MRefCh;  //specifies remaining default values from image properties
		ADStartTP = MStartTP;
		ADEndTP = MEndTP;
		UserDV = ManualDV;
		MaxWeight = 0;
		setBatchMode(true); 
		} 
	
// DV Calculation Menu

	if(analysisType=="Semi-Automatic Z-detection"){  //Semi-Automatic calculation -  skips if manual selection was done
		MenuTitle = "BacLive Image Stablizer: DV Calculation";
		MenuContents = newArray("Autodetect Reference Channel (fluoroescence channel used for detection):","Autodetect Start Timepoint (macro will look here for Z-slice with highest signal):",
		"Autodetect End Timepoint (macro will look here for Z-slice with highest signal):","Max Value Weighting % (increases sensitivity)");
		MenuDefaults = newArray(1,1,TotalTimepoints,0,0);
		ReplyA = userMenu(MenuTitle, MenuContents, MenuDefaults); //menu function with array-specified data and collection fields
		ADRefCh = ReplyA[0]; 
		ADStartTP = ReplyA[1];
		ADEndTP = ReplyA[2];
		MaxWeight = ReplyA[3];
		UserDV = 0;
		}

// Check for user values are valid with respect to source data
	checkValues(ADRefCh, ADStartTP, ADEndTP, UserDV, MaxWeight,TotalChannels,TotalSlices, TotalTimepoints);
	
// Automatic calculation of displacement value unless user as specified a value.
	if (UserDV == 0) {
		TP = ADStartTP; // Sets TP value to Reference Start Timepoint
		print("Using channel "+ADRefCh+", timepoints "+ADStartTP+" and "+ADEndTP+" to detect changes in colony position");
		StartTPMaxSlice = findmaxSlice(ADRefCh, TP, MaxWeight); //uses function findmaxSlice to obtain Start slice with highest value
		TP = ADEndTP;	// Sets TP value to Reference End Timepoint
		EndTPMaxSlice = findmaxSlice(ADRefCh, TP, MaxWeight); //uses function findmaxSlice to obtain End slice with highest value
		DV = (EndTPMaxSlice - StartTPMaxSlice) / (ADEndTP-ADStartTP); 	// Calculate Average Displacement Value per Timepoint
		showMessage("Calculated slice displacement is " + DV + " slices per timepoint");

// If using manual slice chooser or specified a DV value.  
	
	} else {
		DV = UserDV;  // applies user DV value if given
		if (analysisType=="Manual Z-Selection"){
			StartTPMaxSlice = MStartslice; // applies slice chooser start slice if applicable: 
		} else {
			StartTPMaxSlice = 1; // otherwise defaults to 1 
			} }

//Added a "count back" calculation if slices were identified at later TPs so that we can start at TP 1 without losing Z-position syncronisation

	if(ADStartTP>1){
		StartTPMaxSlice = round(StartTPMaxSlice - (ADStartTP*DV));
	}
		
// DV Confirmation and Processing Menu
	MenuTitle = "BacLive Image Stablizer: Processing Menu";
	MenuContents = newArray("Provisional Displacement Value","Processing start timepoint:",
	"Processing end timepoint (will crop outside specified range):","Crop size in Z-slices",
	"OPTIONAL: Z-slice start value:");
	MenuDefaults = newArray(DV,1,TotalTimepoints,30,StartTPMaxSlice);
	ReplyB =  userMenu(MenuTitle, MenuContents, MenuDefaults); //menu function with array-specified data and collection fields
	DV = ReplyB[0];
	ProcessStartTP = ReplyB[1];
	ProcessEndTP = ReplyB[2];
	CropSize = ReplyB[3];
	UserZslice = ReplyB[4];
	
// applies user starting slice position if given 	
	if (UserZslice == 0) { 
		StartSlice = findmaxSlice(ADRefCh, ProcessStartTP, MaxWeight);  // otherwise use slice with max intensity at Processing Start position
	} else {
		StartSlice = UserZslice ;
		}

	print("Stack will be processed using slice "+StartSlice+" as colony starting point");
	
// Stack Processing section

	tempFolder = getDirectory("Chose temp folder for image processing");
	FinalTimepoints = ProcessEndTP - ProcessStartTP + 1;
	print("Processing...");
	for (timepoint=1; timepoint<FinalTimepoints+1; timepoint++) {  //Loop that 1) calculates z start and stop range according to DV 2)  duplicates sub-stack for 1 timepoint  3) deletes that timepoint from the original dataset (to save memory)
		showProgress(timepoint,FinalTimepoints);
		selectImage(SourceID);
//		Stack.setPosition(ADRefCh, 0, timepoint); 
		Zfloat=stackCalc(timepoint,DV,StartSlice,CropSize,TotalSlices,FinalTimepoints); //***New integraled function to determine Float window start and end.  Value is adjusted according to stack limits to avoid problems at early timepoints
		run("Duplicate...", "duplicate slices="+Zfloat[0]+"-"+Zfloat[1]+" frames="+timepoint); // duplicates  frame as determined by timepoint value
		if (timepoint<FinalTimepoints){	
			run("Hyperstack to Stack");
			run("Remove Slice Labels");
			saveAs("tiff", tempFolder+"TP"+(10000+timepoint)+".tif");	//saves substack to temp folder
			close();
			selectImage(SourceID);		
		} else if (timepoint==FinalTimepoints) {  //last timepoint actions
			saveAs("tiff", tempFolder+"TP"+(10000+timepoint)+".tif");	
			close();
			  	}
		}
//  New Stack creation and finishing 
	selectWindow(SourceName);
	close();  // close last timepoint of orignal data
	run("Image Sequence...", "open=["+tempFolder+"TP10001.tif] file=TP sort");
	rename(SourceName);  //restore original name
	setBatchMode("exit and display");  //necessary to exit batch mode for final step to work
	run("Stack to Hyperstack...", "order=xyczt(default) channels="+TotalChannels+" slices="+CropSize+" frames="+FinalTimepoints+" display=Color");  //restore hyperstack with correct values
	
// LUT recovery up to 5 channels
 	lutGiver(TotalChannels, CHluts);
	run("Channels Tool...");
	Stack.setDisplayMode("composite");
	
	for (timepoint=1; timepoint<FinalTimepoints+1; timepoint++){
		File.delete(tempFolder+"TP"+(10000+timepoint)+".tif");
		print("\\Clear");
		}
	showMessage("Temporary files deleted OK \n \nProcessing complete!");

// End of program
// Functions continue below

// Gets LUT color info (up to 5 channels) for final image reconstruction
	function lutGetter(TotalChannels){
		for (i=1; i<(TotalChannels+1); i++) {  //Loops through up to 5 channels assigning their LUTs to three arrays per channel
			Stack.setPosition(i, 1, 1); 	   //if more than 5 channels are present it shouldnt cause an error but LUT info wont be transfered.
	  		getLut(rA, gA, bA);
		    if (i==1){
	    		L = Array.concat(rA,gA,bA); }  // Using a single array to store LUT data
	  		if (i==2){
	   	    	L = Array.concat(L,rA,gA,bA); } //Takes LUT data from previous channel(s) and concatenates with current channel  	
	 		if (i==3){
				L = Array.concat(L,rA,gA,bA); }	 			
	   		if (i==4){
				L = Array.concat(L,rA,gA,bA); }   			
	 		if (i==5){
				L = Array.concat(L,rA,gA,bA); }				
				}
		return L;  //when all channels have been processed, the resulting composite array is sent back as CHluts
		}

// userMenu using three arrays to specify options and default options
	function userMenu(MenuTitle, MenuContents, MenuDefaults){
	Dialog.create(MenuTitle);
 	for (i = 0; i<MenuContents.length; i++) {  //will adjust to the number of menu items
	 		Dialog.addNumber(MenuContents[i], MenuDefaults[i]);
	 		}
	Dialog.show();
	for (i = 0; i<MenuContents.length; i++) {
	 		A = Dialog.getNumber();
	 		if (i==0){
	 			M = A;
	 		} else {
	 			M = Array.concat(M,A);  //collects all the replies into the M array
	 			}
	 		}
	return M;
	}


// Check validity of user values (currently no checking for 2nd menu)
	function checkValues(ADRefCh, ADStartTP, ADEndTP, UserDV, MaxWeight,TotalChannels,TotalSlices, TotalTimepoints){
			if (ADRefCh>TotalChannels || ADRefCh<1) {
				print("Channel out of range");
				exit();	}
			if (ADStartTP>TotalTimepoints || ADEndTP>TotalTimepoints) {
				print("Autodetect Timepoint(s) out of range");
				exit();	}
			if (ADStartTP>TotalTimepoints || ADEndTP>TotalTimepoints ) {
				print("Processing Timepoint(s) out of range");
				exit();	}
			if (MaxWeight<0 || MaxWeight>100)  {
				print("Max Intensity Weighting % should be between 0 and 100");
				exit();	}	
				}
			

	//findmaxSlice function
	function findmaxSlice(ADRefCh, TP, MaxWeight){
			HighestEval = 0;
			selectImage(SourceID);
 			for (i=1; i<=TotalSlices; i++) {
				Stack.setPosition(ADRefCh, i, TP);  //The channel and start frame is determined by the user
				makeRectangle(0, 0, PXWidth, PXHeight);
				run("Measure");
				SliceMean = getResult("Mean"); // grabs the latest result from the Mean Column
				SliceMax = getResult("Max"); // grabs the latest result from the Mean Column
				CE0 = ((SliceMean/255)*(100-MaxWeight))+((SliceMax/255)*MaxWeight); // Provisional combined evaluation using mean (100%-MaxWeight) and max values (MaxWeight).  Max values have been added to help where signal is low.
				run("Clear Results");
				if(i==1){
						CombinedEval = CE0;
						CE1 = CE0;
					} else if (i==2){
						CombinedEval = (CE0+CE1)/2;
						CE2 = CE1;
						CE1 = CE0;												
					} else if (i>3){
						CombinedEval = (CE0+CE1+CE2)/3;
						CE2 = CE1;
						CE1 = CE0;
					}	
				if (CombinedEval > HighestEval){  // this loop identifies the slice with the highest intensity
					HighestEval = CombinedEval;
					SliceChoice = i;
							}
				}
			return SliceChoice;
			}


	//stackCalc function  ***New integrated function
	function stackCalc(timepoint,DV,StartSlice,CropSize,TotalSlices,FinalTimepoints) {
		TPadjustment = (timepoint-1) * DV; // Calculates average change in stack position 
		a = round((StartSlice - (CropSize/2) + TPadjustment)); //Calculates where the Stack should start at a given timepoint; using half crop to place target slice in the centre then obtains nearest integer ; this could be below 1
		if (a < 1) {       //Fixes the first slice as the minimum ceiling
			a = 1;		   // as this is used as the reference for Zfloat end value keeps stack size constant
			}							
		else if ((a + CropSize-1)>TotalSlices){
			a = TotalSlices - (CropSize-1);  //Established maximium ceiling at the bottom of the stack by raising the starting position.
			}
        b = a + (CropSize-1);  //calculates end slice, cropsize-1 to compensate for the first slice; Should be within limits due to previous else if.
        c = newArray(a,b);
        return c;
		}
	
	//  lutGiver Function reappling original LUTs (stored in array CHluts) to each channel
 	function lutGiver(TotalChannels, CHluts) {
			for (activeCh = 1; activeCh<(TotalChannels+1); activeCh++) {  //loop through each channel
				Stack.setPosition(activeCh, 1, 1); 
					for (RGB = 1; RGB < 4; RGB++) {  //Loop through RGB intensities of each LUT
						if (RGB == 1){
							rA = lutChopper(CHluts, activeCh, RGB); } //uses function lutChopper to extract values
						if (RGB == 2){
							gA = lutChopper(CHluts, activeCh, RGB);	 }
						if (RGB == 3){
							bA = lutChopper(CHluts, activeCh, RGB);  }	
								}
					setLut(rA, gA, bA);
					}
 				}

	//  lutChopper extracts individual RGB intensity arrays for each LUT based on activeCh and RGB values
	function lutChopper(CHluts, activeCh, RGB) {
			Lutstart = ((activeCh-1)*768)+(256*(RGB-1));
			Lutend = ((activeCh-1)*768)+(256*(RGB-1))+256;
			A = Array.slice(CHluts,Lutstart,Lutend);
			return A;
		}

	// ScanR - Measure Z slice with max intensity at each timepoint and plot vs timepoint
	function scanR(TotalChannels,SourceID,TotalSlices,TotalTimepoints,PXWidth,PXHeight) {
		containR = newArray(0);
		for (c=1; c<=TotalChannels; c++){  //channel loop
			showProgress(c,TotalChannels);
			selectImage(SourceID);
				for (TP=1; TP<=TotalTimepoints; TP++){ //timepoint loop
					zMax = 0;
					maxMean = 0;
					for (i=1; i<=TotalSlices; i++) {  //Z-slice loop
						Stack.setPosition(c, i, TP);  
						getStatistics(area, mean);
						sliceMean = mean;
						if (sliceMean > maxMean) {
							maxMean = sliceMean;
							zMax = i;
								}
							}							// End Z-slice loop
					containR = Array.concat(containR,zMax);
				       }	// End Timepoint Loop
				 }	 // End Channel Loop	
			return containR;    //This should contain CHx(max slice per timepoints),CHx+1(max slice per timepoints)
			} // Function

				 
	//GraphR - Function to produce multicolor graphs for each channel
	function graphR	(containR,TotalChannels,TotalTimepoints){
				setBatchMode(false);
				XValues=newArray(0);
				ChLegend = "Channel 1";
				for (t=1; t<=TotalTimepoints; t++) {
						XValues = Array.concat(XValues,t);
					}
				lineType = newArray("circle","triangle","cross","dot","line");
			    graphColors = newArray("red","green","blue","yellow","magenta","lightGrey");
				Plot.create("Max Intensity Position", "Timepoint", "Z-slice number");
				for (c=1; c<=TotalChannels; c++){
					Plot.setColor(graphColors[c-1],graphColors[c-1]);  // Add this for color version
					YValues = Array.slice(containR,1+((c-1)*TotalTimepoints),c*TotalTimepoints);
					Plot.add(lineType[c-1], XValues, YValues);
					if (c>1){
						ChLegend = ChLegend+"\nChannel "+c+"";
						}
					}
				Plot.addLegend(ChLegend);
				Plot.show();		
				setBatchMode(true);
				}


				

	