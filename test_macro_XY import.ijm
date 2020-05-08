
//open XY file
run("Results... ", "open=/test.csv"); 

//import XY coordinates as ROIs
for (i = 0; i < nResults; i++) 
   { 
      
      slice = getResult("f", i); 
      x = getResult("x", i); 
      y = getResult("y", i); 

      // Specify ROI size and shape
      run("Specify...", "width=5 height=5 x=&x y=&y slice=&slice centered"); 

      // Add to the ROI manager. 
      roiManager("Add"); 
   } 