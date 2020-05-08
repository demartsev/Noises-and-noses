

function action(input_vid, input_csv, vid_file,  output_dir)
{
run("Movie (FFMPEG)...", "choose=" + input_vid + video_list[j] + " use_virtual_stack first_frame=0 last_frame=-1");

originalName = File.nameWithoutExtension();

//originalNameWithoutExt = replace( originalName , ".avi" , "" ); 
outputName = replace(originalName , "C:/Users/vdemartsev/Documents/Thermal/Processed_FFMPEGS/", "");

//open XY file
run("Clear Results");
run("Results... ", "open=" + input_csv + outputName + ".csv"); 

//import XY coordinates as ROIs
for (i = 0; i < nResults; i++) 
   { 
      
      slice = getResult("f", i); 
      x = getResult("x", i); 
      y = getResult("y", i); 

      // Specify ROI size and shape
      run("Specify...", "width=10 height=5 x=&x y=&y slice=&slice centered"); 

      // Add to the ROI manager. 
      roiManager("Add"); 
     
   } 
      // Clear results table and measure color values 
   	  run("Clear Results");
	  roiManager("Measure");
	  // save results as csv. The measurments needs to be configured within ImageJ
      saveAs("Results", output_dir + outputName + "_out.csv");
      // clear all forms and close all windows
      roiManager("Delete");
      run("Clear Results");
      run("Close");
      close();

}

// input_dir 
input_vid = "C:/Users/vdemartsev/Documents/Thermal/Processed_FFMPEGS/";
input_csv = "C:/Users/vdemartsev/Documents/Thermal/conservative_annotation_095score/ImageJ_input/" ;
// output_dir 
output_dir = "C:/Users/vdemartsev/Documents/Thermal/conservative_annotation_095score/ImageJ_output/";
 
setBatchMode(true)
 
video_list = getFileList("C:/Users/vdemartsev/Documents/Thermal/Processed_FFMPEGS");
csv_list = getFileList("C:/Users/vdemartsev/Documents/Thermal/conservative_annotation_095score/ImageJ_input");
for (j=0; j<video_list.length; j++) 
			action(input_vid, input_csv, video_list[j],  output_dir); 


