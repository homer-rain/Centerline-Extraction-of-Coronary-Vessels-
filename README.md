# Centerline-Extraction-of-Coronary-Vessels, August 2012

The goal of this project was to extract the centerline of coronary vessels in CTAngiography images. The input of the software is 3d images in .raw format and the result is going to be saved
in output.raw file. We need another software like paraview or VolView to visualize the results. 
This software implemented in C++ (windows) using two libraries of ITK3.2 and VTK 5.8. 
I could propose a new method for detection of coronary Ostia (the location in the aorta in that coronary vessels branches starts). The centerline extraction method was implemented based on this (https://www.ncbi.nlm.nih.gov/pubmed/21637981) paper's method.

## Description
This software has three features:
- load CTAngiography images.
- Extract the centerline of Coronary vessels plus a part of aorta.
- Evaluate the precision the extracted vessels.

## Binding to Other libraries
We used 
- ITK3.2
- VTK 5.8



## Further Reading
[Automatic Centerline Extraction of Coronary Arteries in Computed
Tomography Angiography Images](http://confnews.um.ac.ir/images/41/conferences/icee2013/2009_2.pdf)

## install
In order to run this code, you need to 
- Install and build ITK3.2 library.
- Install and build VTK5.8 library.
- Config the code using cmake with the installed libraries. 

and you can use the software.
