cmake_minimum_required(VERSION 2.8)

PROJECT(Vessel)



find_package(ItkVtkGlue REQUIRED)
include(${ItkVtkGlue_USE_FILE})

 
add_executable(vessel vessel.cxx)
target_link_libraries(vessel
  ItkVtkGlue  ${VTK_LIBRARIES} ${ITK_LIBRARIES} vtkVolumeRendering  vtkHybrid  vtkRendering vtkCommon vtkWidgets)
