#pragma once
#include "itkConfidenceConnectedImageFilter.h"
 #include "itkImage.h"
 #include "itkCastImageFilter.h"
 #include "itkCurvatureFlowImageFilter.h"
 #include "itkImageFileReader.h"
 #include "itkImageFileWriter.h"
 #include "itkVTKImageExport.h"
 #include "itkVTKImageImport.h"
 #include "vtkImageImport.h"
 #include "vtkImageExport.h"
 #include "vtkImageReader.h"
 #include "vtkPiecewiseFunction.h"
 #include "vtkVolumeProperty.h"
 #include "vtkVolumeRayCastMapper.h"
 #include "vtkVolumeTextureMapper3D.h"
 #include "vtkFixedPointVolumeRayCastMapper.h"
 #include "vtkVolumeRayCastCompositeFunction.h"
 #include "vtkVolume.h"
 #include "vtkRenderer.h"
 #include "vtkRenderWindow.h"
 #include "vtkRenderWindowInteractor.h"
#include <vtkColorTransferFunction.h>
#include <vtkCamera.h>
#include <vtkImageClip.h>
#include <vtkSphereSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkTesting.h>
#include <vtkStructuredPoints.h>
#include <vtkStructuredPointsReader.h>
#include "vtkcellpicker.h"
#include "vtkPointData.h"
#include <vtkImageData.h>
#include "vtkImagereader.h"
#include "vtkMetaImageReader.h"
#include "vtkGPUVolumeRayCastMapper.h"
#include <vtkAbstractVolumeMapper.h>
#include "vtkMath.h"
#include "vtkPlaneCollection.h"
#include "vtkTransform.h"
#include "vtkVoxel.h"
#include "vtkDoubleArray.h"
#include "vtkBox.h"
#include "vtkDataArray.h"
#include "vtkMetaImageWriter.h"
#include "vtkImageActor.h"
#include "vtkMapper.h"
#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
 
// headers needed for this example
#include <vtkImageData.h>
#include <vtkImageMapper.h>
#include <vtkImageCast.h>
#include <vtkMetaImageWriter.h>
#include <vtkMetaImageReader.h>
#include <vtkImageMandelbrotSource.h>
#include <vtkImageActor.h>
#include "vtkBoxWidget.h"
#include "vtkCommand.h"
#include "vtkPlanes.h"
#include "itkExtractImageFilter.h"

//headers needed for computing correlation

#include "itkNormalizedCorrelationImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkImageKernelOperator.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkMinimumMaximumImageCalculator.h"
#include "QuickView.h"
 #include "itkConnectedComponentImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include <iostream>
#include <string>
#include "ConnectedComponent.h"
class vtkBoxWidgetCallback : public vtkCommand
{
public:
  static vtkBoxWidgetCallback *New()
    { return new vtkBoxWidgetCallback; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      if (this->Mapper)
        {
        vtkPlanes *planes = vtkPlanes::New();
        widget->GetPlanes(planes);
        this->Mapper->SetClippingPlanes(planes);
        planes->Delete();
        }
    }
  void SetMapper(vtkGPUVolumeRayCastMapper* m) 
    { this->Mapper = m; }

protected:
  vtkBoxWidgetCallback() 
    { this->Mapper = 0; }

  vtkGPUVolumeRayCastMapper *Mapper;
};
typedef   float    InputPixelType;
typedef   float    OutputPixelType;
typedef  unsigned char SpherePixelType;
typedef itk::Image< InputPixelType,  3 >   InputImageType;
typedef itk::Image< SpherePixelType,  3 >   SphereImageType;

 typedef   unsigned short    IntPixelType;
typedef itk::Image< IntPixelType,  3 >   IntImageType;
  typedef itk::Image< OutputPixelType, 2 > OutputImageType2D;

class VTKRayCaster
{

 typedef itk::Image< OutputPixelType, 3 > OutputImageType;
  typedef itk::VTKImageExport< OutputImageType > ExportFilterType;
  ExportFilterType::Pointer itkExporter;
public:
	VTKRayCaster(InputImageType::Pointer image,InputImageType::Pointer imageCopy);
	~VTKRayCaster(void);

	  vtkImageImport* vtkImporter;
	  vtkImageData* imageData;
	  vtkAlgorithmOutput* imageDataCopy;
	  vtkAlgorithmOutput* imageout;
	  int Intersect(InputImageType::IndexType m_index,double radius,InputImageType::IndexType* outCenter,InputImageType::Pointer vesselImage,int Rnum,InputImageType::IndexType* queue);
	   vtkTransform *Transform; //use to perform ray transformation
	   int ClipLineWithPlanes(vtkPlaneCollection *planes,const double p1[3], const double p2[3],double &t1, double &t2, int& planeId);
	   double ComputeVolumeOpacity(const int xi[3], const double pcoords[3],vtkImageData *data, vtkDataArray *scalars,vtkPiecewiseFunction *scalarOpacity, vtkPiecewiseFunction *gradientOpacity);
	   vtkDoubleArray *Gradients; //used in volume picking
	   int ClipLineWithExtent(const int extent[6],const double x1[3], const double x2[3],double &t1, double &t2, int &planeId);
	   int writeVTKImage();
	   double Circle_RayCasting(InputImageType::IndexType m_Index,double angle_start,double angle_end,double radius,InputImageType::Pointer vesselImage);
	   void SetTwoDimImage(OutputImageType2D::Pointer image2D);
	   BOOL Find(InputImageType::IndexType index);
		InputImageType::Pointer image;
		InputImageType::Pointer imageCopy;
		void ComputeCorrelation();
		double param_t, intersect[3], paraCoord[3];
		double sourcePnt[3], destinPnt[3], normalVec[3],lastdistinPnt[3];
		double lastMaxOpacity;
		vtkSphereSource *sphere1 ;
		vtkDataArray *sphereNormals;

		vtkDataArray *scalars;
		void *scalarPtr;
		vtkGPUVolumeRayCastMapper* volumeMapper;
		vtkVolumeProperty* volumeProperty;
		vtkPiecewiseFunction* compositeOpacity;
		vtkColorTransferFunction* color;
		vtkPiecewiseFunction *scalarOpacity;
		vtkPiecewiseFunction *gradientOpacity;

		InputImageType::IndexType* Centerline;
		int numOfCenterlinePoint;
		OutputImageType2D::Pointer image2D;
		void writingFile();
		BOOL FindSpherePoint(InputImageType::IndexType index,int index1);
		InputImageType::IndexType spherePoint[1000];
		double FindNeighborIntensity(InputImageType::IndexType index);
		int vesselNum;
		InputImageType::IndexType Circle_RayCasting_Loop(InputImageType::IndexType m_Index,double angle_start,double angle_end,double radius,InputImageType::Pointer vesselImage,InputImageType::IndexType* queue);
		int ConnectedComponentLabeling(InputImageType::IndexType m_Index);
		SphereImageType::Pointer sphereImage;
		int  maximumValue(double *array,float th,InputImageType::IndexType* queue,int index);
		int addtoqueue(double *array,float th,InputImageType::IndexType* queue,int index,int maxIndex);
			int Indices[1000];
			double Relation[1000] ;
			double rayOpacity[1000] ;
			int N_j[1000];
			double F_j[1000];
			int index_queue ;
			IntImageType::Pointer ima;

};
