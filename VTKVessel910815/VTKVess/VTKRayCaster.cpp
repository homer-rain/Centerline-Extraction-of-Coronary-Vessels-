#include "VTKRayCaster.h"



  double sphereCent[3];
int  maximumValue_2(double *array,int index)
{
     int length = index;  // establish size of array
     double max = array[0];       // start with max = first element
	int ind = 0;

     for(int i = 1; i<length; i++)
     {
          if(array[i] > max)
		  {
                max = array[i];
				ind = i;}
     }
	     return ind; 
}

int VTKRayCaster:: addtoqueue(double *array,float th,InputImageType::IndexType* queue,int index,int maxIndex)
{
	InputImageType::IndexType m;
	double max = Relation[maxIndex];
	int length  = index;
	int k = 0;
	 for(int i = 1; i<length; i++)
     {
		 if(abs(array[i] - max)<th&&max>th&&i!=maxIndex&&k++<4)
		  {
			sphere1->GetOutput()->GetPoint(Indices[i], destinPnt);
			m[0] = destinPnt[0]+sphereCent[0];
			m[1] = destinPnt[1]+sphereCent[1];
			m[2] = destinPnt[2]+sphereCent[2];
			//  spherePoint
			if(!Find(m)&&imageCopy->GetPixel(m)>10)
			{
			queue[index_queue][0] = vtkMath::Round( destinPnt[0] )+sphereCent[0];
			queue[index_queue][1] = vtkMath::Round( destinPnt[1]) +sphereCent[1];
			queue[index_queue++][2] = vtkMath::Round(destinPnt[2]) + sphereCent[2];
			//image->SetPixel(queue[index_queue-1],255);
			}

			

		  }
                
     }
	 return 1;
}
int VTKRayCaster:: maximumValue(double *array,float th,InputImageType::IndexType* queue,int index)
{
     int length = index;  // establish size of array
     double max = array[0];       // start with max = first element
	int ind = 0;

     for(int i = 1; i<length; i++)
     {
          if(array[i] > max)
		  {
                max = array[i];
				ind = i;}
     }
	 InputImageType::IndexType m;
	
     return ind;                // return highest value in array
}


int minimumValue(double *array,int index)
{
     int length = index;  // establish size of array
     double max = array[0];       // start with max = first element
	int ind = 0;

     for(int i = 1; i<length; i++)
     {
          if(array[i] > max)
		  {
                max = array[i];
				ind = i;}
     }
     return ind;                // return highest value in array
}
// This function will connect the given itk::VTKImageExport filter to  the given vtkImageImport filter.
 template <typename ITK_Exporter, typename VTK_Importer>

void ConnectPipelines(ITK_Exporter exporter, VTK_Importer* importer)
 {
 importer->SetUpdateInformationCallback(exporter->GetUpdateInformationCallback());
 importer->SetPipelineModifiedCallback(exporter->GetPipelineModifiedCallback());
 importer->SetWholeExtentCallback(exporter->GetWholeExtentCallback());
 importer->SetSpacingCallback(exporter->GetSpacingCallback());
 importer->SetOriginCallback(exporter->GetOriginCallback());
 importer->SetScalarTypeCallback(exporter->GetScalarTypeCallback());
 importer->SetNumberOfComponentsCallback(exporter->GetNumberOfComponentsCallback());
 importer->SetPropagateUpdateExtentCallback(exporter->GetPropagateUpdateExtentCallback());
 importer->SetUpdateDataCallback(exporter->GetUpdateDataCallback());
 importer->SetDataExtentCallback(exporter->GetDataExtentCallback());
 importer->SetBufferPointerCallback(exporter->GetBufferPointerCallback());
 importer->SetCallbackUserData(exporter->GetCallbackUserData());
 }
VTKRayCaster::VTKRayCaster(InputImageType::Pointer image,InputImageType::Pointer imageCopy)
{
index_queue = 0;
 for(int i=0;i<1000;i++)
	 rayOpacity[i] = 0;
  for(int i=0;i<1000;i++)
	  Relation[i] = 0;

for (int i = 0;i<1000;i++)
	N_j[i] = 0;


for (int i = 0;i<1000;i++)
	F_j[i] = 0;

for (int i = 0;i<1000;i++)
	Indices[i] = 0;

	InputImageType::RegionType outputSphere;
	InputImageType::RegionType::IndexType outputStart;
	outputStart[0] = 0;
	outputStart[1] = 0;
	outputStart[2] = 0;
	outputSphere.SetSize( image->GetLargestPossibleRegion().GetSize() );
	outputSphere.SetIndex( outputStart );

	
	 sphereImage = SphereImageType::New();
	sphereImage->SetRegions( outputSphere );

	const InputImageType::SpacingType& spacing = image->GetSpacing();
	const InputImageType::PointType& inputOrigin = image->GetOrigin();
	double outputOrigin[ 3 ];
	for(unsigned int i=0; i< 3; i++)
	{
	outputOrigin[i] = inputOrigin[i] ;
	}
	sphereImage->SetSpacing( spacing );
	sphereImage->SetOrigin( outputOrigin );
	sphereImage->Allocate();
	
	sphereImage->FillBuffer(0);



	/*this->Cell = vtkGenericCell::New();
  this->PointIds = vtkIdList::New();*/
	  this->Gradients = vtkDoubleArray::New();
	   this->Gradients->SetNumberOfComponents(3);
		this->Gradients->SetNumberOfTuples(8);
		 itkExporter = ExportFilterType::New();
		 imageData = vtkImageData::New();
		 this->image = image;
		 typedef itk::ImageFileReader<InputImageType> FileReaderType;
		  FileReaderType::Pointer imageReader = FileReaderType::New();
		  //imageReader->SetFileName("g:\\RawReader\\output_en2.mhd");
		  //imageReader->SetFileName("g:\\CTImage\\output_en.mhd");
		   imageReader->SetFileName("F:\\CTImage\\training\\dataset02\\output_best1.mhd");
  try
  { 
    imageReader->Update();
  }
  catch (itk::ExceptionObject &ex)
  { 
    std::cout << ex << std::endl;
    //return EXIT_FAILURE;
  }
  this->imageCopy = imageReader->GetOutput();
		// rayOpacity = new double(1000);
		 //angle_max = new double(1000);
	lastMaxOpacity = 0;
	itkExporter->SetInput(image);
	/*vtkImageImport**/ vtkImporter = vtkImageImport::New();
	itkExporter->Update();
	ConnectPipelines(itkExporter, vtkImporter);
	vtkImporter->Update();
	imageData = vtkImageData::New();
	imageData = vtkImporter->GetOutput();
	  sphere1 = vtkSphereSource::New();
  sphere1->SetThetaResolution(100);
  sphere1->SetPhiResolution(100);
  sphere1->SetRadius(2);
  
  sphere1->Update();

   // the normals obtained from the outer sphere
  sphereNormals = sphere1->GetOutput()->GetPointData()->GetNormals();

    volumeMapper = vtkGPUVolumeRayCastMapper::New();
  volumeMapper->SetBlendModeToComposite(); // composite first
#if VTK_MAJOR_VERSION <= 5
  volumeMapper->SetInputConnection(imageData->GetProducerPort());
#else
  volumeMapper->SetInputData(imageData);
#endif  
   volumeProperty = vtkVolumeProperty::New();
  volumeProperty->ShadeOff();
  volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);

 compositeOpacity = vtkPiecewiseFunction::New();
	  compositeOpacity->AddPoint(-3024, 0, 0.5, 0.0 );
      compositeOpacity->AddPoint(-155, 0, 0.5, 0.92 );
      compositeOpacity->AddPoint(217, .68, 0.33, 0.45 );
      compositeOpacity->AddPoint(420,.83, 0.5, 0.0);
      compositeOpacity->AddPoint(3071, .80, 0.5, 0.0);
	 volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.
 
	 color = vtkColorTransferFunction::New();
	color->AddRGBPoint( -3024, 0, 0, 0, 0.5, 0.0 );
	color->AddRGBPoint( -155, .55, .25, .15, 0.5, .92 );
	color->AddRGBPoint( 217, .88, .60, .29, 0.33, 0.45 );
	color->AddRGBPoint( 420, 1, .94, .95, 0.5, 0.0 );
	color->AddRGBPoint( 3071, .83, .66, 1, 0.5, 0.0 );
	volumeProperty->SetColor(color);
	volumeMapper->SetBlendModeToComposite();
	volumeProperty->ShadeOn();
	volumeProperty->SetAmbient(0.1);
	volumeProperty->SetDiffuse(0.9);
	volumeProperty->SetSpecular(0.2);
	volumeProperty->SetSpecularPower(10.0);
	volumeProperty->SetScalarOpacityUnitDistance(0.8919);
	vtkVolume* volume = vtkVolume::New();
	volume->SetMapper(volumeMapper);
	volume->SetProperty(volumeProperty);

	////create image
	//		/////////////////////////////////////////////////////////
	//			ima = IntImageType::New();
	//			IntImageType::RegionType region;
	//			IntImageType::IndexType start;
	//			start[0] = 0;
	//			start[1] = 0;
	//			start[2] = 0;
	//			IntImageType::SizeType size;
	//		
	//			size[0] = sphere1->GetRadius()+1;
	//			size[1] = sphere1->GetRadius()+1;
	//			size[2] = sphere1->GetRadius()+1;

	//			region.SetSize(size);
	//			region.SetIndex(start);

	//			ima->SetRegions(region);
	//			ima->Allocate();

	//			/*for(unsigned int r = 0; r < 2*sphere1->GetRadius()+1; r++)
	//			{
	//			for(unsigned int c = 0; c < 2*sphere1->GetRadius()+1; c++)
	//			{
	//				for(unsigned int k = 0; k < 2*sphere1->GetRadius()+1; k++)
	//			{
	//			IntImageType::IndexType pixelIndex;
	//			pixelIndex[0] = r;
	//			pixelIndex[1] = c;
	//			pixelIndex[2] = k;

	//			ima->SetPixel(pixelIndex, 0);
	//			}
	//			}
	//			}*/



  int numComponents = imageData->GetNumberOfScalarComponents();

  scalars = vtkDataArray::CreateDataArray(imageData->GetScalarType());
  scalars->SetNumberOfComponents(numComponents);
  vtkIdType scalarArraySize = numComponents*imageData->GetNumberOfPoints();
  int scalarSize = imageData->GetScalarSize();
  scalarPtr = imageData->GetScalarPointer();


  int independentComponents = volumeProperty->GetIndependentComponents();
  int numIndependentComponents = 1;
  if (independentComponents)
    {
    numIndependentComponents = numComponents;
    }

  // Create a scalar array, it will be needed later


  // Go through each volume component separately
  double tMin = VTK_DOUBLE_MAX;
  for (int component = 0; component < numIndependentComponents; component++)
    { 
    vtkPiecewiseFunction *scalarOpacity = volumeProperty->GetScalarOpacity(component);
    int disableGradientOpacity = volumeProperty->GetDisableGradientOpacity(component);
    vtkPiecewiseFunction *gradientOpacity = 0;
   // if (!disableGradientOpacity && this->UseVolumeGradientOpacity)
      {
		  gradientOpacity = volumeProperty->GetGradientOpacity(component);
      }

    // This is the component used to compute the opacity
    int oComponent = component;
    if (!independentComponents)
      {
      oComponent = numComponents - 1;
      }

    // Make a new array, shifted to the desired component
    scalars->SetVoidArray(static_cast<void *>(static_cast<float *>(scalarPtr)
                                              + scalarSize*oComponent),
                          scalarArraySize, 1);
  }


  Centerline = new InputImageType::IndexType[7000];
  numOfCenterlinePoint = 0;




}

VTKRayCaster::~VTKRayCaster(void)
{
	  
}

int VTKRayCaster:: writeVTKImage()
{

   // adapt path !
	std::string filePath = "g:\\CTImage\\vtkOutput.mhd";
	std::string filePathRaw = "g:\\CTImage\\vtkOutput.raw";

 
  vtkImageCast* castFilter = 
      vtkImageCast::New();
   castFilter->SetOutputScalarTypeToFloat();
   castFilter->SetInputConnection(imageDataCopy);
   castFilter->Update();
 
  vtkMetaImageWriter* writer = vtkMetaImageWriter::New();
  writer->SetInputConnection(castFilter->GetOutputPort());
   writer->SetFileName(filePath.c_str());
   writer->SetRAWFileName(filePathRaw.c_str());
   writer->Write();
 
   // Read and display file for verification that it was written correctly
   vtkMetaImageReader* reader = vtkMetaImageReader::New();
   reader->SetFileName(filePath.c_str());
   reader->Update();
 
 /* vtkImageActor* actor = vtkImageActor::New();
   actor->GetMapper()->SetInputConnection(imageData->GetProducerPort());*/

 
 /* vtkRenderer* renderer = vtkRenderer::New();
  vtkRenderWindow* renderWindow = vtkRenderWindow::New();
   renderWindow->AddRenderer(renderer);
  vtkRenderWindowInteractor* renderWindowInteractor = vtkRenderWindowInteractor::New();
   renderWindowInteractor->SetRenderWindow(renderWindow);
 
   renderer->AddActor(actor);
   renderer->SetBackground(.2, .3, .4);
 
   renderWindow->Render();
   renderWindowInteractor->Start();*/
 
   return EXIT_SUCCESS;
}


int VTKRayCaster::ClipLineWithPlanes(vtkPlaneCollection *planes, const double p1[3], const double p2[3],  double &t1, double &t2, int& planeId)
{
  // The minPlaneId is the index of the plane that t1 lies on
  planeId = -1;
  t1 = 0.0;
  t2 = 1.0;

  vtkCollectionSimpleIterator iter;
  planes->InitTraversal(iter);
  vtkPlane *plane;
  for (int i = 0; (plane = planes->GetNextPlane(iter)); i++)
    {
    // This uses EvaluateFunction instead of FunctionValue because,
    // like the mapper, we want to ignore any transform on the planes.
    double d1 = plane->EvaluateFunction(const_cast<double *>(p1));
    double d2 = plane->EvaluateFunction(const_cast<double *>(p2));

    // If both distances are negative, both points are outside
    if (d1 < 0 && d2 < 0)
      {
      return 0;
      }
    // If only one of the distances is negative, the line crosses the plane
    else if (d1 < 0 || d2 < 0)
      {
      // Compute fractional distance "t" of the crossing between p1 & p2
      double t = 0.0;

      // The "if" here just avoids an expensive division when possible
      if (d1 != 0)
        {
        // We will never have d1==d2 since they have different signs
        t = d1/(d1 - d2);
        }

      // If point p1 was clipped, adjust t1
      if (d1 < 0)
        {
        if (t >= t1)
          {
          t1 = t;
          planeId = i;
          }
        }
      // else point p2 was clipped, so adjust t2
      else
        {
        if (t <= t2)
          {
          t2 = t;
          } 
        }

      // If this happens, there's no line left
      if (t1 > t2)
        {
        return 0;
        }
      }
    }

  return 1;
}


int VTKRayCaster::ClipLineWithExtent(const int extent[6], const double x1[3], const double x2[3],  double &t1, double &t2, int &planeId)
{
  double bounds[6];
  bounds[0] = extent[0]; bounds[1] = extent[1]; bounds[2] = extent[2];
  bounds[3] = extent[3]; bounds[4] = extent[4]; bounds[5] = extent[5];

  int p2;
  return vtkBox::IntersectWithLine(bounds, x1, x2, t1, t2, 0, 0, planeId, p2);
}
double VTKRayCaster::ComputeVolumeOpacity(  const int xi[3], const double pcoords[3],  vtkImageData *data, vtkDataArray *scalars,  vtkPiecewiseFunction *scalarOpacity, vtkPiecewiseFunction *gradientOpacity)
{
  double opacity = 1.0;

  // Get interpolation weights from the pcoords
  double weights[8];
  vtkVoxel::InterpolationFunctions(const_cast<double *>(pcoords), weights);

  // Get the volume extent to avoid out-of-bounds
  int extent[6];
  data->GetExtent(extent);
  int scalarType = data->GetScalarType();

  // Compute the increments for the three directions, checking the bounds
  vtkIdType xInc = 1;
  vtkIdType yInc = extent[1] - extent[0] + 1;
  vtkIdType zInc = yInc*(extent[3] - extent[2] + 1);
  if (xi[0] == extent[1]) { xInc = 0; }
  if (xi[1] == extent[3]) { yInc = 0; }
  if (xi[2] == extent[5]) { zInc = 0; }

  // Use the increments and weights to interpolate the data
  vtkIdType ptId = data->ComputePointId(const_cast<int *>(xi));
  double val = 0.0;
  for (int j = 0; j < 8; j++)
    {
    vtkIdType ptInc = (j & 1)*xInc + ((j>>1) & 1)*yInc + ((j>>2) & 1)*zInc;
    val += weights[j]*scalars->GetComponent(ptId + ptInc, 0);
    }

  // Compute the ScalarOpacity
  if (scalarOpacity)
    {
    opacity *= scalarOpacity->GetValue(val);
    }
  else if (scalarType == VTK_FLOAT || scalarType == VTK_DOUBLE)
    {
    opacity *= val;
    }
  else
    {
    // Assume unsigned char
    opacity *= val/255.0;
    }

  // Compute gradient and GradientOpacity
  if (gradientOpacity)
    {
    data->GetVoxelGradient(xi[0], xi[1], xi[2], scalars, this->Gradients);
    double v[3]; v[0] = v[1] = v[2] = 0.0;
    for (int k = 0; k < 8; k++)
      {
      double *pg = this->Gradients->GetTuple(k);
      v[0] += pg[0]*weights[k]; 
      v[1] += pg[1]*weights[k]; 
      v[2] += pg[2]*weights[k]; 
      }
    double grad = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    opacity *= gradientOpacity->GetValue(grad);
    }


  return opacity;
}




double VTKRayCaster:: FindNeighborIntensity(InputImageType::IndexType m_index)
{
	InputImageType::IndexType temp;
	temp = m_index;
	double sum = imageCopy->GetPixel(temp);
	temp.m_Index[0] = m_index.m_Index[0]+1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp = m_index;
	temp.m_Index[0] = m_index.m_Index[0]-1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+1;
	
	sum =sum+abs( imageCopy->GetPixel(temp));

	
	temp.m_Index[1] = m_index.m_Index[1]-1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	if(vesselNum == 0)

		return sum/9;
	temp = m_index;
	
	temp.m_Index[0] = m_index.m_Index[0]+2;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-2;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+2;
	sum =sum+abs( imageCopy->GetPixel(temp));

	temp = m_index;
	temp.m_Index[0] = m_index.m_Index[0]-2;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-1;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+1;
	
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-2;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+2;
	
	temp = m_index;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]-2;
	sum =sum+abs( imageCopy->GetPixel(temp));
	temp.m_Index[1] = m_index.m_Index[1]+2;
	sum =sum+abs( imageCopy->GetPixel(temp));
	sum = sum/21;




	return sum;




}





double teta_max[500];
double phi_max[500];

 std::ofstream out( "textfile.txt" );

#define VTKCELLPICKER_VOXEL_TOL 1e-6
 int VTKRayCaster::Intersect(InputImageType::IndexType spherecenter,double  radius,InputImageType::IndexType* outCenter,InputImageType::Pointer vesselImage,int Rnum,InputImageType::IndexType* queue )
{

int index1 = 0;

if(!Find(spherecenter))
	return 0;

// InputImageType::IndexType* SphereIndex = new InputImageType::IndexType(500);
 //InputImageType::Pointer imageCopy = image;
 //itkExporter->SetInput(image);
 ///*vtkImageImport**/ vtkImporter = vtkImageImport::New();
 //itkExporter->Update();
 //ConnectPipelines(itkExporter, vtkImporter);
 //vtkImporter->Update();
 //
 //imageData = vtkImporter->GetOutput();

 //double spherecenter[3];
//spherecenter[0] = 66;
//spherecenter[1] = 108;
//spherecenter[2] = 30;



  // init the counter and ray length
  int numIntersected = 0;
  double rayLen = radius; // = 1 - 0.8 + error tolerance
  int sub_id;
  vtkIdType cell_id;
  
//  vtkGenericCell *cell = vtkGenericCell::New();

// Clip the ray with the extent, results go in s1 and s2
  int planeId;
 
 

double max_responce[3];

//////////////////////////////////////////////////////////////////////////
if(Rnum == 2)
{
 typedef itk::ImageFileReader<InputImageType> FileReaderType;
		  FileReaderType::Pointer imageReader = FileReaderType::New();
		  //imageReader->SetFileName("g:\\RawReader\\output_en2.mhd");
		  //imageReader->SetFileName("g:\\CTImage\\FastMarchingFilterOutput2.mhd");
		  imageReader->SetFileName("F:\\CTImage\\training\\dataset02\\output_best1.mhd");
  try
  { 
    imageReader->Update();
  }
  catch (itk::ExceptionObject &ex)
  { 
    std::cout << ex << std::endl;
    //return EXIT_FAILURE;
  }
  this->imageCopy = imageReader->GetOutput();
}	
else if(Rnum == 0||Rnum == 1)
{
	typedef itk::ImageFileReader<InputImageType> FileReaderType;
		  FileReaderType::Pointer imageReader = FileReaderType::New();
		// imageReader->SetFileName("g:\\RawReader\\output_en2.mhd");
		  imageReader->SetFileName( "F:\\CTImage\\training\\dataset02\\vesselness.mhd");
		   //imageReader->SetFileName("g:\\CTImage\\vesselness.mhd");
		 // imageReader->SetFileName("g:\\CTImage\\FastMarchingFilterOutput2.mhd");
  try
  { 
    imageReader->Update();
  }
  catch (itk::ExceptionObject &ex)
  { 
    std::cout << ex << std::endl;
    //return EXIT_FAILURE;
  }
  this->imageCopy = imageReader->GetOutput();
}
		//this->imageCopy = vesselImage;
///////////////////////////////////////////////////////////////////////////
sphere1->SetRadius(radius);
//sphere1->SetStartTheta  (-vnl_math::pi/2);
//sphere1->SetEndTheta( vnl_math::pi/2);
sphere1->Update();


   // the normals obtained from the outer sphere
  sphereNormals = sphere1->GetOutput()->GetPointData()->GetNormals();
double phi = sphere1->GetStartPhi();
  phi = sphere1->GetEndPhi();
  phi = sphere1->GetStartTheta();
  phi = sphere1->GetEndTheta();

  sphereCent[0] = spherecenter.m_Index[0];
  sphereCent[1] = spherecenter.m_Index[1];
  sphereCent[2] = spherecenter.m_Index[2];
  
	double frameRate = 10.0;

	int clip = 1;


  //create image
			/////////////////////////////////////////////////////////
				ima = IntImageType::New();
				InputImageType::RegionType region;
				InputImageType::IndexType start;
				start[0] = 0;
				start[1] = 0;
				start[2] = 0;
				InputImageType::SizeType size;
			
				size[0] = 2*sphere1->GetRadius()+1;
				size[1] = 2*sphere1->GetRadius()+1;
				size[2] = 2*sphere1->GetRadius()+1;

				region.SetSize(size);
				region.SetIndex(start);

				ima->SetRegions(region);
				ima->Allocate();

				/*for(unsigned int r = 0; r < 2*sphere1->GetRadius()+1; r++)
				{
				for(unsigned int c = 0; c < 2*sphere1->GetRadius()+1; c++)
				{
					for(unsigned int k = 0; k < 2*sphere1->GetRadius()+1; k++)
				{
				IntImageType::IndexType pixelIndex;
				pixelIndex[0] = r;
				pixelIndex[1] = c;
				pixelIndex[2] = k;

				ima->SetPixel(pixelIndex, 0);
				}
				}
				}*/

////////////////////////////////////////////////////////////////


    double q1[3], q2[3];
//    this->Transform->TransformPoint(sourcePnt, q1);
   // this->Transform->TransformPoint(destinPnt, q2);

	//sphere1->SetEndPhi(vnl_math::pi);
	//sphere1->SetEndTheta(vnl_math::pi);

int i=1;
index1 = 0;
InputImageType::IndexType temp,temp2;

  //sphere1->SetStartPhi(0);
  //sphere1->SetEndPhi(45);
  //lastdistinPnt[2] = 108;
sourcePnt[0] = sphereCent[0];///*sourcePnt[0]+*/spherecenter.m_Index[0];
sourcePnt[1] = sphereCent[1];///*sourcePnt[1]+*/spherecenter.m_Index[1];
sourcePnt[2] = sphereCent[2];///*sourcePnt[2]+*/spherecenter.m_Index[2];
int threshold = radius;


	
			/////////////////////////////////////////////////////////////

if(Rnum == 1)
{
	
	threshold = (int)sourcePnt[2]+radius;
}

else if(Rnum == 0||Rnum==2)
	threshold = (int)sourcePnt[2];
	vesselNum  = Rnum;
 for ( int i = 0; i < sphere1->GetOutput()->GetNumberOfPoints(); i ++ )

    {
		
		sphere1->GetOutput()->GetPoint(i, destinPnt);
		sphereNormals->GetTuple(i, normalVec);
	
	
    // cast a ray in the negative direction toward sphere1
    //sourcePnt[0] = destinPnt[0] - rayLen * normalVec[0];
    //sourcePnt[1] = destinPnt[1] - rayLen * normalVec[1];
    //sourcePnt[2] = destinPnt[2] - rayLen * normalVec[2];
 
		
	temp.m_Index[0] = vtkMath::Round(destinPnt[0])+sphere1->GetRadius();
	temp.m_Index[1] = vtkMath::Round(destinPnt[1])+sphere1->GetRadius();
	temp.m_Index[2] = vtkMath::Round(destinPnt[2])+sphere1->GetRadius();
	destinPnt[0] = vtkMath::Round( destinPnt[0] )+spherecenter.m_Index[0];
	destinPnt[1] = vtkMath::Round( destinPnt[1]) +spherecenter.m_Index[1];
	destinPnt[2] = vtkMath::Round(destinPnt[2]) + spherecenter.m_Index[2];
	temp2.m_Index[0] = destinPnt[0];
	temp2.m_Index[1] = destinPnt[1];
	temp2.m_Index[2] = destinPnt[2];
	/*if(((int)destinPnt[0]!=lastdistinPnt[0]||(int)destinPnt[1]!=lastdistinPnt[1]||(int)destinPnt[2]!=lastdistinPnt[2]))
	{
		if(destinPnt[2]<=sourcePnt[2]+2)
		{
	lastdistinPnt[0] = (int)destinPnt[0];
	lastdistinPnt[1] = (int)destinPnt[1];
	lastdistinPnt[2] = (int)destinPnt[2];*/


	 if((lastdistinPnt[0]!=vtkMath::Round(destinPnt[0]))||(lastdistinPnt[1]!=vtkMath::Round( destinPnt[1]))||(lastdistinPnt[2]!=vtkMath::Round( destinPnt[2])))
	if(!FindSpherePoint(temp2,index1))
	{
		if(vtkMath::Round(destinPnt[2])<=threshold)
		{
			lastdistinPnt[0] = vtkMath::Round(destinPnt[0]);
			lastdistinPnt[1] = vtkMath::Round(destinPnt[1]);
			lastdistinPnt[2] = vtkMath::Round(destinPnt[2]);
			spherePoint[index1].m_Index[0] = lastdistinPnt[0];
			spherePoint[index1].m_Index[1] = lastdistinPnt[1];
			spherePoint[index1].m_Index[2] = lastdistinPnt[2];
		
		
		
			
			ima->SetPixel(temp,imageCopy->GetPixel(spherePoint[index1]));
 			std::cout<<temp<<(int)ima->GetPixel(temp)<<"\n";

	// imageData = vtkImageData::New();
	
	
  
 //vtkMetaImageReader *reader = vtkMetaImageReader::New();
 //reader->SetFileName("g:\\CTImage\\downsampledimage1.mhd");
 //   reader->Update();
 //  // imageData=reader->GetOutput();
	//imageData = vtkImporter->GetOutput();
	////imageDataCopy = reader->GetOutputPort();
	//imageDataCopy = vtkImporter->GetOutputPort();
   
  double spacing[3], origin[3];
  int extent[6];
  int boxExtent[6];
  boxExtent[0] = sourcePnt[0];
  boxExtent[1] = sourcePnt[1];
  boxExtent[2] = sourcePnt[2];
  boxExtent[3] = destinPnt[0];
  boxExtent[4] = destinPnt[1];
  boxExtent[5] = destinPnt[2];
  imageData->GetSpacing(spacing);
  imageData->GetOrigin(origin);
  imageData->GetExtent(extent);

  
    

 double t1 = 0;
	double t2 = 1;



  double x1[3], x2[3];
  for (int i = 0; i < 3; i++)
    {
		x1[i] = (sourcePnt[i]); /*- origin[i])/spacing[i];*/
		x2[i] = (destinPnt[i]);/*- origin[i])/spacing[i];*/
    }
 /* if (!this->ClipLineWithPlanes(planes, sourcePnt, destinPnt, t1, t2, clippingPlaneId))
      {
      return VTK_DOUBLE_MAX;
      }*/
   


  //volumeMapper->GetClippingPlanes();
 
double s1, s2;

//boxExtent[0] = x1[0];
//boxExtent[1] = x1[1];
//boxExtent[2] = x1[2];
//boxExtent[3] = x2[0];
//boxExtent[4] = x2[1];
//boxExtent[5] = x2[2];
//if (!this->ClipLineWithExtent(extent, x1, x2, s1, s2, planeId))
//    {
//    return VTK_DOUBLE_MAX;
//    }
//  if (s1 >= t1) { t1 = s1; }
//  if (s2 <= t2) { t2 = s2; }
//
//  // Sanity check
//  if (t2 < t1)
//    {
//    return VTK_DOUBLE_MAX;
//    }
      
  // Compute the length of the line intersecting the volume
  double rayLength = sqrt(vtkMath::Distance2BetweenPoints(x1, x2))*(t2 - t1); 

  // This is the minimum increment that will be allowed
  double tTol = VTKCELLPICKER_VOXEL_TOL/rayLength*(t2 - t1);

  // Find out whether there are multiple components in the volume
  
     
  
  //else if ( (volumeMapper = vtkAbstractVolumeMapper::SafeDownCast(m)) )
  //  {
  //  tMin = this->IntersectVolumeWithLine(p1, p2, t1, t2, prop, volumeMapper);
  //  }


	 // Do a ray cast with linear interpolation.
    double opacity = 0.0;
    double lastOpacity = 0.0;
    double lastT = t1;
    double x[3];
    double pcoords[3];
    int xi[3];

    // Ray cast loop
    double t = t1;
    while (t < t2)
      {
      for (int j = 0; j < 3; j++)
        {
        // "t" is the fractional distance between endpoints x1 and x2
        x[j] = x1[j]*(1.0 - t) + x2[j]*t;

        // Paranoia bounds check
        if (x[j] < extent[2*j]) { x[j] = extent[2*j]; }
        else if (x[j] > extent[2*j + 1]) { x[j] = extent[2*j+1]; }
	/*	 if (x[j] > x2[j])
		 {
			{ x[j] = x2[j]; }
			t=2;
		}*/
       
        xi[j] = int(floor(x[j]));
        pcoords[j] = x[j] - xi[j];
        }

	  opacity = this->ComputeVolumeOpacity(xi, pcoords, imageData, scalars,scalarOpacity, gradientOpacity);
	//  teta_max[index] = sphere1->GetThetaResolution();
	//  phi_max[index] = sphere1->GetPhiResolution();;
	  rayOpacity[index1] += opacity;
	  Indices[index1] = i;
	//  InputImageType::IndexType m;
	//  m.m_Index[0] = xi[0];
	//  m.m_Index[1] = xi[1];
	//  m.m_Index[2] = xi[2];
	////  imageCopy->SetPixel(m,RGB(255,255,255));
	//    imageCopy->SetPixel(m,255);
      // If the ray has crossed the isosurface, then terminate the loop
     /* if (opacity > opacityThreshold)
        {
        break;
        }*/

      lastT = t;
      lastOpacity = opacity;

      // Compute the next "t" value that crosses a voxel boundary
      t = 1.0;
      for (int k = 0; k < 3; k++)
        {
        // Skip dimension "k" if it is perpendicular to ray
        if (fabs(x2[k] - x1[k]) > VTKCELLPICKER_VOXEL_TOL*rayLength)
          {
          // Compute the previous coord along dimension "k"
          double lastX = x1[k]*(1.0 - lastT) + x2[k]*lastT;

          // Increment to next slice boundary along dimension "k",
          // including a tolerance value for stabilityin cases
          // where lastX is just less than an integer value.
          double nextX = 0;
          if (x2[k] > x1[k])
            {
            nextX = floor(lastX + VTKCELLPICKER_VOXEL_TOL) + 1;
            }
          else
            {
            nextX = ceil(lastX - VTKCELLPICKER_VOXEL_TOL) - 1;
            }

          // Compute the "t" value for this slice boundary
          double ttry = lastT + (nextX - lastX)/(x2[k] - x1[k]);
          if (ttry > lastT + tTol && ttry < t)
            {
            t = ttry;
            }
          }
      //  } 
      } // End of "while (t <= t2)"
	}
	index1++;
	}
	}
	}

//	 InputImageType::IndexType m;
//	   m.m_Index[0] = sourcePnt[0];
//	  m.m_Index[1] = sourcePnt[1];
//	  m.m_Index[2] = sourcePnt[2];
//	  image->SetPixel(m,255);
//	
//	int max1 = maximumValue(rayOpacity);
////	if(rayOpacity[max1]>lastMaxOpacity)
//	{
//		
//	
//	rayOpacity[max1] = 0;
//	int max2 = maximumValue(rayOpacity);
//	rayOpacity[max2] = 0;
//	int max3 = maximumValue(rayOpacity);
//	sphere1->GetOutput()->GetPoint(Indices[max1], destinPnt);
//
//	 destinPnt[0] = spherecenter[0] + destinPnt[0];
//     destinPnt[1] = spherecenter[1] + destinPnt[1];
//	 destinPnt[2] = spherecenter[2] + destinPnt[2];
//	 
//	
//	   m.m_Index[0] = destinPnt[0];
//	  m.m_Index[1] = destinPnt[1];
//	   m.m_Index[2] = destinPnt[2];
//	  image->SetPixel(m,255);
//	  outCenter[0].m_Index[0] =  m.m_Index[0];
//	 outCenter[0].m_Index[1] =  m.m_Index[1];
//	  outCenter[0].m_Index[2] =  m.m_Index[2];
//	  sphere1->GetOutput()->GetPoint(Indices[max2], destinPnt);
//
//	 destinPnt[0] = spherecenter[0] + destinPnt[0];
//     destinPnt[1] = spherecenter[1] + destinPnt[1];
//	 destinPnt[2] = spherecenter[2] + destinPnt[2];
//	
//	  m.m_Index[0] = destinPnt[0];
//	  m.m_Index[1] = destinPnt[1];
//	  m.m_Index[2] = destinPnt[2];
//	  image->SetPixel(m,255);
//	 // lastMaxOpacity = rayOpacity[max1];
//
//	}


	ConnectedComponent* ccl = new ConnectedComponent();
//	RGBImageType::Pointer intimage =  ccl->filterImage(ima);
	RGBImageType::Pointer intimage =  ccl->itkRelabelComponentImageFilterTest(ima);


	for(unsigned int r = 0; r < sphere1->GetRadius()+1; r++)
	{
	for(unsigned int c = 0; c < sphere1->GetRadius()+1; c++)
	{
	for(unsigned int k = 0; k < sphere1->GetRadius()+1; k++)
	{
		RGBImageType::IndexType pixelIndex;
		pixelIndex[0] = r;
		pixelIndex[1] = c;
		pixelIndex[2] = k;
		RGBImageType::SizeType size = intimage->GetLargestPossibleRegion().GetSize();
		//if((intimage->GetPixel(pixelIndex))>0)
		InputImageType::IndexType p;
		std::cout<<r<<c<<(int)intimage->GetPixel(pixelIndex).GetBlue()<<intimage->GetPixel(pixelIndex)<<"\n";
		p[0]=pixelIndex[0]+spherecenter[0]-sphere1->GetRadius();
		p[1]=pixelIndex[1]+spherecenter[1]-sphere1->GetRadius();
		p[2]=pixelIndex[2]+spherecenter[2]-sphere1->GetRadius();
		image->SetPixel(p,intimage->GetPixel(pixelIndex).GetBlue());
	}
	}
	}
	
	InputImageType::IndexType s,m;
	   s.m_Index[0] = sourcePnt[0];
	  s.m_Index[1] = sourcePnt[1];
	  s.m_Index[2] = sourcePnt[2];

	 
   if( !out )
   {
	   std::cout << "Couldn't open file."  << std::endl;
      return 1;
   }
	out<<"slice "<<spherecenter<<" ";
	double ro = 40;double landa = 0.1;double r_j = 0;
	for(int i=0;i<index1;i++)
	{
		sphere1->GetOutput()->GetPoint(Indices[i], destinPnt);

	destinPnt[0] = vtkMath::Round( destinPnt[0] )+spherecenter.m_Index[0];
	destinPnt[1] = vtkMath::Round( destinPnt[1]) +spherecenter.m_Index[1];
	destinPnt[2] = vtkMath::Round(destinPnt[2]) + spherecenter.m_Index[2];
		m.m_Index[0] = destinPnt[0];
		m.m_Index[1] = destinPnt[1];
		m.m_Index[2] = destinPnt[2];
		if(imageCopy->GetPixel(m)==0)
		{
			/*Relation[i]=0;
			image->SetPixel(m,255);
			outCenter[0].m_Index[0] =  m.m_Index[0];
			outCenter[0].m_Index[1] =  m.m_Index[1];
			outCenter[0].m_Index[2] =  m.m_Index[2];
			return 1;*/
			Relation[i] = 0;
		}
		/*else if(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)==0&&imageCopy->GetPixel(m)==0)
			Relation[i]=0;*/
		else
		{
			std::cout<<"angle = "<</*(angle_max[i]/vnl_math::pi)*180<<*/" "<<rayOpacity[i]<<" "<<imageCopy->GetPixel(m)-imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s))<<"\n";
			Relation[i] = FindNeighborIntensity(m);//abs( imageCopy->GetPixel(m));//abs( rayOpacity[i]*(imageCopy->GetPixel(m)));
			//sphereImage->SetPixel(m,(int)Relation[i]*rayOpacity[i]);
			N_j[i] = ConnectedComponentLabeling(m);
			 r_j =(float) N_j[i]/index1;
			 if(r_j!=0&&N_j[i]>2)
				F_j[i] = 1-(1/(1+exp(-ro*(r_j-landa))));	
			 else
				 F_j[i] = 0;
		}
	}
	out<<"\n";


	///////////////////////////////////////////////////////////////////////////////////////////////

	//typedef itk::ConnectedComponentImageFilter <SphereImageType, SphereImageType >
	//ConnectedComponentImageFilterType;

	//ConnectedComponentImageFilterType::Pointer labelFilter
	//= ConnectedComponentImageFilterType::New ();
	//labelFilter->SetInput(sphereImage);
	//labelFilter->Update();

	//typedef itk::RescaleIntensityImageFilter< SphereImageType, SphereImageType > RescaleFilterType;
	//RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
	//rescaleFilter->SetOutputMinimum(0);
	//rescaleFilter->SetOutputMaximum(255);
	//rescaleFilter->SetInput(labelFilter->GetOutput());
 //typedef itk::ImageRegionConstIterator< SphereImageType > ConstIteratorType;
	//ConstIteratorType inputIt( rescaleFilter->GetOutput(), rescaleFilter->GetOutput()->GetLargestPossibleRegion() );
	////IteratorType outputIt( outputImage, outputRegion );
	//for ( inputIt.GoToBegin(); !inputIt.IsAtEnd();
	//++inputIt)
	//{
	//	out<<inputIt.Get()<<"\n";
	//}

	//typedef  itk::ImageFileWriter<  SphereImageType  > WriterType2;
	// WriterType2::Pointer writer2 = WriterType2::New();
	// writer2->SetFileName( "g:\\CTImage\\component.mhd" );
	// writer2->SetInput(rescaleFilter->GetOutput());
	 
 /* try
    {
   writer2->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }*/

	//////////////////////////////////////////////////////////////////////////////////////////////

	//  image->SetPixel(m,255);
	
	//int max1 = maximumValue(rayOpacity);
int  max1 = maximumValue(Relation,2,queue,index1);
//	if(rayOpacity[max1]>lastMaxOpacity)
	//{
if(Relation[max1] == 0)
	return 0;
	int max2;
	int max3;
	int max4;

	/*rayOpacity[max3] = 0;
	max4 = maximumValue(rayOpacity);*/
	double min = 0;
	sphere1->GetOutput()->GetPoint(Indices[max1], destinPnt);

	destinPnt[0] = vtkMath::Round( destinPnt[0] )+spherecenter.m_Index[0];
	destinPnt[1] = vtkMath::Round( destinPnt[1]) +spherecenter.m_Index[1];
	destinPnt[2] = vtkMath::Round(destinPnt[2]) + spherecenter.m_Index[2];
	
	InputImageType::IndexType NewCenter;
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	 // image->SetPixel(m,255);
	  BOOL flag1=false;
	  BOOL flag2 = false;
	  flag1 = Find(m);
	if(flag1)
	 {
		 rayOpacity[max1] = 0;
		 Relation[max1] = 0;
		//max2 = maximumValue(rayOpacity);
		 max2 = maximumValue(Relation,2,queue,index1);
		 max1=max2;

	
	 }
	else{ 
		 if(/*abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&*/abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
	 {
		 std::cout<</*(angle_max[max1]/vnl_math::pi)*180<<*/"  x="<<destinPnt[0]<<"y = "<<destinPnt[1]<<" z="<<destinPnt[2]<<"\n";
		 	 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			   std::cout<<"I(p)-I(x)0"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			  image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  writingFile();
			 outCenter[0].m_Index[0] =  m.m_Index[0];
			outCenter[0].m_Index[1] =  m.m_Index[1];
			outCenter[0].m_Index[2] =  m.m_Index[2];
			
			while(imageCopy->GetPixel(m)>10)
											
											
											{
												std ::cout<<imageCopy->GetPixel(m);
												NewCenter = Circle_RayCasting_Loop(m,-1*vnl_math::pi/2,1*vnl_math::pi/2,2,image,queue);
												{
												m.m_Index[0] = NewCenter.m_Index[0];
												m.m_Index[1] = NewCenter.m_Index[1];
												m.m_Index[2] = NewCenter.m_Index[2];
												//rad =2;
												
												}
												


											}
			out<<"1 "<<Relation[max1]<<" ";
			//Circle_RayCasting_Loop
			
			// if(imageCopy->GetPixel(m)>10)
			 {
				//addtoqueue(Relation,2,queue,index1,max1);
				 queue[index_queue++] = m;
				return 1;//angle_max[max2];
			}
			//else
			//	return 0;
			 
	 }
	 else
	 {
		
		// Centerline[numOfCenterlinePoint].m_Index[0]++;
		//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
		// return 0;
	  }
	}

	//rayOpacity[max1] = 0;
	//Relation[max1] = 0;
	//max2 = maximumValue(rayOpacity);
	max2 = maximumValue(Relation,2,queue,index1);
	max1=max2;
	sphere1->GetOutput()->GetPoint(Indices[max1], destinPnt);

	destinPnt[0] = vtkMath::Round( destinPnt[0] )+spherecenter.m_Index[0];
	destinPnt[1] = vtkMath::Round( destinPnt[1]) +spherecenter.m_Index[1];
	destinPnt[2] = vtkMath::Round(destinPnt[2]) + spherecenter.m_Index[2];
	 
	
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	  //image->SetPixel(m,255);

	
	
		
	if(flag1 &&Find(m))
	{
		std::cout<</*(angle_max[max1]/vnl_math::pi)*180<<*/" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max2]<<"\n";
		rayOpacity[max2] = 0;
		Relation[max2] = 0;
		 max3 = maximumValue(Relation,2,queue,index1);
		 max1=max3;
	
	}
	else if(flag1&&!Find(m))
	{
		//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
		{
			std::cout<</*(angle_max[max1]/vnl_math::pi)*180<<*/"  x="<<destinPnt[0]<<"y = "<<destinPnt[1]<<" z="<<destinPnt[2]<<"\n";
			 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			 std::cout<<"I(p)-I(x)1"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 writingFile();
			 outCenter[0].m_Index[0] =  m.m_Index[0];
			outCenter[0].m_Index[1] =  m.m_Index[1];
			outCenter[0].m_Index[2] =  m.m_Index[2];
			out<<"2 "<<Relation[max1]<<" ";
			while(imageCopy->GetPixel(m)>10)
											
											
											{
												std ::cout<<imageCopy->GetPixel(m);
												NewCenter = Circle_RayCasting_Loop(m,-1*vnl_math::pi/2,1*vnl_math::pi/2,2,image,queue);
												{
												m.m_Index[0] = NewCenter.m_Index[0];
												m.m_Index[1] = NewCenter.m_Index[1];
												m.m_Index[2] = NewCenter.m_Index[2];
												//rad =2;
												
												}
												


											}
		/*	if(imageCopy->GetPixel(m)>10)*/
			{
			//	addtoqueue(Relation,2,queue,index1,max1);
				 queue[index_queue++] = m;
				return 1;//angle_max[max2];
			}
			/*else
				return 0;*/
		}
//	else
		{
			//Centerline[numOfCenterlinePoint].m_Index[0]++;
			//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
			//return 0;
		}
			 //return angle_max[max1];
	}
	//rayOpacity[max2] = 0;
	//Relation[max2] = 0;
	//max3 = maximumValue(rayOpacity);
	max3 = maximumValue(Relation,2,queue,index1);
	max1=max3;
	 sphere1->GetOutput()->GetPoint(Indices[max1], destinPnt);

	destinPnt[0] = vtkMath::Round( destinPnt[0] )+spherecenter.m_Index[0];
	destinPnt[1] = vtkMath::Round( destinPnt[1]) +spherecenter.m_Index[1];
	destinPnt[2] = vtkMath::Round(destinPnt[2]) + spherecenter.m_Index[2];

	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	 
	  flag2 = Find(m);
	  if(flag1 &&flag2)
	{
		rayOpacity[max3] = 0;
		Relation[max3] = 0;
	  // max4 = maximumValue(rayOpacity);
		max4 = maximumValue(Relation,2,queue,index1);
		 max1=max4;
		 for(int i=0;i<index1;i++)
		 {
			sphere1->GetOutput()->GetPoint(Indices[max1], destinPnt);

			destinPnt[0] = vtkMath::Round( destinPnt[0] )+spherecenter.m_Index[0];
			destinPnt[1] = vtkMath::Round( destinPnt[1]) +spherecenter.m_Index[1];
			destinPnt[2] = vtkMath::Round(destinPnt[2]) + spherecenter.m_Index[2];

			m.m_Index[0] = destinPnt[0];
			m.m_Index[1] = destinPnt[1];
			m.m_Index[2] = destinPnt[2];
			if(!Find(m))
			{
				//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
				{
					std::cout<</*(angle_max[max1]/vnl_math::pi)*180<<*/"  x="<<destinPnt[0]<<"y = "<<destinPnt[1]<<" z="<<destinPnt[2]<<"\n";
				Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
				Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
				Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
				 std::cout<<"I(p)-I(x) 2"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
				 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				 //vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				 writingFile();
				 outCenter[0].m_Index[0] =  m.m_Index[0];
				outCenter[0].m_Index[1] =  m.m_Index[1];
				outCenter[0].m_Index[2] =  m.m_Index[2];
				out<<"4 "<<Relation[max1]<<" "<<i<<" ";
					while(imageCopy->GetPixel(m)>10)
											
											//while(center.m_Index[2]>50)
											{
												std ::cout<<imageCopy->GetPixel(m);
												NewCenter = Circle_RayCasting_Loop(m,-1*vnl_math::pi/2,1*vnl_math::pi/2,2,image,queue);
												{
												m.m_Index[0] = NewCenter.m_Index[0];
												m.m_Index[1] = NewCenter.m_Index[1];
												m.m_Index[2] = NewCenter.m_Index[2];
												//rad =2;
												
												}
												


											}
			/*	if(imageCopy->GetPixel(m)>10)*/
				{
			//	addtoqueue(Relation,2,queue,index1,max1);
					 queue[index_queue++] = m;
				return 1;//angle_max[max2];
				}
			/*else
				return 0;*/
				 
				// break;
				}
//			else
				{
				//Centerline[numOfCenterlinePoint].m_Index[0]++;
				//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
				//return 0;
				//min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				//int ind = 0;
				//for(int j=0;j<index;j++)
				//{
				//std::cout<<angle_max[max4]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max4]<<"\n";
				//m.m_Index[0] = spherecenter[0] + (radius *vcl_cos(angle_max[j]));
				//m.m_Index[1] = spherecenter[1] +(radius *vcl_sin(angle_max[j]));
				//m.m_Index[2] = spherecenter[2];
				//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<min)
				//{
				//min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				//ind = j;
				////break;
				//}


				//}
				//return angle_max[ind];
				}
				 //return angle_max[max1];
			}
			else
			{
				rayOpacity[max1] = 0;
				//max1 = maximumValue(rayOpacity);
				Relation[max1] = 0;
				max1 = maximumValue(Relation,2,queue,index1);
			}
	  }
		 return 0;

	}
	else if(flag1&&!flag2)
	{
		//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
		{
			std::cout<</*(angle_max[max1]/vnl_math::pi)*180<<*/"  x="<<destinPnt[0]<<"y = "<<destinPnt[1]<<" z="<<destinPnt[2]<<"\n";
			 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			 std::cout<<"I(p)-I(x) 3"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			 writingFile();
			 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  writingFile();
			 outCenter[0].m_Index[0] =  m.m_Index[0];
			outCenter[0].m_Index[1] =  m.m_Index[1];
			outCenter[0].m_Index[2] =  m.m_Index[2];
			out<<"3 "<<Relation[max1]<<" ";
			while(imageCopy->GetPixel(m)>10)
											
											//while(center.m_Index[2]>50)
											{
												std ::cout<<imageCopy->GetPixel(m);
												NewCenter = Circle_RayCasting_Loop(m,-1*vnl_math::pi/2,1*vnl_math::pi/2,2,image,queue);
												{
												m.m_Index[0] = NewCenter.m_Index[0];
												m.m_Index[1] = NewCenter.m_Index[1];
												m.m_Index[2] = NewCenter.m_Index[2];
												//rad =2;
												
												}
												


											}
			//if(imageCopy->GetPixel(m)>10)
			{
				//addtoqueue(Relation,2,queue,index1,max1);
				 queue[index_queue++] = m;
				return 1;//angle_max[max2];
			}
			/*else
				return 0;*/
		}
		/*else
		{
		*/		/*for(int i=0;i<index;i++)
				{
				m.m_Index[0] = spherecenter[0] + (radius *vcl_cos(angle_max[i]));
				m.m_Index[1] = spherecenter[1] +(radius *vcl_sin(angle_max[i]));
				m.m_Index[2] = spherecenter[2];
				if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<min)
				{
				min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				break;
				}*/


				//}
			//Centerline[numOfCenterlinePoint].m_Index[0]++;
			//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
			//return 0;
		//}
			  
	}

	  //lastMaxOpacity = rayOpacity[max1];

	//}

	 typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;
	 WriterType::Pointer writer = WriterType::New();
	 writer->SetFileName( "g:\\CTImage\\Sphere.mhd" );
	 writer->SetInput(image);
	 
  try
    {
   writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
	return 1;
}



void VTKRayCaster:: writingFile()
{
typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	 WriterType::Pointer writer = WriterType::New();
	 writer->SetFileName( "g:\\CTImage\\Sphere.mhd" );
	 writer->SetInput(image);
	 
  try
    {
   writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
 

}


double VTKRayCaster::Circle_RayCasting(InputImageType::IndexType m_Index,double angle_start,double angle_end,double radius,InputImageType::Pointer vesselImage)
{

	//itkExporter->SetInput(vesselImage);
	///*vtkImageImport**/ vtkImporter = vtkImageImport::New();
	//itkExporter->Update();
	//ConnectPipelines(itkExporter, vtkImporter);
	//vtkImporter->Update();
	//imageData = vtkImageData::New();
	//imageData = vtkImporter->GetOutput();
	 InputImageType::IndexType m;

double angle_max[100];


for(int i=0;i<100;i++)
	angle_max[i] = 0;

 //InputImageType::Pointer imageCopy = image;

 //typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	// WriterType::Pointer writer = WriterType::New();
	// writer->SetFileName( "g:\\CTImage\\Sphere.mhd" );
	// writer->SetInput(image);
	// 
 // try
 //   {
 //  writer->Update();
 //   }
 // catch( itk::ExceptionObject & excep )
 //   {
 //   std::cerr << "Exception caught !" << std::endl;
 //   std::cerr << excep << std::endl;
 //   }


 

 double spherecenter[3];
 spherecenter[0] = m_Index[0];
 spherecenter[1] = m_Index[1];
 spherecenter[2] = m_Index[2];
  

  // init the counter and ray length
  int numIntersected = 0;
  double rayLen = 10.0000001; // = 1 - 0.8 + error tolerance
  int sub_id;
  vtkIdType cell_id;
  double param_t, intersect[3], paraCoord[3];
  double sourcePnt[3], destinPnt[3], normalVec[3];
//  vtkGenericCell *cell = vtkGenericCell::New();

// Clip the ray with the extent, results go in s1 and s2
  int planeId;
  



//  vtkGPUVolumeRayCastMapper* volumeMapper = vtkGPUVolumeRayCastMapper::New();
//  volumeMapper->SetBlendModeToComposite(); // composite first
//#if VTK_MAJOR_VERSION <= 5
//  volumeMapper->SetInputConnection(imageData->GetProducerPort());
//#else
//  volumeMapper->SetInputData(imageData);
//#endif  
//  vtkVolumeProperty* volumeProperty = vtkVolumeProperty::New();
//  volumeProperty->ShadeOff();
//  volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);
//
// vtkPiecewiseFunction* compositeOpacity = vtkPiecewiseFunction::New();
//	  compositeOpacity->AddPoint(-3024, 0, 0.5, 0.0 );
//      compositeOpacity->AddPoint(-155, 0, 0.5, 0.92 );
//      compositeOpacity->AddPoint(217, .68, 0.33, 0.45 );
//      compositeOpacity->AddPoint(420,.83, 0.5, 0.0);
//      compositeOpacity->AddPoint(3071, .80, 0.5, 0.0);
//	 volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.
// 
//	vtkColorTransferFunction* color = vtkColorTransferFunction::New();
//	color->AddRGBPoint( -3024, 0, 0, 0, 0.5, 0.0 );
//	color->AddRGBPoint( -155, .55, .25, .15, 0.5, .92 );
//	color->AddRGBPoint( 217, .88, .60, .29, 0.33, 0.45 );
//	color->AddRGBPoint( 420, 1, .94, .95, 0.5, 0.0 );
//	color->AddRGBPoint( 3071, .83, .66, 1, 0.5, 0.0 );
//	volumeProperty->SetColor(color);
//	volumeMapper->SetBlendModeToComposite();
//	volumeProperty->ShadeOn();
//	volumeProperty->SetAmbient(0.1);
//	volumeProperty->SetDiffuse(0.9);
//	volumeProperty->SetSpecular(0.2);
//	volumeProperty->SetSpecularPower(10.0);
//	volumeProperty->SetScalarOpacityUnitDistance(0.8919);
//	vtkVolume* volume = vtkVolume::New();
//	volume->SetMapper(volumeMapper);
//	volume->SetProperty(volumeProperty);
//	
//
//
//  
//	double frameRate = 10.0;
//
//	int clip = 1;
//
//	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
//   vtkRenderer *renderer = vtkRenderer::New();
  // Add a box widget if the clip option was selected
//  vtkBoxWidget *box = vtkBoxWidget::New();
//  if (clip)
//    {
//    box->SetInteractor(iren);
//    box->SetPlaceFactor(1.01);
//   // if ( reductionFactor < 1.0 )
//  //    {      
//		  box->SetInput(imageData);
//   //   }
//   // else
//   /*   {
//		  box->SetInput(imageData);
//      }*/
//    
//    box->SetDefaultRenderer(renderer);
//    box->InsideOutOn();
//    box->PlaceWidget();
//    vtkBoxWidgetCallback *callback = vtkBoxWidgetCallback::New();
//	callback->SetMapper(volumeMapper);
//    box->AddObserver(vtkCommand::InteractionEvent, callback);
//    callback->Delete();
//   // box->EnabledOn();
////    box->GetSelectedFaceProperty()->SetOpacity(0.0);
//  }
//
//  // Create the renderer, render window and interactor
// 
//
//  vtkRenderWindow *renWin = vtkRenderWindow::New();
//  renWin->AddRenderer(renderer);
//
//
//  iren->SetRenderWindow(renWin);
//  iren->SetDesiredUpdateRate(frameRate / (1+clip) );
//  
//  iren->GetInteractorStyle()->SetDefaultRenderer(renderer);



 /* renWin->SetSize(600,600);
  renWin->Render();*/

 /* if ( !volumeMapper->IsRenderSupported(renWin, volumeProperty) )
    {
      cout << "This mapper is unsupported on this platform" << endl;
      exit(EXIT_FAILURE);
    }*/
  

  // Add the volume to the scene
 // renderer->AddVolume( volume );

 // renderer->ResetCamera();

 // // interact with data
 //renWin->Render();

 // iren->Start();
  

		/*vtkPlanes *Planes = vtkPlanes::New();
        box->GetPlanes(Planes);
        volumeMapper->SetClippingPlanes(Planes);*/

//writeVTKImage();


	//double radius = 10;
  vtkPlaneCollection *planes = 0;
   int clippingPlaneId = -1;
    planes = volumeMapper->GetClippingPlanes();
 /*  if (volumeMapper && (planes = volumeMapper->GetClippingPlanes())
      && (planes->GetNumberOfItems() > 0))
    {*/
    // This is a bit ugly: need to transform back to world coordinates
    double q1[3], q2[3];
//    this->Transform->TransformPoint(sourcePnt, q1);
   // this->Transform->TransformPoint(destinPnt, q2);

	destinPnt[0] = 0;
	destinPnt[1] = 0;
	destinPnt[2] = 0;
	sourcePnt[0] = spherecenter[0];
	sourcePnt[1] = spherecenter[1];
	sourcePnt[2] = spherecenter[2];
	int index2 = 0;
	 double spacing[3], origin[3];
	 int extent[6];
  imageData->GetSpacing(spacing);
  imageData->GetOrigin(origin);
  imageData->GetExtent(extent);
  for(double angle = angle_start;angle <= angle_end/*1*vnl_math::pi*/; angle += vnl_math::pi/180.0 )
      {

      destinPnt[0] = (radius *vcl_cos(angle));
      destinPnt[1] = (radius *vcl_sin(angle));
	  destinPnt[2] = 0;
	  if((lastdistinPnt[0]!=vtkMath::Round(destinPnt[0]))||(lastdistinPnt[1]!=vtkMath::Round( destinPnt[1]))||(lastdistinPnt[2]!=vtkMath::Round( destinPnt[2])))
	{
	lastdistinPnt[0] = vtkMath::Round(destinPnt[0]);
	lastdistinPnt[1] = vtkMath::Round(destinPnt[1]);
	lastdistinPnt[2] = vtkMath::Round(destinPnt[2]);
	
 
    // cast a ray in the negative direction toward sphere1
   /* sourcePnt[0] = destinPnt[0] - rayLen * normalVec[0];
    sourcePnt[1] = destinPnt[1] - rayLen * normalVec[1];
    sourcePnt[2] = destinPnt[2] - rayLen * normalVec[2];*/
 
	

	destinPnt[0] = destinPnt[0] +spherecenter[0];
	destinPnt[1] = destinPnt[1] +spherecenter[1];
	destinPnt[2] = destinPnt[2] +spherecenter[2];
 
	rayOpacity[index2]=0;
	



	// imageData = vtkImageData::New();
	
	
//  
// /*vtkMetaImageReader *reader = vtkMetaImageReader::New();
// reader->SetFileName("g:\\CTImage\\downsampledimage1.mhd");
//    reader->Update();*/
//   // imageData=reader->GetOutput();
//	imageData = vtkImporter->GetOutput();
//	//imageDataCopy = reader->GetOutputPort();
//	imageDataCopy = vtkImporter->GetOutputPort();
//   

//  int boxExtent[6];
//  boxExtent[0] = sourcePnt[0];
//  boxExtent[1] = sourcePnt[1];
//  boxExtent[2] = sourcePnt[2];
//  boxExtent[3] = destinPnt[0];
//  boxExtent[4] = destinPnt[1];
//  boxExtent[5] = destinPnt[2];

//
//  
//    
//
 double t1 = 0;
	double t2 = 1;


	/*m.m_Index[0] =(int) destinPnt[0];
	m.m_Index[1] =(int) destinPnt[1];
	m.m_Index[2] = (int)destinPnt[2];
	image->SetPixel(m,255);*/

  double x1[3], x2[3];
  for (int i = 0; i < 3; i++)
    {
		x1[i] = (sourcePnt[i]); /*- origin[i])/spacing[i];*/
		x2[i] = (destinPnt[i]);/*- origin[i])/spacing[i];*/
    }
// /* if (!this->ClipLineWithPlanes(planes, sourcePnt, destinPnt, t1, t2, clippingPlaneId))
//      {
//      return VTK_DOUBLE_MAX;
//      }*/
//   
//
//
//  volumeMapper->GetClippingPlanes();
// 
//double s1, s2;
//
////boxExtent[0] = x1[0];
////boxExtent[1] = x1[1];
////boxExtent[2] = x1[2];
////boxExtent[3] = x2[0];
////boxExtent[4] = x2[1];
////boxExtent[5] = x2[2];
////if (!this->ClipLineWithExtent(extent, x1, x2, s1, s2, planeId))
////    {
////    return VTK_DOUBLE_MAX;
////    }
////  if (s1 >= t1) { t1 = s1; }
////  if (s2 <= t2) { t2 = s2; }
////
////  // Sanity check
////  if (t2 < t1)
////    {
////    return VTK_DOUBLE_MAX;
////    }
//      
  // Compute the length of the line intersecting the volume
  double rayLength = sqrt(vtkMath::Distance2BetweenPoints(x1, x2))*(t2 - t1); 

//  // This is the minimum increment that will be allowed
  double tTol = VTKCELLPICKER_VOXEL_TOL/rayLength*(t2 - t1);

  // Find out whether there are multiple components in the volume
  int numComponents = imageData->GetNumberOfScalarComponents();
  int independentComponents = volumeProperty->GetIndependentComponents();
  int numIndependentComponents = 1;
  if (independentComponents)
    {
    numIndependentComponents = numComponents;
    }

  // Create a scalar array, it will be needed later
  vtkDataArray *scalars = vtkDataArray::CreateDataArray(imageData->GetScalarType());
  scalars->SetNumberOfComponents(numComponents);
  vtkIdType scalarArraySize = numComponents*imageData->GetNumberOfPoints();
  int scalarSize = imageData->GetScalarSize();
  void *scalarPtr = imageData->GetScalarPointer();

  // Go through each volume component separately
  double tMin = VTK_DOUBLE_MAX;
  for (int component = 0; component < numIndependentComponents; component++)
    { 
    scalarOpacity = volumeProperty->GetScalarOpacity(component);
    int disableGradientOpacity = volumeProperty->GetDisableGradientOpacity(component);
    gradientOpacity = 0;
   // if (!disableGradientOpacity && this->UseVolumeGradientOpacity)
      {
		  gradientOpacity = volumeProperty->GetGradientOpacity(component);
      }
//
    // This is the component used to compute the opacity
    int oComponent = component;
    if (!independentComponents)
      {
      oComponent = numComponents - 1;
      }

    // Make a new array, shifted to the desired component
    scalars->SetVoidArray(static_cast<void *>(static_cast<float *>(scalarPtr) + scalarSize*oComponent), scalarArraySize, 1);

     
  
  //else if ( (volumeMapper = vtkAbstractVolumeMapper::SafeDownCast(m)) )
  //  {
  //  tMin = this->IntersectVolumeWithLine(p1, p2, t1, t2, prop, volumeMapper);
  //  }

//
	 // Do a ray cast with linear interpolation.
    double opacity = 0.0;
    double lastOpacity = 0.0;
    double lastT = t1;
    double x[3];
    double pcoords[3];
    int xi[3];
	int index3 = 0;

    // Ray cast loop
    double t = t1;
	//rayOpacity[index] = 0;
    while (t < t2)
      {
      for (int j = 0; j < 3; j++)
        {
        // "t" is the fractional distance between endpoints x1 and x2
        x[j] = x1[j]*(1.0 - t) + x2[j]*t;

        // Paranoia bounds check
        if (x[j] < extent[2*j]) { x[j] = extent[2*j]; }
        else if (x[j] > extent[2*j + 1]) { x[j] = extent[2*j+1]; 
		}
	/*	 if (x[j] > x2[j])
		 {
			{ x[j] = x2[j]; }
			t=2;
		}*/
//       
        xi[j] = int(floor(x[j]));
        pcoords[j] = x[j] - xi[j];
        }

	  opacity = this->ComputeVolumeOpacity(xi, pcoords, imageData, scalars,
                                           scalarOpacity, gradientOpacity);
//	  std::cout<<"opacity \n"<<opacity;
	  rayOpacity[index2]+=opacity;
	  angle_max[index2] = angle;

	  
//	  
//	 
//      // If the ray has crossed the isosurface, then terminate the loop
//     /* if (opacity > opacityThreshold)
//        {
//        break;
//        }*/
//
      lastT = t;
//      lastOpacity = opacity;
//
//      // Compute the next "t" value that crosses a voxel boundary
      t = 1.0;
      for (int k = 0; k < 3; k++)
        {
        // Skip dimension "k" if it is perpendicular to ray
        if (fabs(x2[k] - x1[k]) > VTKCELLPICKER_VOXEL_TOL*rayLength)
          {
          // Compute the previous coord along dimension "k"
          double lastX = x1[k]*(1.0 - lastT) + x2[k]*lastT;

          // Increment to next slice boundary along dimension "k",
          // including a tolerance value for stabilityin cases
          // where lastX is just less than an integer value.
          double nextX = 0;
          if (x2[k] > x1[k])
            {
            nextX = floor(lastX + VTKCELLPICKER_VOXEL_TOL) + 1;
            }
          else
            {
            nextX = ceil(lastX - VTKCELLPICKER_VOXEL_TOL) - 1;
            }

          // Compute the "t" value for this slice boundary
          double ttry = lastT + (nextX - lastX)/(x2[k] - x1[k]);
          if (ttry > lastT + tTol && ttry < t)
            {
            t = ttry;
            }
          }
        } 
      } // End of "while (t <= t2)"
	
	}
	index2++;
	}
	}
	InputImageType::IndexType s;
	   s.m_Index[0] = sourcePnt[0];
	  s.m_Index[1] = sourcePnt[1];
	  s.m_Index[2] = sourcePnt[2];

	for(int i=0;i<index2;i++)
	{
		destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[i]));
		destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[i]));
		destinPnt[2] = spherecenter[2];
		m.m_Index[0] = destinPnt[0];
		m.m_Index[1] = destinPnt[1];
		m.m_Index[2] = destinPnt[2];
		if(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)==0&&imageCopy->GetPixel(m)>0)
		{
			Relation[i]=0;
		//	return 0;
			if(!Find(m))
			{
				std::cout<<imageCopy->GetPixel(m)<<"\n";
				Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
				Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
				Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
				image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				writingFile();
				return angle_max[i];
			}
		}
		else if(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)==0&&imageCopy->GetPixel(m)==0)
			Relation[i]=0;
		else
		{
		std::cout<</*"angle = "<<(angle_max[i]/vnl_math::pi)*180<<*/"ray Opacity = "<<rayOpacity[i]<<" m = "<<imageCopy->GetPixel(m)<<" s = "<<imageCopy->GetPixel(s)<<"difference "<<imageCopy->GetPixel(m)-imageCopy->GetPixel(s)<<" "<<rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s))<<"\n";
		Relation[i] =abs( rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)));
		}
		/*std::cout<<"x = "<<m.m_Index[0]<<" y =  "<<m.m_Index[1]<<"\n";
		std::cout<<"angle = "<<angle_max[i]<<" "<<rayOpacity[i]<<" "<<imageCopy->GetPixel(m)-imageCopy->GetPixel(s)<<" "<<rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s))<<"\n";
		Relation[i] =abs( rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)));*/
	}
	

	//  image->SetPixel(m,255);
	
	//int max1 = maximumValue(rayOpacity);
	int max1 = maximumValue_2(Relation,index2);
//	if(rayOpacity[max1]>lastMaxOpacity)
	//{
		
	int max2;
	int max3;
	int max4;

	/*rayOpacity[max3] = 0;
	max4 = maximumValue(rayOpacity);*/
	double min = 0;
	 destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max1]));
     destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max1]));
	 destinPnt[2] = spherecenter[2];
	 
	
	
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	  image->SetPixel(m,255);
	  BOOL flag1=false;
	  BOOL flag2 = false;
	  flag1 = Find(m);
	if(flag1)
	 {
		 rayOpacity[max1] = 0;
		 Relation[max1] = 0;
		//max2 = maximumValue(rayOpacity);
		 max2 = maximumValue_2(Relation,index2);
		 max1=max2;

	
	 }
	else{ 
		 if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
	 {
		 		std::cout<<angle_max[max1]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max1]<<"\n";
		 	 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			   std::cout<<"I(p)-I(x)0"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			  image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			   vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  writingFile();
			  return angle_max[max1];
			 
	 }
	 else
	 {
		
		// Centerline[numOfCenterlinePoint].m_Index[0]++;
		//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
		// return 0;
	  }
	}

	rayOpacity[max1] = 0;
	Relation[max1] = 0;
	//max2 = maximumValue(rayOpacity);
	max2 = maximumValue_2(Relation,index2);
	max1=max2;
	 destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max2]));
     destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max2]));
	 destinPnt[2] = spherecenter[2];
	 
	
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	  //image->SetPixel(m,255);

	
	
		
	if(flag1 &&Find(m))
	{
		std::cout<<angle_max[max2]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max2]<<"\n";
		rayOpacity[max2] = 0;
		Relation[max2] = 0;
		 max3 = maximumValue_2(rayOpacity,index2);
		 max1=max3;
	
	}
	else if(flag1&&!Find(m))
	{
		if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0){
			std::cout<<angle_max[max2]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max2]<<"\n";
			 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			 std::cout<<"I(p)-I(x)1"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 writingFile();
		
			return angle_max[max2];
		}
		else
		{
			//Centerline[numOfCenterlinePoint].m_Index[0]++;
			//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
			//return 0;
		}
			 //return angle_max[max1];
	}
	rayOpacity[max2] = 0;
	Relation[max2] = 0;
	//max3 = maximumValue(rayOpacity);
	max3 = maximumValue_2(Relation,index2);
	max1=max3;
	 destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max3]));
     destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max3]));
	 destinPnt[2] = spherecenter[2];
	 
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	 
	  flag2 = Find(m);
	  if(flag1 &&flag2)
	{
		rayOpacity[max3] = 0;
		Relation[max3] = 0;
	  // max4 = maximumValue(rayOpacity);
		max4 = maximumValue_2(Relation,index2);
		 max1=max4;
		 for(int i=0;i<index2;i++)
		 {
			destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max1]));
			destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max1]));
			destinPnt[2] = spherecenter[2];

			m.m_Index[0] = destinPnt[0];
			m.m_Index[1] = destinPnt[1];
			m.m_Index[2] = destinPnt[2];
			if(!Find(m))
			{
				if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0){
					std::cout<<angle_max[max4]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max4]<<"\n";
				Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
				Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
				Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
				 std::cout<<"I(p)-I(x) 2"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
				 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				 writingFile();
			
				return angle_max[max1];
				 
				// break;
				}
				else
				{
				//Centerline[numOfCenterlinePoint].m_Index[0]++;
				//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
				//return 0;
				min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				//int ind = 0;
				//for(int j=0;j<index;j++)
				//{
				//std::cout<<angle_max[max4]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max4]<<"\n";
				//m.m_Index[0] = spherecenter[0] + (radius *vcl_cos(angle_max[j]));
				//m.m_Index[1] = spherecenter[1] +(radius *vcl_sin(angle_max[j]));
				//m.m_Index[2] = spherecenter[2];
				//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<min)
				//{
				//min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				//ind = j;
				////break;
				//}


				//}
				//return angle_max[ind];
				}
				 //return angle_max[max1];
			}
			else
			{
				rayOpacity[max1] = 0;
				//max1 = maximumValue(rayOpacity);
				Relation[max1] = 0;
				max1 = maximumValue_2(Relation,index2);
			}
	  }

	}
	else if(flag1&&!flag2)
	{
		if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
		{
			std::cout<<angle_max[max3]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max3]<<"\n";
			 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			 std::cout<<"I(p)-I(x) 3"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			 writingFile();
			 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  writingFile();
			 
			 return angle_max[max3];
		}
		else
		{
				for(int i=0;i<index2;i++)
				{
				m.m_Index[0] = spherecenter[0] + (radius *vcl_cos(angle_max[i]));
				m.m_Index[1] = spherecenter[1] +(radius *vcl_sin(angle_max[i]));
				m.m_Index[2] = spherecenter[2];
				if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<min)
				{
				min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				break;
				}


				}
			//Centerline[numOfCenterlinePoint].m_Index[0]++;
			//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
			//return 0;
		//}
			  
	}

	  //lastMaxOpacity = rayOpacity[max1];

	}

	typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	 WriterType::Pointer writer = WriterType::New();
	 writer->SetFileName( "g:\\CTImage\\Sphere.mhd" );
	 writer->SetInput(image);
	 
  try
    {
   writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  return angle_max[max1];
}

InputImageType::IndexType VTKRayCaster::Circle_RayCasting_Loop(InputImageType::IndexType m_Index,double angle_start,double angle_end,double radius,InputImageType::Pointer vesselImage,InputImageType::IndexType* queue)
{

	//itkExporter->SetInput(vesselImage);
	///*vtkImageImport**/ vtkImporter = vtkImageImport::New();
	//itkExporter->Update();
	//ConnectPipelines(itkExporter, vtkImporter);
	//vtkImporter->Update();
	//imageData = vtkImageData::New();
	//imageData = vtkImporter->GetOutput();
	 InputImageType::IndexType m;

double angle_max[100];

for(int i=0;i<100;i++)
	angle_max[i] = 0;



	typedef itk::ImageFileReader<InputImageType> FileReaderType;
		  FileReaderType::Pointer imageReader = FileReaderType::New();
		// imageReader->SetFileName("g:\\RawReader\\output_en2.mhd");
		  imageReader->SetFileName( "F:\\CTImage\\training\\dataset02\\vesselness.mhd");
		   //imageReader->SetFileName("g:\\CTImage\\vesselness.mhd");
		 // imageReader->SetFileName("g:\\CTImage\\FastMarchingFilterOutput2.mhd");
  try
  { 
    imageReader->Update();
  }
  catch (itk::ExceptionObject &ex)
  { 
    std::cout << ex << std::endl;
    //return EXIT_FAILURE;
  }
  this->imageCopy = imageReader->GetOutput();

 //InputImageType::Pointer imageCopy = image;

 //typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	// WriterType::Pointer writer = WriterType::New();
	// writer->SetFileName( "g:\\CTImage\\Sphere.mhd" );
	// writer->SetInput(image);
	// 
 // try
 //   {
 //  writer->Update();
 //   }
 // catch( itk::ExceptionObject & excep )
 //   {
 //   std::cerr << "Exception caught !" << std::endl;
 //   std::cerr << excep << std::endl;
 //   }


 

 double spherecenter[3];
 spherecenter[0] = m_Index[0];
 spherecenter[1] = m_Index[1];
 spherecenter[2] = m_Index[2];
  

  // init the counter and ray length
  int numIntersected = 0;
  double rayLen = 10.0000001; // = 1 - 0.8 + error tolerance
  int sub_id;
  vtkIdType cell_id;
  double param_t, intersect[3], paraCoord[3];
  double sourcePnt[3], destinPnt[3], normalVec[3];
//  vtkGenericCell *cell = vtkGenericCell::New();

// Clip the ray with the extent, results go in s1 and s2
  int planeId;
  



//  vtkGPUVolumeRayCastMapper* volumeMapper = vtkGPUVolumeRayCastMapper::New();
//  volumeMapper->SetBlendModeToComposite(); // composite first
//#if VTK_MAJOR_VERSION <= 5
//  volumeMapper->SetInputConnection(imageData->GetProducerPort());
//#else
//  volumeMapper->SetInputData(imageData);
//#endif  
//  vtkVolumeProperty* volumeProperty = vtkVolumeProperty::New();
//  volumeProperty->ShadeOff();
//  volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);
//
// vtkPiecewiseFunction* compositeOpacity = vtkPiecewiseFunction::New();
//	  compositeOpacity->AddPoint(-3024, 0, 0.5, 0.0 );
//      compositeOpacity->AddPoint(-155, 0, 0.5, 0.92 );
//      compositeOpacity->AddPoint(217, .68, 0.33, 0.45 );
//      compositeOpacity->AddPoint(420,.83, 0.5, 0.0);
//      compositeOpacity->AddPoint(3071, .80, 0.5, 0.0);
//	 volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.
// 
//	vtkColorTransferFunction* color = vtkColorTransferFunction::New();
//	color->AddRGBPoint( -3024, 0, 0, 0, 0.5, 0.0 );
//	color->AddRGBPoint( -155, .55, .25, .15, 0.5, .92 );
//	color->AddRGBPoint( 217, .88, .60, .29, 0.33, 0.45 );
//	color->AddRGBPoint( 420, 1, .94, .95, 0.5, 0.0 );
//	color->AddRGBPoint( 3071, .83, .66, 1, 0.5, 0.0 );
//	volumeProperty->SetColor(color);
//	volumeMapper->SetBlendModeToComposite();
//	volumeProperty->ShadeOn();
//	volumeProperty->SetAmbient(0.1);
//	volumeProperty->SetDiffuse(0.9);
//	volumeProperty->SetSpecular(0.2);
//	volumeProperty->SetSpecularPower(10.0);
//	volumeProperty->SetScalarOpacityUnitDistance(0.8919);
//	vtkVolume* volume = vtkVolume::New();
//	volume->SetMapper(volumeMapper);
//	volume->SetProperty(volumeProperty);
//	
//
//
//  
//	double frameRate = 10.0;
//
//	int clip = 1;
//
//	vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
//   vtkRenderer *renderer = vtkRenderer::New();
  // Add a box widget if the clip option was selected
//  vtkBoxWidget *box = vtkBoxWidget::New();
//  if (clip)
//    {
//    box->SetInteractor(iren);
//    box->SetPlaceFactor(1.01);
//   // if ( reductionFactor < 1.0 )
//  //    {      
//		  box->SetInput(imageData);
//   //   }
//   // else
//   /*   {
//		  box->SetInput(imageData);
//      }*/
//    
//    box->SetDefaultRenderer(renderer);
//    box->InsideOutOn();
//    box->PlaceWidget();
//    vtkBoxWidgetCallback *callback = vtkBoxWidgetCallback::New();
//	callback->SetMapper(volumeMapper);
//    box->AddObserver(vtkCommand::InteractionEvent, callback);
//    callback->Delete();
//   // box->EnabledOn();
////    box->GetSelectedFaceProperty()->SetOpacity(0.0);
//  }
//
//  // Create the renderer, render window and interactor
// 
//
//  vtkRenderWindow *renWin = vtkRenderWindow::New();
//  renWin->AddRenderer(renderer);
//
//
//  iren->SetRenderWindow(renWin);
//  iren->SetDesiredUpdateRate(frameRate / (1+clip) );
//  
//  iren->GetInteractorStyle()->SetDefaultRenderer(renderer);



 /* renWin->SetSize(600,600);
  renWin->Render();*/

 /* if ( !volumeMapper->IsRenderSupported(renWin, volumeProperty) )
    {
      cout << "This mapper is unsupported on this platform" << endl;
      exit(EXIT_FAILURE);
    }*/
  

  // Add the volume to the scene
 // renderer->AddVolume( volume );

 // renderer->ResetCamera();

 // // interact with data
 //renWin->Render();

 // iren->Start();
  

		/*vtkPlanes *Planes = vtkPlanes::New();
        box->GetPlanes(Planes);
        volumeMapper->SetClippingPlanes(Planes);*/

//writeVTKImage();


	//double radius = 10;
  vtkPlaneCollection *planes = 0;
   int clippingPlaneId = -1;
    planes = volumeMapper->GetClippingPlanes();
 /*  if (volumeMapper && (planes = volumeMapper->GetClippingPlanes())
      && (planes->GetNumberOfItems() > 0))
    {*/
    // This is a bit ugly: need to transform back to world coordinates
    double q1[3], q2[3];
//    this->Transform->TransformPoint(sourcePnt, q1);
   // this->Transform->TransformPoint(destinPnt, q2);

	destinPnt[0] = 0;
	destinPnt[1] = 0;
	destinPnt[2] = 0;
	sourcePnt[0] = spherecenter[0];
	sourcePnt[1] = spherecenter[1];
	sourcePnt[2] = spherecenter[2];
int index3 = 0;
	 double spacing[3], origin[3];
	 int extent[6];
  imageData->GetSpacing(spacing);
  imageData->GetOrigin(origin);
  imageData->GetExtent(extent);
  for(double angle = angle_start;angle <= angle_end/*1*vnl_math::pi*/; angle += vnl_math::pi/180.0 )
      {

      destinPnt[0] = (radius *vcl_cos(angle));
      destinPnt[1] = (radius *vcl_sin(angle));
	  destinPnt[2] = 0;
	  if((lastdistinPnt[0]!=vtkMath::Round(destinPnt[0]))||(lastdistinPnt[1]!=vtkMath::Round( destinPnt[1]))||(lastdistinPnt[2]!=vtkMath::Round( destinPnt[2])))
	{
	lastdistinPnt[0] = vtkMath::Round(destinPnt[0]);
	lastdistinPnt[1] = vtkMath::Round(destinPnt[1]);
	lastdistinPnt[2] = vtkMath::Round(destinPnt[2]);
	
 
    // cast a ray in the negative direction toward sphere1
   /* sourcePnt[0] = destinPnt[0] - rayLen * normalVec[0];
    sourcePnt[1] = destinPnt[1] - rayLen * normalVec[1];
    sourcePnt[2] = destinPnt[2] - rayLen * normalVec[2];*/
 
	

	destinPnt[0] = destinPnt[0] +spherecenter[0];
	destinPnt[1] = destinPnt[1] +spherecenter[1];
	destinPnt[2] = destinPnt[2] +spherecenter[2];
 
	rayOpacity[index3]=0;
	



	// imageData = vtkImageData::New();
	
	
//  
// /*vtkMetaImageReader *reader = vtkMetaImageReader::New();
// reader->SetFileName("g:\\CTImage\\downsampledimage1.mhd");
//    reader->Update();*/
//   // imageData=reader->GetOutput();
//	imageData = vtkImporter->GetOutput();
//	//imageDataCopy = reader->GetOutputPort();
//	imageDataCopy = vtkImporter->GetOutputPort();
//   

//  int boxExtent[6];
//  boxExtent[0] = sourcePnt[0];
//  boxExtent[1] = sourcePnt[1];
//  boxExtent[2] = sourcePnt[2];
//  boxExtent[3] = destinPnt[0];
//  boxExtent[4] = destinPnt[1];
//  boxExtent[5] = destinPnt[2];

//
//  
//    
//
 double t1 = 0;
	double t2 = 1;


	/*m.m_Index[0] =(int) destinPnt[0];
	m.m_Index[1] =(int) destinPnt[1];
	m.m_Index[2] = (int)destinPnt[2];
	image->SetPixel(m,255);*/

  double x1[3], x2[3];
  for (int i = 0; i < 3; i++)
    {
		x1[i] = (sourcePnt[i]); /*- origin[i])/spacing[i];*/
		x2[i] = (destinPnt[i]);/*- origin[i])/spacing[i];*/
    }
// /* if (!this->ClipLineWithPlanes(planes, sourcePnt, destinPnt, t1, t2, clippingPlaneId))
//      {
//      return VTK_DOUBLE_MAX;
//      }*/
//   
//
//
//  volumeMapper->GetClippingPlanes();
// 
//double s1, s2;
//
////boxExtent[0] = x1[0];
////boxExtent[1] = x1[1];
////boxExtent[2] = x1[2];
////boxExtent[3] = x2[0];
////boxExtent[4] = x2[1];
////boxExtent[5] = x2[2];
////if (!this->ClipLineWithExtent(extent, x1, x2, s1, s2, planeId))
////    {
////    return VTK_DOUBLE_MAX;
////    }
////  if (s1 >= t1) { t1 = s1; }
////  if (s2 <= t2) { t2 = s2; }
////
////  // Sanity check
////  if (t2 < t1)
////    {
////    return VTK_DOUBLE_MAX;
////    }
//      
  // Compute the length of the line intersecting the volume
  double rayLength = sqrt(vtkMath::Distance2BetweenPoints(x1, x2))*(t2 - t1); 

//  // This is the minimum increment that will be allowed
  double tTol = VTKCELLPICKER_VOXEL_TOL/rayLength*(t2 - t1);

  // Find out whether there are multiple components in the volume
  int numComponents = imageData->GetNumberOfScalarComponents();
  int independentComponents = volumeProperty->GetIndependentComponents();
  int numIndependentComponents = 1;
  if (independentComponents)
    {
    numIndependentComponents = numComponents;
    }

  // Create a scalar array, it will be needed later
  vtkDataArray *scalars = vtkDataArray::CreateDataArray(imageData->GetScalarType());
  scalars->SetNumberOfComponents(numComponents);
  vtkIdType scalarArraySize = numComponents*imageData->GetNumberOfPoints();
  int scalarSize = imageData->GetScalarSize();
  void *scalarPtr = imageData->GetScalarPointer();

  // Go through each volume component separately
  double tMin = VTK_DOUBLE_MAX;
  for (int component = 0; component < numIndependentComponents; component++)
    { 
    scalarOpacity = volumeProperty->GetScalarOpacity(component);
    int disableGradientOpacity = volumeProperty->GetDisableGradientOpacity(component);
    gradientOpacity = 0;
   // if (!disableGradientOpacity && this->UseVolumeGradientOpacity)
      {
		  gradientOpacity = volumeProperty->GetGradientOpacity(component);
      }
//
    // This is the component used to compute the opacity
    int oComponent = component;
    if (!independentComponents)
      {
      oComponent = numComponents - 1;
      }

    // Make a new array, shifted to the desired component
    scalars->SetVoidArray(static_cast<void *>(static_cast<float *>(scalarPtr) + scalarSize*oComponent), scalarArraySize, 1);

     
  
  //else if ( (volumeMapper = vtkAbstractVolumeMapper::SafeDownCast(m)) )
  //  {
  //  tMin = this->IntersectVolumeWithLine(p1, p2, t1, t2, prop, volumeMapper);
  //  }

//
	 // Do a ray cast with linear interpolation.
    double opacity = 0.0;
    double lastOpacity = 0.0;
    double lastT = t1;
    double x[3];
    double pcoords[3];
    int xi[3];

    // Ray cast loop
    double t = t1;
	//rayOpacity[index] = 0;
    while (t < t2)
      {
      for (int j = 0; j < 3; j++)
        {
        // "t" is the fractional distance between endpoints x1 and x2
        x[j] = x1[j]*(1.0 - t) + x2[j]*t;

        // Paranoia bounds check
        if (x[j] < extent[2*j]) { x[j] = extent[2*j]; }
        else if (x[j] > extent[2*j + 1]) { x[j] = extent[2*j+1]; 
		}
	/*	 if (x[j] > x2[j])
		 {
			{ x[j] = x2[j]; }
			t=2;
		}*/
//       
        xi[j] = int(floor(x[j]));
        pcoords[j] = x[j] - xi[j];
        }

	  opacity = this->ComputeVolumeOpacity(xi, pcoords, imageData, scalars,
                                           scalarOpacity, gradientOpacity);
//	  std::cout<<"opacity \n"<<opacity;
	  rayOpacity[index3]+=opacity;
	  angle_max[index3] = angle;

	  
//	  
//	 
//      // If the ray has crossed the isosurface, then terminate the loop
//     /* if (opacity > opacityThreshold)
//        {
//        break;
//        }*/
//
      lastT = t;
//      lastOpacity = opacity;
//
//      // Compute the next "t" value that crosses a voxel boundary
      t = 1.0;
      for (int k = 0; k < 3; k++)
        {
        // Skip dimension "k" if it is perpendicular to ray
        if (fabs(x2[k] - x1[k]) > VTKCELLPICKER_VOXEL_TOL*rayLength)
          {
          // Compute the previous coord along dimension "k"
          double lastX = x1[k]*(1.0 - lastT) + x2[k]*lastT;

          // Increment to next slice boundary along dimension "k",
          // including a tolerance value for stabilityin cases
          // where lastX is just less than an integer value.
          double nextX = 0;
          if (x2[k] > x1[k])
            {
            nextX = floor(lastX + VTKCELLPICKER_VOXEL_TOL) + 1;
            }
          else
            {
            nextX = ceil(lastX - VTKCELLPICKER_VOXEL_TOL) - 1;
            }

          // Compute the "t" value for this slice boundary
          double ttry = lastT + (nextX - lastX)/(x2[k] - x1[k]);
          if (ttry > lastT + tTol && ttry < t)
            {
            t = ttry;
            }
          }
        } 
      } // End of "while (t <= t2)"
	
	}
	index3++;
	}
	}
	InputImageType::IndexType s;
	   s.m_Index[0] = sourcePnt[0];
	  s.m_Index[1] = sourcePnt[1];
	  s.m_Index[2] = sourcePnt[2];

	for(int i=0;i<index3;i++)
	{
		destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[i]));
		destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[i]));
		destinPnt[2] = spherecenter[2];
		m.m_Index[0] = destinPnt[0];
		m.m_Index[1] = destinPnt[1];
		m.m_Index[2] = destinPnt[2];
		if(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)==0&&imageCopy->GetPixel(m)>0)
		{
			Relation[i]=0;
		//	return 0;
			if(!Find(m))
			{
				std::cout<<imageCopy->GetPixel(m)<<"\n";
				Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
				Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
				Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
				image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				writingFile();
			
				return Centerline[numOfCenterlinePoint-1];
			}
		}
		else if(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)==0&&imageCopy->GetPixel(m)==0)
			Relation[i]=0;
		else
		{
		std::cout<</*"angle = "<<(angle_max[i]/vnl_math::pi)*180<<*/"ray Opacity = "<<rayOpacity[i]<<" m = "<<imageCopy->GetPixel(m)<<" s = "<<imageCopy->GetPixel(s)<<"difference "<<imageCopy->GetPixel(m)-imageCopy->GetPixel(s)<<" "<<rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s))<<"\n";
		Relation[i] =abs( rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)));
		}
		/*std::cout<<"x = "<<m.m_Index[0]<<" y =  "<<m.m_Index[1]<<"\n";
		std::cout<<"angle = "<<angle_max[i]<<" "<<rayOpacity[i]<<" "<<imageCopy->GetPixel(m)-imageCopy->GetPixel(s)<<" "<<rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s))<<"\n";
		Relation[i] =abs( rayOpacity[i]/(imageCopy->GetPixel(m)-imageCopy->GetPixel(s)));*/
	}
	

	//  image->SetPixel(m,255);
	
	//int max1 = maximumValue(rayOpacity);
	int max1 = maximumValue_2(Relation,index3);
//	if(rayOpacity[max1]>lastMaxOpacity)
	//{
		
	int max2;
	int max3;
	int max4;

	/*rayOpacity[max3] = 0;
	max4 = maximumValue(rayOpacity);*/
	double min = 0;
	 destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max1]));
     destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max1]));
	 destinPnt[2] = spherecenter[2];
	 
	
	
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	  image->SetPixel(m,255);
	  BOOL flag1=false;
	  BOOL flag2 = false;
	  flag1 = Find(m);
	if(flag1)
	 {
		 rayOpacity[max1] = 0;
		 Relation[max1] = 0;
		//max2 = maximumValue(rayOpacity);
		 max2 = maximumValue(Relation,4,queue,index3);
		 max1=max2;

	
	 }
	else{ 
		// if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
	 {
		 		std::cout<<angle_max[max1]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max1]<<"\n";
		 	 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			   std::cout<<"I(p)-I(x)0"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			  image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			   vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  writingFile();
			 
			  return Centerline[numOfCenterlinePoint-1];
			 
	 }
// else
	 {
		
		// Centerline[numOfCenterlinePoint].m_Index[0]++;
		//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
		// return 0;
	  }
	}

	rayOpacity[max1] = 0;
	Relation[max1] = 0;
	//max2 = maximumValue(rayOpacity);
	max2 = maximumValue(Relation,4,queue,index3);
	max1=max2;
	 destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max2]));
     destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max2]));
	 destinPnt[2] = spherecenter[2];
	 
	
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	  //image->SetPixel(m,255);

	
	
		
	if(flag1 &&Find(m))
	{
		std::cout<<angle_max[max2]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max2]<<"\n";
		rayOpacity[max2] = 0;
		Relation[max2] = 0;
		 max3 = maximumValue(rayOpacity,4,queue,index3);
		 max1=max3;
	
	}
	else if(flag1&&!Find(m))
	{
		//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
		{
			std::cout<<angle_max[max2]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max2]<<"\n";
			 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			 std::cout<<"I(p)-I(x)1"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 writingFile();
			
			return Centerline[numOfCenterlinePoint-1];
		}
//	else
		{
			//Centerline[numOfCenterlinePoint].m_Index[0]++;
			//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
			//return 0;
		}
			 //return angle_max[max1];
	}
	rayOpacity[max2] = 0;
	Relation[max2] = 0;
	//max3 = maximumValue(rayOpacity);
	max3 = maximumValue(Relation,4,queue,index3);
	max1=max3;
	 destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max3]));
     destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max3]));
	 destinPnt[2] = spherecenter[2];
	 
	  m.m_Index[0] = destinPnt[0];
	  m.m_Index[1] = destinPnt[1];
	  m.m_Index[2] = destinPnt[2];
	 
	  flag2 = Find(m);
	  if(flag1 &&flag2)
	{
		rayOpacity[max3] = 0;
		Relation[max3] = 0;
	  // max4 = maximumValue(rayOpacity);
		max4 = maximumValue(Relation,4,queue,index3);
		 max1=max4;
		 for(int i=0;i<index3;i++)
		 {
			destinPnt[0] = spherecenter[0] + (radius *vcl_cos(angle_max[max1]));
			destinPnt[1] = spherecenter[1] + (radius *vcl_sin(angle_max[max1]));
			destinPnt[2] = spherecenter[2];

			m.m_Index[0] = destinPnt[0];
			m.m_Index[1] = destinPnt[1];
			m.m_Index[2] = destinPnt[2];
			if(!Find(m))
			{
			//	if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
				{
					std::cout<<angle_max[max4]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max4]<<"\n";
				Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
				Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
				Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
				 std::cout<<"I(p)-I(x) 2"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
				 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
				 writingFile();
			
				return Centerline[numOfCenterlinePoint-1];
				 
				// break;
				}
//			else
				{
				//Centerline[numOfCenterlinePoint].m_Index[0]++;
				//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
				//return 0;
				min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				//int ind = 0;
				//for(int j=0;j<index;j++)
				//{
				//std::cout<<angle_max[max4]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max4]<<"\n";
				//m.m_Index[0] = spherecenter[0] + (radius *vcl_cos(angle_max[j]));
				//m.m_Index[1] = spherecenter[1] +(radius *vcl_sin(angle_max[j]));
				//m.m_Index[2] = spherecenter[2];
				//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<min)
				//{
				//min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				//ind = j;
				////break;
				//}


				//}
				//return angle_max[ind];
				}
				 //return angle_max[max1];
			}
			else
			{
				rayOpacity[max1] = 0;
				//max1 = maximumValue(rayOpacity);
				Relation[max1] = 0;
				max1 = maximumValue(Relation,4,queue,index3);
			}
	  }

	}
	else if(flag1&&!flag2)
	{
		//if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<0.3&&abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))>0)
		{
			std::cout<<angle_max[max3]<<" "<<destinPnt[0]<<" "<<destinPnt[1]<<" "<<rayOpacity[max3]<<"\n";
			 Centerline[numOfCenterlinePoint].m_Index[0] = destinPnt[0];
			 Centerline[numOfCenterlinePoint].m_Index[1] = destinPnt[1];
			 Centerline[numOfCenterlinePoint++].m_Index[2] = destinPnt[2];
			 std::cout<<"I(p)-I(x) 3"<<imageCopy->GetPixel(s)<<" "<<imageCopy->GetPixel(m)<<" "<<imageCopy->GetPixel(s)-imageCopy->GetPixel(m)<<"\n";
			 writingFile();
			 image->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			 vesselImage->SetPixel(Centerline[numOfCenterlinePoint-1],255);
			  writingFile();
		
			 return Centerline[numOfCenterlinePoint-1];
		}
	//else
		{
				for(int i=0;i<index3;i++)
				{
				m.m_Index[0] = spherecenter[0] + (radius *vcl_cos(angle_max[i]));
				m.m_Index[1] = spherecenter[1] +(radius *vcl_sin(angle_max[i]));
				m.m_Index[2] = spherecenter[2];
				if(abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m))<min)
				{
				min = abs(imageCopy->GetPixel(s)-imageCopy->GetPixel(m));
				break;
				}


				}
			//Centerline[numOfCenterlinePoint].m_Index[0]++;
			//Circle_RayCasting(Centerline[numOfCenterlinePoint],-1*vnl_math::pi/6,vnl_math::pi,2,imageCopy);
			//return 0;
		//}
			  
	}

	  //lastMaxOpacity = rayOpacity[max1];

	}

	typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	 WriterType::Pointer writer = WriterType::New();
	 writer->SetFileName( "g:\\CTImage\\Sphere.mhd" );
	 writer->SetInput(image);
	 
  try
    {
   writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  return Centerline[numOfCenterlinePoint-2];
}

BOOL VTKRayCaster::Find(InputImageType::IndexType index)
{
	for(int i=0;i<numOfCenterlinePoint;i++)
	{
		if(Centerline[i].m_Index[0]==index.m_Index[0]&&Centerline[i].m_Index[1]==index.m_Index[1]&&Centerline[i].m_Index[2]==index.m_Index[2])
		{
			return true;
		}
		else if(vcl_sqrt(pow((float)( Centerline[i].m_Index[0] - index[0]),2)+pow((float)( Centerline[i].m_Index[1] - index[1]),2)+pow((float)( Centerline[i].m_Index[2] - index[2]),2))<1)
			return true;
	}
	return false;
}

BOOL VTKRayCaster::FindSpherePoint(InputImageType::IndexType Index,int index)
{
	for(int i=0;i<index;i++)
	{
		if(spherePoint[i].m_Index[0]==Index.m_Index[0]&&spherePoint[i].m_Index[1]==Index.m_Index[1]&&spherePoint[i].m_Index[2]==Index.m_Index[2])
		{
			return true;
		}
	}
	return false;
}



void VTKRayCaster::SetTwoDimImage(OutputImageType2D::Pointer image2D)
{
	this->image2D = image2D;
}

int VTKRayCaster::ConnectedComponentLabeling(InputImageType::IndexType m_Index)
{
	int NumOfConnectedComponent = 0;
	InputImageType::IndexType temp;
	double centerIntensity = imageCopy->GetPixel(m_Index);
	temp = m_Index;
	for(int i=m_Index[2]-1;i<m_Index[2]+1;i++)
	{
		temp[2] = i;
		double sum = imageCopy->GetPixel(temp);
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp.m_Index[0] = m_Index.m_Index[0]+1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp.m_Index[1] = m_Index.m_Index[1]-1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp.m_Index[1] = m_Index.m_Index[1]+1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp = m_Index;
		temp[2] = i;
		temp.m_Index[0] = m_Index.m_Index[0]-1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp.m_Index[1] = m_Index.m_Index[1]-1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp.m_Index[1] = m_Index.m_Index[1]+1;
		
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		
		temp.m_Index[1] = m_Index.m_Index[1]-1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
		temp.m_Index[1] = m_Index.m_Index[1]+1;
		sum =abs( imageCopy->GetPixel(temp));
		if(abs(sum-centerIntensity)<=1)
			NumOfConnectedComponent++;
	}
	return NumOfConnectedComponent;
}
void VTKRayCaster::ComputeCorrelation()
{
	 

	


	typedef itk::RegionOfInterestImageFilter< OutputImageType2D,OutputImageType2D > ExtractFilterType;
 
  ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();
 
  OutputImageType2D::IndexType start;
  start[0] = 105;
  start[1] = 109;
  //start[2] = 109;
 // start.Fill(50);
 
  OutputImageType2D::SizeType patchSize;
  patchSize[0] = 3;
  patchSize[1] = 3;
  //patchSize[2] = 0;
 // patchSize.Fill(51);
 
  OutputImageType2D::RegionType desiredRegion(start,patchSize);
 
  extractFilter->SetRegionOfInterest(desiredRegion);
  extractFilter->SetInput(image2D);
  extractFilter->Update();

  typedef itk::NormalizedCorrelationImageFilter<OutputImageType2D, OutputImageType2D, OutputImageType2D> CorrelationFilterType;
 
  itk::ImageKernelOperator<float> kernelOperator;
  kernelOperator.SetImageKernel(extractFilter->GetOutput());
 
  // The radius of the kernel must be the radius of the patch, NOT the size of the patch
  itk::Size<2> radius = extractFilter->GetOutput()->GetLargestPossibleRegion().GetSize();
  radius[0] = (radius[0]-1) / 2;
  radius[1] = (radius[1]-1) / 2;
  //radius[2] = (radius[2]-1) / 2;
 
  kernelOperator.CreateToRadius(radius);
 
  CorrelationFilterType::Pointer correlationFilter = CorrelationFilterType::New();
  correlationFilter->SetInput(image2D);
  correlationFilter->SetTemplate(kernelOperator);
  correlationFilter->Update();
 
  typedef itk::MinimumMaximumImageCalculator <OutputImageType2D>    MinimumMaximumImageCalculatorType;

  MinimumMaximumImageCalculatorType::Pointer minimumMaximumImageCalculatorFilter
          = MinimumMaximumImageCalculatorType::New ();
  minimumMaximumImageCalculatorFilter->SetImage(correlationFilter->GetOutput());
  minimumMaximumImageCalculatorFilter->Compute();
 
  itk::Index<2> maximumCorrelationPatchCenter = minimumMaximumImageCalculatorFilter->GetIndexOfMaximum();
  std::cout << "Maximum: " << maximumCorrelationPatchCenter << std::endl;
 
  // Note that the best correlation is at the center of the patch we extracted (ie. (75, 75) rather than the corner (50,50)
 
  typedef itk::RescaleIntensityImageFilter< OutputImageType2D, OutputImageType2D > RescaleFilterType;
  typedef itk::ImageFileWriter<OutputImageType2D> WriterType;
  
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(correlationFilter->GetOutput());
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

 


   QuickView viewer;
 // viewer.AddImage(reader->GetOutput());
  //viewer.AddImage(extractFilter->GetOutput());
 // viewer.AddImage(correlationFilter->GetOutput());
 // viewer.AddImage(bestPatchExtractFilter->GetOutput());
  //viewer.Visualize();


}