#include "RegionGrowing.h"

const     unsigned int    Dimension = 2;

// Software Guide : EndCodeSnippet





RegionGrowing::RegionGrowing(ImageType::Pointer image)
{

	this->image = image;
	//image = ImageType::New();
	//CreateImage();

}

RegionGrowing::~RegionGrowing(void)
{
}

void RegionGrowing:: CreateImage ()
{
  // Create an image with 2 connected components
	ImageType::RegionType region ;
  ImageType::IndexType start;
  start.m_Index[0]=0;
  start.m_Index[1]=0;
  start.m_Index[2]=0;
  
 
  ImageType::SizeType size;
  unsigned int NumRows = 256;
  unsigned int NumCols = 256;
  size[0] = NumRows;
  size[1] = NumCols;
  size[2] = 136;
 
  region.SetSize(size);
  region.SetIndex(start);
 
  image->SetRegions(region);
  image->Allocate();
}


ImageType::Pointer RegionGrowing::Opening()
{
	  ErodeFilterType::Pointer  grayscaleErode  = ErodeFilterType::New();
	 StructuringElementType  structuringElement;
	
	  structuringElement.SetRadius( 6 );  // 3x3 structuring element
	
	  structuringElement.CreateStructuringElement();
	
	 grayscaleErode->SetKernel(  structuringElement );
	
	 grayscaleErode->SetInput(image);
	
	// grayscaleErode->Update();
	 //return grayscaleErode->GetOutput();
	 return image;
}
OutputImageType::Pointer RegionGrowing:: AortaSegmentation(OutputImageType::Pointer image2d,OutputImageType::IndexType radiusIndex,int axialnum)

{
 typedef   float           InternalPixelType;
  const     unsigned int    Dimension = 2;
  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  // Software Guide : EndCodeSnippet


  typedef float                            OutputPixelType;
  typedef itk::Image< OutputPixelType, Dimension > OutputImageType;

  typedef itk::CastImageFilter< InternalImageType, OutputImageType >    CastingFilterType;


 
  CastingFilterType::Pointer caster = CastingFilterType::New();



  // We instantiate reader and writer types
  //
   typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	 WriterType::Pointer writer = WriterType::New();
	 writer->SetFileName( "g:\\CTImage\\SegmentedAorta.mhd" );
  
  typedef   itk::CurvatureFlowImageFilter<InternalImageType, InternalImageType>
    CurvatureFlowImageFilterType;
  
  CurvatureFlowImageFilterType::Pointer smoothing =
                         CurvatureFlowImageFilterType::New();
  
  typedef itk::NeighborhoodConnectedImageFilter<InternalImageType,  InternalImageType > ConnectedFilterType;
 
  ConnectedFilterType::Pointer neighborhoodConnected = ConnectedFilterType::New();
 
  smoothing->SetInput( image2d );
  neighborhoodConnected->SetInput( smoothing->GetOutput() );
  caster->SetInput( neighborhoodConnected->GetOutput() );
 writer->SetInput( caster->GetOutput() );
 
  smoothing->SetNumberOfIterations( 5 );
  smoothing->SetTimeStep( 0.125 );
  // Software Guide : EndCodeSnippet


  

  const InternalPixelType lowerThreshold = atof( "1150" );
  const InternalPixelType upperThreshold = atof( "1420" );

  // Software Guide : BeginCodeSnippet
  neighborhoodConnected->SetLower(  lowerThreshold  );
  neighborhoodConnected->SetUpper(  upperThreshold  );
 
  InternalImageType::SizeType   radius;

  radius[0] = 2;   // two pixels along X
  radius[1] = 2;   // two pixels along Y

  neighborhoodConnected->SetRadius( radius );
  

  InternalImageType::IndexType  index;

  index[0] = radiusIndex[0];
  index[1] = radiusIndex[1];


  // Software Guide : BeginCodeSnippet
  neighborhoodConnected->SetSeed( index );
  neighborhoodConnected->SetReplaceValue( 255 );
 /* QuickView viewer;
  viewer.AddImage<OutputImageType>(caster->GetOutput());
  viewer.Visualize();*/
  /*pasteFilter->SetSourceImage(caster->GetOutput());
  pasteFilter->SetSourceRegion(caster->GetOutput()->GetLargestPossibleRegion());*/
 ImageType::IndexType index3d;

 index3d.m_Index[0] = 0;
 index3d.m_Index[1]=0;
 index3d.m_Index[2] = axialnum;
 OutputImageType::IndexType index2d;
 //pasteFilter->SetDestinationIndex(index2);

 
  try
    {
   writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  for(int i=0;i<256;i++)
	 for(int j=0;j<256;j++)
	 {
		 index2d.m_Index[0]=j;
		 index2d.m_Index[1] = i;
		 index3d.m_Index[0]=j;
		 index3d.m_Index[1] = i;
		 index3d.m_Index[2] = axialnum;
		 if(caster->GetOutput()->GetPixel(index2d)==255)
		 image->SetPixel(index3d,caster->GetOutput()->GetPixel(index2d));
	 }
  return caster->GetOutput();


}







ImageType::Pointer RegionGrowing:: VesselSegmentation(ImageType::Pointer image3d,ImageType::IndexType radiusIndex,int axialnum)

{
 typedef   float           InternalPixelType;
  const     unsigned int    Dimension = 3;
  typedef itk::Image< InternalPixelType, 3 >  InternalImageType;
  // Software Guide : EndCodeSnippet


  typedef float                            OutputPixelType;
  typedef itk::Image< OutputPixelType, 3 > OutputImageType;

  typedef itk::CastImageFilter< ImageType, ImageType >    CastingFilterType;

 typedef itk::ImageFileReader<ImageType> FileReaderType;
		  FileReaderType::Pointer imageReader = FileReaderType::New();
		  //imageReader->SetFileName("g:\\RawReader\\output_en2.mhd");
		  imageReader->SetFileName("g:\\Rawreader\\output_en2.mhd");
  try
  { 
    imageReader->Update();
  }
  catch (itk::ExceptionObject &ex)
  { 
    std::cout << ex << std::endl;
    //return EXIT_FAILURE;
  }
  image3d = imageReader->GetOutput();


   typedef itk::RescaleIntensityImageFilter<
               ImageType, ImageType >  RescaleFilterType;
     RescaleFilterType::Pointer    rescaleFilter    = RescaleFilterType::New();
  CastingFilterType::Pointer caster = CastingFilterType::New();
  rescaleFilter->SetInput(    image3d );
  rescaleFilter->SetOutputMinimum(  0 );
  rescaleFilter->SetOutputMaximum( 255 );
    rescaleFilter->Update();
  // We instantiate reader and writer types
  //
   typedef  itk::ImageFileWriter<  OutputImageType  > WriterType;

	 WriterType::Pointer writer = WriterType::New();
	 writer->SetFileName( "g:\\CTImage\\SegmentedAorta3d.mhd" );
  
  typedef   itk::CurvatureFlowImageFilter<ImageType, ImageType>    CurvatureFlowImageFilterType;
  
  CurvatureFlowImageFilterType::Pointer smoothing =  CurvatureFlowImageFilterType::New();
  
  typedef itk::NeighborhoodConnectedImageFilter<InternalImageType,  InternalImageType > ConnectedFilterType;
 
  ConnectedFilterType::Pointer neighborhoodConnected = ConnectedFilterType::New();
 
  smoothing->SetInput( rescaleFilter->GetOutput() );
  neighborhoodConnected->SetInput( smoothing->GetOutput() );
  caster->SetInput( neighborhoodConnected->GetOutput() );
  writer->SetInput( caster->GetOutput() );
//writer->SetInput( rescaleFilter->GetOutput() );
  //writer->SetInput(image3d);
  smoothing->SetNumberOfIterations( 5 );
  smoothing->SetTimeStep( 0.125 );
  // Software Guide : EndCodeSnippet


  std::cout<<rescaleFilter->GetOutput()->GetPixel(radiusIndex)<<"-------------------- \n";

  const InternalPixelType lowerThreshold = atof( "150" );
  const InternalPixelType upperThreshold = atof( "240" );

  // Software Guide : BeginCodeSnippet
  neighborhoodConnected->SetLower(  lowerThreshold  );
  neighborhoodConnected->SetUpper(  upperThreshold  );
 
  InternalImageType::SizeType   radius;

  radius[0] = 2;   // two pixels along X
  radius[1] = 2;   // two pixels along Y
  radius[2] = 2;   // two pixels along Z

  neighborhoodConnected->SetRadius( radius );
  

  InternalImageType::IndexType  index;

  index[0] = radiusIndex[0];
  index[1] = radiusIndex[1];
  index[2] = radiusIndex[2];


  // Software Guide : BeginCodeSnippet
  neighborhoodConnected->SetSeed( index );
  neighborhoodConnected->SetReplaceValue( 255 );
 /* QuickView viewer;
  viewer.AddImage<OutputImageType>(caster->GetOutput());
  viewer.Visualize();*/
  /*pasteFilter->SetSourceImage(caster->GetOutput());
  pasteFilter->SetSourceRegion(caster->GetOutput()->GetLargestPossibleRegion());*/
 ImageType::IndexType index3d;

 index3d.m_Index[0] = 0;
 index3d.m_Index[1] = 0;
 index3d.m_Index[2] = axialnum;
 OutputImageType::IndexType index2d;
 //pasteFilter->SetDestinationIndex(index2);

 
  try
    {
   writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  //for(int i=0;i<256;i++)
	 //for(int j=0;j<256;j++)
	 //{
		// index2d.m_Index[0]=j;
		// index2d.m_Index[1] = i;
		// index3d.m_Index[0]=j;
		// index3d.m_Index[1] = i;
		// index3d.m_Index[2] = axialnum;
		// 
		// image->SetPixel(index3d,caster->GetOutput()->GetPixel(index2d));
	 //}
  return caster->GetOutput();


}


