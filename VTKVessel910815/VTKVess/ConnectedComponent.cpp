#include "ConnectedComponent.h"
  const     unsigned int    Dimension = 3;
ConnectedComponent::ConnectedComponent(void)
{
}

ConnectedComponent::~ConnectedComponent(void)
{
}


#include "itkScalarConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkImageFileReader.h"
#include "itkVTKImageIO.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkSimpleFilterWatcher.h"
#include "itkImageRegionIterator.h"
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkImageRegionIteratorWithIndex.h"
 
RGBImageType::Pointer ConnectedComponent::filterImage(InternalImageType1::Pointer image)
{
   const     unsigned int    Dimension = 3;
  typedef itk::Image< MaskPixelType, Dimension >  MaskImageType;
  typedef itk::Image<unsigned short,Dimension> OutputImageType;


 

  typedef itk::ImageFileReader< InternalImageType1 > ReaderType;
  typedef itk::ImageFileWriter<  RGBImageType  > WriterType;

  
  typedef itk::ScalarConnectedComponentImageFilter< InternalImageType1, OutputImageType, MaskImageType > FilterType;
  typedef itk::RelabelComponentImageFilter< OutputImageType, OutputImageType > RelabelType;
  
 
  WriterType::Pointer writer = WriterType::New();

  FilterType::Pointer filter = FilterType::New();
  RelabelType::Pointer relabel = RelabelType::New();
  
  itk::SimpleFilterWatcher watcher(filter);
  watcher.QuietOn();


  
  // create a mask containing the upper left hand corner and
  // a chunk out of the middle
  MaskImageType::Pointer mask = MaskImageType::New();
  mask->SetRegions(image->GetLargestPossibleRegion());
  mask->CopyInformation(image);
  mask->Allocate();
  mask->FillBuffer(itk::NumericTraits<MaskPixelType>::Zero);

  MaskImageType::RegionType maskRegion = mask->GetLargestPossibleRegion();
  MaskImageType::SizeType maskSize = maskRegion.GetSize();

  MaskImageType::RegionType region;
  MaskImageType::SizeType size;
  MaskImageType::IndexType index;
  
  // use upper left corner
  index.Fill(0);
  for (unsigned int i=0; i<MaskImageType::ImageDimension; i++)
    {
    size[i] = static_cast<unsigned long> (0.5 * maskSize[i]);
    }
  region.SetIndex (index);
  region.SetSize (size);

  itk::ImageRegionIterator<MaskImageType> mit(mask,region);
  while (!mit.IsAtEnd())
    {
    mit.Set(itk::NumericTraits<MaskPixelType>::max());
    ++mit;
    }

  // use middle section
  for (unsigned int i=0; i<MaskImageType::ImageDimension; i++)
    {
    index[i] = static_cast<long> (0.375 * maskSize[i]);
    size[i] = static_cast<unsigned long> (0.25 * maskSize[i]);
    }
  region.SetIndex (index);
  region.SetSize (size);

  itk::ImageRegionIterator<MaskImageType> mit2(mask,region);
  while (!mit2.IsAtEnd())
    {
    mit2.Set(itk::NumericTraits<MaskPixelType>::max());
    ++mit2;
    }

  filter->SetInput (image);
  filter->SetMaskImage (mask);
  filter->SetDistanceThreshold( 5 );
  filter->SetFunctor(filter->GetFunctor());

  //if (argc > 4)
  //  {
  //  int fullyConnected = atoi( argv[4] );
  //  filter->SetFullyConnected( fullyConnected );
  //  }
  relabel->SetInput( filter->GetOutput() );
  //if (argc > 5)
  //  {
  //  int minSize = atoi( argv[5] );
  //  relabel->SetMinimumObjectSize( minSize );
  //  std::cerr << "minSize: " << minSize << std::endl;
  //  }
  
  try
    {
    relabel->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Relabel: exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

  // Remap the labels to viewable colors
  RGBImageType::Pointer colored = RGBImageType::New();
  colored->SetRegions( filter->GetOutput()->GetBufferedRegion() );
  colored->Allocate();

  unsigned short numObjects = relabel->GetNumberOfObjects();
  
  std::vector<RGBPixelType> colormap;
  RGBPixelType px;
  colormap.resize( numObjects+1 );
  itk::Statistics::MersenneTwisterRandomVariateGenerator::GetInstance()->SetSeed(1031571);
  itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer rvgen = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
  for (unsigned short i=0; i < colormap.size(); ++i)
    {
    px.SetRed(
      static_cast<unsigned char>(255*rvgen->GetUniformVariate( 0.3333, 1.0 ) ));
    px.SetGreen(
      static_cast<unsigned char>(255*rvgen->GetUniformVariate( 0.3333, 1.0 ) ));
    px.SetBlue(
      static_cast<unsigned char>(255*rvgen->GetUniformVariate( 0.3333, 1.0 ) ));
    colormap[i] = px;
    }
  
  itk::ImageRegionIteratorWithIndex<OutputImageType>
    it(relabel->GetOutput(), relabel->GetOutput()->GetBufferedRegion());
  itk::ImageRegionIteratorWithIndex<RGBImageType> cit(colored,
                                             colored->GetBufferedRegion());
 int* arrayOut =new int[  colormap.size()];
 for(int i=0;i< colormap.size();i++)
	 arrayOut[i] = 0;
  while( !it.IsAtEnd() )
    {
    if (it.Get() == 0)
      {
      cit.Set(RGBPixelType(static_cast<unsigned char>(0)));
      }
    else
      {
      cit.Set( colormap[it.Get()] );
	  arrayOut[it.Get()]++;
	  std::cout<<colormap[it.Get()]<<" "<<it.GetIndex()<<std::endl;
      }
    ++it;
    ++cit;
    }
  
  try
    {
    writer->SetInput (colored);
    writer->SetFileName( "output.mhd" );
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  return colored;
  }

  RGBImageType::Pointer ConnectedComponent:: itkRelabelComponentImageFilterTest(InternalImageType1::Pointer image )
{
  /*if( argc < 5 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0];
    std::cerr << " inputImage  outputImage threshold_low threshold_hi" << std::endl;
    return EXIT_FAILURE;
    }*/

  typedef   unsigned short  InternalPixelType;
  typedef   unsigned long   LabelPixelType;
  typedef   unsigned char   WritePixelType;
  const     unsigned int    Dimension = 3;

  typedef itk::Image< InternalPixelType, Dimension >  InternalImageType;
  typedef itk::Image< LabelPixelType, Dimension>   LabelImageType;
  typedef itk::Image<WritePixelType, Dimension> WriteImageType;

  typedef itk::ImageFileReader< InternalImageType > ReaderType;
  typedef itk::ImageFileWriter<  WriteImageType  > WriterType;


  typedef itk::ChangeInformationImageFilter<InternalImageType> ChangeFilterType;
  typedef itk::BinaryThresholdImageFilter< InternalImageType, InternalImageType > ThresholdFilterType;
  typedef itk::ConnectedComponentImageFilter< InternalImageType, LabelImageType > ConnectedComponentType;
  typedef itk::RelabelComponentImageFilter< LabelImageType, LabelImageType > RelabelComponentType;
  typedef itk::BinaryThresholdImageFilter<LabelImageType, WriteImageType> FinalThresholdFilterType;
  typedef itk::LabelStatisticsImageFilter< InternalImageType, LabelImageType> StatisticsFilterType;

  typedef itk::NumericTraits<InternalPixelType>::RealType RealType;

  typedef itk::Statistics::Histogram<RealType> HistogramType;

  typedef HistogramType::IndexType HIndexType;
  int NumBins = 13;
  RealType LowerBound = 51.0;
  RealType UpperBound = 252.0;

  ReaderType::Pointer reader = ReaderType::New();
  WriterType::Pointer writer = WriterType::New();
  ChangeFilterType::Pointer change = ChangeFilterType::New();
  ThresholdFilterType::Pointer threshold = ThresholdFilterType::New();
  ConnectedComponentType::Pointer connected = ConnectedComponentType::New();
  RelabelComponentType::Pointer relabel = RelabelComponentType::New();
  FinalThresholdFilterType::Pointer finalThreshold = FinalThresholdFilterType::New();
  StatisticsFilterType::Pointer statistics = StatisticsFilterType::New();

  FilterWatcher watcher(relabel);
  FilterWatcher statswatcher(statistics);

  //reader->SetFileName( argv[1] );

  // try changing the spacing on the output to test the sorting on
  // physical size
  ChangeFilterType::SpacingType changeSpacing;
  changeSpacing[0] = 1;
  changeSpacing[1] = 0.5;
  change->SetInput( image );
  change->SetOutputSpacing(changeSpacing);
  change->ChangeSpacingOn();

  // Create a binary input image to label
  InternalPixelType threshold_low, threshold_hi;
  threshold_low = atoi( "48");
  threshold_hi = atoi( "80");

  threshold->SetInput (change->GetOutput());
  threshold->SetInsideValue(itk::NumericTraits<InternalPixelType>::One);
  threshold->SetOutsideValue(itk::NumericTraits<InternalPixelType>::Zero);
  threshold->SetLowerThreshold(threshold_low);
  threshold->SetUpperThreshold(threshold_hi);
  threshold->Update();

  // Label the components in the image and relabel them so that object
  // numbers increase as the size of the objects decrease.
  connected->SetInput (threshold->GetOutput());
  relabel->SetInput( connected->GetOutput() );
  relabel->SetNumberOfObjectsToPrint( 5 );
  std::cout << "Modified time of relabel's output = " << relabel->GetOutput()->GetMTime() << std::endl;
  relabel->Update();
  std::cout << "NumberOfObjects: " << relabel->GetNumberOfObjects() << " OriginalNumberOfObjects: " <<
    relabel->GetOriginalNumberOfObjects() << " MinimumObjectSize: " << relabel->GetMinimumObjectSize() << std::endl;

  // pull out the largest object
  finalThreshold->SetInput( relabel->GetOutput() );
  finalThreshold->SetLowerThreshold( 1 ); // object #1
  finalThreshold->SetUpperThreshold( 1 ); // object #1
  finalThreshold->SetInsideValue(255);
  finalThreshold->SetOutsideValue(itk::NumericTraits<WritePixelType>::Zero);

  try
    {
    writer->SetInput (finalThreshold->GetOutput());
    writer->SetFileName( "output.mhd" );
    writer->Update();
    std::cout << "Modified time of relabel's output = " << relabel->GetOutput()->GetMTime() << std::endl;
    writer->Update();
    std::cout << "Modified time of relabel's output = " << relabel->GetOutput()->GetMTime() << std::endl;
    relabel->Modified();
    relabel->Update();
    std::cout << "Modified time of relabel's output = " << relabel->GetOutput()->GetMTime() << std::endl;
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }


  try
    {
    statistics->SetInput( change->GetOutput() );
    statistics->SetLabelInput( relabel->GetOutput() );
    statistics->SetHistogramParameters(NumBins, LowerBound, UpperBound);
    statistics->UseHistogramsOn();
    statistics->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught during statistics calculation!"
              << std::endl;
    std::cerr << excep << std::endl;
    }
  try
    {
    HistogramType::Pointer histogram;
    unsigned long printNum = statistics->GetNumberOfLabels();
    if (printNum > 10)
      {
      printNum = 10;
      }
    for (unsigned int ii=0; ii < printNum; ++ii)
      {
      std::cout << "Label " << ii << ": " << (statistics->HasLabel(ii) ? "Exists" : "Does not exist") << std::endl;
      std::cout << "\tCount = " << statistics->GetCount(ii) << std::endl;
      std::cout << "\tMinimum = " << statistics->GetMinimum(ii) << std::endl;
      std::cout << "\tMaximum = " << statistics->GetMaximum(ii) << std::endl;
      std::cout << "\tMean = " << statistics->GetMean(ii) << std::endl;
      std::cout << "\tSigma = " << statistics->GetSigma(ii) << std::endl;
      std::cout << "\tVariance = " << statistics->GetVariance(ii) << std::endl;
      std::cout << "\tSum = " << statistics->GetSum(ii) << std::endl;
      std::cout << "\tMedian = " << statistics->GetMedian(ii) << std::endl;
      std::cout << "\tRegion = " << statistics->GetRegion(ii) << std::endl;
      const StatisticsFilterType::BoundingBoxType bbox =
        statistics->GetBoundingBox(ii);

      std::cout << "\tBounding box = ";
      for ( unsigned int jj = 0; jj <  bbox.size(); jj++)
        {
        std::cout << bbox[jj] << " ";
        }
      std::cout << std::endl;
      if (statistics->HasLabel(ii))
        {
        std::cout << "\tHistogram Frequencies:" << std::endl;
        histogram = statistics->GetHistogram(ii);
        for (int jj=0;jj<=NumBins;jj++)
          {
          std::cout << histogram->GetFrequency(jj) << ", ";
          }
        std::cout <<  std::endl;
        }
      }

    printNum = 2;
    for (unsigned int ii=statistics->GetNumberOfObjects();
           ii < statistics->GetNumberOfObjects()+printNum; ++ii)
      {
      std::cout << "Label " << ii << ": " << (statistics->HasLabel(ii) ? "Exists" : "Does not exist") << std::endl;
      std::cout << "\tCount = " << statistics->GetCount(ii) << std::endl;
      std::cout << "\tMinimum = " << statistics->GetMinimum(ii) << std::endl;
      std::cout << "\tMaximum = " << statistics->GetMaximum(ii) << std::endl;
      std::cout << "\tMean = " << statistics->GetMean(ii) << std::endl;
      std::cout << "\tSigma = " << statistics->GetSigma(ii) << std::endl;
      std::cout << "\tVariance = " << statistics->GetVariance(ii) << std::endl;
      std::cout << "\tSum = " << statistics->GetSum(ii) << std::endl;
      std::cout << "\tMedian = " << statistics->GetMedian(ii) << std::endl;
      if (statistics->HasLabel(ii))
        {
        std::cout << "\tEvery tenth Histogram Frequencies:" << std::endl;
        histogram = statistics->GetHistogram(ii);
        for (int jj=0;jj<=NumBins;jj++)
          {
          std::cout << histogram->GetFrequency(jj) << ", ";
          }
        std::cout <<  std::endl;
        }
      }

    }
  catch (...)
    {
    std::cerr << "Exception caught while printing statistics" << std::endl;
    }

  return EXIT_SUCCESS;
}

//IntImageType::Pointer ConnectedComponent::filterImage(IntImageType::Pointer image)
//{
//
//
// 
//  typedef itk::ScalarConnectedComponentImageFilter<IntImageType, IntImageType> ScalarConnectedComponentImageFilterType;
//  ScalarConnectedComponentImageFilterType::Pointer scalarConnectedComponentFilter = ScalarConnectedComponentImageFilterType::New();
//  scalarConnectedComponentFilter->SetInput(image);
//  scalarConnectedComponentFilter->Update();
// 
//  typedef  itk::ImageFileWriter< IntImageType  > WriterType;
//  WriterType::Pointer writer = WriterType::New();
//  writer->SetFileName("output.mhd");
//  writer->SetInput(scalarConnectedComponentFilter->GetOutput());
//  writer->Update();
//
// 
//  //QuickView viewer;
//  //viewer.AddImage<ImageType>( image, false );  
//  //viewer.AddImage<ImageType>( rescaleFilter->GetOutput(), false );  
//  //viewer.Visualize();
//  return scalarConnectedComponentFilter->GetOutput();
//}