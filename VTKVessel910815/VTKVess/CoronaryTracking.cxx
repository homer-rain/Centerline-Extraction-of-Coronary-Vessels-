#include "CoronaryTracking.h"

CoronaryTracking::CoronaryTracking(void)
{
}

CoronaryTracking::~CoronaryTracking(void)
{
}
void CoronaryTracking::track(InputImageType::IndexType *m_CircleIndex,InputImageType::Pointer image,int indexCounter)
{
	  const unsigned int Dimension = 3;
	typedef itk::ImageRegionConstIterator< InputImageType > ConstIteratorType;
	typedef itk::ImageRegionIterator< InputImageType>       IteratorType;
	typedef itk::ImageFileWriter< InputImageType > WriterType;
	const InputImageType::SpacingType& spacing = image->GetSpacing();
	const InputImageType::PointType& inputOrigin = image->GetOrigin();
	InputImageType::RegionType::IndexType outputStart;
	InputImageType::RegionType outputRegion;

  outputStart[0] = 0;
  outputStart[1] = 0;
  InputImageType::RegionType::SizeType  size;
  size = image->GetBufferedRegion().GetSize();
  outputRegion.SetSize( size );
  outputRegion.SetIndex( outputStart );


  InputImageType::Pointer outputImage = InputImageType::New();
    outputImage->SetRegions( outputRegion );

	double   outputOrigin[ Dimension ];

	for(unsigned int i=0; i< Dimension; i++)
	{
	outputOrigin[i] = inputOrigin[i] + spacing[i] /** inputStart[i]*/;
	}
	outputImage->SetSpacing( spacing );
	outputImage->SetOrigin(  outputOrigin );
	outputImage->Allocate();


	double maxValue1 = image->GetPixel(m_CircleIndex[0]);
	 double maxValue2 = image->GetPixel(m_CircleIndex[0]);
	 InputImageType::IndexType indexOfMaxValue1;
	 InputImageType::IndexType indexOfMaxValue2;
	 int index = 0;
	 for(int i=1;i<indexCounter;i++)
	  {
		  if(maxValue1<image->GetPixel(m_CircleIndex[i]))
		  {
			  maxValue1 = image->GetPixel(m_CircleIndex[i]);
			  indexOfMaxValue1 = m_CircleIndex[i];
			  index =i;

		  }
		  image->SetPixel(m_CircleIndex[i],(0,255,0));
		//  outputImage->SetPixel(m_CircleIndex[i],(0,255,0));

		 
	  }
 //std::cerr <<image->GetPixel(indexOfMaxValue1) << std::endl;;

 for(int i=1;i<indexCounter;i++)
  {
	  if(maxValue2<image->GetPixel(m_CircleIndex[i])&&i!=index)
		 // if((m_CircleIndex[i].m_Index[0]!=indexOfMaxValue1.m_Index[0])&&(m_CircleIndex[i].m_Index[1]!=indexOfMaxValue1.m_Index[1])&&(m_CircleIndex[i].m_Index[2]!=indexOfMaxValue1.m_Index[2]))
	  {
		  maxValue2 = image->GetPixel(m_CircleIndex[i]);
		  indexOfMaxValue2 = m_CircleIndex[i];

	  }

//	  std::cout <<image->GetPixel(m_CircleIndex[i]) << std::endl;;
  }
 // std::cerr <<image->GetPixel(indexOfMaxValue2) << std::endl;;
WriterType::Pointer writer = WriterType::New();
  writer->SetFileName( "aorta.mhd" );
  writer->SetInput( outputImage );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
//    return -1;   
    }

	 ConstIteratorType inputIt(   image, image->GetRequestedRegion()  );
    IteratorType      outputIt(  outputImage,         outputRegion );

  inputIt.GoToBegin();
 // outputIt.GoToBegin();

  while( !inputIt.IsAtEnd() )
    {
		if(inputIt.Get()>80)
			outputIt.Set(  inputIt.Get()  );
	//	std::cout<<inputIt.GetIndex()<<" "<<	inputIt.Get()<<"\n";

		++inputIt;
    ++outputIt;
    }




  writer->SetFileName( "vessel80.mhd" );
  writer->SetInput( outputImage );

  try
    {
    writer->Update();
    }
  catch ( itk::ExceptionObject &err)
    {
    std::cerr << "ExceptionObject caught !" << std::endl; 
    std::cerr << err << std::endl; 
//    return -1;   
    }
}
