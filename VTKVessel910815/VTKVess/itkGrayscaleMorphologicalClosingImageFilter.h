 /*=========================================================================
 
   Program:   Insight Segmentation & Registration Toolkit
   Module:    $RCSfile: itkGrayscaleMorphologicalClosingImageFilter.h,v $
   Language:  C++
   Date:      $Date: -- :: $
   Version:   $Revision: . $
 
   Copyright (c) Insight Software Consortium. All rights reserved.
   See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.
 
      This software is distributed WITHOUT ANY WARRANTY; without even
      the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
      PURPOSE.  See the above copyright notices for more information.
 
 =========================================================================*/
 #ifndef __itkGrayscaleMorphologicalClosingImageFilter_h
 #define __itkGrayscaleMorphologicalClosingImageFilter_h
 
 // First make sure that the configuration is available.
 // This line can be removed once the optimized versions
 // gets integrated into the main directories.
 #include "itkConfigure.h"
 
 #ifdef ITK_USE_CONSOLIDATED_MORPHOLOGY
 #include "itkOptGrayscaleMorphologicalClosingImageFilter.h"
 #else
 
 
 #include "itkImageToImageFilter.h"
 
 namespace itk {
 
 template<class TInputImage, class TOutputImage, class TKernel>
 class ITK_EXPORT GrayscaleMorphologicalClosingImageFilter : 
     public ImageToImageFilter<TInputImage, TOutputImage>
 {
 public:
   typedef GrayscaleMorphologicalClosingImageFilter      Self;
   typedef ImageToImageFilter<TInputImage, TOutputImage> Superclass;
   typedef SmartPointer<Self>                            Pointer;
   typedef SmartPointer<const Self>                      ConstPointer;
 
   itkNewMacro(Self);  
 
   itkTypeMacro(GrayscaleMorphologicalClosingImageFilter, 
                ImageToImageFilter);
 
   typedef TInputImage                              InputImageType;
   typedef TOutputImage                             OutputImageType;
   typedef typename InputImageType::Pointer         InputImagePointer;
   typedef typename OutputImageType::RegionType     OutputImageRegionType;
 
   typedef typename TInputImage::PixelType PixelType;
 
   typedef TKernel KernelType;
 
   itkSetMacro(Kernel, KernelType);
 
   itkGetConstReferenceMacro(Kernel, KernelType);
 
   itkStaticConstMacro(InputImageDimension, unsigned int,
                       TInputImage::ImageDimension);
   itkStaticConstMacro(OutputImageDimension, unsigned int,
                       TOutputImage::ImageDimension);
   itkStaticConstMacro(KernelDimension, unsigned int,
                       TKernel::NeighborhoodDimension);
 
 #ifdef ITK_USE_CONCEPT_CHECKING
 
   itkConceptMacro(SameTypeCheck,
     (Concept::SameType<PixelType, typename TOutputImage::PixelType>));
   itkConceptMacro(SameDimensionCheck,
     (Concept::SameDimension<InputImageDimension, OutputImageDimension>));
   itkConceptMacro(SameDimensionCheck,
     (Concept::SameDimension<InputImageDimension, KernelDimension>));
   itkConceptMacro(InputLessThanComparableCheck,
     (Concept::LessThanComparable<PixelType>));
   itkConceptMacro(InputGreaterThanComparableCheck,
     (Concept::GreaterThanComparable<PixelType>));
   itkConceptMacro(KernelGreaterThanIntCheck,
     (Concept::GreaterThanComparable<typename TKernel::PixelType, int>));
 
 #endif
 
 protected:
   GrayscaleMorphologicalClosingImageFilter();
   ~GrayscaleMorphologicalClosingImageFilter() {};
   void PrintSelf(std::ostream& os, Indent indent) const;
 
   void GenerateInputRequestedRegion();
 
   void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));
 
   void  GenerateData ();
 
 private:
   GrayscaleMorphologicalClosingImageFilter(const Self&); //purposely not implemented
   void operator=(const Self&); //purposely not implemented
 
   KernelType m_Kernel;
 }; // end of class
 
 } // end namespace itk
   
 #ifndef ITK_MANUAL_INSTANTIATION
 #include "itkGrayscaleMorphologicalClosingImageFilter.txx"
 #endif
 
 #endif
 
 #endif
 
