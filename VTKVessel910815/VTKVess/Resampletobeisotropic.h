#pragma once
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

// Software Guide : BeginCodeSnippet
#include "itkResampleImageFilter.h"
#include "itkRecursiveGaussianImageFilter.h"
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// The resampling filter will need a Transform in order to map point
// coordinates and will need an interpolator in order to compute intensity
// values for the new resampled image. In this particular case we use the
// \doxygen{IdentityTransform} because the image is going to be resampled by
// preserving the physical extent of the sampled region. The Linear
// interpolator is used as a common trade-off, although arguably we should use
// one type of interpolator for the in-plane subsampling process and another
// one for the inter-slice supersampling, but again, one should wonder why to
// enter into technical sophistication here, when what we are doing is to
// cover-up for an improper acquisition of medical data, and we are just trying
// to make it look as if it was correctly acquired.
//
// Software Guide : EndLatex 


// Software Guide : BeginCodeSnippet
#include "itkIdentityTransform.h"
#include "itkLinearInterpolateImageFunction.h"
// Software Guide : EndCodeSnippet




// Software Guide : BeginLatex
//
// Note that, as part of the preprocessing of the image, in this example we are
// also rescaling the range of intensities. This operation has already been
// described as Intensity Windowing. In a real clinical application, this step
// requires careful consideration of the range of intensities that contain
// information about the anatomical structures that are of interest for the
// current clinical application. It practice you may want to remove this step
// of intensity rescaling.
//
// Software Guide : EndLatex 

// Software Guide : BeginCodeSnippet
#include "itkIntensityWindowingImageFilter.h"
// Software Guide : EndCodeSnippet
class Resampletobeisotropic
{
public:
	Resampletobeisotropic(void);
	~Resampletobeisotropic(void);
	int ResampleVolumeToBeIsotropic();

};

