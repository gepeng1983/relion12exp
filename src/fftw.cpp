/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/***************************************************************************
 *
 * Authors:    Roberto Marabini                 (roberto@cnb.csic.es)
 *             Carlos Oscar S. Sorzano          (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/fftw.h"
#include "src/args.h"
#include <string.h>
#include <pthread.h>

static pthread_mutex_t fftw_plan_mutex = PTHREAD_MUTEX_INITIALIZER;

//#define DEBUG_PLANS

// Constructors and destructors --------------------------------------------
FourierTransformer::FourierTransformer()
{
    init();

#ifdef DEBUG_PLANS
    std::cerr << "INIT this= "<<this<< std::endl;
#endif
}

FourierTransformer::~FourierTransformer()
{
    clear();
#ifdef DEBUG_PLANS
    std::cerr << "CLEARED this= "<<this<< std::endl;
#endif
}

FourierTransformer::FourierTransformer(const FourierTransformer& op)
{
	// Clear current object
	clear();
	// New object is an extact copy of op
    *this = op;
}

void FourierTransformer::init()
{
    fReal            = NULL;
    fComplex         = NULL;
    fPlanForward     = NULL;
    fPlanBackward    = NULL;
    dataPtr          = NULL;
    complexDataPtr   = NULL;
    threadsSetOn=false;
    nthreads = 1;
}

void FourierTransformer::clear()
{
    fFourier.clear();
    // Clean-up all other FFTW-allocated memory
    destroyPlans();
    // Initialise all pointers to NULL, set nthreads to 1 and threadsSetOn to false
    init();

}

void FourierTransformer::cleanup()
{
	// First clear object and destroy plans
    clear();
    // Then clean up all the junk fftw keeps lying around
    // SOMEHOW THE FOLLOWING IS NOT ALLOWED WHEN USING MULTPLE TRANSFORMER OBJECTS....
	if (threadsSetOn)
    	fftw_cleanup_threads();
    else
    	fftw_cleanup();
#ifdef DEBUG_PLANS
    std::cerr << "CLEANED-UP this= "<<this<< std::endl;
#endif

}

void FourierTransformer::destroyPlans()
{

    // Anything to do with plans has to be protected for threads!
    pthread_mutex_lock(&fftw_plan_mutex);
    if (fPlanForward !=NULL)
    {
#ifdef DEBUG_PLANS
    	std::cerr << " DESTROY fPlanForward= " << fPlanForward  <<" this= "<<this<<std::endl;
#endif
    	fftw_destroy_plan(fPlanForward);
    }
    if (fPlanBackward!=NULL)
    {
#ifdef DEBUG_PLANS
    	std::cerr << " DESTROY fPlanBackward= " << fPlanBackward  <<" this= "<<this<< std::endl;
#endif
   	fftw_destroy_plan(fPlanBackward);
    }
    pthread_mutex_unlock(&fftw_plan_mutex);

}

void FourierTransformer::setThreadsNumber(int tNumber)
{
    if (tNumber!=1)
    {
        threadsSetOn=true;
        nthreads = tNumber;
        pthread_mutex_lock(&fftw_plan_mutex);
        if(fftw_init_threads()==0)
            REPORT_ERROR("FFTW cannot init threads (setThreadsNumber)");
        fftw_plan_with_nthreads(nthreads);
        pthread_mutex_unlock(&fftw_plan_mutex);
    }
}

// Initialization ----------------------------------------------------------
const MultidimArray<double> &FourierTransformer::getReal() const
{
    return (*fReal);
}

const MultidimArray<std::complex<double> > &FourierTransformer::getComplex() const
{
    return (*fComplex);
}


void FourierTransformer::setReal(MultidimArray<double> &input)
{
    bool recomputePlan=false;
    if (fReal==NULL)
        recomputePlan=true;
    else if (dataPtr!=MULTIDIM_ARRAY(input))
        recomputePlan=true;
    else
        recomputePlan=!(fReal->sameShape(input));
    fFourier.resize(ZSIZE(input),YSIZE(input),XSIZE(input)/2+1);
    fReal=&input;

    if (recomputePlan)
    {
        int ndim=3;
        if (ZSIZE(input)==1)
        {
            ndim=2;
            if (YSIZE(input)==1)
                ndim=1;
        }
        int *N = new int[ndim];
        switch (ndim)
        {
        case 1:
            N[0]=XSIZE(input);
            break;
        case 2:
            N[0]=YSIZE(input);
            N[1]=XSIZE(input);
            break;
        case 3:
            N[0]=ZSIZE(input);
            N[1]=YSIZE(input);
            N[2]=XSIZE(input);
            break;
        }

        // Destroy both forward and backward plans if they already exist
        destroyPlans();
        // Make new plans
        pthread_mutex_lock(&fftw_plan_mutex);
        fPlanForward = fftw_plan_dft_r2c(ndim, N, MULTIDIM_ARRAY(*fReal),
                                         (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_ESTIMATE);
        fPlanBackward = fftw_plan_dft_c2r(ndim, N,
                                          (fftw_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal),
                                          FFTW_ESTIMATE);
        pthread_mutex_unlock(&fftw_plan_mutex);

        if (fPlanForward == NULL || fPlanBackward == NULL)
            REPORT_ERROR("FFTW plans cannot be created");

#ifdef DEBUG_PLANS
        std::cerr << " SETREAL fPlanForward= " << fPlanForward << " fPlanBackward= " << fPlanBackward  <<" this= "<<this<< std::endl;
#endif

        delete [] N;
        dataPtr=MULTIDIM_ARRAY(*fReal);
    }
}

void FourierTransformer::setReal(MultidimArray<std::complex<double> > &input)
{
    bool recomputePlan=false;
    if (fComplex==NULL)
        recomputePlan=true;
    else if (complexDataPtr!=MULTIDIM_ARRAY(input))
        recomputePlan=true;
    else
        recomputePlan=!(fComplex->sameShape(input));
    fFourier.resize(input);
    fComplex=&input;

    if (recomputePlan)
    {
        int ndim=3;
        if (ZSIZE(input)==1)
        {
            ndim=2;
            if (YSIZE(input)==1)
                ndim=1;
        }
        int *N = new int[ndim];
        switch (ndim)
        {
        case 1:
            N[0]=XSIZE(input);
            break;
        case 2:
            N[0]=YSIZE(input);
            N[1]=XSIZE(input);
            break;
        case 3:
            N[0]=ZSIZE(input);
            N[1]=YSIZE(input);
            N[2]=XSIZE(input);
            break;
        }

        pthread_mutex_lock(&fftw_plan_mutex);
        if (fPlanForward!=NULL)
            fftw_destroy_plan(fPlanForward);
        fPlanForward=NULL;
        fPlanForward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(*fComplex),
                                     (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_FORWARD, FFTW_ESTIMATE);
        if (fPlanBackward!=NULL)
            fftw_destroy_plan(fPlanBackward);
        fPlanBackward=NULL;
        fPlanBackward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(fFourier),
                                      (fftw_complex*) MULTIDIM_ARRAY(*fComplex), FFTW_BACKWARD, FFTW_ESTIMATE);
        if (fPlanForward == NULL || fPlanBackward == NULL)
            REPORT_ERROR("FFTW plans cannot be created");
        delete [] N;
        complexDataPtr=MULTIDIM_ARRAY(*fComplex);
        pthread_mutex_unlock(&fftw_plan_mutex);
    }
}

void FourierTransformer::setFourier(MultidimArray<std::complex<double> > &inputFourier)
{
    memcpy(MULTIDIM_ARRAY(fFourier),MULTIDIM_ARRAY(inputFourier),
           MULTIDIM_SIZE(inputFourier)*2*sizeof(double));
}

// Transform ---------------------------------------------------------------
void FourierTransformer::Transform(int sign)
{
    if (sign == FFTW_FORWARD)
    {
        fftw_execute(fPlanForward);

        // Normalisation of the transform
        unsigned long int size=0;
        if(fReal!=NULL)
            size = MULTIDIM_SIZE(*fReal);
        else if (fComplex!= NULL)
            size = MULTIDIM_SIZE(*fComplex);
        else
            REPORT_ERROR("No complex nor real data defined");

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fFourier)
            DIRECT_MULTIDIM_ELEM(fFourier,n) /= size;
    }
    else if (sign == FFTW_BACKWARD)
        fftw_execute(fPlanBackward);

}

void FourierTransformer::FourierTransform()
{
    Transform(FFTW_FORWARD);
}

void FourierTransformer::inverseFourierTransform()
{
    Transform(FFTW_BACKWARD);
}

// Inforce Hermitian symmetry ---------------------------------------------
void FourierTransformer::enforceHermitianSymmetry()
{
    int ndim=3;
    if (ZSIZE(*fReal)==1)
    {
        ndim=2;
        if (YSIZE(*fReal)==1)
            ndim=1;
    }
    long int yHalf=YSIZE(*fReal)/2;
    if (YSIZE(*fReal)%2==0)
        yHalf--;
    long int zHalf=ZSIZE(*fReal)/2;
    if (ZSIZE(*fReal)%2==0)
        zHalf--;
    switch (ndim)
    {
    case 2:
        for (long int i=1; i<=yHalf; i++)
        {
        	long int isym=intWRAP(-i,0,YSIZE(*fReal)-1);
            std::complex<double> mean=0.5*(
                                          DIRECT_A2D_ELEM(fFourier,i,0)+
                                          conj(DIRECT_A2D_ELEM(fFourier,isym,0)));
            DIRECT_A2D_ELEM(fFourier,i,0)=mean;
            DIRECT_A2D_ELEM(fFourier,isym,0)=conj(mean);
        }
        break;
    case 3:
        for (long int k=0; k<ZSIZE(*fReal); k++)
        {
        	long int ksym=intWRAP(-k,0,ZSIZE(*fReal)-1);
            for (long int i=1; i<=yHalf; i++)
            {
            	long int isym=intWRAP(-i,0,YSIZE(*fReal)-1);
                std::complex<double> mean=0.5*(
                                              DIRECT_A3D_ELEM(fFourier,k,i,0)+
                                              conj(DIRECT_A3D_ELEM(fFourier,ksym,isym,0)));
                DIRECT_A3D_ELEM(fFourier,k,i,0)=mean;
                DIRECT_A3D_ELEM(fFourier,ksym,isym,0)=conj(mean);
            }
        }
        for (long int k=1; k<=zHalf; k++)
        {
        	long int ksym=intWRAP(-k,0,ZSIZE(*fReal)-1);
            std::complex<double> mean=0.5*(
                                          DIRECT_A3D_ELEM(fFourier,k,0,0)+
                                          conj(DIRECT_A3D_ELEM(fFourier,ksym,0,0)));
            DIRECT_A3D_ELEM(fFourier,k,0,0)=mean;
            DIRECT_A3D_ELEM(fFourier,ksym,0,0)=conj(mean);
        }
        break;
    }
}

/* FFT Magnitude  ------------------------------------------------------- */
void FFT_magnitude(const MultidimArray< std::complex<double> > &v,
                   MultidimArray<double> &mag)
{
    mag.resize(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
        DIRECT_MULTIDIM_ELEM(mag, n) = abs(DIRECT_MULTIDIM_ELEM(v, n));
}

/* FFT Phase ------------------------------------------------------- */
void FFT_phase(const MultidimArray< std::complex<double> > &v,
               MultidimArray<double> &phase)
{
    phase.resize(v);
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
        DIRECT_MULTIDIM_ELEM(phase, n) = atan2(DIRECT_MULTIDIM_ELEM(v,n).imag(), DIRECT_MULTIDIM_ELEM(v, n).real());
}

void randomizePhasesBeyond(MultidimArray<double> &v, int index)
{
    MultidimArray< std::complex< double > > FT;
    FourierTransformer transformer;

    transformer.FourierTransform(v, FT, false);

    int index2 = index*index;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
    {
    	if (kp*kp + ip*ip + jp*jp >= index2)
    	{
    		double mag = abs(DIRECT_A3D_ELEM(FT, k, i, j));
    		double phas = rnd_unif(0., 2.*PI);
    		double realval = mag * cos(phas);
    		double imagval = mag * sin(phas);
    		DIRECT_A3D_ELEM(FT, k, i, j) = std::complex<double>(realval, imagval);
    	}
    }

    // Inverse transform
    transformer.inverseFourierTransform();

}




// Fourier ring correlation -----------------------------------------------
//#define SAVE_REAL_PART
void frc_dpr(MultidimArray< double > &m1,
             MultidimArray< double > &m2,
             MultidimArray< double > &frc,
             double sampling_rate,
             MultidimArray< double > &freq,
             MultidimArray< double > &frc_noise,
             MultidimArray< double > &dpr,
             MultidimArray< double > &error_l2,
             bool skip_all_but_frc)
{
    if (!m1.sameShape(m2))
        REPORT_ERROR("MultidimArrays have different shapes!");

    MultidimArray< std::complex< double > > FT1;
    FourierTransformer transformer1;
    transformer1.FourierTransform(m1, FT1, false);

    MultidimArray< std::complex< double > > FT2;
    FourierTransformer transformer2;
    transformer2.FourierTransform(m2, FT2, false);

    MultidimArray< int > radial_count(XSIZE(m1)/2+1);
    MultidimArray<double> num, den1, den2;
    Matrix1D<double> f(3);
    num.initZeros(radial_count);
    den1.initZeros(radial_count);
    den2.initZeros(radial_count);

    frc.initZeros(radial_count);
    if (!skip_all_but_frc)
    {
    	dpr.initZeros(radial_count);
		frc_noise.initZeros(radial_count);
    	freq.initZeros(radial_count);
		error_l2.initZeros(radial_count);
    }

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT1)
    {
    	double R = sqrt(kp*kp + ip*ip + jp*jp);
        if (R>0.5*XSIZE(m1))
            continue;
        int idx=ROUND(R);
        std::complex<double> z1=DIRECT_A3D_ELEM(FT1, k, i, j);
        std::complex<double> z2=DIRECT_A3D_ELEM(FT2, k, i, j);
        double absz1=abs(z1);
        double absz2=abs(z2);
        num(idx)+=real(conj(z1) * z2);
        den1(idx)+= absz1*absz1;
        den2(idx)+= absz2*absz2;
        if (!skip_all_but_frc) //this takes to long for a huge volume
        {
        	error_l2(idx)+=norm(z1-z2);
            double phaseDiff=realWRAP(RAD2DEG((atan2(z1.imag(), z1.real())) -
                                              (atan2(z2.imag(), z2.real()))),-180, 180);
            dpr(idx)+=sqrt((absz1+absz2)*phaseDiff*phaseDiff/(absz1+absz2));
        }
        radial_count(idx)++;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D(frc)
    {
        frc(i) = num(i)/sqrt(den1(i)*den2(i));
        if (!skip_all_but_frc)
        {
        	frc_noise(i) = 2 / sqrt((double) radial_count(i));
        	freq(i) = (double) i / (XSIZE(m1) * sampling_rate);
        	error_l2(i)/=radial_count(i);
            dpr(i)/=radial_count(i);
        }

    }
}
void getFSC(MultidimArray< double > & m1,
		    MultidimArray< double > & m2,
		    MultidimArray< double >& fsc)
{
	MultidimArray< double > dummy;
	frc_dpr(m1, m2, fsc, 1., dummy, dummy, dummy, dummy, true);
}

void selfScaleToSizeFourier(long int Ydim, long int Xdim, MultidimArray<double>& Mpmem, int nThreads)
{

    //Mmem = *this
    //memory for fourier transform output
    MultidimArray<std::complex<double> > MmemFourier;
    // Perform the Fourier transform
    FourierTransformer transformerM;
    transformerM.setThreadsNumber(nThreads);
    transformerM.FourierTransform(Mpmem, MmemFourier, false);

    // Create space for the downsampled image and its Fourier transform
    Mpmem.resize(Ydim, Xdim);
    MultidimArray<std::complex<double> > MpmemFourier;
    FourierTransformer transformerMp;
    transformerMp.setReal(Mpmem);
    transformerMp.getFourierAlias(MpmemFourier);
    long int ihalf = XMIPP_MIN((YSIZE(MpmemFourier)/2+1),(YSIZE(MmemFourier)/2+1));
    long int xsize = XMIPP_MIN((XSIZE(MmemFourier)),(XSIZE(MpmemFourier)));
    //Init with zero
    MpmemFourier.initZeros();
    for (long int i=0; i<ihalf; i++)
        for (long int j=0; j<xsize; j++)
            MpmemFourier(i,j)=MmemFourier(i,j);
    for (long int i=YSIZE(MpmemFourier)-1, n=1; n < ihalf-1; i--, n++)
    {
    	long int ip = YSIZE(MmemFourier) - n;
    	for (long int j=0; j<xsize; j++)
            MpmemFourier(i,j)=MmemFourier(ip,j);
    }

    // Transform data
    transformerMp.inverseFourierTransform();
}


// Shift an image through phase-shifts in its Fourier Transform
void shiftImageInFourierTransform(MultidimArray<std::complex<double> > &in,
		                          MultidimArray<std::complex<double> > &out,
								  TabSine &tab_sin, TabCosine &tab_cos,
								  double oridim, double xshift, double yshift, double zshift)
{
	// Only work on non-zero shifts
	if (ABS(xshift) < XMIPP_EQUAL_ACCURACY && ABS(yshift) < XMIPP_EQUAL_ACCURACY && ABS(zshift) < XMIPP_EQUAL_ACCURACY)
	{
		out = in;
	}
	else
	{
		out.resize(in);
		xshift /= -oridim;
		yshift /= -oridim;
		zshift /= -oridim;
		double dotp, a, b, c, d, ac, bd, ab_cd, x, y;
		switch (in.getDim())
		{
		case 1:
			for (long int j = 0; j < XSIZE(in); j++)
			{
				x = j;
				dotp = 2 * PI * (x * xshift);
				a = tab_cos(dotp);
				b = tab_sin(dotp);
				c = DIRECT_A1D_ELEM(in, j).real();
				d = DIRECT_A1D_ELEM(in, j).imag();
				ac = a * c;
				bd = b * d;
				ab_cd = (a + b) * (c + d); // (ab_cd-ac-bd = ad+bc : but needs 4 multiplications)
				DIRECT_A1D_ELEM(out, j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
			}
			break;
		case 2:
			for (long int i=0; i<XSIZE(in); i++)
				for (long int j=0; j<XSIZE(in); j++)
				{
					x = j;
					y = i;
					dotp = 2 * PI * (x * xshift + y * yshift);
					a = tab_cos(dotp);
					b = tab_sin(dotp);
					c = DIRECT_A2D_ELEM(in, i, j).real();
					d = DIRECT_A2D_ELEM(in, i, j).imag();
					ac = a * c;
					bd = b * d;
					ab_cd = (a + b) * (c + d);
					DIRECT_A2D_ELEM(out, i, j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
				}
			for (long int i=YSIZE(in)-1; i>=XSIZE(in); i--)
			{
				y = i - YSIZE(in);
				for (long int j=0; j<XSIZE(in); j++)
				{
					x = j;
					dotp = 2 * PI * (x * xshift + y * yshift);
					a = tab_cos(dotp);
					b = tab_sin(dotp);
					c = DIRECT_A2D_ELEM(in, i, j).real();
					d = DIRECT_A2D_ELEM(in, i, j).imag();
					ac = a * c;
					bd = b * d;
					ab_cd = (a + b) * (c + d);
					DIRECT_A2D_ELEM(out, i, j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
				}
			}
			break;
		case 3:
			REPORT_ERROR("shiftImageInFourierTransform ERROR: Dim=3 not implemented yet...");
			break;
		default:
			REPORT_ERROR("shiftImageInFourierTransform ERROR: dimension should be 1, 2 or 3!");
		}
	}
}
// Shift an image through phase-shifts in its Fourier Transform (without pretabulated sine and cosine)
void shiftImageInFourierTransform(MultidimArray<std::complex<double> > &in,
		                          MultidimArray<std::complex<double> > &out,
		                          double oridim, double xshift, double yshift, double zshift)
{
	// Only work on non-zero shifts
	if (ABS(xshift) < XMIPP_EQUAL_ACCURACY && ABS(yshift) < XMIPP_EQUAL_ACCURACY && ABS(zshift) < XMIPP_EQUAL_ACCURACY)
	{
		out = in;
	}
	else
	{
		out.resize(in);
		xshift /= -oridim;
		yshift /= -oridim;
		zshift /= -oridim;
		double dotp, a, b, c, d, ac, bd, ab_cd, x, y;
		switch (in.getDim())
		{
		case 1:
			for (long int j = 0; j < XSIZE(in); j++)
			{
				x = j;
				dotp = 2 * PI * (x * xshift);
				a = cos(dotp);
				b = sin(dotp);
				c = DIRECT_A1D_ELEM(in, j).real();
				d = DIRECT_A1D_ELEM(in, j).imag();
				ac = a * c;
				bd = b * d;
				ab_cd = (a + b) * (c + d); // (ab_cd-ac-bd = ad+bc : but needs 4 multiplications)
				DIRECT_A1D_ELEM(out, j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
			}
			break;
		case 2:
			for (long int i=0; i<XSIZE(in); i++)
				for (long int j=0; j<XSIZE(in); j++)
				{
					x = j;
					y = i;
					dotp = 2 * PI * (x * xshift + y * yshift);
					a = cos(dotp);
					b = sin(dotp);
					c = DIRECT_A2D_ELEM(in, i, j).real();
					d = DIRECT_A2D_ELEM(in, i, j).imag();
					ac = a * c;
					bd = b * d;
					ab_cd = (a + b) * (c + d);
					DIRECT_A2D_ELEM(out, i, j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
				}
			for (long int i=YSIZE(in)-1; i>=XSIZE(in); i--)
			{
				y = i - YSIZE(in);
				for (long int j=0; j<XSIZE(in); j++)
				{
					x = j;
					dotp = 2 * PI * (x * xshift + y * yshift);
					a = cos(dotp);
					b = sin(dotp);
					c = DIRECT_A2D_ELEM(in, i, j).real();
					d = DIRECT_A2D_ELEM(in, i, j).imag();
					ac = a * c;
					bd = b * d;
					ab_cd = (a + b) * (c + d);
					DIRECT_A2D_ELEM(out, i, j) = std::complex<double>(ac - bd, ab_cd - ac - bd);
				}
			}
			break;
		case 3:
			REPORT_ERROR("shiftImageInFourierTransform ERROR: Dim=3 not implemented yet...");
			break;
		default:
			REPORT_ERROR("shiftImageInFourierTransform ERROR: dimension should be 1, 2 or 3!");
		}
	}
}

void getSpectrum(MultidimArray<double> &Min,
                 MultidimArray<double> &spectrum,
                 int spectrum_type)
{

    MultidimArray<std::complex<double> > Faux;
    int xsize = XSIZE(Min);
    Matrix1D<double> f(3);
    MultidimArray<double> count(xsize);
    FourierTransformer transformer;

    spectrum.initZeros(xsize);
    count.initZeros();
    transformer.FourierTransform(Min, Faux, false);
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
    {
    	long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
        if (spectrum_type == AMPLITUDE_SPECTRUM)
            spectrum(idx) += abs(dAkij(Faux, k, i, j));
        else
            spectrum(idx) += norm(dAkij(Faux, k, i, j));
        count(idx) += 1.;
    }

    for (long int i = 0; i < xsize; i++)
        if (count(i) > 0.)
            spectrum(i) /= count(i);

}

void divideBySpectrum(MultidimArray<double> &Min,
                      MultidimArray<double> &spectrum,
                      bool leave_origin_intact)
{

    MultidimArray<double> div_spec(spectrum);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrum)
    {
        if (ABS(dAi(spectrum,i)) > 0.)
            dAi(div_spec,i) = 1./dAi(spectrum,i);
        else
            dAi(div_spec,i) = 1.;
    }
    multiplyBySpectrum(Min,div_spec,leave_origin_intact);
}

void multiplyBySpectrum(MultidimArray<double> &Min,
                        MultidimArray<double> &spectrum,
                        bool leave_origin_intact)
{

    MultidimArray<std::complex<double> > Faux;
    Matrix1D<double> f(3);
    MultidimArray<double> lspectrum;
    FourierTransformer transformer;
    double dim3 = XSIZE(Min)*YSIZE(Min)*ZSIZE(Min);

    transformer.FourierTransform(Min, Faux, false);
    lspectrum=spectrum;
    if (leave_origin_intact)
        lspectrum(0)=1.;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
    {
    	long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
        dAkij(Faux, k, i, j) *=  lspectrum(idx) * dim3;
    }
    transformer.inverseFourierTransform();

}


void whitenSpectrum(MultidimArray<double> &Min,
                    MultidimArray<double> &Mout,
                    int spectrum_type,
                    bool leave_origin_intact)
{

    MultidimArray<double> spectrum;
    getSpectrum(Min,spectrum,spectrum_type);
    Mout=Min;
    divideBySpectrum(Mout,spectrum,leave_origin_intact);

}

void adaptSpectrum(MultidimArray<double> &Min,
                   MultidimArray<double> &Mout,
                   const MultidimArray<double> &spectrum_ref,
                   int spectrum_type,
                   bool leave_origin_intact)
{

    MultidimArray<double> spectrum;
    getSpectrum(Min,spectrum,spectrum_type);
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrum)
    {
        dAi(spectrum, i) = (dAi(spectrum, i) > 0.) ? dAi(spectrum_ref,i)/ dAi(spectrum, i) : 1.;
    }
    Mout=Min;
    multiplyBySpectrum(Mout,spectrum,leave_origin_intact);

}

/** Kullback-Leibner divergence */
double getKullbackLeibnerDivergence(MultidimArray<std::complex<double> > &Fimg,
		MultidimArray<std::complex<double> > &Fref, MultidimArray<double> &sigma2,
		MultidimArray<double> &p_i, MultidimArray<double> &q_i, int highshell, int lowshell )
{
	// First check dimensions are OK
	if (!Fimg.sameShape(Fref))
		REPORT_ERROR("getKullbackLeibnerDivergence ERROR: Fimg and Fref are not of the same shape.");

	if (highshell < 0)
		highshell = XSIZE(Fimg) - 1;
	if (lowshell < 0)
		lowshell = 0;

	if (highshell > XSIZE(sigma2))
		REPORT_ERROR("getKullbackLeibnerDivergence ERROR: highshell is larger than size of sigma2 array.");

	if (highshell < lowshell)
		REPORT_ERROR("getKullbackLeibnerDivergence ERROR: highshell is smaller than lowshell.");

	// Initialize the histogram
	MultidimArray<int> histogram;
	int histogram_size = 101;
	int histogram_origin = histogram_size / 2;
	double sigma_max = 10.;
	double histogram_factor = histogram_origin / sigma_max;
	histogram.initZeros(histogram_size);

	// This way this will work in both 2D and 3D
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fimg)
	{
		int ires = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
		if (ires >= lowshell && ires <= highshell)
		{
			// Use FT of masked image for noise estimation!
			double diff_real = (DIRECT_A3D_ELEM(Fref, k, i, j)).real() - (DIRECT_A3D_ELEM(Fimg, k, i, j)).real();
			double diff_imag = (DIRECT_A3D_ELEM(Fref, k, i, j)).imag() - (DIRECT_A3D_ELEM(Fimg, k, i, j)).imag();
			double sigma = sqrt(DIRECT_A1D_ELEM(sigma2, ires));

			// Divide by standard deviation to normalise all the difference
			diff_real /= sigma;
			diff_imag /= sigma;

			// Histogram runs from -10 sigma to +10 sigma
			diff_real += sigma_max;
			diff_imag += sigma_max;

			// Make histogram on-the-fly;
			// Real part
			int ihis = ROUND(diff_real * histogram_factor);
			if (ihis < 0)
				ihis = 0;
			else if (ihis >= histogram_size)
				ihis = histogram_size - 1;
			histogram(ihis)++;
			// Imaginary part
			ihis = ROUND(diff_imag * histogram_factor);
			if (ihis < 0)
				ihis = 0;
			else if (ihis > histogram_size)
				ihis = histogram_size;
			histogram(ihis)++;

		}
	}

	// Normalise the histogram and the discretised analytical Gaussian
	double norm = (double)histogram.sum();
	double gaussnorm = 0.;
	for (int i = 0; i < histogram_size; i++)
	{
		double x = (double)i / histogram_factor;
		gaussnorm += gaussian1D(x - sigma_max, 1. , 0.);
	}

	// Now calculate the actual Kullback-Leibner divergence
	double kl_divergence = 0.;
	p_i.resize(histogram_size);
	q_i.resize(histogram_size);
	for (int i = 0; i < histogram_size; i++)
	{
		// Data distribution
		p_i(i) = (double)histogram(i) / norm;
		// Theoretical distribution
		double x = (double)i / histogram_factor;
		q_i(i) = gaussian1D(x - sigma_max, 1. , 0.) / gaussnorm;

		if (p_i(i) > 0.)
			kl_divergence += p_i(i) * log (p_i(i) / q_i(i));
	}
	kl_divergence /= (double)histogram_size;

	return kl_divergence;



}
