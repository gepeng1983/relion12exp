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

#ifndef __XmippFFTW_H
#define __XmippFFTW_H

#include <complex>
#include <fftw3.h>
#include "src/multidim_array.h"
#include "src/funcs.h"
#include "src/tabfuncs.h"

/** @defgroup FourierW FFTW Fourier transforms
  * @ingroup DataLibrary
  */

/** For all direct elements in the complex array in FFTW format.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * physical definition. It also defines 'kp', 'ip' and 'jp', which are the logical coordinates
 * It also works for 1D or 2D FFTW transforms
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(V)
 * {
 *     int r2 = jp*jp + ip*ip + kp*kp;
 *
 *     std::cout << "element at physical coords: "<< i<<" "<<j<<" "<<k<<" has value: "<<DIRECT_A3D_ELEM(m, k, i, j) << std::endl;
 *     std::cout << "its logical coords are: "<< ip<<" "<<jp<<" "<<kp<<std::endl;
 *     std::cout << "its distance from the origin = "<<sqrt(r2)<<std::endl;
 *
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(V) \
    for (long int k = 0, kp = 0; k<ZSIZE(V); k++, kp = (k < XSIZE(V)) ? k : k - ZSIZE(V)) \
    	for (long int i = 0, ip = 0 ; i<YSIZE(V); i++, ip = (i < XSIZE(V)) ? i : i - YSIZE(V)) \
    		for (long int j = 0, jp = 0; j<XSIZE(V); j++, jp = j)

/** For all direct elements in the complex array in FFTW format.
 *  The same as above, but now only for 2D images (this saves some time as k is not sampled
 */
#define FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(V) \
	for (long int i = 0, ip = 0 ; i<YSIZE(V); i++, ip = (i < XSIZE(V)) ? i : i - YSIZE(V)) \
		for (long int j = 0, jp = 0; j<XSIZE(V); j++, jp = j)

/** FFTW volume element: Logical access.
 *
 * @code
 *
 * FFTW_ELEM(V, -1, -2, 1) = 1;
 * val = FFTW_ELEM(V, -1, -2, 1);
 * @endcode
 */
#define FFTW_ELEM(V, kp, ip, jp) \
    DIRECT_A3D_ELEM((V),((kp < 0) ? (kp + ZSIZE(V)) : (kp)), ((ip < 0) ? (ip + YSIZE(V)) : (ip)), (jp))

/** FFTW 2D image element: Logical access.
 *
 * @code
 *
 * FFTW2D_ELEM(V, --2, 1) = 1;
 * val = FFTW2D_ELEM(V, -2, 1);
 * @endcode
 */
#define FFTW2D_ELEM(V, ip, jp) \
    DIRECT_A2D_ELEM((V), ((ip < 0) ? (ip + YSIZE(V)) : (ip)), (jp))

/** Fourier Transformer class.
 * @ingroup FourierW
 *
 * The memory for the Fourier transform is handled by this object.
 * However, the memory for the real space image is handled externally
 * and this object only has a pointer to it.
 *
 * Here you have an example of use
 * @code
 * FourierTransformer transformer;
 * MultidimArray< std::complex<double> > Vfft;
 * transformer.FourierTransform(V(),Vfft,false);
 * MultidimArray<double> Vmag;
 * Vmag.resize(Vfft);
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(Vmag)
 *     Vmag(k,i,j)=20*log10(abs(Vfft(k,i,j)));
 * @endcode
 */
class FourierTransformer
{
public:
    /** Real array, in fact a pointer to the user array is stored. */
    MultidimArray<double> *fReal;

     /** Complex array, in fact a pointer to the user array is stored. */
    MultidimArray<std::complex<double> > *fComplex;

    /** Fourier array  */
    MultidimArray< std::complex<double> > fFourier;

    /* fftw Forawrd plan */
    fftw_plan fPlanForward;

    /* fftw Backward plan */
    fftw_plan fPlanBackward;

    /* number of threads*/
    int nthreads;

    /* Threads has been used in this program*/
    bool threadsSetOn;

// Public methods
public:
    /** Default constructor */
    FourierTransformer();

    /** Destructor */
    ~FourierTransformer();

    /** Copy constructor
     *
     * The created FourierTransformer is a perfect copy of the input array but with a
     * different memory assignment.
     *
     */
    FourierTransformer(const FourierTransformer& op);

    /** Set Number of threads
     * This function, which should be called once, performs any
     * one-time initialization required to use threads on your
     * system.
     *
     *  The nthreads argument indicates the number of threads you
     *  want FFTW to use (or actually, the maximum number). All
     *  plans subsequently created with any planner routine will use
     *  that many threads. You can call fftw_plan_with_nthreads,
     *  create some plans, call fftw_plan_with_nthreads again with a
     *  different argument, and create some more plans for a new
     *  number of threads. Plans already created before a call to
     *  fftw_plan_with_nthreads are unaffected. If you pass an
     *  nthreads argument of 1 (the default), threads are
     *  disabled for subsequent plans. */
    void setThreadsNumber(int tNumber);

    /** Compute the Fourier transform of a MultidimArray, 2D and 3D.
        If getCopy is false, an alias to the transformed data is returned.
        This is a faster option since a copy of all the data is avoided,
        but you need to be careful that an inverse Fourier transform may
        change the data.
        */
    template <typename T, typename T1>
        void FourierTransform(T& v, T1& V, bool getCopy=true)
        {
            setReal(v);
            Transform(FFTW_FORWARD);
            if (getCopy) getFourierCopy(V);
            else         getFourierAlias(V);
        }

    /** Compute the Fourier transform.
        The data is taken from the matrix with which the object was
        created. */
    void FourierTransform();

    /** Inforce Hermitian symmetry.
        If the Fourier transform risks of losing Hermitian symmetry,
        use this function to renforce it. */
    void enforceHermitianSymmetry();

    /** Compute the inverse Fourier transform.
        The result is stored in the same real data that was passed for
        the forward transform. The Fourier coefficients are taken from
        the internal Fourier coefficients */
    void inverseFourierTransform();

    /** Compute the inverse Fourier transform.
        New data is provided for the Fourier coefficients and the output
        can be any matrix1D, 2D or 3D. It is important that the output
        matrix is already resized to the right size before entering
        in this function. */
    template <typename T, typename T1>
        void inverseFourierTransform(T& V, T1& v)
        {
            setReal(v);
            setFourier(V);
            Transform(FFTW_BACKWARD);
        }

    /** Get Fourier coefficients. */
    template <typename T>
        void getFourierAlias(T& V) {V.alias(fFourier); return;}

    /** Get Fourier coefficients. */
    template <typename T>
        void getFourierCopy(T& V) {
            V.resize(fFourier);
            memcpy(MULTIDIM_ARRAY(V),MULTIDIM_ARRAY(fFourier),
                MULTIDIM_SIZE(fFourier)*2*sizeof(double));
        }

    /** Return a complete Fourier transform (two halves).
    */
    template <typename T>
        void getCompleteFourier(T& V) {
            V.resize(*fReal);
            int ndim=3;
            if (ZSIZE(*fReal)==1)
            {
                ndim=2;
                if (YSIZE(*fReal)==1)
                    ndim=1;
            }
            switch (ndim)
            {
                case 1:
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(V)
                        if (i<XSIZE(fFourier))
                            DIRECT_A1D_ELEM(V,i)=DIRECT_A1D_ELEM(fFourier,i);
                        else
                            DIRECT_A1D_ELEM(V,i)=
                                conj(DIRECT_A1D_ELEM(fFourier,
                                    XSIZE(*fReal)-i));
                    break;
                case 2:
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(V)
                        if (j<XSIZE(fFourier))
                            DIRECT_A2D_ELEM(V,i,j)=
                                DIRECT_A2D_ELEM(fFourier,i,j);
                        else
                            DIRECT_A2D_ELEM(V,i,j)=
                                conj(DIRECT_A2D_ELEM(fFourier,
                                    (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                     XSIZE(*fReal)-j));
                    break;
                case 3:
                    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V)
                        if (j<XSIZE(fFourier))
                            DIRECT_A3D_ELEM(V,k,i,j)=
                                DIRECT_A3D_ELEM(fFourier,k,i,j);
                        else
                            DIRECT_A3D_ELEM(V,k,i,j)=
                                conj(DIRECT_A3D_ELEM(fFourier,
                                    (ZSIZE(*fReal)-k)%ZSIZE(*fReal),
                                    (YSIZE(*fReal)-i)%YSIZE(*fReal),
                                     XSIZE(*fReal)-j));
                    break;
            }
        }

    /** Set one half of the FT in fFourier from the input complete Fourier transform (two halves).
        The fReal and fFourier already should have the right sizes
    */
    template <typename T>
        void setFromCompleteFourier(T& V) {
        int ndim=3;
        if (ZSIZE(*fReal)==1)
        {
            ndim=2;
            if (YSIZE(*fReal)==1)
                ndim=1;
        }
        switch (ndim)
        {
        case 1:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fFourier)
                DIRECT_A1D_ELEM(fFourier,i)=DIRECT_A1D_ELEM(V,i);
            break;
        case 2:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(fFourier)
                DIRECT_A2D_ELEM(fFourier,i,j) = DIRECT_A2D_ELEM(V,i,j);
            break;
        case 3:
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(fFourier)
                DIRECT_A3D_ELEM(fFourier,k,i,j) = DIRECT_A3D_ELEM(V,k,i,j);
            break;
        }
    }

// Internal methods
public:
    /* Pointer to the array of doubles with which the plan was computed */
    double * dataPtr;

    /* Pointer to the array of complex<double> with which the plan was computed */
    std::complex<double> * complexDataPtr;

    /* Initialise all pointers to NULL */
    void init();

    /** Clear object */
    void clear();

    /** This calls fftw_cleanup.
     * NOTE!! When using multiple threads, only ONE thread can call this function, as it cleans up things that are shared among all threads...
 	 *  Therefore, this cleanup is something that needs to be done manually...
    */
    void cleanup();

    /** Destroy both forward and backward fftw planes (mutex locked */
    void destroyPlans();

    /** Computes the transform, specified in Init() function
        If normalization=true the forward transform is normalized
        (no normalization is made in the inverse transform)
        If normalize=false no normalization is performed and therefore
        the image is scaled by the number of pixels.
    */
    void Transform(int sign);

    /** Get the Multidimarray that is being used as input. */
    const MultidimArray<double> &getReal() const;
    const MultidimArray<std::complex<double> > &getComplex() const;

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<double> &img);

    /** Set a Multidimarray for input.
        The data of img will be the one of fComplex. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<std::complex<double> > &img);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(MultidimArray<std::complex<double> > &imgFourier);
};

/** FFT Magnitude 1D
 * @ingroup FourierOperations
 */
void FFT_magnitude(const MultidimArray< std::complex< double > > & v,
                   MultidimArray< double >& mag);

/** FFT Phase 1D
 * @ingroup FourierOperations
 */
void FFT_phase(const MultidimArray< std::complex< double > > & v,
               MultidimArray< double >& phase);

// Randomize phases beyond the given shell (index)
void randomizePhasesBeyond(MultidimArray<double> &I, int index);

/** Center an array, to have its origin at the origin of the FFTW
 *
 */
template <typename T>
void CenterFFT(MultidimArray< T >& v, bool forward)
{
    if ( v.getDim() == 1 )
    {
        // 1D
        MultidimArray< T > aux;
        int l, shift;

        l = XSIZE(v);
        aux.resize(l);
        shift = (int)(l / 2);

        if (!forward)
            shift = -shift;

        // Shift the input in an auxiliar vector
        for (int i = 0; i < l; i++)
        {
            int ip = i + shift;

            if (ip < 0)
                ip += l;
            else if (ip >= l)
                ip -= l;

            aux(ip) = DIRECT_A1D_ELEM(v, i);
        }

        // Copy the vector
        for (int i = 0; i < l; i++)
            DIRECT_A1D_ELEM(v, i) = DIRECT_A1D_ELEM(aux, i);
    }
    else if ( v.getDim() == 2 )
    {
        // 2D
        MultidimArray< T > aux;
        int l, shift;

        // Shift in the X direction
        l = XSIZE(v);
        aux.resize(l);
        shift = (int)(l / 2);

        if (!forward)
            shift = -shift;

        for (int i = 0; i < YSIZE(v); i++)
        {
            // Shift the input in an auxiliar vector
            for (int j = 0; j < l; j++)
            {
                int jp = j + shift;

                if (jp < 0)
                    jp += l;
                else if (jp >= l)
                    jp -= l;

                aux(jp) = DIRECT_A2D_ELEM(v, i, j);
            }

            // Copy the vector
            for (int j = 0; j < l; j++)
                DIRECT_A2D_ELEM(v, i, j) = DIRECT_A1D_ELEM(aux, j);
        }

        // Shift in the Y direction
        l = YSIZE(v);
        aux.resize(l);
        shift = (int)(l / 2);

        if (!forward)
            shift = -shift;

        for (int j = 0; j < XSIZE(v); j++)
        {
            // Shift the input in an auxiliar vector
            for (int i = 0; i < l; i++)
            {
                int ip = i + shift;

                if (ip < 0)
                    ip += l;
                else if (ip >= l)
                    ip -= l;

                aux(ip) = DIRECT_A2D_ELEM(v, i, j);
            }

            // Copy the vector
            for (int i = 0; i < l; i++)
                DIRECT_A2D_ELEM(v, i, j) = DIRECT_A1D_ELEM(aux, i);
        }
    }
    else if ( v.getDim() == 3 )
    {
        // 3D
        MultidimArray< T > aux;
        int l, shift;

        // Shift in the X direction
        l = XSIZE(v);
        aux.resize(l);
        shift = (int)(l / 2);

        if (!forward)
            shift = -shift;

        for (int k = 0; k < ZSIZE(v); k++)
            for (int i = 0; i < YSIZE(v); i++)
            {
                // Shift the input in an auxiliar vector
                for (int j = 0; j < l; j++)
                {
                    int jp = j + shift;

                    if (jp < 0)
                        jp += l;
                    else if (jp >= l)
                        jp -= l;

                    aux(jp) = DIRECT_A3D_ELEM(v, k, i, j);
                }

                // Copy the vector
                for (int j = 0; j < l; j++)
                    DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, j);
            }

        // Shift in the Y direction
        l = YSIZE(v);
        aux.resize(l);
        shift = (int)(l / 2);

        if (!forward)
            shift = -shift;

        for (int k = 0; k < ZSIZE(v); k++)
            for (int j = 0; j < XSIZE(v); j++)
            {
                // Shift the input in an auxiliar vector
                for (int i = 0; i < l; i++)
                {
                    int ip = i + shift;

                    if (ip < 0)
                        ip += l;
                    else if (ip >= l)
                        ip -= l;

                    aux(ip) = DIRECT_A3D_ELEM(v, k, i, j);
                }

                // Copy the vector
                for (int i = 0; i < l; i++)
                    DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, i);
            }

        // Shift in the Z direction
        l = ZSIZE(v);
        aux.resize(l);
        shift = (int)(l / 2);

        if (!forward)
            shift = -shift;

        for (int i = 0; i < YSIZE(v); i++)
            for (int j = 0; j < XSIZE(v); j++)
            {
                // Shift the input in an auxiliar vector
                for (int k = 0; k < l; k++)
                {
                    int kp = k + shift;
                    if (kp < 0)
                        kp += l;
                    else if (kp >= l)
                        kp -= l;

                    aux(kp) = DIRECT_A3D_ELEM(v, k, i, j);
                }

                // Copy the vector
                for (int k = 0; k < l; k++)
                    DIRECT_A3D_ELEM(v, k, i, j) = DIRECT_A1D_ELEM(aux, k);
            }
    }
    else
    {
    	v.printShape();
    	REPORT_ERROR("CenterFFT ERROR: Dimension should be 1, 2 or 3");
    }
}



// Window an FFTW-centered Fourier-transform to a given size
template<class T>
void windowFourierTransform(MultidimArray<T > &in,
			  			    MultidimArray<T > &out,
			  			    long int newdim)
{
	// Check size of the input array
	if (YSIZE(in) > 1 && YSIZE(in)/2 + 1 != XSIZE(in))
		REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
	long int newhdim = newdim/2 + 1;

	// If same size, just return input
	if (newhdim == XSIZE(in))
	{
		out = in;
		return;
	}

	// Otherwise apply a windowing operation
	// Initialise output array
	switch (in.getDim())
	{
	case 1:
		out.initZeros(newhdim);
		break;
	case 2:
		out.initZeros(newdim, newhdim);
		break;
	case 3:
		out.initZeros(newdim, newdim, newhdim);
		break;
	default:
    	REPORT_ERROR("windowFourierTransform ERROR: dimension should be 1, 2 or 3!");
    }
	if (newhdim > XSIZE(in))
	{
		long int max_r2 = (XSIZE(in) -1) * (XSIZE(in) - 1);
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(in)
		{
			// Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
			if (kp*kp + ip*ip + jp*jp <= max_r2)
				FFTW_ELEM(out, kp, ip, jp) = FFTW_ELEM(in, kp, ip, jp);
		}
	}
	else
	{
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(out)
		{
			FFTW_ELEM(out, kp, ip, jp) = FFTW_ELEM(in, kp, ip, jp);
		}
	}
}

// A resize operation in Fourier-space (i.e. changing the sampling of the Fourier Transform) by windowing in real-space
// If recenter=true, the real-space array will be recentered to have its origin at the origin of the FT
template<class T>
void resizeFourierTransform(MultidimArray<T > &in,
			  			    MultidimArray<T > &out,
			  			    long int newdim, bool do_recenter=true)
{
	// Check size of the input array
	if (YSIZE(in) > 1 && YSIZE(in)/2 + 1 != XSIZE(in))
		REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
	long int newhdim = newdim/2 + 1;
	long int olddim = 2* (XSIZE(in) - 1);

	// If same size, just return input
	if (newhdim == XSIZE(in))
	{
		out = in;
		return;
	}

	// Otherwise apply a windowing operation
	MultidimArray<std::complex<double> > Fin;
	MultidimArray<double> Min;
	FourierTransformer transformer;
	long int x0, y0, z0, xF, yF, zF;
	x0 = y0 = z0 = FIRST_XMIPP_INDEX(newdim);
	xF = yF = zF = LAST_XMIPP_INDEX(newdim);

	// Initialise output array
	switch (in.getDim())
	{
	case 1:
		Min.resize(olddim);
		y0=yF=z0=zF=0;
		break;
	case 2:
		Min.resize(olddim, olddim);
		z0=zF=0;
		break;
	case 3:
		Min.resize(olddim, olddim, olddim);
		break;
	default:
    	REPORT_ERROR("resizeFourierTransform ERROR: dimension should be 1, 2 or 3!");
    }

	// This is to handle double-valued input arrays
	Fin.resize(ZSIZE(in), YSIZE(in), XSIZE(in));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(in)
	{
		DIRECT_MULTIDIM_ELEM(Fin, n) = DIRECT_MULTIDIM_ELEM(in, n);
	}
	transformer.inverseFourierTransform(Fin, Min);
	Min.setXmippOrigin();
	if (do_recenter)
		CenterFFT(Min, false);

	// Now do the actual windowing in real-space
	Min.window(z0, y0, x0, zF, yF, xF);
	Min.setXmippOrigin();

	// If upsizing: mask the corners to prevent aliasing artefacts
	if (newdim > olddim)
	{
		FOR_ALL_ELEMENTS_IN_ARRAY3D(Min)
		{
			if (k*k + i*i + j*j > olddim*olddim/4)
			{
				A3D_ELEM(Min, k, i, j) = 0.;
			}
		}
	}

	// Recenter FFT back again
	if (do_recenter)
		CenterFFT(Min, true);

	// And do the inverse Fourier transform
	transformer.clear();
	transformer.FourierTransform(Min, out);
}


/** Autocorrelation function of a Xmipp vector
 * @ingroup FourierOperations
 *
 * Fast calculation of the autocorrelation vector of a given one using Fast
 * Fourier Transform. (Using the correlation theorem)
 */
template <typename T>
void auto_correlation_vector(const MultidimArray< T > & Img, MultidimArray< double >& R)
{
    Img.checkDimension(1);

    // Compute the Fourier Transform
    MultidimArray< std::complex< double > > FFT1;
    FourierTransformer transformer1;
    R=Img;
    transformer1.FourierTransform(R, FFT1, false);

    // Multiply FFT1 * FFT1'
    double dSize=XSIZE(Img);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(FFT1)
        FFT1(i) *= dSize * conj(FFT1(i));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Autocorrelation function of a Xmipp vector
 * @ingroup FourierOperations
 *
 * Fast calcuation of the correlation matrix on two matrices using Fast Fourier
 * Transform. (Using the correlation theorem). The output matrix must be already
 * resized
 */
template <typename T>
void correlation_vector(const MultidimArray< T > & m1,
                        const MultidimArray< T > & m2,
                        MultidimArray< double >& R)
{
    m1.checkDimension(2);
    m2.checkDimension(2);

    // Compute the Fourier Transforms
    MultidimArray< std::complex< double > > FFT1, FFT2;
    FourierTransformer transformer1, transformer2;
    R=m1;
    transformer1.FourierTransform(R, FFT1, false);
    transformer2.FourierTransform((MultidimArray<T> &)m2, FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=XSIZE(m1);
    FOR_ALL_ELEMENTS_IN_ARRAY1D(FFT1)
        FFT1(i) *= dSize * conj(FFT2(i));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Compute the correlation vector without using Fourier.
 * @ingroup FourierOperations
 *
 *  The results are the same as the previous ones but this function
 *  is threadsafe while the previous one is not.
 *
 *  It is assumed that the two vectors v1, and v2 are of the same size. */
template <class T>
void correlation_vector_no_Fourier(const MultidimArray<T> &v1, const MultidimArray<T> &v2,
    MultidimArray<T> &result)
{
    v1.checkDimension(1);
    v2.checkDimension(1);

    result.initZeros(v1);
    result.setXmippOrigin();
    long int N=XSIZE(v1)-1;
    FOR_ALL_ELEMENTS_IN_ARRAY1D(result)
        for (long int k=0; k<XSIZE(v1); k++)
            result(i)+=DIRECT_A1D_ELEM(v1,intWRAP(k+i,0,N))*
                       DIRECT_A1D_ELEM(v2,k);
    STARTINGX(result)=0;
}


/** Correlation of two nD images
 * @ingroup FourierOperations
 *
 * Fast calcuation of the correlation matrix on two matrices using Fast Fourier
 * Transform. (Using the correlation theorem). The output matrix must be already
 * resized
 */
template <typename T>
void correlation_matrix(const MultidimArray< T > & m1,
                        const MultidimArray< T > & m2,
                        MultidimArray< double >& R)
{
    // Compute the Fourier Transforms
    MultidimArray< std::complex< double > > FFT1, FFT2;
    FourierTransformer transformer1, transformer2;
    R=m1;
    transformer1.FourierTransform(R, FFT1, false);
    transformer2.FourierTransform((MultidimArray<T> &)m2, FFT2, false);

    // Multiply FFT1 * FFT2'
    double dSize=MULTIDIM_SIZE(R);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(FFT1)
        FFT1(k, i, j) *= dSize * conj(FFT2(k, i, j));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}

/** Autocorrelation function of an image
 * @ingroup FourierOperations
 *
 * Fast calcuation of the autocorrelation matrix of a given one using Fast
 * Fourier Transform. (Using the correlation theorem)
 */
template <typename T>
void auto_correlation_matrix(const MultidimArray< T > & Img, MultidimArray< double >& R)
{
    // Compute the Fourier Transform
    MultidimArray< std::complex< double > > FFT1;
    FourierTransformer transformer1;
    R=Img;
    transformer1.FourierTransform(R, FFT1);

    // Multiply FFT1 * FFT1'
    double dSize=MULTIDIM_SIZE(Img);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(FFT1)
        FFT1(k, i, j) *= dSize * conj(FFT1(k, i, j));

    // Invert the product, in order to obtain the correlation image
    transformer1.inverseFourierTransform();

    // Center the resulting image to obtain a centered autocorrelation
    CenterFFT(R, true);
}


/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * @ingroup FourierOperations
 */
void frc_dpr(MultidimArray< double > & m1,
             MultidimArray< double > & m2,
             MultidimArray< double >& frc,
             double sampling_rate,
             MultidimArray< double >& freq,
             MultidimArray< double >& frc_noise,
             MultidimArray< double >& dpr,
             MultidimArray< double >& error_l2,
             bool skip_all_but_frc=false);

/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * @ingroup FourierOperations
 * Simpler I/O than above.
 */
void getFSC(MultidimArray< double > & m1,
		    MultidimArray< double > & m2,
		    MultidimArray< double >& fsc);

/** Scale matrix using Fourier transform
 * @ingroup FourierOperations
 * Ydim and Xdim define the output size, Mpmem is the matrix to scale
 */
void selfScaleToSizeFourier(long int Ydim, long int Xdim, MultidimArray<double>& Mpmem, int nthreads=1);

// Shift an image through phase-shifts in its Fourier Transform
// Note that in and out may be the same array, in that case in is overwritten with the result
// if oridim is in pixels, xshift, yshift and zshift should be in pixels as well!
// or both can be in Angstroms
void shiftImageInFourierTransform(MultidimArray<std::complex<double> > &in,
								  MultidimArray<std::complex<double> > &out,
								  TabSine &tab_sin, TabCosine &tab_cos,
								  double oridim, double xshift, double yshift = 0., double zshift = 0.);

// Shift an image through phase-shifts in its Fourier Transform (without tabulated sine and cosine)
// Note that in and out may be the same array, in that case in is overwritten with the result
// if oridim is in pixels, xshift, yshift and zshift should be in pixels as well!
// or both can be in Angstroms
void shiftImageInFourierTransform(MultidimArray<std::complex<double> > &in,
						          MultidimArray<std::complex<double> > &out,
								  double oridim, double xshift, double yshift = 0., double zshift = 0.);

#define POWER_SPECTRUM 0
#define AMPLITUDE_SPECTRUM 1

/** Get the amplitude or power_class spectrum of the map in Fourier space.
 * @ingroup FourierOperations
    i.e. the radial average of the (squared) amplitudes of all Fourier components
*/
void getSpectrum(MultidimArray<double> &Min,
                 MultidimArray<double> &spectrum,
                 int spectrum_type=POWER_SPECTRUM);

/** Divide the input map in Fourier-space by the spectrum provided.
 * @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void divideBySpectrum(MultidimArray<double> &Min,
                      MultidimArray<double> &spectrum,
                      bool leave_origin_intact=false);

/** Multiply the input map in Fourier-space by the spectrum provided.
 * @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void multiplyBySpectrum(MultidimArray<double> &Min,
                        MultidimArray<double> &spectrum,
                        bool leave_origin_intact=false);

/** Perform a whitening of the amplitude/power_class spectrum of a 3D map
 * @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void whitenSpectrum(MultidimArray<double> &Min,
                    MultidimArray<double> &Mout,
                    int spectrum_type=AMPLITUDE_SPECTRUM,
                    bool leave_origin_intact=false);

/** Adapts Min to have the same spectrum as spectrum_ref
 * @ingroup FourierOperations
    If only_amplitudes==true, the amplitude rather than the power_class spectrum will be equalized
*/
void adaptSpectrum(MultidimArray<double> &Min,
                   MultidimArray<double> &Mout,
                   const MultidimArray<double> &spectrum_ref,
                   int spectrum_type=AMPLITUDE_SPECTRUM,
                   bool leave_origin_intact=false);

/** Kullback-Leibner divergence */
double getKullbackLeibnerDivergence(MultidimArray<std::complex<double> > &Fimg,
		MultidimArray<std::complex<double> > &Fref, MultidimArray<double> &sigma2,
		MultidimArray<double> &p_i, MultidimArray<double> &q_i,
		int highshell = -1, int lowshell = -1);


#endif
