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
#include "src/mask.h"

// Mask out corners outside sphere (replace by average value)
// Apply a soft mask (raised cosine with cosine_width pixels width)
void softMaskOutsideMap(MultidimArray<double> &vol, double radius, double cosine_width, MultidimArray<double> *Mnoise)
{

	vol.setXmippOrigin();
	double r, radius_p, raisedcos, sum_bg = 0., sum = 0.;
	if (radius < 0)
		radius = (double)XSIZE(vol)/2.;
	radius_p = radius + cosine_width;


	if (Mnoise == NULL)
	{
		// Calculate average background value
		FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
		{
			r = sqrt((double)(k*k + i*i + j*j));
			if (r < radius)
				continue;
			else if (r > radius_p)
			{
				sum    += 1.;
				sum_bg += A3D_ELEM(vol, k, i, j);
			}
			else
			{
				raisedcos = 0.5 + 0.5 * cos(PI * (radius_p - r) / cosine_width );
				sum += raisedcos;
				sum_bg += raisedcos * A3D_ELEM(vol, k, i, j);
			}
		}
		sum_bg /= sum;
	}

	// Apply noisy or average background value
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		r = sqrt((double)(k*k + i*i + j*j));
		if (r < radius)
		{
			continue;
		}
		else if (r > radius_p)
		{
			A3D_ELEM(vol, k, i, j) = (Mnoise == NULL) ? sum_bg : A3D_ELEM(*Mnoise, k, i, j);
		}
		else
		{
			raisedcos = 0.5 + 0.5 * cos(PI * (radius_p - r) / cosine_width );
			double add = (Mnoise == NULL) ?  sum_bg : A3D_ELEM(*Mnoise, k, i, j);
			A3D_ELEM(vol, k, i, j) = (1 - raisedcos) * A3D_ELEM(vol, k, i, j) + raisedcos * add;
		}
	}

}

void softMaskOutsideMap(MultidimArray<double> &vol, MultidimArray<double> &msk, bool invert_mask)
{

	if (msk.computeMax() > 1. || msk.computeMin() < 0.)
	{
		std::cerr << " msk.computeMax()= " << msk.computeMax() << " msk.computeMin()= " << msk.computeMin() << std::endl;
		REPORT_ERROR("ERROR: Values in the solvent mask should be between zero and one.");
	}
	if (!(msk.sameShape(vol)))
		REPORT_ERROR("ERROR: Solvent mask does not have the same size as the reference vol.");

	// Replace solvent by the average value in the solvent region
	double sum = 0.;
	double sum_bg = 0.;
	double solv;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(msk)
	{
		solv = (invert_mask) ? DIRECT_A3D_ELEM(msk, k, i, j) : 1. - DIRECT_A3D_ELEM(msk, k, i, j);
		sum    += solv;
		sum_bg += solv * DIRECT_A3D_ELEM(vol, k, i, j);
	}
	sum_bg /= sum;

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(msk)
	{
		solv = (invert_mask) ? DIRECT_A3D_ELEM(msk, k, i, j) : 1. - DIRECT_A3D_ELEM(msk, k, i, j);
		DIRECT_A3D_ELEM(vol, k, i, j) = ( 1. - solv) * DIRECT_A3D_ELEM(vol, k, i, j) + solv * sum_bg;
	}


}

