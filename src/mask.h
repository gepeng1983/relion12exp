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

#ifndef MASK_H_
#define MASK_H_

#include "src/multidim_array.h"
#include "src/fftw.h"

// Mask out corners outside sphere (replace by average value)
// Apply a soft mask (raised cosine with cosine_width pixels width)
void softMaskOutsideMap(MultidimArray<double> &vol, double radius = -1., double cosine_width = 3, MultidimArray<double> *Mnoise = NULL);

// Apply a soft mask and set density outside the mask at the average value of those pixels in the original map
void softMaskOutsideMap(MultidimArray<double> &vol, MultidimArray<double> &msk, bool invert_mask = false);


#endif /* MASK_H_ */
