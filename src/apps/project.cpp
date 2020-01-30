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

#include <src/projector.h>
#include <src/backprojector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/ctf.h>
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/metadata_table.h>
#include <src/ml_model.h>
#include <src/exp_model.h>
#include <src/healpix_sampling.h>
class project_parameters
{
public:

	FileName fn_map, fn_ang, fn_out, fn_img, fn_model, fn_sym;
	double rot, tilt, psi, xoff, yoff, angpix, maxres, stddev_white_noise, particle_diameter, ana_prob_range, ana_prob_step;
	int padding_factor;
	int r_max, r_min_nn, interpolator, highres_kl;
    bool do_only_one, do_ctf, ctf_phase_flipped, do_timing, do_add_noise, do_kl_divergence, do_subtract_exp, do_zero_mask, dont_mask, do_ana_prob, do_image_scale, do_nr_outliers;
	// I/O Parser
	IOParser parser;
	MlModel model;
	Experiment data;

// The following is for doing a speed test of the projector and backprojector
//#define TIMING
#ifdef TIMING
    Timer timer;
	int TIMING_PROJ, TIMING_BACKPROJ, TIMING_PREPARE, TIMING_RECONSTRUCT, TIMING_FFTW, TIMING_INVFFTW;
#endif

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_map = parser.getOption("--i", "Input map to be projected");
		fn_out = parser.getOption("--o", "Rootname for output projections", "proj");
       	do_ctf = parser.checkOption("--ctf", "Apply CTF to reference projections");
       	ctf_phase_flipped = parser.checkOption("--ctf_phase_flip", "Flip phases of the CTF in the output projections");
       	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1"));
       	fn_ang = parser.getOption("--ang", "STAR file with orientations for multiple projections (if None, assume single projection)","None");
       	rot = textToFloat(parser.getOption("--rot", "First Euler angle (for a single projection)", "0"));
       	tilt = textToFloat(parser.getOption("--tilt", "Second Euler angle (for a single projection)", "0"));
       	psi = textToFloat(parser.getOption("--psi", "Third Euler angle (for a single projection)", "0"));
       	xoff = textToFloat(parser.getOption("--xoff", "Origin X-offsets (in pixels) (for a single projection)", "0"));
       	yoff = textToFloat(parser.getOption("--yoff", "Origin Y-offsets (in pixels) (for a single projection)", "0"));
       	do_add_noise = parser.checkOption("--add_noise", "Add noise to the output projections (only with --ang)");
       	stddev_white_noise = textToFloat(parser.getOption("--white_noise", "Standard deviation of added white Gaussian noise", "0"));
       	fn_model = parser.getOption("--model_noise", "Model STAR file with power spectra for coloured Gaussian noise", "");
       	do_subtract_exp = parser.checkOption("--subtract_exp", "Subtract experimental image (in --ang) from the projection");
       	do_only_one = (fn_ang == "None");

       	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
       	padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));

       	if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation"))
       		interpolator = NEAREST_NEIGHBOUR;
       	else
       		interpolator = TRILINEAR;

       	int kl_section = parser.addSection("Kullback-Leibner divergence calculations");
       	do_kl_divergence = parser.checkOption("--kl_divergence", "Calculate Kullback-Leibner divergence (with --ang and --model_noise)");
       	highres_kl = textToInteger(parser.getOption("--kl_highest_shell", "Calculate Kullback-Leibner divergence until this resolution shell", "10"));
       	particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images for KL divergence (in Angstroms)", "-1"));
       	do_zero_mask = parser.checkOption("--zero_mask","Mask surrounding background in particles to zero for KL divergence (by default the solvent area is filled with random noise)");
       	dont_mask = parser.checkOption("--dont_mask","Dont mask surrounding background in experimental particles when subtracting those )");
       	do_ana_prob = parser.checkOption("--ana_prob","Analyse angular probability distributions");
       	ana_prob_range  = textToFloat(parser.getOption("--ana_prob_range", "Range (in degrees) to sample rot and tilt for analysis of angular probability distributions", "5."));
       	ana_prob_step  = textToFloat(parser.getOption("--ana_prob_step", "Range (in degrees) to sample rot and tilt for analysis of angular probability distributions", "0.5"));
       	fn_sym = parser.getOption("--sym", "Symmetry group (for analysis of angular probability distributions)", "c1");
       	do_image_scale = parser.checkOption("--image_scale","Analyse per-image scale correction");
       	do_nr_outliers = parser.checkOption("--nr_outliers","Analyse nr of outlier Fourier components");

       	// Hidden
       	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

       	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	}

	void project()
	{

    	MetaDataTable DFo;
    	Matrix2D<double> A3D;
    	FileName fn_expimg;
    	double kl_divergence;

    	MultidimArray<std::complex<double> > F3D, F2D, Fexpimg;
    	MultidimArray<double> Fctf, dummy;
    	Image<double> vol, img, expimg;
    	FourierTransformer transformer, transformer_expimg;

#ifdef TIMING
    	//DIFFF = timer.setNew("difff");
    	TIMING_PROJ = timer.setNew("projector");
    	TIMING_BACKPROJ = timer.setNew("back-projector");
    	TIMING_PREPARE = timer.setNew("prepare-projectpr");
    	TIMING_RECONSTRUCT = timer.setNew("reconstruct");
    	TIMING_FFTW = timer.setNew("fftw-2D");
    	TIMING_INVFFTW = timer.setNew("inverse fftw-2D");
#endif
		std::cerr << " Reading map: " << fn_map << std::endl;
    	vol.read(fn_map);
    	std::cerr << " Done reading map!" << std::endl;

    	if (!do_only_one)
    	{
    		std::cerr << " Reading STAR file with all angles " << fn_ang << std::endl;
    		data.read(fn_ang);
    		std::cerr << " Done reading STAR file!" << std::endl;
    	}

    	// Now that we have the size of the volume, check r_max
   		if (maxres < 0.)
   			r_max = XSIZE(vol());
   		else
   			r_max = CEIL(XSIZE(vol()) * angpix / maxres);

    	// Set right size of F2D and initialize to zero
    	img().resize(YSIZE(vol()), XSIZE(vol()));
    	transformer.setReal(img());
    	transformer.getFourierAlias(F2D);

#ifdef TIMING
		timer.tic(TIMING_PREPARE);
		// Set up the back-projector
		BackProjector backprojector((int)XSIZE(vol()), 3, "C1", interpolator, padding_factor, r_min_nn);
		backprojector.initZeros();
#endif
    	// Set up the projector
    	Projector projector((int)XSIZE(vol()), interpolator, padding_factor, r_min_nn);
    	projector.computeFourierTransformMap(vol(), dummy, 2* r_max);

#ifdef TIMING
		timer.toc(TIMING_PREPARE);
#endif

    	if (do_only_one)
    	{
    		Euler_rotation3DMatrix(rot, tilt, psi, A3D);
    		F2D.initZeros();
    		projector.get2DFourierTransform(F2D, A3D, IS_NOT_INV);
            if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001)
            	shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), -xoff , -yoff );
        	transformer.inverseFourierTransform();
        	// Shift the image back to the center...
        	CenterFFT(img(), false);
        	img.write(fn_out);
        	std::cerr<<" Done writing "<<fn_out<<std::endl;
    	}
    	else
    	{
            init_progress_bar(data.numberOfParticles());
            DFo.clear();
            rot = tilt = psi = xoff = yoff = 0.;

            // KL-divergence preps
            if (do_kl_divergence)
            {
            	if (fn_model != "")
            		model.read(fn_model);
            	else
            		REPORT_ERROR("ERROR: When calculating KL-divergences provide --model_noise");
            }


            // Can only add noise to multiple images
            if (do_add_noise)
            {
            	if (fn_model != "")
            		model.read(fn_model);
            	else if (stddev_white_noise > 0.)
            		stddev_white_noise /= XSIZE(vol()) * sqrt(2); // fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
            	else
            		REPORT_ERROR("ERROR: When adding noise provide either --model_noise or --white_noise");
            }


            long int imgno = 0;
            for (long int part_id = 0; part_id < data.numberOfParticles(); part_id++)
            {

    			for (int iseries = 0; iseries < data.getNrImagesInSeries(part_id); iseries++, imgno++)
    			{

    				MetaDataTable MDimg, MDmic;
    				// Extract the relevant MetaDataTable row from MDimg
    				MDimg = data.getMetaDataImage(part_id, iseries);
    				int mic_id = data.getMicrographId(part_id, iseries);

    				MDimg.getValue(EMDL_ORIENT_ROT, rot);
					MDimg.getValue(EMDL_ORIENT_TILT, tilt);
					MDimg.getValue(EMDL_ORIENT_PSI, psi);
					MDimg.getValue(EMDL_ORIENT_ORIGIN_X, xoff);
					MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, yoff);
					//std::cerr << " rot= " << rot << " tilt= " << tilt << " psi= " << psi << " xoff= " << xoff << " yoff= " << yoff<< std::endl;

					Euler_rotation3DMatrix(rot, tilt, psi, A3D);
					F2D.initZeros();
#ifdef TIMING
		timer.tic(TIMING_PROJ);
#endif
					projector.project(F2D, A3D, IS_NOT_INV);
#ifdef TIMING
		timer.toc(TIMING_PROJ);
#endif

					if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001)
						shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), -xoff , -yoff );

					// Apply CTF if necessary
					CTF ctf;
					if (do_ctf)
					{
						MDmic = data.getMetaDataMicrograph(part_id, iseries);
						ctf.read(MDmic, MDimg);
						Fctf.resize(F2D);
						ctf.getFftwImage(Fctf, XSIZE(vol()), XSIZE(vol()), angpix, ctf_phase_flipped, false, false, true);
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
						{
							DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						}
					}

					if (do_subtract_exp || do_kl_divergence || do_ana_prob)
					{
						MDimg.getValue(EMDL_IMAGE_NAME, fn_expimg);
						expimg.read(fn_expimg);
						expimg().setXmippOrigin();

						/// Mask the image as in ml_optimiser.cpp
						Matrix1D<double> shift(2);
						int my_old_xoff = ROUND(xoff);
						int my_old_yoff = ROUND(yoff);
						XX(shift) = my_old_xoff;
						YY(shift) = my_old_yoff;
						// Shift the image to its center
						selfTranslate(expimg(), shift, DONT_WRAP);
						// Then apply soft circular mask

                                                if (!dont_mask)
                                                {
						if (!do_zero_mask)
						{
							// Create noisy image for outside the mask
							MultidimArray<std::complex<double> > Fnoise;
							MultidimArray<double> Mnoise;
							Mnoise.resize(img());
							transformer.setReal(Mnoise);
							transformer.getFourierAlias(Fnoise);
							// Fill Fnoise with random numbers, use power spectrum of the noise for its variance
							FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fnoise)
							{
								int ires = ROUND( sqrt( (double)(ip * ip + jp * jp) ) );
								if (ires >= 0 && ires < XSIZE(Fnoise))
								{
									double sigma = sqrt(DIRECT_A1D_ELEM(model.sigma2_noise[mic_id], ires));
									DIRECT_A2D_ELEM(Fnoise, i, j).real() = rnd_gaus(0., sigma);
									DIRECT_A2D_ELEM(Fnoise, i, j).imag() = rnd_gaus(0., sigma);
								}
								else
								{
									DIRECT_A2D_ELEM(Fnoise, i, j) = 0.;
								}
							}
							// Back to real space Mnoise
							transformer.inverseFourierTransform();
							Mnoise.setXmippOrigin();
							softMaskOutsideMap(img(), particle_diameter / (2. * angpix), 5, &Mnoise);
						}
						else
						{
							softMaskOutsideMap(expimg(), particle_diameter / (2. * angpix), 5);
						}
                                                }

						//Then shift back
						selfTranslate(expimg(), -shift, DONT_WRAP);
						// Then calculate the FT
						CenterFFT(expimg(), true);
						transformer_expimg.FourierTransform(expimg(), Fexpimg);

						if (do_kl_divergence)
						{

						    //Write an output ASCII file for each image
							std::ofstream  fh;
							fn_img.compose(fn_out,imgno+1,"kldat");
							fh.open((fn_img).c_str(), std::ios::out);
						    if (!fh)
						        REPORT_ERROR( (std::string)"MetaDataTable::write Cannot write to file: " + fn_img);

						    // Calculate KL-divergence per resolution shell first
						    MultidimArray<double> p_i, q_i;
						    for (int ires = 1; ires <= highres_kl; ires++)
						    {
						    	kl_divergence = getKullbackLeibnerDivergence(Fexpimg, F2D, model.sigma2_noise[mic_id], p_i, q_i, ires, ires - 1);
						    	fh << "# ires= "<< ires <<" KL= "<< kl_divergence << std::endl;
						    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p_i)
						    	{
						    		fh  << (n - 50.)/5. << " " << DIRECT_MULTIDIM_ELEM(p_i, n) << " " << DIRECT_MULTIDIM_ELEM(q_i, n) << std::endl;
						    	}
						    }
						    kl_divergence = getKullbackLeibnerDivergence(Fexpimg, F2D, model.sigma2_noise[mic_id], p_i, q_i, highres_kl, 0);
					    	fh << "# ALL res  KL= "<< kl_divergence << std::endl;
					    	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(p_i)
					    	{
					    		fh <<(n - 50.)/5. << " " << DIRECT_MULTIDIM_ELEM(p_i, n) << " " << DIRECT_MULTIDIM_ELEM(q_i, n) << std::endl;
					    	}
							MDimg.setValue(EMDL_PARTICLE_KL_DIVERGENCE, kl_divergence);

							fh.close();

						}
						else if (do_nr_outliers)
						{

							// This way this will work in both 2D and 3D
							double test_fom = 0.;
							FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
							{
								int ires = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
								if (ires <= highres_kl)
								{
									// Use FT of masked image for noise estimation!
									double diff_real = (DIRECT_A3D_ELEM(F2D, k, i, j)).real() - (DIRECT_A3D_ELEM(Fexpimg, k, i, j)).real();
									double diff_imag = (DIRECT_A3D_ELEM(F2D, k, i, j)).imag() - (DIRECT_A3D_ELEM(Fexpimg, k, i, j)).imag();
									double sigma = sqrt(DIRECT_A1D_ELEM(model.sigma2_noise[mic_id], ires));
									// Divide by standard deviation to normalise all the difference
									diff_real /= sigma;
									diff_imag /= sigma;

									if (diff_real > 3.)
										test_fom += 1.;
									if (diff_imag > 3.)
										test_fom += 1.;
								}
							}
							MDimg.setValue(EMDL_PARTICLE_FOM, test_fom);

						}
						else if (do_image_scale)
						{
							// For per-image scale correction
							MultidimArray<double> scale_correction_sumXA, scale_correction_sumA2;
							scale_correction_sumXA.initZeros(highres_kl + 1);
							scale_correction_sumA2.initZeros(highres_kl + 1);
							double tot_sumXA = 0., tot_sumA2 = 0.;
							FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F2D)
							{
								int ires = ROUND( sqrt( (double)(ip*ip + jp*jp) ) );
								if (ires <= highres_kl)
								{
									double sumXA = (DIRECT_A2D_ELEM(F2D, i, j)).real() * (DIRECT_A2D_ELEM(Fexpimg, i, j)).real();
									sumXA += (DIRECT_A2D_ELEM(F2D, i, j)).imag() * (DIRECT_A2D_ELEM(Fexpimg, i, j)).imag();
									double sumA2 = (DIRECT_A2D_ELEM(F2D, i, j)).real() * (DIRECT_A2D_ELEM(F2D, i, j)).real();
									sumA2 += (DIRECT_A2D_ELEM(F2D, i, j)).imag() * (DIRECT_A2D_ELEM(F2D, i, j)).imag();
									DIRECT_A1D_ELEM(scale_correction_sumXA, ires) += sumXA;
									DIRECT_A1D_ELEM(scale_correction_sumA2, ires) += sumA2;
									tot_sumXA += sumXA;
									tot_sumA2 += sumA2;
								}
							}

							//Write an output ASCII file for each image
							std::ofstream  fh;
							fn_img.compose(fn_out,imgno+1,"imagescale");
							fh.open((fn_img).c_str(), std::ios::out);
							if (!fh)
								REPORT_ERROR( (std::string)"MetaDataTable::write Cannot write to file: " + fn_img);

							// Divide sumXA by sumA2
							for (int ires = 0; ires <= highres_kl; ires++)
							{
								DIRECT_A1D_ELEM(scale_correction_sumXA, ires) /= DIRECT_A1D_ELEM(scale_correction_sumA2, ires);
								fh << ires << " " << DIRECT_A1D_ELEM(scale_correction_sumXA, ires) << std::endl;
							}
							tot_sumXA /= tot_sumA2;
							fh << " # "<<tot_sumXA<<std::endl;
							fh.close();
							MDimg.setValue(EMDL_PARTICLE_FOM, tot_sumXA);

						}
						else if (do_ana_prob)
						{

							Image<double> P_sim, P_exp;
							int psize = 2 * CEIL(ana_prob_range / ana_prob_step) + 1;
							P_sim().initZeros(psize, psize);
							P_exp().initZeros(psize, psize);
							P_sim().setXmippOrigin();
							P_exp().setXmippOrigin();

							// Loop over all rot, tilt and psi angles
							MultidimArray<std::complex<double> > lF2D = F2D;
							double min_exp_diff2 = 99.e99, min_sim_diff2 = 99.e99;
							int max_k, max_i, max_j;
						    for (long int kk=STARTINGZ(P_exp()); kk<=FINISHINGZ(P_exp()); kk++) \
						        for (long int ii=STARTINGY(P_exp()); ii<=FINISHINGY(P_exp()); ii++) \
						            for (long int jj=STARTINGX(P_exp()); jj<=FINISHINGX(P_exp()); jj++)
							{
								// Modify angles
						        double lrot = rot + (kk * ana_prob_step);
								double ltilt = tilt + (ii * ana_prob_step);
								double lpsi = psi + (jj * ana_prob_step);

								// Make the reference projection
								Euler_rotation3DMatrix(lrot, ltilt, lpsi, A3D);
								lF2D.initZeros();
								projector.project(lF2D, A3D, IS_NOT_INV);
								if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001)
									shiftImageInFourierTransform(lF2D, lF2D, XSIZE(vol()), -xoff , -yoff );
								if (do_ctf)
								{
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(lF2D)
									{
										DIRECT_MULTIDIM_ELEM(lF2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
									}
								}

								// Calculate exp and sim probability
								double exp_diff2 = 0., sim_diff2 = 0.;
								FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F2D)
								{
									int ires = ROUND( sqrt( (double)(ip*ip + jp*jp) ) );
									if (ires <= highres_kl)
									{
										double exp_diff_real = (DIRECT_A2D_ELEM(lF2D, i, j)).real() - (DIRECT_A2D_ELEM(Fexpimg, i, j)).real();
										double exp_diff_imag = (DIRECT_A2D_ELEM(lF2D, i, j)).imag() - (DIRECT_A2D_ELEM(Fexpimg, i, j)).imag();
										exp_diff2 += (exp_diff_real * exp_diff_real + exp_diff_imag * exp_diff_imag) / (2. * DIRECT_A1D_ELEM(model.sigma2_noise[mic_id], ires));
										double sim_diff_real = (DIRECT_A2D_ELEM(lF2D, i, j)).real() - (DIRECT_A2D_ELEM(F2D, i, j)).real();
										double sim_diff_imag = (DIRECT_A2D_ELEM(lF2D, i, j)).imag() - (DIRECT_A2D_ELEM(F2D, i, j)).imag();
										sim_diff2 += (sim_diff_real * sim_diff_real + sim_diff_imag * sim_diff_imag) / (2. * DIRECT_A1D_ELEM(model.sigma2_noise[mic_id], ires));
									}
								}
								A3D_ELEM(P_exp(), kk, ii, jj) = exp_diff2;
								A3D_ELEM(P_sim(), kk, ii, jj) = sim_diff2;

								// Keep track of minimum diff2 values
								if (exp_diff2 < min_exp_diff2)
								{
									min_exp_diff2 = exp_diff2;
									max_k = kk;
									max_i = ii;
									max_j = jj;
								}
								if (sim_diff2 < min_sim_diff2)
									min_sim_diff2 = sim_diff2;

							}

							// Translate P_exp with its maximum to the center: this allows for small mis-alignments:
						    // it's the shape of the peak we're most interested in, not the position
						    Matrix1D<double> shift_exp(2);
						    XX(shift_exp) = (double)max_j;
						    YY(shift_exp) = (double)max_i;
						    //ZZ(shift_exp) = max_k;
						    selfTranslate(P_sim(), shift_exp, DONT_WRAP, 99.e99);


						    // Now exponentiate, and calculate correlation
				            double sum_exp2 = 0., sum_sim2 = 0., sum_corr = 0.;
				            FOR_ALL_ELEMENTS_IN_ARRAY3D(P_sim())
				            {
								double exp_diff2 = A3D_ELEM(P_exp(), k, i, j) - min_exp_diff2;
								if (exp_diff2 > 700.)
									A3D_ELEM(P_exp(), k, i, j) = 0.;
								else
									A3D_ELEM(P_exp(), k, i, j) = exp(-exp_diff2);

								double sim_diff2 = A3D_ELEM(P_sim(), k, i, j) - min_sim_diff2;
								if (sim_diff2 > 700.)
									A3D_ELEM(P_sim(), k, i, j) = 0.;
								else
									A3D_ELEM(P_sim(), k, i, j) = exp(-sim_diff2);

								sum_corr += A3D_ELEM(P_exp(), k, i, j) * A3D_ELEM(P_sim(), k, i, j);
								sum_exp2 += A3D_ELEM(P_exp(), k, i, j) * A3D_ELEM(P_exp(), k, i, j);
								sum_sim2 += A3D_ELEM(P_sim(), k, i, j) * A3D_ELEM(P_sim(), k, i, j);

				            }

							double ccf = sum_corr/(sqrt(sum_exp2)*sqrt(sum_sim2));
							fn_img.compose(fn_out,imgno+1,"anaprob_exp.spi");
							P_exp.write(fn_img);
							fn_img.compose(fn_out,imgno+1,"anaprob_sim.spi");
							P_sim.write(fn_img);
							std::cerr << " fn_expimg= " << fn_expimg << " imgno= "<<imgno+1<< " ccf= " <<ccf << std::endl;

				            // Now calculate cross-correlation coefficient
				            MDimg.setValue(EMDL_PARTICLE_FOM, ccf);
						}


						if (do_subtract_exp)
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
							{
								DIRECT_MULTIDIM_ELEM(F2D, n) = DIRECT_MULTIDIM_ELEM(Fexpimg, n) - DIRECT_MULTIDIM_ELEM(F2D, n);
							}
						}

					}

					// Apply Gaussian noise
					if (do_add_noise)
					{
						if (fn_model !="")
						{
							// Add coloured noise
							FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F2D)
							{
								int ires = ROUND( sqrt( (double)(ip*ip + jp*jp) ) );
								ires = XMIPP_MIN(ires, model.ori_size/2); // at freqs higher than Nyquist: use last sigma2 value
								double sigma = sqrt(DIRECT_A1D_ELEM(model.sigma2_noise[mic_id], ires));
								DIRECT_A2D_ELEM(F2D, i, j).real() += rnd_gaus(0., sigma);
								DIRECT_A2D_ELEM(F2D, i, j).imag() += rnd_gaus(0., sigma);
							}
						}
						else
						{
							// Add white noise
							FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(F2D)
							{
								DIRECT_A2D_ELEM(F2D, i, j).real() += rnd_gaus(0., stddev_white_noise);
								DIRECT_A2D_ELEM(F2D, i, j).imag() += rnd_gaus(0., stddev_white_noise);
							}
						}
					}

#ifdef TIMING
					else
					{
						Fctf.resize(F2D);
						Fctf.initConstant(1.);
					}

					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}

					timer.tic(TIMING_BACKPROJ);
					// Also time the corresponding backprojection operation
					backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
					timer.toc(TIMING_BACKPROJ);
					timer.tic(TIMING_INVFFTW);
					transformer.inverseFourierTransform();
					timer.toc(TIMING_INVFFTW);
					timer.tic(TIMING_FFTW);
					transformer.FourierTransform();
					timer.toc(TIMING_FFTW);

#else
					DFo.addObject();
					DFo.setObject(MDimg.getObject());
					if (!do_kl_divergence)
					{
						transformer.inverseFourierTransform();
						// Shift the image back to the center...
						CenterFFT(img(), false);
						fn_img.compose(fn_out,imgno+1,"spi");
						img.write(fn_img);
						// Set the image name to the output projection
						DFo.setValue(EMDL_IMAGE_NAME,fn_img);
					}
#endif
    			} // end loop iseries

    			if (part_id%60==0) progress_bar(imgno);

            }
            progress_bar(data.numberOfParticles());

            // Write out STAR file with all information
            fn_img = fn_out + ".star";
            DFo.write(fn_img);
            std::cout<<" Done writing "<<imgno<<" images in "<<fn_img<<std::endl;

#ifdef TIMING
            timer.tic(TIMING_RECONSTRUCT);
            backprojector.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy, false, 1, -1);
            timer.toc(TIMING_RECONSTRUCT);
            vol.write("timing_reconstruction.spi");
            std::cerr << "Written timing_reconstruction.spi" << std::endl;
            timer.printTimes(false);

            timer.toc(TIMING_RECONSTRUCT);
#endif

    	} // end else do_only_one

	}// end project function

};

int main(int argc, char *argv[])
{
	time_config();
	project_parameters prm;

	try
    {
		prm.read(argc, argv);

		prm.project();
    }

    catch (RelionError XE)
    {
        prm.usage();
        std::cout << XE;
        exit(1);
    }

    return 0;

}
