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


#ifndef PREPROCESSING_H_
#define PREPROCESSING_H_
#include  <glob.h>
#include  <vector>
#include  <string>
#include  <stdlib.h>
#include  <stdio.h>
#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include <src/fftw.h>
#include <src/time.h>

class Preprocessing
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Output rootname
	FileName fn_in, fn_out;

	////////////////////////////////////// CTF estimation
	// Perform CTF estimation (using CTFFIND3)?
	bool do_ctffind;

	// Filenames (may include wildcards) for all MRC-format micrographs to be processed by CTFFIND3
	FileName fns_ctffind_in;

	// Filenames of all the micrographs to estimate the CTF from
	std::vector<FileName> fn_micrographs;

	// Dimension of squared area of the micrograph to use for CTF estimation
	int ctf_win;

	// CTFFIND3 executable
	FileName fn_ctffind3_exe;

	////// CTFFIND3 parameters
	// Size of the box to calculate FFTw
	double box_size;

	// Minimum and maximum resolution (in A) to be taken into account
	double resol_min, resol_max;

	// Defocus search parameters (in A, positive is underfocus)
	double min_defocus, max_defocus, step_defocus;

	// Voltage (kV)
	double Voltage;

	// Spherical aberration
	double Cs;

	// Amplitude contrast (e.g. 0.07)
	double AmplitudeConstrast;

	// Magnification
	double Magnification;

	// Detector pixel size (um)
	double PixelSize;

	////////////////////////////////////// Extract particles from the micrographs
	// Perform particle extraction?
	bool do_extract;

	// Extract particles from movies instead of single micrographs
	bool do_movie_extract;

	// First frame to extract from movies
	int movie_first_frame;

	// Last frame to extract from movies
	int movie_last_frame;

	// Number of individual movie frames to average over
	int avg_n_frames;

	// Rootname to identify movies, e.g. mic001_movie.mrcs will be the movie of mic001.mrc if fn_movie="movie"
	FileName fn_movie;

	// Filenames (may include wildcards) for all coordinate files to be used for particle extraction
	FileName fns_coords_in;

	// Format of the coordinate files: ximdisp, boxer or xmipp
	FileName fn_coord_format;

	// Filenames of all the coordinate files to use for particle extraction
	std::vector<FileName> fn_coords;

	// Box size to extract the particles in
	int extract_size;

	////////////////////////////////////// Post-extraction image modifications
	// Perform re-scaling of extracted images
	bool do_rescale;
	int scale;

	// Perform re-windowing of extracted images
	bool do_rewindow;
	int window;

	// Perform normalization of the extract images
	bool do_normalise;

	// Perform contrast inversion of the extracted images
	bool do_invert_contrast;

	// Standard deviations to remove black and white dust
	double white_dust_stddev, black_dust_stddev;

	// Radius of a circle in the extracted images outside of which one calculates background mean and stddev
	int bg_radius;

	// Use input stack to perform the image modifications
	FileName fn_operate_in;

	// Name of output stack (only when fn_operate in is given)
	FileName fn_operate_out;

	//////////////////////////////////// Output STAR file
	bool do_join_starfile;

public:
	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// General Running
	void run();

	// Run CTFFIND3 to get CTF parameters
	void runCtffind();

	// Execute CTFFIND3 for a single micrograph
	void executeCtffind3(FileName fn_mic);

	// Get micrograph metadata
	bool getCtffind3Results(FileName fn_mic, double &defU, double &defV, double &defAng, double &CC,
			double &HT, double &CS, double &AmpCnst, double &XMAG, double &DStep);

	// join all STAR files into one
	// This is done separate from runExtractParticles to allow particle extraction to be done in parallel...
	void joinAllStarFiles();

	// Extract particles from the micrographs
	void runExtractParticles();

	// Read coordinates from text files
	void readCoordinates(FileName fn_coord, std::vector< Matrix1D<int> > &all_pos);

	// For the given coordinate file, read the micrograph and/or movie and extract all particles
	void extractParticlesFromFieldOfView(FileName fn_coord);

	// Actually extract particles. This can be from one (average) micrgraph or from a single frame from a movie
	void extractParticlesFromOneFrame(std::vector< Matrix1D<int> > &pos,
			FileName fn_mic, int iframe, int n_frames, FileName fn_stack, MetaDataTable &MD, long int &my_current_nr_images, long int my_total_nr_images,
			double &all_avg, double &all_stddev, double &all_minval, double &all_maxval);

	// Perform per-image operations (e.g. normalise, rescaling, rewindowing and inverting contrast) on an input stack (or STAR file)
	void runOperateOnInputFile(FileName fn_perimage_in);

	// Here normalisation, windowing etc is performed on an individual image and it is written to disc
	void performPerImageOperations(Image<double> &Ipart, FileName fn_stack, int nframes, long int image_nr, long int nr_of_images,
			double &all_avg, double &all_stddev, double &all_minval, double &all_maxval);

	// For image normalization
	void normalise(Image<double> &I);

	void calculateBackgroundAvgStddev(Image<double> &I, double &avg, double &stddev);

	// For dust removal
	void removeDust(Image<double> &I, bool is_white, double thresh, double avg, double stddev);

	// for contrast inversion
	void invert_contrast(Image<double> &I);

	// for image re-scaling
	void rescale(Image<double> &I, int mysize);

	// for image re-windowing
	void rewindow(Image<double> &I, int mysize);

	// Get micrograph name from the rootname
	// The rootname may have an additional string after the uniqye micrograph name
	// That way, multiple "families" of distinct particle types may be extracted from the same micrographs
	FileName getMicrographNameFromRootName(FileName fn_root);

	// The inverse of the function above
	FileName getRootNameFromMicrographName(FileName fn_mic);

};

#endif /* PREPROCESSING_H_ */
