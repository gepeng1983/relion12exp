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
#include "src/preprocessing.h"

void Preprocessing::read(int argc, char **argv, int rank)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_out = parser.getOption("--o", "Output rootname", "particles");
	// Dont allow for directories here!
	fn_out = fn_out.removeDirectories();

	int ctf_section = parser.addSection("CTF estimation");
	fns_ctffind_in = parser.getOption("--ctffind","All micrographs for which to run Niko Grigorieff's CTFFIND3 (may contain wildcards e.g. \"mics/*.mrc\")","");
        // Use a smaller squared part of the micrograph to estimate CTF (e.g. to avoid film labels...)
        ctf_win =  textToInteger(parser.getOption("--ctfWin", "Size (in pixels) of a centered, squared window to use for CTF-estimation", "-1"));

	fn_ctffind3_exe = parser.getOption("--ctffind3_exe","Location of ctffind3 executable (or through RLN_CTFFIND3_EXECUTABLE environment variable)","");
	// First parameter line in CTFFIND
	Cs = textToFloat(parser.getOption("--CS", "Spherical Aberration (mm) ","2.0"));
	Voltage = textToFloat(parser.getOption("--HT", "Voltage (kV)","300"));
	AmplitudeConstrast = textToFloat(parser.getOption("--AmpCnst", "Amplitude constrast", "0.1"));
	Magnification = textToFloat(parser.getOption("--XMAG", "Magnification", "60000"));
	PixelSize = textToFloat(parser.getOption("--DStep", "Detector pixel size (um)", "14"));
	// Second parameter line in CTFFIND
	box_size = textToFloat(parser.getOption("--Box", "Size of the boxes to calculate FFTs", "512"));
	resol_min = textToFloat(parser.getOption("--ResMin", "Minimum resolution (in A) to include in calculations", "100"));
	resol_max = textToFloat(parser.getOption("--ResMax", "Maximum resolution (in A) to include in calculations", "7"));
	min_defocus = textToFloat(parser.getOption("--dFMin", "Minimum defocus value (in A) to search", "10000"));
	max_defocus = textToFloat(parser.getOption("--dFMax", "Maximum defocus value (in A) to search", "50000"));
	step_defocus = textToFloat(parser.getOption("--FStep", "defocus step size (in A) for search", "250"));

	int particle_section = parser.addSection("Particle selection");
	fns_coords_in = parser.getOption("--coord_files","The coordinate files for all particles to be output (may contain wildcards e.g. \"mics/*.box\")","");

	int extract_section = parser.addSection("Particle extraction");
	do_extract = parser.checkOption("--extract", "Extract all particles from the micrographs");
	fn_coord_format = parser.getOption("--coord_format","Format of the coordinate files, choose from: ximdisp, boxer and xmipp2","ximdisp");
	extract_size = textToInteger(parser.getOption("--extract_size", "Size of the box to extract the particles in (in pixels)", "-1"));
	do_movie_extract = parser.checkOption("--extract_movies", "Also extract particles from movie stacks (e.g. from DDDs)");
	avg_n_frames = textToInteger(parser.getOption("--avg_movie_frames", "Average over this number of individual movie frames", "1"));
	movie_first_frame = textToInteger(parser.getOption("--first_movie_frame", "Extract from this movie frame onwards", "1"));
	movie_first_frame--; // (start counting at 0, not 1)
	movie_last_frame = textToInteger(parser.getOption("--last_movie_frame", "Extract until this movie frame (default=all movie frames)", "0"));
	movie_last_frame--; // (start counting at 0, not 1)
	fn_movie = parser.getOption("--movie_rootname", "Common name to relate movies to the single micrographs (e.g. mic001_movie.mrcs related to mic001.mrc)", "movie");

	int perpart_section = parser.addSection("Particle operations");
	scale  = textToInteger(parser.getOption("--scale", "Re-scale the particles to this size (in pixels)", "-1"));
	window  = textToInteger(parser.getOption("--window", "Re-window the particles to this size (in pixels)", "-1"));
	do_normalise = parser.checkOption("--norm", "Normalise the background to average zero and stddev one");
	bg_radius = textToInteger(parser.getOption("--bg_radius", "Radius of the circular mask that will be used to define the background area (in pixels)", "-1"));
	white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
	black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));
	do_invert_contrast = parser.checkOption("--invert_contrast", "Invert the contrast in the input images");
	fn_operate_in = parser.getOption("--operate_on", "The operations above are applied to all extracted particles. Use this option to operate on an input stack/STAR file", "");

	// Initialise verb for non-parallel execution
	verb = 1;

	if (!( checkParameter(argc, argv, "--o") || checkParameter(argc, argv, "--operate_on") ))
		REPORT_ERROR("Provide either --o or --operate_on");

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void Preprocessing::usage()
{
	parser.writeUsage(std::cerr);
}

void Preprocessing::initialise()
{

	if (fns_ctffind_in == "" && fns_coords_in == "" && fn_operate_in == "" && !do_join_starfile)
		REPORT_ERROR("Provide either --ctffind, --extract, --join_starfile or --operate_on");

	// Set up which micrographs to estimate CTFs from
	if (fns_ctffind_in != "")
	{
		do_ctffind = true;
		if (fn_ctffind3_exe == "")
		{
			char * penv;
			penv = getenv ("RLN_CTFFIND3_EXECUTABLE");
			if (penv!=NULL)
				fn_ctffind3_exe = (std::string)penv;
		}
		if (verb > 0)
		{
			std::cout << " Using CTFFINDs executable in: " << fn_ctffind3_exe << std::endl;
			std::cout << " to estimate CTF parameters for the following micrographs: " << std::endl;
		}

		// Get the filenames of all micrographs to be processed by CTFFIND
		fns_ctffind_in.globFiles(fn_micrographs);
		if (verb > 0)
			for(unsigned  int  i = 0; i < fn_micrographs.size(); ++i)
				std::cout << "  * " << fn_micrographs[i] << std::endl;
	}
	else
		do_ctffind = false;

	// Set up which coordinate files to extract particles from (or to join STAR file for)
	do_join_starfile = false;
	if (fns_coords_in != "")
	{
		do_join_starfile = true;
		if (do_extract && verb > 0)
		{
			if (!(fn_coord_format == "boxer" || fn_coord_format == "xmipp2" || fn_coord_format == "ximdisp"))
				REPORT_ERROR("Preprocessing::initialise ERROR: unrecognised coordinate file format: " + fn_coord_format);

			if (extract_size < 0)
				REPORT_ERROR("Preprocessing::initialise ERROR: please provide the size of the box to extract particle using --extract_size ");

			std::cout << " Extract particles based on the following (" << fn_coord_format << ") coordinate files: " << std::endl;
		}
		else if (!do_extract && verb > 0)
		{
			std::cout << " Creating output STAR file for particles based on the following coordinate files: " << std::endl;
		}

		// Get the filenames of all micrographs to be processed by CTFFIND
		fns_coords_in.globFiles(fn_coords);
		if (verb > 0)
			for(unsigned  int  i = 0; i < fn_coords.size(); ++i)
				std::cout << "  * " << fn_coords[i] << std::endl;
	}


	if (do_extract || fn_operate_in != "")
	{
		// Check whether to do re-scaling
		do_rescale = (scale > 0);
		if (do_rescale && scale%2 != 0)
			REPORT_ERROR("ERROR: only re-scaling to even-sized images is allowed in RELION...");

		// Check whether to do re-windowing
		do_rewindow = (window > 0);
		if (do_rewindow && window%2 != 0)
			REPORT_ERROR("ERROR: only re-windowing to even-sized images is allowed in RELION...");

		// Check for bg_radius in case of normalisation
		if (do_normalise && bg_radius < 0)
			REPORT_ERROR("ERROR: please provide a radius for a circle that defines the background area when normalising...");
	}

}

void Preprocessing::run()
{

	if (do_ctffind)
		runCtffind();

	if (do_extract)
		runExtractParticles();

	if (do_join_starfile)
		joinAllStarFiles();

	if (fn_operate_in != "")
		runOperateOnInputFile(fn_operate_in);

	if (verb > 0)
		std::cout << " Done!" <<std::endl;
}

void Preprocessing::runCtffind()
{

	int barstep;
	if (verb > 0)
	{
		std::cout << " Estimating CTF parameters using CTFFIND3 ..." << std::endl;
		init_progress_bar(fn_micrographs.size());
		barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
	}

	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
    {
    	if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

    	executeCtffind3(fn_micrographs[imic]);
	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

}

void Preprocessing::executeCtffind3(FileName fn_mic)
{

	FileName fn_root = fn_mic.withoutExtension();
	FileName fn_script = fn_root + "_ctffind3.com";
	FileName fn_log = fn_root + "_ctffind3.log";
	FileName fn_ctf = fn_root + ".ctf";
        FileName fn_mic_win;

	std::ofstream  fh;
	fh.open((fn_script).c_str(), std::ios::out);
	if (!fh)
	 REPORT_ERROR( (std::string)"Preprocessing::execute_ctffind3 cannot create file: " + fn_script);

        // If given, then put a square window of ctf_win on the micrograph for CTF estimation
        if (ctf_win > 0)
        {
            // Window micrograph to a smaller, squared sub-micrograph to estimate CTF on
            fn_mic_win = fn_root + "_win.mrc";
            // Read in micrograph, window and write out again
            Image<double> I;
            I.read(fn_mic);
            I().setXmippOrigin();
            I().window(FIRST_XMIPP_INDEX(ctf_win), FIRST_XMIPP_INDEX(ctf_win), LAST_XMIPP_INDEX(ctf_win), LAST_XMIPP_INDEX(ctf_win));
            // Calculate mean, stddev, min and max
            double avg, stddev, minval, maxval;
            I().computeStats(avg, stddev, minval, maxval);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
            I.write(fn_mic_win);
        }
        else
            fn_mic_win = fn_mic;

	// Write script to run x3dna
	fh << "#!/usr/bin/env csh"<<std::endl;
	fh << fn_ctffind3_exe << " > " << fn_log << " << EOF"<<std::endl;
	fh << fn_mic_win << std::endl;
	fh << fn_ctf << std::endl;
	// CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]
	fh << Cs << ", " << Voltage << ", " << AmplitudeConstrast << ", " << Magnification << ", " << PixelSize<< std::endl;
	// Box, ResMin[A], ResMax[A], dFMin[A], dFMax[A], FStep
	fh << box_size << ", " << resol_min << ", " << resol_max << ", " << min_defocus << ", " << max_defocus << ", " << step_defocus << std::endl;
	fh <<"EOF"<<std::endl;
	fh.close();

	// Execute ctffind3
	if (!system(NULL))
	 REPORT_ERROR("There is a problem with the system call to run ctffind3");
	FileName fn_cmnd = "csh "+ fn_script;
	system ( fn_cmnd.c_str() );

        // Remove windowed file again
        if (ctf_win > 0)
        {
            if( remove( fn_mic_win.c_str() ) != 0 )
                REPORT_ERROR( "Error deleting windowed micrograph file..." );
        }

}

bool Preprocessing::getCtffind3Results(FileName fn_microot, double &defU, double &defV, double &defAng, double &CC,
		double &HT, double &CS, double &AmpCnst, double &XMAG, double &DStep)
{

	FileName fn_log = fn_microot + "_ctffind3.log";
	std::ifstream in(fn_log.data(), std::ios_base::in);
    if (in.fail())
    	return false;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Proceed until the next "Final values" statement
    // The loop statement may be necessary for data blocks that have a list AND a table inside them
    bool Final_is_found = false;
    bool Cs_is_found = false;
    std::string line;
    std::vector<std::string> words;
    while (getline(in, line, '\n'))
    {
        // Find data_ lines

    	 if (line.find("CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]") != std::string::npos)
    	 {
    		 Cs_is_found = true;
    		 getline(in, line, '\n');
    		 tokenize(line, words);
    		 if (words.size() < 5)
    			 REPORT_ERROR("ERROR: Unexpected number of words on data line with CS[mm], HT[kV], etc in " + fn_log);
    		 CS = textToFloat(words[0]);
    		 HT = textToFloat(words[1]);
    		 AmpCnst = textToFloat(words[2]);
    		 XMAG = textToFloat(words[3]);
    		 DStep = textToFloat(words[4]);
    	 }

    	if (line.find("Final Values") != std::string::npos)
        {
        	Final_is_found = true;
            tokenize(line, words);
            if (words.size() < 6)
            	REPORT_ERROR("ERROR: Unexpected number of words on Final values line in " + fn_log);
            defU = textToFloat(words[0]);
            defV = textToFloat(words[1]);
            defAng = textToFloat(words[2]);
            CC = textToFloat(words[3]);
        }
    }
    if (!Cs_is_found)
    	REPORT_ERROR("ERROR: cannot find line with Cs[mm], HT[kV], etc values in " + fn_log);
    if (!Final_is_found)
    	REPORT_ERROR("ERROR: cannot find line with Final values in " + fn_log);

    in.close();

    return true;

}


void Preprocessing::joinAllStarFiles()
{

	MetaDataTable MDout, MDonestack;

	// Fix order of the labels in the output file
	MDout.addLabel(EMDL_MICROGRAPH_NAME);
	MDout.addLabel(EMDL_IMAGE_COORD_X);
	MDout.addLabel(EMDL_IMAGE_COORD_Y);
	MDout.addLabel(EMDL_IMAGE_NAME);

	std::cout << " Joining all metadata in one STAR file..." << std::endl;
	bool has_other_ctfs, has_this_ctf;
	double defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep;
	has_other_ctfs = false;
	for (long int ipos = 0; ipos < fn_coords.size(); ipos++)
    {
		FileName fn_mic = getMicrographNameFromRootName((fn_coords[ipos]).withoutExtension());
		// If the micrograph did not exist, particles were not extracted: just continue with the next one
		if (fn_mic == "")
			continue;
		FileName fn_microot = getRootNameFromMicrographName(fn_mic);

		// Gather the results from ctffind
		has_this_ctf = getCtffind3Results(fn_microot, defU, defV, defAng, CC,
				HT, CS, AmpCnst, XMAG, DStep);

			// Re-scaled detector pixel size
		if (has_this_ctf && do_rescale)
                    DStep *= (double)extract_size/(double)scale;

		if (ipos == 0 && has_this_ctf)
		{
			// Set has_other_ctfs to true if the first micrograph has a logfile
			// In that case add all CTF labels to the output MDtable
			has_other_ctfs = true;
			MDout.addLabel(EMDL_CTF_DEFOCUSU);
			MDout.addLabel(EMDL_CTF_DEFOCUSV);
			MDout.addLabel(EMDL_CTF_DEFOCUS_ANGLE);
			MDout.addLabel(EMDL_CTF_VOLTAGE);
			MDout.addLabel(EMDL_CTF_CS);
			MDout.addLabel(EMDL_CTF_Q0);
			MDout.addLabel(EMDL_CTF_MAGNIFICATION);
			MDout.addLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE);
			MDout.addLabel(EMDL_CTF_FOM);
		}

		if (!has_this_ctf && has_other_ctfs)
			REPORT_ERROR("joinAllStarFiles%ERROR: Exiting because of missing CTFFIND3 logfiles for micrograph " + fn_mic);

		if (has_this_ctf && !has_other_ctfs)
			REPORT_ERROR("joinAllStarFiles%ERROR: Exiting because of missing CTFFIND3 logfiles ...");


		FileName fn_star = "Particles/" + fn_mic.withoutExtension() + "_" + fn_out + ".star";
		MDonestack.read(fn_star);

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDonestack)
		{
			// This was double, right?
			//MDonestack.setValue(EMDL_MICROGRAPH_NAME, fn_mic);

			if (has_this_ctf)
			{
				MDonestack.setValue(EMDL_CTF_DEFOCUSU, defU);
				MDonestack.setValue(EMDL_CTF_DEFOCUSV, defV);
				MDonestack.setValue(EMDL_CTF_DEFOCUS_ANGLE, defAng);
				MDonestack.setValue(EMDL_CTF_VOLTAGE, HT);
				MDonestack.setValue(EMDL_CTF_Q0, AmpCnst);
				MDonestack.setValue(EMDL_CTF_CS, CS);
				MDonestack.setValue(EMDL_CTF_MAGNIFICATION, XMAG);
				MDonestack.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, DStep);
				MDonestack.setValue(EMDL_CTF_FOM, CC);
			}

			// Add the entire line to the joined STAR-file
			MDout.addObject(MDonestack.getObject());

			// Remove the individual star file to clean up
			//remove(fn_star.c_str());
		}

    }

	// Write out the joined star file
	FileName fn_tmp;
	if (do_movie_extract)
		fn_tmp = fn_out + "_" + fn_movie + ".star";
	else
		fn_tmp = fn_out + ".star";
	MDout.write(fn_tmp);
	std::cout << " Written out STAR file with all particles in " << fn_tmp << std::endl;

}

void Preprocessing::runExtractParticles()
{

	int barstep;
	if (verb > 0)
	{
		std::cout << " Extracting particles from the micrographs ..." << std::endl;
		init_progress_bar(fn_coords.size());
		barstep = XMIPP_MAX(1, fn_coords.size() / 60);
	}

	FileName fn_olddir = "";
	for (long int ipos = 0; ipos < fn_coords.size(); ipos++)
    {

		FileName fn_dir = "Particles/" + fn_coords[ipos].beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}

		if (verb > 0 && ipos % barstep == 0)
			progress_bar(ipos);

    	extractParticlesFromFieldOfView(fn_coords[ipos]);
	}

	if (verb > 0)
		progress_bar(fn_coords.size());

}


void Preprocessing::readCoordinates(FileName fn_coord, std::vector< Matrix1D<int> > &all_pos)
{
    all_pos.clear();
    std::ifstream in(fn_coord.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "Preprocessing::readCoordinates ERROR: File " + fn_coord + " does not exists" );

    // Start reading the ifstream at the top
    in.seekg(0);
    std::string line;
    Matrix1D<int> onepos(2);
    int n = 0;
    while (getline(in, line, '\n'))
    {
    	std::vector<std::string> words;
    	tokenize(line, words);

    	if (fn_coord_format=="ximdisp" || fn_coord_format=="xmipp2")
    	{
    		if (n==0)
    		{
    			// first ximdisp line should be: "x         y      density"
    			// first xmipp2  line should be: "# <X position> <Y position>"
    			if (words.size() < 3)
    	           	REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on first line of " + fn_coord);

    			if (fn_coord_format=="ximdisp" && (words[0] != "x" || words[1] != "y"  || words[2] != "density"))
    				std::cout << " Warning first line of ximdisp coordinates file " << fn_coord << " is not: x         y      density" << std::endl;

    			if (fn_coord_format=="xmipp2" &&  (words[0] != "#" || words[1] != "<X" || words[2] != "position>"))
    				std::cout << " Warning first line of xmipp2 coordinates file " << fn_coord << " is not: # <X position> <Y position>" << std::endl;
    		}
    		else
    		{
    			// All other lines contain x, y as first two entries (ignore anything else)
    			if (words.size() < 2)
    				REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on data line of " + fn_coord);

    			XX(onepos) = textToInteger(words[0]);
    			YY(onepos) = textToInteger(words[1]);
    			all_pos.push_back(onepos);
    		}
    	}
    	else if (fn_coord_format=="boxer")
    	{
   			if (words.size() < 4)
    			REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on data line of " + fn_coord);

   			XX(onepos) = textToInteger(words[0]) + textToInteger(words[2]) / 2; // boxer stores corner and particle together with box size
    		YY(onepos) = textToInteger(words[1]) + textToInteger(words[3]) / 2;
    		all_pos.push_back(onepos);
    	}
        n++;

    }
    in.close();

}

void Preprocessing::extractParticlesFromFieldOfView(FileName fn_coord)
{
	std::vector< Matrix1D<int> > pos;

	// Read in the coordinates file
	readCoordinates(fn_coord, pos);
	// Warn for small groups
	if (pos.size() < 10)
	{
		std:: cout << "WARNING: there are only " << pos.size() << " particles in " << fn_coord <<". Consider joining multiple micrographs into one group. "<< std::endl;
	}

	// Check the micrograph exists
	FileName fn_mic;
	fn_mic = getMicrographNameFromRootName(fn_coord.withoutExtension());
	// Return if the micrograph does not exist
	if (fn_mic == "")
	{
		std::cout << "WARNING: cannot find micrograph for coordinate file " << fn_coord << " with " << pos.size() << " particles" << std::endl;
		return;
	}

	// Name of the output stack
	// Add the same root as the output STAR file (that way one could extract two "families" of different particle stacks)
	FileName fn_stack = "Particles/" + fn_mic.withoutExtension() + "_" + fn_out + ".mrcs";
	// Name of this micrographs STAR file
	FileName fn_star = "Particles/" + fn_mic.withoutExtension() + "_" + fn_out + ".star";

	// Read the header of the micrograph to see how many frames there are.
	Image<double> Imic;
	Imic.read(fn_mic, false, -1, false, true); // readData = false, select_image = -1, mapData= false, is_2D = true);

	int xdim, ydim, zdim;
	long int ndim;
	Imic.getDimensions(xdim, ydim, zdim, ndim);

	// Just to be sure...
	if (do_movie_extract && ndim < 2)
		std::cout << "WARNING: movie " << fn_mic << " does not have multiple frames..." << std::endl;

	MetaDataTable MD;
	long int my_current_nr_images = 0;
	long int my_total_nr_images;
	double all_avg = 0;
	double all_stddev = 0;
	double all_minval = 99.e99;
	double all_maxval = -99.e99;

	// The total number of images to be extracted
	my_total_nr_images = pos.size() * ndim;

	// To deal with default movie_last_frame value
	if (movie_last_frame < 0)
		movie_last_frame = ndim - 1;

	for (long int iframe = movie_first_frame; iframe <= movie_last_frame; iframe += avg_n_frames)
	{
		extractParticlesFromOneFrame(pos, fn_mic, iframe, ndim, fn_stack, MD, my_current_nr_images, my_total_nr_images,
				all_avg, all_stddev, all_minval, all_maxval);

		// Keep track of total number of images extracted thus far
		my_current_nr_images += pos.size();

	}

	MD.setName("images");
	MD.write(fn_star);


}

// Actually extract particles. This can be from one (average) micrograph or from a single movie frame
void Preprocessing::extractParticlesFromOneFrame(std::vector< Matrix1D<int> > &pos,
		FileName fn_mic, int iframe, int n_frames,
		FileName fn_stack, MetaDataTable &MD, long int &my_current_nr_images, long int my_total_nr_images,
		double &all_avg, double &all_stddev, double &all_minval, double &all_maxval)
{

	Image<double> Ipart, Imic, Itmp;


	FileName fn_frame;
	// If movies, then average over avg_n_frames
	if (n_frames > 1)
	{
		for (int ii =0; ii < avg_n_frames; ii++)
		{
			int iiframe = iframe + ii;
			// If we run over the size of the movie, then discard these frames
			if (iiframe >= n_frames)
				return;
			fn_frame.compose(iiframe + 1, fn_mic);
			if (ii==0)
			{
				Imic.read(fn_frame, true, -1, false, true); // readData = true, select_image = -1, mapData= false, is_2D = true
			}
			else
			{
				Itmp.read(fn_frame, true, -1, false, true); // readData = true, select_image = -1, mapData= false, is_2D = true
				Imic() += Itmp();
				Itmp.clear();
			}
		}
	}
	else
	{
		fn_frame = fn_mic;
		Imic.read(fn_frame);
	}

	// Now window all particles from the micrograph
	for (long int ipos = 0; ipos < pos.size(); ipos++)
	{
		long int xpos = XX(pos[ipos]);
		long int ypos = YY(pos[ipos]);
		long int x0, xF, y0, yF;
		x0 = xpos + FIRST_XMIPP_INDEX(extract_size);
		xF = xpos + LAST_XMIPP_INDEX(extract_size);
		y0 = ypos + FIRST_XMIPP_INDEX(extract_size);
		yF = ypos + LAST_XMIPP_INDEX(extract_size);

		// extract one particle in Ipart
		Imic().window(Ipart(), y0, x0, yF, xF);

		// Discard particles that are completely outside the micrograph and print a warning
		if (yF < 0 || y0 >= YSIZE(Imic()) || xF < 0 || x0 >= XSIZE(Imic()))
		{
			std::cout << " Warning! ignoring particle " << ipos + 1 << " on micrograph " << fn_frame << " that lies completely outside the micrograph..." << std::endl;
		}
		else
		{
			// Check boundaries: fill pixels outside the boundary with the nearest ones inside
			// This will create lines at the edges, rather than zeros
			Ipart().setXmippOrigin();

			// X-boundaries
			if (x0 < 0 || xF >= XSIZE(Imic()) )
			{
				FOR_ALL_ELEMENTS_IN_ARRAY2D(Ipart())
				{
					if (j + xpos < 0)
						A2D_ELEM(Ipart(), i, j) = A2D_ELEM(Ipart(), i, -xpos);
					else if (j + xpos >= XSIZE(Imic()))
						A2D_ELEM(Ipart(), i, j) = A2D_ELEM(Ipart(), i, XSIZE(Imic()) - xpos - 1);
				}
			}

			// Y-boundaries
			if (y0 < 0 || yF >= YSIZE(Imic()))
			{
				FOR_ALL_ELEMENTS_IN_ARRAY2D(Ipart())
				{
					if (i + ypos < 0)
						A2D_ELEM(Ipart(), i, j) = A2D_ELEM(Ipart(), -ypos, j);
					else if (i + ypos >= YSIZE(Imic()))
						A2D_ELEM(Ipart(), i, j) = A2D_ELEM(Ipart(), YSIZE(Imic()) - ypos - 1, j);
				}
			}

			// performPerImageOperations will also append the particle to the output stack in fn_stack
			performPerImageOperations(Ipart, fn_stack, n_frames, my_current_nr_images + ipos, my_total_nr_images, all_avg, all_stddev, all_minval, all_maxval);

			// Also store all the particles information in the STAR file
			MD.addObject();
			FileName fn_img;
			//TODO: check this!
			fn_img.compose(my_current_nr_images + ipos + 1, fn_stack); // start image counting in stacks at 1!
			if (do_movie_extract)
			{
				FileName fn_part, fn_mic;
				long int dum;
				fn_frame.decompose(dum, fn_mic);
				fn_part.compose(ipos + 1,  "Particles/" + getRootNameFromMicrographName(fn_mic)); // start image counting in stacks at 1!
				// for automated re-alignment of particles in relion_refine: have rlnParticleName equal to rlnImageName in non-movie star file
				fn_part += "_" + fn_out + ".mrcs";
				MD.setValue(EMDL_PARTICLE_NAME, fn_part);
			}
			MD.setValue(EMDL_IMAGE_NAME, fn_img);
			MD.setValue(EMDL_MICROGRAPH_NAME, fn_frame);
			MD.setValue(EMDL_IMAGE_COORD_X, (double)XX(pos[ipos]));
			MD.setValue(EMDL_IMAGE_COORD_Y, (double)YY(pos[ipos]));
		}
	}


}


void Preprocessing::runOperateOnInputFile(FileName fn_operate_on)
{
	Image<double> Ipart, Iout;
	MetaDataTable MD;
	long int Nimg;

	FileName fn_stack = fn_out+".mrcs";
	FileName fn_star = fn_out+".star";

	if (fn_operate_on.isStarFile())
	{
		// Readt STAR file and get total number of images
		MD.read(fn_operate_on);
		Nimg = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			Nimg++;
		}
		MD.firstObject(); // reset pointer to the first object in the table
	}
	else
	{
		// Read the header of the stack to see how many images there
		Iout.read(fn_operate_on, false);
		Nimg = NSIZE(Iout());
	}

	double all_avg = 0;
	double all_stddev = 0;
	double all_minval = 99.e99;
	double all_maxval = -99.e99;
	init_progress_bar(Nimg);
	int barstep = XMIPP_MAX(1, Nimg / 120);
	for (long int i = 0; i < Nimg; i++)
	{
		FileName fn_tmp;

		// Read in individual miages from the stack
		Ipart.clear();
		if (fn_operate_on.isStarFile())
		{
			MD.getValue(EMDL_IMAGE_NAME, fn_tmp);
			Ipart.read(fn_tmp);

			// Set the new name at this point in the MDtable, e.g. as 000001@out.mrcs
			fn_tmp.compose(i+1,fn_stack);
			MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
			if (i < (Nimg - 1))
				MD.nextObject();
		}
		else
		{
			Ipart.read(fn_operate_on, true, i);
			// Set the new name at this point in the MDtable, e.g. as 000001@out.mrcs
			fn_tmp.compose(i+1,fn_stack);
			MD.addObject();
			MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
		}

		performPerImageOperations(Ipart, fn_stack, 1, i, Nimg, all_avg, all_stddev, all_minval, all_maxval);

		// progress bar
		if (i % barstep == 0) progress_bar(i);

	}
	progress_bar(Nimg);

	std::cout << " Done writing to " << fn_stack << std::endl;
	MD.setName("images");
	MD.write(fn_star);
	std::cout << " Also written a STAR file with the image names as " << fn_star << std::endl;

}


void Preprocessing::performPerImageOperations(Image<double> &Ipart, FileName fn_stack, int nframes, long int image_nr, long int nr_of_images,
		double &all_avg, double &all_stddev, double &all_minval, double &all_maxval)
{

	Ipart().setXmippOrigin();

	if (do_rescale) rescale(Ipart, scale);

	if (do_rewindow) rewindow(Ipart, window);

	if (do_normalise) normalise(Ipart);

	if (do_invert_contrast) invert_contrast(Ipart);

	// For movies: multiple the image intensities by sqrt(nframes) so the stddev in the average of the normalised frames is again 1
	if (nframes > 1)
		Ipart() *= sqrt((double)nframes/(double)avg_n_frames);

	// Calculate mean, stddev, min and max
	double avg, stddev, minval, maxval;
	Ipart().computeStats(avg, stddev, minval, maxval);

	// Keep track of overall statistics
	all_minval = XMIPP_MIN(minval, all_minval);
	all_maxval = XMIPP_MAX(maxval, all_maxval);
	all_avg	+= avg;
	all_stddev += stddev*stddev;

	// Last particle: reset the min, max, avg and stddev values in the main header
	if (image_nr == (nr_of_images - 1))
	{
		all_avg /= nr_of_images;
		all_stddev = sqrt(all_stddev/nr_of_images);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, all_minval);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, all_maxval);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, all_avg);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, all_stddev);
	}

	// Write this particle to the stack on disc
	// First particle: write stack in overwrite mode, from then on just append to it
	if (image_nr == 0)
		Ipart.write(fn_stack, -1, (nr_of_images > 1), WRITE_OVERWRITE);
	else
		Ipart.write(fn_stack, -1, false, WRITE_APPEND);


}

void Preprocessing::normalise(Image<double> &I)
{
	int bg_radius2 = bg_radius * bg_radius;
	double avg, stddev;

	if (2*bg_radius > XSIZE(I()))
		REPORT_ERROR("Preprocessing::normalise ERROR: 2*bg_radius is larger than image size!");

	// Calculate initial avg and stddev values
	calculateBackgroundAvgStddev(I, avg, stddev);

	// Remove white and black noise
	if (white_dust_stddev > 0.)
		removeDust(I, true, white_dust_stddev, avg, stddev);
	if (black_dust_stddev > 0.)
		removeDust(I, false, black_dust_stddev, avg, stddev);

	// If some dust was removed: recalculate avg and stddev
	if (white_dust_stddev > 0. || black_dust_stddev > 0.)
		calculateBackgroundAvgStddev(I, avg, stddev);

	// Subtract avg and divide by stddev for all pixels
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I())
		DIRECT_MULTIDIM_ELEM(I(), n) = (DIRECT_MULTIDIM_ELEM(I(), n) - avg) / stddev;
}

void Preprocessing::calculateBackgroundAvgStddev(Image<double> &I, double &avg, double &stddev)
{
	int bg_radius2 = bg_radius * bg_radius;
	double n = 0.;
	avg = 0.;
	stddev = 0.;

	// Calculate avg in the background pixels
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		if (k*k + i*i + j*j > bg_radius2)
		{
			avg += A3D_ELEM(I(), k, i, j);
			n += 1.;
		}
	}
	avg /= n;

	// Calculate stddev in the background pixels
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		if (k*k + i*i + j*j > bg_radius2)
		{
			double aux = A3D_ELEM(I(), k, i, j) - avg;
			stddev += aux * aux;
		}
	}
	stddev = sqrt(stddev/n);
}

void Preprocessing::removeDust(Image<double> &I, bool is_white, double thresh, double avg, double stddev)
{
	FOR_ALL_ELEMENTS_IN_ARRAY3D(I())
	{
		double aux =  A3D_ELEM(I(), k, i, j);
		if (is_white && aux - avg > thresh * stddev)
			A3D_ELEM(I(), k, i, j) = rnd_gaus(avg, stddev);
		else if (!is_white && aux - avg < -thresh * stddev)
			A3D_ELEM(I(), k, i, j) = rnd_gaus(avg, stddev);
	}
}

void Preprocessing::invert_contrast(Image<double> &I)
{
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I())
	{
		DIRECT_MULTIDIM_ELEM(I(), n) *= -1;
	}
}

void Preprocessing::rescale(Image<double> &I, int mysize)
{
	// 3D scaling not implemented for now...
	if (I().getDim() == 3)
		REPORT_ERROR("Rescaling of 3D images has not been implemented yet...");

	int olddim = XSIZE(I());
	selfScaleToSizeFourier(mysize, mysize, I());

	// Also modify the scale in the MDmainheader (if present)
	double oldscale, newscale;
    if (I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, oldscale))
    {
    	newscale = oldscale * (double)olddim / (double)mysize;
    	I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, newscale);
    }
    if (I.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, oldscale))
    {
    	newscale = oldscale * (double)olddim / (double)mysize;
    	I.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, newscale);
    }
    // Don't change sampling in Z because only 2D images implemented anyway...

}

void Preprocessing::rewindow(Image<double> &I, int mysize)
{
	// Check 2D or 3D dimensionality
	if (I().getDim() == 2)
	{
		I().window(FIRST_XMIPP_INDEX(mysize), FIRST_XMIPP_INDEX(mysize),
				   LAST_XMIPP_INDEX(mysize),  LAST_XMIPP_INDEX(mysize));
	}
	else if (I().getDim() == 3)
	{
		I().window(FIRST_XMIPP_INDEX(mysize), FIRST_XMIPP_INDEX(mysize), FIRST_XMIPP_INDEX(mysize),
				   LAST_XMIPP_INDEX(mysize),  LAST_XMIPP_INDEX(mysize),  LAST_XMIPP_INDEX(mysize));
	}

}

FileName Preprocessing::getMicrographNameFromRootName(FileName fn_root)
{
	FileName fn_mic, fn_mic_nos;
	//Search the unique part of the micrograph name (i.e. coord name may have additional text after the micrograph name without extension...
	// e.g. fn_root="mic001_all.pos" may correspond to mic001.mrc
	for (int i = 0; i < fn_root.length(); i++)
	{
		if (do_movie_extract)
		{
			fn_mic = fn_root.substr(0, fn_root.length() - i) + "_" + fn_movie + ".mrcs";
			// For movies also allow name without the s of .mrcs
			fn_mic_nos = fn_root.substr(0, fn_root.length() - i) + "_" + fn_movie + ".mrc";
		}
		else
			fn_mic = fn_root.substr(0, fn_root.length() - i) + ".mrc";
		std::vector<FileName> fn_mics;
		if (fn_mic.globFiles(fn_mics) == 1)
		{
			fn_mic = fn_mics[0];
			break;
		}
		else if (do_movie_extract && fn_mic_nos.globFiles(fn_mics) == 1)
		{
			// For movies also allow name without the s of .mrcs
			fn_mic = fn_mics[0];
			break;
		}
		if (i == fn_root.length() - 1)
		{
			return "";
		}
	}

	return fn_mic;

}

FileName Preprocessing::getRootNameFromMicrographName(FileName fn_mic)
{
	if (do_movie_extract)
	{
		if (fn_mic.contains(".mrcs"))
			return fn_mic.without("_" + fn_movie + ".mrcs");
		else if (fn_mic.contains(".mrc"))
			return fn_mic.without("_" + fn_movie + ".mrc");
		else
			REPORT_ERROR("Preprocessing::getRootNameFromMicrographName ERROR: movie name does not end in \"_\" + movie-identifier + \".mrc\" or \".mrcs\"");
	}
	else
		return fn_mic.without(".mrc");


}
