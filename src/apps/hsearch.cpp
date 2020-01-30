
#include "src/helix.h"
#include "src/image.h"
#include "src/args.h"

FileName fnin, fnout;

double hrise, hturn, hinner, houter, apix, htake, hturninc, hriseinc;

long int hsuper;

int main(int argc, char ** argv)
{
	IOParser parser;

	try {
	
		parser.setCommandLine(argc, argv);
	
		parser.addSection("Input / output");
		fnin=parser.getOption("--i", "Input volume", "");
//		fnout=parser.getOption("--o", "Output volume", "");
		parser.addSection("Helical parameters");
		hsuper=textToInteger(parser.getOption("--hsuper", "Helicising super-sampling factor in integer", "1"));
		hrise=textToFloat(parser.getOption("--hrise", "Helical rise per subunit in Angstrom", "10.0"));
		hturn=textToFloat(parser.getOption("--hturn", "Helical turn per subunit in Angstrom", "10.0"));
		hinner=textToFloat(parser.getOption("--hinner", "Helical inner radius in Angstrom", "10.0"));
		houter=textToFloat(parser.getOption("--houter", "Helical outer radius in Angstrom", "100.0"));
		htake=textToFloat(parser.getOption("--htake", "Fraction to use from center of input volume ", "0.50"));
		apix=textToFloat(parser.getOption("--angpix", "Pixel size in Angstrom per pixel", "1.00"));
		hriseinc=textToFloat(parser.getOption("--hrisestep", "Step size for search of helical rise per subunit in Angstrom", "0.1"));
		hturninc=textToFloat(parser.getOption("--hrisestep", "Step size for search of helical turn per subunit in Angstrom", "0.1"));


		Image<double> im;

		im.read(fnin);

		Heliciser<double> h(apix, hrise, hturn, hinner, houter, hsuper, htake);
	
		h.hsearch(im.data, hriseinc, hturninc);

	}	

	catch (RelionError XE)
	{
		parser.writeUsage(std::cerr);
		std::cerr << XE;
		exit(1);
	}

	exit(0);

}


