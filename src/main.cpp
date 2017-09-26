/** Model of a chiral topological superconductor embedded in a surrounding
 *  superconductor with Rashba spin-orbit interaction (or with a p-wave order
 *  parameter component). The spin-polarized LDOS is either calculated along a
 *  1D cut crossing the island if the flag cut1D is true (=1), or otherwise
 *  accross the full 2D surface. Parameters are read in from the file
 *  "Parameters".
 *
 *  TBTK version: v0.9.4.
 *
 *  @author Kristofer Björnson
 */

#include "ChebyshevSolver.h"
#include "CPropertyExtractor.h"
#include "FileParser.h"
#include "FileWriter.h"
#include "Model.h"
#include "ParameterSet.h"
#include "SpinPolarizedLDOS.h"

#include <complex>

#include "Timer.h"

using namespace std;
using namespace TBTK;

const complex<double> i(0, 1);

int main(int argc, char **argv){
	Timer::tick();

	//Read parameters from the file "Parameters".
	ParameterSet *parameterSet = FileParser::readParameterSet(
		"Parameters"
	);

	//Lattice size.
	const int SIZE_X = parameterSet->getInt("SIZE_X");;
	const int SIZE_Y = parameterSet->getInt("SIZE_Y");
	const double RADIUS = parameterSet->getDouble("RADIUS");
	const double BOUNDARY_WIDTH = parameterSet->getDouble(
		"BOUNDARY_WIDTH"
	);

	//Chebyshev parameters.
	const int NUM_COEFFICIENTS = parameterSet->getInt("NUM_COEFFICIENTS");
	const int ENERGY_RESOLUTION = parameterSet->getInt(
		"ENERGY_RESOLUTION"
	);
	const double SCALE_FACTOR = parameterSet->getDouble("SCALE_FACTOR");
	const double LOWER_BOUND = parameterSet->getDouble("LOWER_BOUND");
	const double UPPER_BOUND = parameterSet->getDouble("UPPER_BOUND");

	//Model parameters.
	complex<double> mu = parameterSet->getComplex("mu");
	complex<double> t = parameterSet->getComplex("t");
	complex<double> D_s = parameterSet->getComplex("D_s");
	complex<double> D_t = parameterSet->getComplex("D_t");
	complex<double> alpha = parameterSet->getComplex("alpha");
	complex<double> V_z = parameterSet->getComplex("V_z");

	//Flag indicating whether to restrict the calculation to a 1D cut. If
	//false, the calculation will be performed over the full 2D surface.
	bool cut1D = parameterSet->getBool("cut1D");

	//Setup the strength of the Zeeman term.
	complex<double> **magnetization = new complex<double>*[SIZE_X];
	for(int x = 0; x < SIZE_X; x++){
		magnetization[x] = new complex<double>[SIZE_Y];
		for(int y = 0; y < SIZE_Y; y++){
			int rx = x - SIZE_X/2;
			int ry = y - SIZE_Y/2;
			double r = sqrt(rx*rx + ry*ry);
			magnetization[x][y] = V_z*(
				M_PI/2. - atan((r-RADIUS)/BOUNDARY_WIDTH)
			)/(M_PI);
		}
	}

	//Setup the order parameter.
	complex<double> orderparameterS[SIZE_X][SIZE_Y];
	complex<double> orderparameterP[SIZE_X][SIZE_Y];
	for(int x = 0; x < SIZE_X; x++){
		for(int y = 0; y < SIZE_Y; y++){
			int rx = x - SIZE_X/2;
			int ry = y - SIZE_Y/2;
			double r = sqrt(rx*rx + ry*ry);
			orderparameterS[x][y] = D_s;
			orderparameterP[x][y] = D_t;
		}
	}

	//Create model and set up hopping parameters.
	Model model;
	for(int x = 0; x < SIZE_X; x++){
		for(int y = 0; y < SIZE_Y; y++){
			for(int s = 0; s < 2; s++){
				//Add hopping amplitudes corresponding to the
				//chemical potential
				model << HoppingAmplitude(
					-mu,
					{x, y, s},
					{x, y, s}
				);
				model << HoppingAmplitude(
					mu,
					{x, y, s+2},
					{x, y, s+2}
				);

				//Add hopping amplitudes corresponding to the
				//Zeeman term
				model << HoppingAmplitude(
					2.*magnetization[x][y]*(s-1/2.),
					{x, y, s},
					{x, y, s}
				);
				model << HoppingAmplitude(
					-2.*magnetization[x][y]*(s-1/2.),
					{x, y, s+2},
					{x, y, s+2}
				);

				//Add hopping amplitudes corresponding to t and
				//the Rashba spin-orbit interaction
				if(x+1 < SIZE_X){
					model << HoppingAmplitude(
						-t,
						{(x+1)%SIZE_X, y, s},
						{x, y, s}
					) + HC;
					model << HoppingAmplitude(
						t,
						{(x+1)%SIZE_X, y, s+2},
						{x, y, s+2}
					) + HC;
					model << HoppingAmplitude(
						-alpha*2.*(1/2. - s),
						{(x+1)%SIZE_X, y, (s+1)%2},
						{x, y, s}
					) + HC;
					model << HoppingAmplitude(
						alpha*2.*(1/2. - s),
						{(x+1)%SIZE_X, y, (s+1)%2+2},
						{x, y, s+2}
					) + HC;
					model << HoppingAmplitude(
						orderparameterP[x][y]*2.*(1/2. - s),
						{(x+1)%SIZE_X, y, s+2},
						{x, y, s}
					) + HC;
					model << HoppingAmplitude(
						-orderparameterP[x][y]*2.*(1/2. - s),
						{(x+1)%SIZE_X, y, s},
						{x, y, s+2}
					) + HC;
				}
				if(y+1 < SIZE_Y){
					model << HoppingAmplitude(
						-t,
						{x, (y+1)%SIZE_Y, s},
						{x, y, s}
					) + HC;
					model << HoppingAmplitude(
						t,
						{x, (y+1)%SIZE_Y, s+2},
						{x, y, s+2}
					) + HC;
					model << HoppingAmplitude(
						-i*alpha,
						{x, (y+1)%SIZE_Y, (s+1)%2},
						{x, y, s}
					) + HC;
					model << HoppingAmplitude(
						-i*alpha,
						{x, (y+1)%SIZE_Y, (s+1)%2+2},
						{x, y, s+2}
					) + HC;
					model << HoppingAmplitude(
						i*orderparameterP[x][y],
						{x, (y+1)%SIZE_Y, s+2},
						{x, y, s}
					) + HC;
					model << HoppingAmplitude(
						i*orderparameterP[x][y],
						{x, (y+1)%SIZE_Y, s},
						{x, y, s+2}
					) + HC;
				}

				//Add hopping amplitudes corresponding to
				//superconducting order parameter
				model << HoppingAmplitude(
					2.*orderparameterS[x][y]*(s - 1/2.),
					{x, y, 3-s},
					{x, y, s}
				) + HC;
			}
		}
	}

	//Construct model. (The second call is needed to use the Chebyshev
	//solver).
	model.construct();
	model.constructCOO();

	//Setup the ChebyshevSolver.
	ChebyshevSolver cSolver;
	cSolver.setModel(model);
	cSolver.setScaleFactor(SCALE_FACTOR);

	//Set the filename and remove any file already in the folder.
	FileWriter::setFileName("TBTKResults.h5");
	FileWriter::clear();

	//Create PropertyExtractor. The parameter are in order: The
	//ChebyshevSolver, number of expansion coefficients used in the
	//Cebyshev expansion, energy resolution with which the Green's function
	// is evaluated, whether calculate expansion functions using a GPU or
	//not, whether to evaluate the Green's function using a GPU or not,
	//whether to use a lookup table for the Green's function or not
	//(required if the Green's function is evaluated on a GPU), and the
	//lower and upper bound between which the Green's function is evaluated
	//(has to be inside the interval [-SCALE_FACTOR, SCALE_FACTOR]).
	CPropertyExtractor pe(
		cSolver,
		NUM_COEFFICIENTS,
		true,
		false,
		true
	);
	pe.setEnergyWindow(LOWER_BOUND, UPPER_BOUND, ENERGY_RESOLUTION);

	//Calculate and save the spin-polarized LDOS
	if(cut1D){
		//Calculate the spin-polarized LDOS along the cut
		//x = [0, SIZE_X-1], y = SIZE_Y/2.
		Property::SpinPolarizedLDOS spinPolarizedLDOS
			= pe.calculateSpinPolarizedLDOS(
				{IDX_X,		SIZE_Y/2,	IDX_SPIN},
				{SIZE_X,	1,		2}
			);
		//Write the spin-polarized local density of states to file.
		FileWriter::writeSpinPolarizedLDOS(spinPolarizedLDOS);
	}
	else{
		//Calculate the spin-polarized LDOS over the full surface.
		Property::SpinPolarizedLDOS spinPolarizedLDOS
			= pe.calculateSpinPolarizedLDOS(
				{IDX_X,		IDX_Y,		IDX_SPIN},
				{SIZE_X,	SIZE_Y,		2}
			);
		//Write the spin-polarized local density of states to file.
		FileWriter::writeSpinPolarizedLDOS(spinPolarizedLDOS);
	}

	Timer::tock();

	return 0;
}
