// Translate this file with
//
// g++ -O3 --std=c++11 assignment-2019.c -o assignment
//
// Run it with
//
// ./assignment
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2019 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>
#include <time.h>

clock_t start, end;
double cpu_time_used;


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Set up scenario from the command line.
 *
 * This operation is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-3) / 7;

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



/**
 * This is the only operation you are allowed to change in the assignment.
 */
void updateBody() {
  maxV = 0.0;
  minDx = std::numeric_limits < double > ::max();

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double * force0 = new double[NumberOfBodies]();
  double * force1 = new double[NumberOfBodies]();
  double * force2 = new double[NumberOfBodies]();

  for (int i = 0; i < NumberOfBodies; i++) {
    for (int j = i + 1; j < NumberOfBodies; j++) {
      const double distance = sqrt(
        (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
        (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
        (x[j][2] - x[i][2]) * (x[j][2] - x[i][2])
      );

      minDx = std::min( minDx,distance );

      const double m2perD3 = mass[j] * mass[i] / distance / distance / distance;
      double forceHolder = (x[j][0] - x[i][0]) * m2perD3;
      force0[i] += forceHolder;
      force0[j] -= forceHolder;

      forceHolder = (x[j][1] - x[i][1]) * m2perD3;
      force1[i] += forceHolder;
      force1[j] -= forceHolder;

      forceHolder = (x[j][2] - x[i][2]) * m2perD3;
      force2[i] += forceHolder;
      force2[j] -= forceHolder;
    }

    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];

    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

    // Check for a collision
    if (minDx <= 1e-2) {
      std::cout << "minDX < 1e-2" << std::endl;
      for (int i=0; i<NumberOfBodies; i++) {
        for (int j=i+1; j<NumberOfBodies; j++) {
          const double distanceBetweenIJ = 
            (x[j][0] - x[i][0]) * (x[j][0] - x[i][0]) +
            (x[j][1] - x[i][1]) * (x[j][1] - x[i][1]) +
            (x[j][2] - x[i][2]) * (x[j][2] - x[i][2]);
          std::cout << distanceBetweenIJ << std::endl;
          if (distanceBetweenIJ <= 1e-4) {
            int min, max;
            if (i < j) {
              min = i;
              max = j;
            } else {
              min = j;
              max = i;
            }
            std::cout << "Collision between min: " << min << " and max : " << max << std::endl;
            NumberOfBodies--;
            std::cout << "There remains: " << NumberOfBodies << std::endl;
            const double oneOverMass = 1 / (mass[i] + mass[j]);
            const double weightedMassMin = mass[min] * oneOverMass;
            const double weightedMassMax = mass[max] * oneOverMass;

            x[min][0] = weightedMassMin * x[min][0] + weightedMassMax * x[max][0];
            x[min][1] = weightedMassMin * x[min][1] + weightedMassMax * x[max][1];
            x[min][2] = weightedMassMin * x[min][2] + weightedMassMax * x[max][2];

            x[max][0] = x[NumberOfBodies][0];
            x[max][1] = x[NumberOfBodies][1];
            x[max][2] = x[NumberOfBodies][2];

            if (NumberOfBodies == 1) {
              std::cout << "Position of final object: " << x[0][0] << ", " <<
               x[0][1] << ", " << x[0][1] << "." << std::endl;
              std::_Exit(0);
            }

            v[min][0] = weightedMassMin * v[min][0] + weightedMassMax * v[max][0];
            v[min][1] = weightedMassMin * v[min][1] + weightedMassMax * v[max][1];
            v[min][2] = weightedMassMin * v[min][2] + weightedMassMax * v[max][2];

            v[max][0] = v[NumberOfBodies][0];
            v[max][1] = v[NumberOfBodies][1];
            v[max][2] = v[NumberOfBodies][2];

            mass[min] = mass[i] + mass[j];
            mass[max] = mass[NumberOfBodies];
          }
        }
      }
    }

    maxV = std::max(maxV, v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]);
  }

  t += timeStepSize;
  maxV = std::sqrt(maxV);

  delete[] force0;
  delete[] force1;
  delete[] force2;
}


/**
 * Main routine.
 *
 * Not to be changed in assignment.
 */
int main(int argc, char** argv) {
  start = clock();
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0   0 0 0 1.0   0   0 1.0 0 1.0 0 1.0 0   0 1.0 \t One spiralling around the other one" << std::endl
              << "0.01  100.0 3.0 0 0   0 1.0   0 0.4 0   0 0   0 0   0 0.2 2.0 0 0 0 0 0 1.0 \t Three body setup from first lecture" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
/*
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;*/  

      tPlot += tPlotDelta;
    }
  }

  closeParaviewVideoFile();
  
  end = clock();

  cpu_time_used = ((double) end - start) / CLOCKS_PER_SEC;

  std::cout << cpu_time_used << std::endl;

  return 0;
}
