#ifndef INC_AMBER_NETCDF_H
#define INC_AMBER_NETCDF_H
/* 
 * Daniel R. Roe 
 * 2010-12-07
 * A C implementation of routines for reading and writing the Amber Netcdf 
 * trajectory format.
 * Based on Cpptraj implementation.
 * Original implementation of netcdf in Amber by Jon Mongan.
 */
struct AmberNetcdf {
  int ncid;           // Netcdf ID of the file when open
  int frameID;        // ID of frame variable
  int ncframe;        // Number of frames in the file
  int currentFrame;   // Current frame
  int atomID;         // ID of atom variable
  int ncatom;         // # of atoms
  int coordID;        // ID of coords variable
  int cellAngleID;    // ID of box angle variable
  int cellLengthID;   // ID of box length variable
  float temp0;        // Temperature of current frame (if TempVarID!=-1)

  int spatialID;
  int labelID;
  int cell_spatialID;
  int cell_angularID;
  int spatialVID;
  int timeVID;
  int cell_spatialVID;
  int cell_angularVID;
  int TempVarID;
};

#endif
