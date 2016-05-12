/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "tools/Communicator.h"
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR TEMPLATE
///*
//This file provides a template for if you want to introduce a new CV.
//
//<!-----You should add a description of your CV here---->
//This code is used to calculate the density of atoms from GROUP within the union of spherical shells
//
//\par Examples
//
//The following plumed.dat file tells plumed to calculate the N/Ntw (defined in GROUP) values in the union spherical shells
//================================
//STATICUNIONSPHSH ...
//GROUPB=223-53763:3
//NVOL=1
//RLOW=-0.5 RHIGH=0.6 
//X=0.2 Y=0.2 Z=0.2
//SIGMA=0.01 CUTOFF=0.02
//LABEL=sphsh
//... STATICUNIONSPHSH
//================================
//
//Here we add in a restraint for the calculated Ntw
//================================
//RESTRAINT ARG=sphsh.Ntw AT=XXX KAPPA=KKK SLOPE=LLL LABEL=restraint
//================================
//
//and finally, print out the quantity we want
//================================
//PRINT STRIDE=50 ARG=sphsh.N,sphsh.Ntw,restraint.bias FILE=SPHSH
//================================
//
//<!---You should put an example of how to use your CV here--->
//
//*/
////+ENDPLUMEDOC

class StaticUnionSphSh: public Colvar {
  bool pbc;
  bool serial;

private:
// private variables only for this routine
  unsigned nvols;
  unsigned natoms;
  vector<double> rlo, rhi; 
  vector<Vector> vol_coor;
  double sig, rcut;
  double sig2, rcut2;
  double normconst, preerf, prelinear;

public:
  StaticUnionSphSh(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(StaticUnionSphSh,"STATICUNIONSPHSH")

void StaticUnionSphSh::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("PAIR",false,"Pair only 1st element of the 1st group with 1st element in the second, etc");
  keys.add("atoms","GROUP","List of atoms which we want to measure the N and Ntw values in side the union of volumes");
  keys.add("compulsory", "NVOL", "1", "Number of union volumes");
  keys.add("numbered", "RLOW", "the inner radius of the spherical shell");
  keys.add("numbered", "RHIGH", "the outter radius of the spherical shell");
  keys.add("numbered", "X", "the x-coordinate of the spherical shell");
  keys.add("numbered", "Y", "the y-coordinate of the spherical shell");
  keys.add("numbered", "Z", "the z-coordinate of the spherical shell");
  keys.add("compulsory", "SIGMA", "0.01", "the width of the coarse-graining Gaussian");
  keys.add("compulsory", "CUTOFF", "0.02", "the cutoff distance of the coarse-graining Gaussian");
  
  componentsAreNotOptional(keys);
  keys.addOutputComponent("N","default","number of target atom inside the volume");
  keys.addOutputComponent("Ntw","default","coarse-grained number of target atom inside the volume");
}

StaticUnionSphSh::StaticUnionSphSh(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{
  // serial for debug only 
  parseFlag("SERIAL",serial);

  // read in the list of atoms 
  vector<AtomNumber> full_lista;
  parseAtomList("GROUP",full_lista);
  natoms = full_lista.size();  

  // pbc condition
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  // read in INDUS parameters
  parse("NVOL",nvols);
  unsigned nsize=0;
  unsigned npos=0;
  rlo.resize( nvols );
  rhi.resize( nvols );
  vol_coor.resize( nvols );
  for(unsigned i=0;i<nvols;++i) 
  {
    rlo[i]=0.0;
    rhi[i]=0.0;
    for(unsigned j=0;j<3;++j)
      vol_coor[i][j]=0.0;
  }
  // read in volume sizes and positions
  for(unsigned i=0;i<nvols;++i)
  {
      if( !parseNumbered( "RLOW", i+1, rlo[i] ) ) break;
      parseNumbered( "RHIGH", i+1, rhi[i] );
      nsize ++;
  }
  for(unsigned i=0;i<nvols;++i)
  {
      if( !parseNumbered( "X", i+1, vol_coor[i][0] ) ) break;
      parseNumbered( "Y", i+1, vol_coor[i][1] );
      parseNumbered( "Z", i+1, vol_coor[i][2] );
      npos++;
  }
  // if only one size or only one position is defined, set all volumes to that size or position
  if( nsize==0 )
  {
     parse("RLOW",rlo[0]);
     parse("RHIGH",rhi[0]);
     for(unsigned i=1;i<nvols;++i)
     {
       rlo[i]=rlo[0];
       rhi[i]=rhi[0];
     }
     nsize = nvols;
  }
  if ( npos==0 )
  {
     parse("X",vol_coor[0][0]);
     parse("Y",vol_coor[0][1]);
     parse("Z",vol_coor[0][2]);
     for (unsigned i=1;i<nvols;++i)
        vol_coor[i]=vol_coor[0];
     npos = nvols;
  }

  if ( nsize != nvols ) error("missing size parameters for union volumes!");
  if ( npos != nvols ) error("missing positions for union volumes!");

  // read in smearing length and cutoff
  parse("SIGMA",sig);
  parse("CUTOFF",rcut);

  // Precalculate a few constants
  rcut2 = rcut*rcut;
  sig2 = 2.0*sig*sig;
  // Normalisation for a gaussian truncated at x = rcut,
  // shifted and scaled to be continuous, has units of length
  normconst = sqrt( M_PI * sig2 ) * erf( rcut / (sqrt(2.0)*sig) )
                - 2*rcut*exp( - rcut2 / sig2 );
  preerf = sqrt( 0.5 * M_PI * sig * sig ) / normconst; // dimensionless
  prelinear = exp( - rcut2 / sig2 ) / normconst; // units of 1/L

  // the two components need to be the same type: "ComponentWithDerivatives", otherwise "resize derivatives" error will occur  
  addComponentWithDerivatives("N"); componentIsNotPeriodic("N");
  addComponentWithDerivatives("Ntw"); componentIsNotPeriodic("Ntw");

  // return fullatomlist
  requestAtoms(full_lista);

  // print out the informations read from the input file
  log.printf("  between %u static volumes and %u atoms\n",static_cast<unsigned>(nvols),static_cast<unsigned>(full_lista.size()));
  log.printf("  group of static volumes):\n");
  for(unsigned int i=0;i<nvols;++i){
   log.printf("  %f %f %f %f %f\n", rhi[i], rlo[i], vol_coor[i][0], vol_coor[i][1], vol_coor[i][2]);
  }
  log.printf("  \n group of atoms:\n");
  for(unsigned int i=0;i<full_lista.size();++i){
   if ( (i+1) % 25 == 0 ) log.printf("  \n");
   log.printf("  %d", full_lista[i].serial());
  }
  log.printf("  \n");
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");
  //if(dopair) log.printf("  with PAIR option\n");
}

// calculator
void StaticUnionSphSh::calculate()
{
  double totN=0.0;
  double totNtw=0.0;

  vector<Vector> derivatives(natoms);
  Tensor virial;

  // Paralellize the computation to each processor
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }else{
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }

  // each processor selects equal amount of atoms and calculates the contribution to N or/and Ntw
  // note that each atom needs to run throught all the volumes to decide the union contribtuion
  // WARNING: this paralellization step is tested to be quite slow in the case of large number of union volumes, need to improve
  // SUGGESTION: find a way to utilize the neighbourlist effectively
  // and provide tables for the erf functions
  for(unsigned i=rank;i<natoms;i+=stride) 
  {
    // Initialize per-atom quantities
    bool bInV, bwellInV, bwellOutV;
    bInV = false;
    bwellInV = false;
    bwellOutV = true;

    double NTwiddle_i=1.0;
    //Vector gradNTwiddle;
    vector<double> hTwiddle(nvols);
    vector<Vector> gradhTwiddle(nvols);

    // loop over over all volumes
    for (unsigned j=0;j<nvols;j++)
    {
      // Calculating rsphi
      double rsphi;
      Vector dr;

      if(pbc){
        dr = pbcDistance(vol_coor[j],getPosition(i)); // r_atm(i) - r_vol(j)
      } else {
        dr = delta(vol_coor[j],getPosition(i));
      }
      rsphi = dr.modulo();
      
      double thisNTwiddle = 0.0;
      double thisGradNTwiddle = 0.0; 

      double lo = rlo[j];
      double hi = rhi[j];
      double loMinus = lo - rcut;
      double loPlus = lo + rcut;
      double hiMinus = hi - rcut;
      double hiPlus = hi + rcut;

      // Could atom i possibly contribute to the NTwiddle of this volume j?
      if( (rsphi >= loMinus) && (rsphi <= hiPlus) )
      {
        bwellOutV = false;
        // Could it contribute to N?
        if( (rsphi >= lo) && (rsphi <= hi) && (!bInV) )
        {
          totN += 1.0;
          bInV = true;
        }

        // Calculate contribution to NTilde and dNTilde/dr_i
        thisNTwiddle += 1.0;
        if( rsphi < loPlus )  // Boundary region
        {
          thisNTwiddle -= 0.5 + preerf * erf( sqrt(0.5) * (lo - rsphi)/sig )
                              - prelinear * (lo - rsphi);
          thisGradNTwiddle -= ( exp( -0.5*((rsphi - lo)/sig)*((rsphi - lo)/sig) )
                              - exp( -rcut2/sig2 ) ) / normconst;
        }
        if( rsphi > hiMinus )  // Boundary region
        {
          thisNTwiddle -= 0.5 + preerf * erf( sqrt(0.5) * (rsphi - hi)/sig )
                              - prelinear * (rsphi - hi);
          thisGradNTwiddle += ( exp( -0.5*((rsphi - hi)/sig)*((rsphi - hi)/sig) )
                              - exp( -rcut2 / sig2 ) ) / normconst;
        }
   
        // Loading hTwiddle and gradhTwiddle
        hTwiddle[j] = thisNTwiddle;
        gradhTwiddle[j] = -1.0*thisGradNTwiddle*dr/rsphi;

        // Accumulate NTwiddle of intersections of complements
        NTwiddle_i *= ( 1.0 - thisNTwiddle );

        // Is it well within V?
        if( (rsphi >= loPlus) && (rsphi <= hiMinus) )
        {
          bwellInV = true;
          NTwiddle_i = 0.0;
          break; // break out of loop over all pairs
        }
      } // if atom i contributes 
      else continue;
    } // loop over all pairs

    // Accumulate local NTwiddle 
    totNtw += ( 1.0 - NTwiddle_i ); 

    // Assign gradNTwiddle
    if( (! bwellInV) && (! bwellOutV) )
    {
      for( unsigned k = 0; k < nvols; k++)
      {
        Vector cellgradh;
        cellgradh = gradhTwiddle[k];

        for( unsigned m = 0; m < k; m++)
          cellgradh *= (1.0 - hTwiddle[m]);
        for( unsigned m = k+1; m < nvols; m++)
          cellgradh *= (1.0 - hTwiddle[m]);
 
        //the force exerted by cell k on atom i
        derivatives[i] += cellgradh;
      }
    } // if gradNTwiddle exists (boundary region)

  } // loop over all atoms in group B at local processor
 
  // Sum over all the processors
  if(!serial){
    comm.Sum(totN);
    comm.Sum(totNtw);
    comm.Sum(derivatives);
  }

  // Distribute computed values to the CV components
  Value* valueN=getPntrToComponent("N");
  Value* valueNtw=getPntrToComponent("Ntw");
  valueN->set(totN);
  valueNtw->set(totNtw); 

  // update the forces on each atoms and the virial contributions on the simulation box
  for (int i=0; i<natoms; i++)
  {
    setAtomsDerivatives(valueNtw,i,derivatives[i]);
    virial += -0.5*Tensor(getPosition(i), derivatives[i]);
  }
  setBoxDerivatives(virial);

} // end of calculate()
}
}
