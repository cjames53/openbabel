/**********************************************************************
forcefielduff.cpp - UFF force field.
 
Copyright (C) 2007-2008 by Geoffrey Hutchison
Some portions Copyright (C) 2006-2008 by Tim Vandermeersch

This file is part of the Open Babel project.
For more information, see <http://openbabel.sourceforge.net/>
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/mol.h>
#include <openbabel/locale.h>

#include "forcefielduff.h"


using namespace std;

// This implementation was created based on open code and reference websites:
// http://towhee.sourceforge.net/forcefields/uff.html
// http://rdkit.org/
// http://franklin.chm.colostate.edu/mmac/uff.html
// (for the last, use the Wayback Machine: http://www.archive.org/

// As well, the main UFF paper:
// Rappe, A. K., et. al.; J. Am. Chem. Soc. (1992) 114(25) p. 10024-10035.

namespace OpenBabel
{
  template<bool gradients>
  void OBFFBondCalculationUFF::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }
 
    vector3 vab, da, db;
    double delta2, dE;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    delta = rab - r0; // we pre-compute the r0 below
    delta2 = delta * delta;
    energy = kb * delta2; // we fold the 1/2 into kb below

    if (gradients) {
      dE = 2.0 * kb * delta;
      da *= dE;
      db *= dE;
      da.Get(force_a);
      db.Get(force_b);
    }
  }
  
  template<bool gradients>
  double OBForceFieldUFF::E_Bond()
  {
    vector<OBFFBondCalculationUFF>::iterator i;
    double energy = 0.0;
        
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nB O N D   S T R E T C H I N G\n\n");
      OBFFLog("ATOM TYPES  BOND    BOND       IDEAL       FORCE\n");
      OBFFLog(" I      J   TYPE   LENGTH     LENGTH     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("------------------------------------------------------------------------\n");
    }
 
    for (i = _bondcalculations.begin(); i != _bondcalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;
      
      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
      }
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%-5s %-5s  %4.2f%8.3f   %8.3f     %8.3f   %8.3f   %8.3f\n",
                (*i).a->GetType(), (*i).b->GetType(), 
                (*i).bt, (*i).rab, (*i).r0, (*i).kb, (*i).delta, (*i).energy);
        OBFFLog(_logbuf);
      }
    }
    
    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL BOND STRETCHING ENERGY = %8.3f %s\n",  energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }
    return energy;
  }
  
  template<bool gradients>
  void OBFFAngleCalculationUFF::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c)) {
      energy = 0.0;
      return;
    }
 
    vector3 da, db, dc;
		double dE;
    
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      theta = OBForceField::VectorAngleDerivative(da, db, dc) * DEG_TO_RAD;
    } else {
      theta = a->GetAngle(b, c) * DEG_TO_RAD;
		}

    if (!isfinite(theta))
      theta = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us
    
    double cosT = cos(theta);
    switch (coord) {
    case 1: // sp -- linear case, minima at 180 degrees, max (amplitude 2*ka) at 0, 360
      // Fixed typo from Rappe paper (i.e., it's NOT 1 - cosT)
      energy = ka*(1.0 + cosT);
      break;
    case 2: // sp2 -- trigonal planar, min at 120, 240, max at 0, 360 (amplitude 2*ka)
      // Rappe form: (1 - cos 3*theta) -- minima at 0, 360 (bad...)
      energy = (ka/4.5) * (1.0 + (1.0 + cosT)*(4.0*cosT));
      break;
    case 4: // square planar // min at 90, 180, 270, max at 0, 360 (amplitude 2*ka)
    case 6: // octahedral
      // Rappe form: (1 - cos 4*theta) -- minima at 0, 360 (bad...)
      energy = ka * (1.0 + cosT)*cosT*cosT;
      break;
    default: // general (sp3) coordination
      energy = ka*(c0 + c1*cosT + c2*(2.0*cosT*cosT - 1.0)); // use cos 2t = (2cos^2 - 1)
    }  
        
    if (gradients) {
      double sinT = sin(theta);
      
      switch (coord) {
      case 1: // sp -- linear case
        dE = -ka * sinT;
        break;
      case 2: // sp2 -- trigonal planar
        dE = -(ka*4.0/4.5) * (sinT + sin(2.0*theta));
        break;
      case 4: // square planar
      case 6: // octahedral
        dE = -ka * cosT * (2.0 + 3.0 * cosT) * sinT;
        break;
      default: // general (sp3) coordination
        dE = -ka * (c1*sinT + 2.0 * c2*sin(2.0 * theta));
      }
      
      
      da *= dE; // da = dTheta/dx * dE/dTheta
      db *= dE; // da = dTheta/dx * dE/dTheta
      dc *= dE; // da = dTheta/dx * dE/dTheta
      da.Get(force_a);
      db.Get(force_b);
      dc.Get(force_c);
    }
  }
  
  template<bool gradients>
  double OBForceFieldUFF::E_Angle()
  {
    vector<OBFFAngleCalculationUFF>::iterator i;
    double energy = 0.0;

    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nA N G L E   B E N D I N G\n\n");
      OBFFLog("ATOM TYPES       VALENCE     IDEAL      FORCE\n");
      OBFFLog(" I    J    K      ANGLE      ANGLE     CONSTANT      DELTA      ENERGY\n");
      OBFFLog("-----------------------------------------------------------------------------\n");
    }
       
    for (i = _anglecalculations.begin(); i != _anglecalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;
      
      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
        AddGradient((*i).force_c, (*i).idx_c);
      }

      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%-5s %-5s %-5s%8.3f  %8.3f     %8.3f   %8.3f   %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).theta * RAD_TO_DEG, (*i).theta0, (*i).ka, (*i).delta, (*i).energy);
        OBFFLog(_logbuf);
      }
    }
 
    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ANGLE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }
    return energy;
  }
  
  template<bool gradients>
  void OBFFTorsionCalculationUFF::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c, idx_d)) {
      energy = 0.0;
      return;
    }
 
    vector3 da, db, dc, dd;
    double cosine;
    double dE;
    
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      dd = d->GetVector();
      tor = OBForceField::VectorTorsionDerivative(da, db, dc, dd);
      if (!isfinite(tor))
        tor = 1.0e-3;
      tor *= DEG_TO_RAD;
    } else {
      vector3 vab, vbc, vcd, abbc, bccd;
      vab = a->GetVector() - b->GetVector();
      vbc = b->GetVector() - c->GetVector();
      vcd = c->GetVector() - d->GetVector();
      abbc = cross(vab, vbc);
      bccd = cross(vbc, vcd);

      double dotAbbcBccd = dot(abbc,bccd);
      tor = acos(dotAbbcBccd / (abbc.length() * bccd.length()));
      if (IsNearZero(dotAbbcBccd) || !isfinite(tor)) { // stop any NaN or infinity
        tor = 1.0e-3; // rather than NaN
      }
      else if (dotAbbcBccd > 0.0) {
        tor = -tor;
      }
    }

    cosine = cos(tor * n);
    energy = V * (1.0 - cosNPhi0*cosine);
    
    if (gradients) {
      dE = -(V * n * cosNPhi0 * sin(n * tor));
      da *= dE;
      db *= dE;
      dc *= dE;
      dd *= dE;
      da.Get(force_a);
      db.Get(force_b);
      dc.Get(force_c);
      dd.Get(force_d);
    }
  }
  
  template<bool gradients>
  double OBForceFieldUFF::E_Torsion() 
  {
    vector<OBFFTorsionCalculationUFF>::iterator i;
    double energy = 0.0;
 
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nT O R S I O N A L\n\n");
      OBFFLog("----ATOM TYPES-----    FORCE         TORSION\n");
      OBFFLog(" I    J    K    L     CONSTANT        ANGLE         ENERGY\n");
      OBFFLog("----------------------------------------------------------------\n");
    }
    
    for (i = _torsioncalculations.begin(); i != _torsioncalculations.end(); ++i) {

      i->template Compute<gradients>();
      energy += i->energy;
      
      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
        AddGradient((*i).force_c, (*i).idx_c);
        AddGradient((*i).force_d, (*i).idx_d);
      }
 
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%-5s %-5s %-5s %-5s%6.3f       %8.3f     %8.3f\n",
                (*i).a->GetType(), (*i).b->GetType(), 
                (*i).c->GetType(), (*i).d->GetType(), (*i).V, 
                (*i).tor * RAD_TO_DEG, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL TORSIONAL ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  /*
  //  a
  //   \
  //    b---d      plane = a-b-c
  //   /  
  //  c
  */
  template<bool gradients>
  void OBFFOOPCalculationUFF::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b, idx_c, idx_d)) {
      energy = 0.0;
      return;
    }
 
    vector3 da, db, dc, dd;
    double dE;
    
    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      dc = c->GetVector();
      dd = d->GetVector();
      angle = OBForceField::VectorOOPDerivative(da, db, dc, dd) * DEG_TO_RAD;  

	    if (!isfinite(angle))
	      angle = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;

      // somehow we already get the -1 from the OOPDeriv -- so we'll omit it here
      dE = koop * (c1*sin(angle) + 2.0 * c2 * sin(2.0*angle));
      da *= dE;
      db *= dE;
      dc *= dE;
      dd *= dE;
      da.Get(force_a);
      db.Get(force_b);
      dc.Get(force_c);
      dd.Get(force_d);
    } else {
      angle = DEG_TO_RAD*Point2PlaneAngle(d->GetVector(), a->GetVector(), b->GetVector(), c->GetVector());
      if (!isfinite(angle))
        angle = 0.0; // doesn't explain why GetAngle is returning NaN but solves it for us;
    }
      
    energy = koop * (c0 + c1 * cos(angle) + c2 * cos(2.0*angle));
  }

  template<bool gradients>
  double OBForceFieldUFF::E_OOP() 
  {
     vector<OBFFOOPCalculationUFF>::iterator i;
     double energy = 0.0;
    
     IF_OBFF_LOGLVL_HIGH {
       OBFFLog("\nO U T - O F - P L A N E   B E N D I N G\n\n");
       OBFFLog("ATOM TYPES                 OOP     FORCE \n");
       OBFFLog(" I    J     K     L       ANGLE   CONSTANT     ENERGY\n");
       OBFFLog("----------------------------------------------------------\n");
     }

     for (i = _oopcalculations.begin(); i != _oopcalculations.end(); ++i) {
      i->template Compute<gradients>();
      energy += i->energy;
      
      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
        AddGradient((*i).force_c, (*i).idx_c);
        AddGradient((*i).force_d, (*i).idx_d);
      }
 
      IF_OBFF_LOGLVL_HIGH {
         snprintf(_logbuf, BUFF_SIZE, "%-5s %-5s %-5s %-5s%8.3f   %8.3f     %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), (*i).c->GetType(), (*i).d->GetType(), 
                 (*i).angle * RAD_TO_DEG, (*i).koop, (*i).energy);
         OBFFLog(_logbuf);
       }
     }

     IF_OBFF_LOGLVL_HIGH {
       snprintf(_logbuf, BUFF_SIZE, "     TOTAL OUT-OF-PLANE BENDING ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
       OBFFLog(_logbuf);
     }
     return energy;
  }

  template<bool gradients>
  void OBFFVDWCalculationUFF::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }
 
    vector3 da, db;
    double term6, term12, dE, term7, term13;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);
    
    if (IsNearZero(rab, 1.0e-3))
      rab = 1.0e-3;
    
    term7 = term13 = term6 = ka / rab;

    term6 = term6 * term6 * term6; // ^3
    term6 = term6 * term6; // ^6
    term12 = term6 * term6; // ^12
   
    energy = kab * ((term12) - (2.0 * term6));
    
    if (gradients) { 
      term13 = term13 * term12; // ^13
      term7 = term7 * term6; // ^7
      dE = kab * 12.0 * (term7/ka - term13/ka);
      da *= dE;
      db *= dE;
      da.Get(force_a);
      db.Get(force_b);
    }
  }
  
  template<bool gradients>
  double OBForceFieldUFF::E_VDW()
  {
    vector<OBFFVDWCalculationUFF>::iterator i;
    double energy = 0.0;
     
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nV A N   D E R   W A A L S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J        Rij       kij       ENERGY\n");
      OBFFLog("-----------------------------------------\n");
      //          XX   XX     -000.000  -000.000  -000.000  -000.000
    }
    
    unsigned int j = 0; 
    for (i = _vdwcalculations.begin(); i != _vdwcalculations.end(); ++i, ++j) {
      // Cut-off check
      if (_cutoff)
        if (!_vdwpairs.BitIsSet(j)) 
          continue;
     
      i->template Compute<gradients>();
      energy += i->energy;
      
      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
      }
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%-5s %-5s %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).kab, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL VAN DER WAALS ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  template<bool gradients>
  void OBFFElectrostaticCalculationUFF::Compute()
  {
    if (OBForceField::IgnoreCalculation(idx_a, idx_b)) {
      energy = 0.0;
      return;
    }
 
    vector3 da, db;
    double dE, rab2;

    if (gradients) {
      da = a->GetVector();
      db = b->GetVector();
      rab = OBForceField::VectorLengthDerivative(da, db);
    } else
      rab = a->GetDistance(b);

    if (IsNearZero(rab, 1.0e-3))
      rab = 1.0e-3;

    energy = qq / rab;

    if (gradients) {
      rab2 = rab * rab;
      dE = -qq / rab2;
      da *= dE;
      db *= dE;
      da.Get(force_a);
      db.Get(force_b);
    } 
  }
  
  template<bool gradients>
  double OBForceFieldUFF::E_Electrostatic()
  {
    vector<OBFFElectrostaticCalculationUFF>::iterator i;
    double energy = 0.0;
     
    IF_OBFF_LOGLVL_HIGH {
      OBFFLog("\nE L E C T R O S T A T I C   I N T E R A C T I O N S\n\n");
      OBFFLog("ATOM TYPES\n");
      OBFFLog(" I    J           Rij   332.17*QiQj  ENERGY\n");
      OBFFLog("-------------------------------------------\n");
      //            XX   XX     -000.000  -000.000  -000.000  
    }

    unsigned int j = 0; 
    for (i = _electrostaticcalculations.begin(); i != _electrostaticcalculations.end(); ++i, ++j) {
      // Cut-off check
      if (_cutoff)
        if (!_elepairs.BitIsSet(j)) 
          continue;
      
      i->template Compute<gradients>();
      energy += i->energy;
      
      if (gradients) {
        AddGradient((*i).force_a, (*i).idx_a);
        AddGradient((*i).force_b, (*i).idx_b);
      }
      
      IF_OBFF_LOGLVL_HIGH {
        snprintf(_logbuf, BUFF_SIZE, "%-5s %-5s   %8.3f  %8.3f  %8.3f\n", (*i).a->GetType(), (*i).b->GetType(), 
                (*i).rab, (*i).qq, (*i).energy);
        OBFFLog(_logbuf);
      }
    }

    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "     TOTAL ELECTROSTATIC ENERGY = %8.3f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }

  //***********************************************
  //Make a global instance
  OBForceFieldUFF theForceFieldUFF("UFF", true);
  //***********************************************

  OBForceFieldUFF::~OBForceFieldUFF()
  {
  }

  OBForceFieldUFF &OBForceFieldUFF::operator=(OBForceFieldUFF &src)
  {
    _mol = src._mol;

    _ffparams    = src._ffparams;

    _bondcalculations          = src._bondcalculations;
    _anglecalculations         = src._anglecalculations;
    _torsioncalculations       = src._torsioncalculations;
    _oopcalculations           = src._oopcalculations;
    _vdwcalculations           = src._vdwcalculations;
    _electrostaticcalculations = src._electrostaticcalculations;
    _init                      = src._init;

    return *this;
  }

  double CalculateBondDistance(OBFFParameter *i, OBFFParameter *j, double bondorder)
  {
    double ri, rj;
    double chiI, chiJ;
    double rbo, ren;
    ri = i->_dpar[0];
    rj = j->_dpar[0];
    chiI = i->_dpar[8];
    chiJ = j->_dpar[8];

    // precompute the equilibrium geometry
    // From equation 3
    rbo = -0.1332*(ri+rj)*log(bondorder);
    // From equation 4
    ren = ri*rj*(pow((sqrt(chiI) - sqrt(chiJ)),2.0)) / (chiI*ri + chiJ*rj);
    // From equation 2
    // NOTE: See http://towhee.sourceforge.net/forcefields/uff.html
    // There is a typo in the published paper
    return(ri + rj + rbo - ren);
  }

  bool OBForceFieldUFF::SetupCalculations()
  {
    OBFFParameter *parameterA, *parameterB, *parameterC;
    OBAtom *a, *b, *c, *d;
    
    IF_OBFF_LOGLVL_LOW
      OBFFLog("\nS E T T I N G   U P   C A L C U L A T I O N S\n\n");
 
    // 
    // Bond Calculations
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP BOND CALCULATIONS...\n");
    
    OBFFBondCalculationUFF bondcalc;
    double bondorder;

    _bondcalculations.clear();
    
    FOR_BONDS_OF_MOL(bond, _mol) {
      a = bond->GetBeginAtom();
      b = bond->GetEndAtom();

      // skip this bond if the atoms are ignored 
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
        continue;

      // if there are any groups specified, check if the two bond atoms are in a single intraGroup
      if (HasGroups()) {
        bool validBond = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()))
            validBond = true;
        }
        if (!validBond)
          continue;
      }
 
      bondorder = bond->GetBondOrder(); 
      if (bond->IsAromatic())
        bondorder = 1.5;
      if (bond->IsAmide())
        bondorder = 1.41;
      
      bondcalc.a = a;
      bondcalc.b = b;
      bondcalc.bt = bondorder;

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);
			
      if (parameterA == NULL || parameterB == NULL) { 
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR BOND %d-%d (IDX)...\n", 
              a->GetIdx(), b->GetIdx());
          OBFFLog(_logbuf);
        }
        return false;
      }
 
      bondcalc.r0 = CalculateBondDistance(parameterA, parameterB, bondorder);

      // here we fold the 1/2 into the kij from equation 1a
      // Otherwise, this is equation 6 from the UFF paper.
      bondcalc.kb = (0.5 * KCAL_TO_KJ * 664.12 
        * parameterA->_dpar[5] * parameterB->_dpar[5])
        / (bondcalc.r0 * bondcalc.r0 * bondcalc.r0);

      bondcalc.SetupPointers();
      _bondcalculations.push_back(bondcalc);
    }

    //
    // Angle Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ANGLE CALCULATIONS...\n");
 
    OBFFAngleCalculationUFF anglecalc;
 
    _anglecalculations.clear();
    
    double sinT0;
		double rab, rbc, rac;
		OBBond *bondPtr;
    FOR_ANGLES_OF_MOL(angle, _mol) {
      b = _mol.GetAtom((*angle)[0] + 1);
      a = _mol.GetAtom((*angle)[1] + 1);
      c = _mol.GetAtom((*angle)[2] + 1);
      
      // skip this angle if the atoms are ignored 
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) || _constraints.IsIgnored(c->GetIdx()) ) 
        continue;
 
      // if there are any groups specified, check if the three angle atoms are in a single intraGroup
      if (HasGroups()) {
        bool validAngle = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()) && 
              _intraGroup[i].BitIsOn(c->GetIdx()))
            validAngle = true;
        }
        if (!validAngle)
          continue;
      }
 
      anglecalc.a = a;
      anglecalc.b = b;
      anglecalc.c = c;

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);
      parameterC = GetParameterUFF(c->GetType(), _ffparams);

      if (parameterA == NULL || parameterB == NULL || parameterC == NULL) { 
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR ANGLE %d-%d-%d (IDX)...\n", 
              a->GetIdx(), b->GetIdx(), c->GetIdx());
          OBFFLog(_logbuf);
        }
        return false;
      }
 

      anglecalc.coord = parameterB->_ipar[0]; // coordination of central atom

      anglecalc.zi = parameterA->_dpar[5];
      anglecalc.zk = parameterC->_dpar[5];
      anglecalc.theta0 = parameterB->_dpar[1];
      anglecalc.cosT0 = cos(anglecalc.theta0 * DEG_TO_RAD);
      sinT0 = sin(anglecalc.theta0 * DEG_TO_RAD);
      anglecalc.c2 = 1.0 / (4.0 * sinT0 * sinT0);
      anglecalc.c1 = -4.0 * anglecalc.c2 * anglecalc.cosT0;
      anglecalc.c0 = anglecalc.c2*(2.0*anglecalc.cosT0*anglecalc.cosT0 + 1.0);

			// Precompute the force constant
			bondPtr = _mol.GetBond(a,b);
			bondorder = bondPtr->GetBondOrder(); 
      if (bondPtr->IsAromatic())
        bondorder = 1.5;
      if (bondPtr->IsAmide())
        bondorder = 1.41;      
			rab = CalculateBondDistance(parameterA, parameterB, bondorder);

			bondPtr = _mol.GetBond(b,c);
			bondorder = bondPtr->GetBondOrder(); 
      if (bondPtr->IsAromatic())
        bondorder = 1.5;
      if (bondPtr->IsAmide())
        bondorder = 1.41;
			rbc = CalculateBondDistance(parameterB, parameterC, bondorder);
			rac = sqrt(rab*rab + rbc*rbc - 2.0 * rab*rbc*anglecalc.cosT0);

			// Equation 13 from paper -- corrected by Towhee
			// Note that 1/(rij * rjk) cancels with rij*rjk in eqn. 13
			anglecalc.ka = (644.12 * KCAL_TO_KJ) * (anglecalc.zi * anglecalc.zk / (pow(rac, 5.0)));
			anglecalc.ka *= (3.0*rab*rbc*(1.0 - anglecalc.cosT0*anglecalc.cosT0) - rac*rac*anglecalc.cosT0);
      
      anglecalc.SetupPointers();
      _anglecalculations.push_back(anglecalc);
    }
    
    //
    // Torsion Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP TORSION CALCULATIONS...\n");
 
    OBFFTorsionCalculationUFF torsioncalc;
    double torsiontype;
    double phi0 = 0.0;

    _torsioncalculations.clear();
 
    double vi, vj;
    FOR_TORSIONS_OF_MOL(t, _mol) {
      a = _mol.GetAtom((*t)[0] + 1);
      b = _mol.GetAtom((*t)[1] + 1);
      c = _mol.GetAtom((*t)[2] + 1);
      d = _mol.GetAtom((*t)[3] + 1);

      // skip this torsion if the atoms are ignored 
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
           _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) ) 
        continue;
 
      // if there are any groups specified, check if the four torsion atoms are in a single intraGroup
      if (HasGroups()) {
        bool validTorsion = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()) && 
              _intraGroup[i].BitIsOn(c->GetIdx()) && _intraGroup[i].BitIsOn(d->GetIdx()))
            validTorsion = true;
        }
        if (!validTorsion)
          continue;
      }
 
      OBBond *bc = _mol.GetBond(b, c);
      torsiontype = bc->GetBondOrder(); 
      if (bc->IsAromatic())
        torsiontype = 1.5;
      if (bc->IsAmide())
        torsiontype = 1.41;
      
      torsioncalc.a = a;
      torsioncalc.b = b;
      torsioncalc.c = c;
      torsioncalc.d = d;
      torsioncalc.tt = torsiontype;

      parameterB = GetParameterUFF(b->GetType(), _ffparams);
      parameterC = GetParameterUFF(c->GetType(), _ffparams);

      if (parameterB == NULL || parameterC == NULL) { 
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR TORSION X-%d-%d-X (IDX)...\n", 
              b->GetIdx(), c->GetIdx());
          OBFFLog(_logbuf);
        }
        return false;
      }
 
      if (parameterB->_ipar[0] == 3 && parameterC->_ipar[0] == 3) {
        // two sp3 centers
        phi0 = 60.0;
        torsioncalc.n = 3;
        vi = parameterB->_dpar[6];
        vj = parameterC->_dpar[6];

        // exception for a pair of group 6 sp3 atoms
        switch (b->GetAtomicNum()) {
        case 8:
          vi = 2.0;
          torsioncalc.n = 2;
          phi0 = 90.0;
          break;
        case 16:
        case 34:
        case 52:
        case 84:
          vi = 6.8;
          torsioncalc.n = 2;
          phi0 = 90.0;
        }
        switch (c->GetAtomicNum()) {
        case 8:
          vj = 2.0;
          torsioncalc.n = 2;
          phi0 = 90.0;
          break;
        case 16:
        case 34:
        case 52:
        case 84:
          vj = 6.8;
          torsioncalc.n = 2;
          phi0 = 90.0;
        }

        torsioncalc.V = 0.5 * KCAL_TO_KJ * sqrt(vi * vj);

      } else if (parameterB->_ipar[0] == 2 && parameterC->_ipar[0] == 2) {
        // two sp2 centers
        phi0 = 180.0;
        torsioncalc.n = 2;
        torsioncalc.V = 0.5 * KCAL_TO_KJ * 5.0 *
          sqrt(parameterB->_dpar[7]*parameterC->_dpar[7]) * 
          (1.0 + 4.18 * log(torsiontype));
      } else if ((parameterB->_ipar[0] == 2 && parameterC->_ipar[0] == 3)
                 || (parameterB->_ipar[0] == 3 && parameterC->_ipar[0] == 2)) {
        // one sp3, one sp2
        phi0 = 0.0;
        torsioncalc.n = 6;
        torsioncalc.V = 0.5 * KCAL_TO_KJ * 1.0;

        // exception for group 6 sp3
        if (parameterC->_ipar[0] == 3) {
          switch (c->GetAtomicNum()) {
          case 8:
          case 16:
          case 34:
          case 52:
          case 84:
            torsioncalc.n = 2;
            phi0 = 90.0;
          }
        }
        if (parameterB->_ipar[0] == 3) {
          switch (b->GetAtomicNum()) {
          case 8:
          case 16:
          case 34:
          case 52:
          case 84:
            torsioncalc.n = 2;
            phi0 = 90.0;
          }
        }
      }

      if (IsNearZero(torsioncalc.V)) // don't bother calcuating this torsion
        continue;

      // still need to implement special case of sp2-sp3 with sp2-sp2

      torsioncalc.cosNPhi0 = cos(torsioncalc.n * DEG_TO_RAD * phi0);
      torsioncalc.SetupPointers();
      _torsioncalculations.push_back(torsioncalc);     
    }
    
    //
    // OOP/Inversion Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP OOP CALCULATIONS...\n");
    OBFFOOPCalculationUFF oopcalc;

    _oopcalculations.clear();
 
    double phi;
    // The original Rappe paper in JACS isn't very clear about the parameters
    // The following was adapted from Towhee
    FOR_ATOMS_OF_MOL(atom, _mol) {
      b = (OBAtom*) &*atom;

      switch (b->GetAtomicNum()) {
      case 6: // carbon
      case 7: // nitrogen
      case 8: // oxygen
      case 15: // phos.
      case 33: // as
      case 51: // sb
      case 83: // bi
        break;
      default: // no inversion term for this element
        continue;
      }

      a = NULL;
      c = NULL;
      d = NULL;
      
      if (EQn(b->GetType(), "N_3", 3) ||
          EQn(b->GetType(), "N_2", 3) ||
          EQn(b->GetType(), "N_R", 3) ||
          EQn(b->GetType(), "O_2", 3) ||
          EQn(b->GetType(), "O_R", 3)) {
        oopcalc.c0 = 1.0;
        oopcalc.c1 = -1.0;
        oopcalc.c2 = 0.0;
        oopcalc.koop = 6.0 * KCAL_TO_KJ;
      }
      else if (EQn(b->GetType(), "P_3+3", 5) ||
               EQn(b->GetType(), "As3+3", 5) ||
               EQn(b->GetType(), "Sb3+3", 5) ||
               EQn(b->GetType(), "Bi3+3", 5)) {

        if (EQn(b->GetType(), "P_3+3", 5))
          phi = 84.4339 * DEG_TO_RAD;
        else if (EQn(b->GetType(), "As3+3", 5))
          phi = 86.9735 * DEG_TO_RAD;
        else if (EQn(b->GetType(), "Sb3+3", 5))
          phi = 87.7047 * DEG_TO_RAD;
        else
          phi = 90.0 * DEG_TO_RAD;

        oopcalc.c1 = -4.0 * cos(phi);
        oopcalc.c2 = 1.0;
        oopcalc.c0 = -1.0*oopcalc.c1 * cos(phi) + oopcalc.c2*cos(2.0*phi);
        oopcalc.koop = 22.0 * KCAL_TO_KJ;
      }
      else if (!(EQn(b->GetType(), "C_2", 3) || EQn(b->GetType(), "C_R", 3)))
        continue; // inversion not defined for this atom type

      // C atoms, we should check if we're bonded to O
      FOR_NBORS_OF_ATOM(nbr, b) {
        if (a == NULL)
          a = (OBAtom*) &*nbr;
        else if (c == NULL)
          c = (OBAtom*) &*nbr;
        else
          d = (OBAtom*) &*nbr;
      }
	  
      if ((a == NULL) || (c == NULL) || (d == NULL))
        continue;
 
      // skip this oop if the atoms are ignored 
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) ||
           _constraints.IsIgnored(c->GetIdx()) || _constraints.IsIgnored(d->GetIdx()) ) 
        continue;
 
      // if there are any groups specified, check if the four oop atoms are in a single intraGroup
      if (HasGroups()) {
        bool validOOP = false;
        for (unsigned int i=0; i < _intraGroup.size(); ++i) {
          if (_intraGroup[i].BitIsOn(a->GetIdx()) && _intraGroup[i].BitIsOn(b->GetIdx()) && 
              _intraGroup[i].BitIsOn(c->GetIdx()) && _intraGroup[i].BitIsOn(d->GetIdx()))
            validOOP = true;
        }
        if (!validOOP)
          continue;
      }
 
      if (EQn(b->GetType(), "C_2", 3) || EQn(b->GetType(), "C_R", 3)) {
        oopcalc.c0 = 1.0;
        oopcalc.c1 = -1.0;
        oopcalc.c2 = 0.0;
        oopcalc.koop = 6.0 * KCAL_TO_KJ;
        if (EQn(a->GetType(), "O_2", 3) ||
            EQn(c->GetType(), "O_2", 3) ||
            EQn(d->GetType(), "O_2", 3)) {
          oopcalc.koop = 50.0 * KCAL_TO_KJ;
        }
      }
            
      // A-B-CD || C-B-AD  PLANE = ABC
      oopcalc.a = a;
      oopcalc.b = b;
      oopcalc.c = c;
      oopcalc.d = d;
      oopcalc.koop /= 3.0; // three OOPs to consider
	    
      oopcalc.SetupPointers();
      _oopcalculations.push_back(oopcalc);
      
      // C-B-DA || D-B-CA  PLANE BCD
      oopcalc.a = d;
      oopcalc.d = a;
      
      oopcalc.SetupPointers();
      _oopcalculations.push_back(oopcalc);
      
      // A-B-DC || D-B-AC  PLANE ABD
      oopcalc.a = a;
      oopcalc.c = d;
      oopcalc.d = c;
	    
      oopcalc.SetupPointers();
      _oopcalculations.push_back(oopcalc);
    } // for all atoms

    // 
    // VDW Calculations
    //
    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP VAN DER WAALS CALCULATIONS...\n");
    
    OBFFVDWCalculationUFF vdwcalc;

    _vdwcalculations.clear();
    
    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);

      // skip this vdw if the atoms are ignored 
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
        continue;
  
      // if there are any groups specified, check if the two atoms are in a single _interGroup or if
      // two two atoms are in one of the _interGroups pairs.
      if (HasGroups()) {
        bool validVDW = false;
        for (unsigned int i=0; i < _interGroup.size(); ++i) {
          if (_interGroup[i].BitIsOn(a->GetIdx()) && _interGroup[i].BitIsOn(b->GetIdx())) 
            validVDW = true;
        }
        for (unsigned int i=0; i < _interGroups.size(); ++i) {
          if (_interGroups[i].first.BitIsOn(a->GetIdx()) && _interGroups[i].second.BitIsOn(b->GetIdx())) 
            validVDW = true;
          if (_interGroups[i].first.BitIsOn(b->GetIdx()) && _interGroups[i].second.BitIsOn(a->GetIdx())) 
            validVDW = true;
        }
 
        if (!validVDW)
          continue;
      }
 
      if (a->IsConnected(b)) {
        continue;
      }
      if (a->IsOneThree(b)) {
        continue;
      }

      parameterA = GetParameterUFF(a->GetType(), _ffparams);
      parameterB = GetParameterUFF(b->GetType(), _ffparams);

      if (parameterA == NULL || parameterB == NULL) { 
        IF_OBFF_LOGLVL_LOW {
          snprintf(_logbuf, BUFF_SIZE, "    COULD NOT FIND PARAMETERS FOR VDW INTERACTION %d-%d (IDX)...\n", 
              a->GetIdx(), b->GetIdx());
          OBFFLog(_logbuf);
        }
        return false;
      }
 
      vdwcalc.Ra = parameterA->_dpar[2];
      vdwcalc.ka = parameterA->_dpar[3];
      vdwcalc.Rb = parameterB->_dpar[2];
      vdwcalc.kb = parameterB->_dpar[3];

      vdwcalc.a = &*a;
      vdwcalc.b = &*b;
     
      //this calculations only need to be done once for each pair, 
      //we do them now and save them for later use
      vdwcalc.kab = KCAL_TO_KJ * sqrt(vdwcalc.ka * vdwcalc.kb);
      
      // 1-4 scaling
      // This isn't mentioned in the UFF paper, but is common for other methods
      //       if (a->IsOneFour(b))
      //         vdwcalc.kab *= 0.5;

      // ka now represents the xij in equation 20 -- the expected vdw distance
      vdwcalc.ka = sqrt(vdwcalc.Ra * vdwcalc.Rb);
      
      vdwcalc.SetupPointers();
      _vdwcalculations.push_back(vdwcalc);
    }
    
    // NOTE: No electrostatics are set up
    // If you want electrostatics with UFF (not a good idea), you will need to call SetupElectrostatics
    
    return true;
  }

  bool OBForceFieldUFF::SetupElectrostatics()
  {
    // 
    // Electrostatic Calculations
    //
    OBAtom *a, *b;

    IF_OBFF_LOGLVL_LOW
      OBFFLog("SETTING UP ELECTROSTATIC CALCULATIONS...\n");
 
    OBFFElectrostaticCalculationUFF elecalc;

    _electrostaticcalculations.clear();
    
    // Note that while the UFF paper mentions an electrostatic term,
    // it does not actually use it. Both Towhee and the UFF FAQ
    // discourage the use of electrostatics with UFF.

    FOR_PAIRS_OF_MOL(p, _mol) {
      a = _mol.GetAtom((*p)[0]);
      b = _mol.GetAtom((*p)[1]);
      
      // skip this ele if the atoms are ignored 
      if ( _constraints.IsIgnored(a->GetIdx()) || _constraints.IsIgnored(b->GetIdx()) )
        continue;
  
      // if there are any groups specified, check if the two atoms are in a single _interGroup or if
      // two two atoms are in one of the _interGroups pairs.
      if (HasGroups()) {
        bool validEle = false;
        for (unsigned int i=0; i < _interGroup.size(); ++i) {
          if (_interGroup[i].BitIsOn(a->GetIdx()) && _interGroup[i].BitIsOn(b->GetIdx())) 
            validEle = true;
        }
        for (unsigned int i=0; i < _interGroups.size(); ++i) {
          if (_interGroups[i].first.BitIsOn(a->GetIdx()) && _interGroups[i].second.BitIsOn(b->GetIdx())) 
            validEle = true;
          if (_interGroups[i].first.BitIsOn(b->GetIdx()) && _interGroups[i].second.BitIsOn(a->GetIdx())) 
            validEle = true;
        }
 
        if (!validEle)
          continue;
      }
 
      if (a->IsConnected(b)) {
        continue;
      }
      if (a->IsOneThree(b)) {
        continue;
      }
      
      // Remember that at the moment, this term is not currently used
      // These are also the Gasteiger charges, not the Qeq mentioned in the UFF paper
      elecalc.qq = KCAL_TO_KJ * 332.0637 * a->GetPartialCharge() * b->GetPartialCharge();
      
      if (elecalc.qq) {
        elecalc.a = &*a;
        elecalc.b = &*b;
        
        elecalc.SetupPointers();
        _electrostaticcalculations.push_back(elecalc);
      }
    }
    return true;
  }

  bool OBForceFieldUFF::ParseParamFile()
  {
    vector<string> vs;
    char buffer[BUFF_SIZE];
    
    OBFFParameter parameter;

    // open data/UFF.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "UFF.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open UFF.prm", obError);
      return false;
    }

    // Set the locale for number parsing to avoid locale issues: PR#1785463
    obLocale.SetLocale();

    while (ifs.getline(buffer, BUFF_SIZE)) {
      tokenize(vs, buffer);
      if (vs.size() < 13)
        continue;
      
      if (EQn(buffer, "param", 5)) {
        // set up all parameters from this
        parameter.clear();
        parameter._a = vs[1]; // atom type
        parameter._dpar.push_back(atof(vs[2].c_str())); // r1
        parameter._dpar.push_back(atof(vs[3].c_str())); // theta0
        parameter._dpar.push_back(atof(vs[4].c_str())); // x1
        parameter._dpar.push_back(atof(vs[5].c_str())); // D1
        parameter._dpar.push_back(atof(vs[6].c_str())); // zeta
        parameter._dpar.push_back(atof(vs[7].c_str())); // Z1
        parameter._dpar.push_back(atof(vs[8].c_str())); // Vi
        parameter._dpar.push_back(atof(vs[9].c_str())); // Uj
        parameter._dpar.push_back(atof(vs[10].c_str())); // Xi
        parameter._dpar.push_back(atof(vs[11].c_str())); // Hard
        parameter._dpar.push_back(atof(vs[12].c_str())); // Radius

        char coord = vs[1][2]; // 3rd character of atom type
        switch (coord) {
        case '1': // linear
          parameter._ipar.push_back(1);
          break;
        case '2': // trigonal planar (sp2)
        case 'R': // aromatic (N_R)
          parameter._ipar.push_back(2);
          break;
        case '3': // tetrahedral (sp3)
          parameter._ipar.push_back(3);
          break;
        case '4': // square planar
          parameter._ipar.push_back(4);
          break;
        case '5': // trigonal bipyramidal -- not actually in parameterization
          parameter._ipar.push_back(5);
          break;
        case '6': // octahedral
          parameter._ipar.push_back(6);
          break;
        case '7': // pentagonal bipyramidal -- not actually in parameterization
            parameter._ipar.push_back(7);
            break;
        default: // general case (unknown coordination)
          // These atoms appear to generally be linear coordination like Cl
          parameter._ipar.push_back(1);
        }
        
        _ffparams.push_back(parameter);
      }
    }
	
    if (ifs)
      ifs.close();
 
    // return the locale to the original one
    obLocale.RestoreLocale();
 
    return 0;
  }
  
  bool OBForceFieldUFF::SetTypes()
  {
    vector<vector<int> > _mlist; //!< match list for atom typing
    vector<pair<OBSmartsPattern*,string> > _vexttyp; //!< external atom type rules
    vector<vector<int> >::iterator j;
    vector<pair<OBSmartsPattern*,string> >::iterator i;
    OBSmartsPattern *sp;
    vector<string> vs;
    char buffer[BUFF_SIZE];
 
    _mol.SetAtomTypesPerceived();
    
    // open data/UFF.prm
    ifstream ifs;
    if (OpenDatafile(ifs, "UFF.prm").length() == 0) {
      obErrorLog.ThrowError(__FUNCTION__, "Cannot open UFF.prm", obError);
      return false;
    }

    while (ifs.getline(buffer, BUFF_SIZE)) {
      if (EQn(buffer, "atom", 4)) {
      	tokenize(vs, buffer);

        sp = new OBSmartsPattern;
        if (sp->Init(vs[1])) {
          _vexttyp.push_back(pair<OBSmartsPattern*,string> (sp,vs[2]));
        }
        else {
          delete sp;
          sp = NULL;
          obErrorLog.ThrowError(__FUNCTION__, " Could not parse atom type table from UFF.prm", obInfo);
          return false;
        }    
      }
    }
 
    for (i = _vexttyp.begin();i != _vexttyp.end();++i) {
      if (i->first->Match(_mol)) {
        _mlist = i->first->GetMapList();
        for (j = _mlist.begin();j != _mlist.end();++j) {
          _mol.GetAtom((*j)[0])->SetType(i->second);
        }
      }
    }
    IF_OBFF_LOGLVL_LOW {
      OBFFLog("\nA T O M   T Y P E S\n\n");
      OBFFLog("IDX\tTYPE\n");
      
      FOR_ATOMS_OF_MOL (a, _mol) {
        snprintf(_logbuf, BUFF_SIZE, "%d\t%s\n", a->GetIdx(), a->GetType());
        OBFFLog(_logbuf);
      }

    }
     
    if (ifs)
      ifs.close();

    return true;
  }
  
  double OBForceFieldUFF::Energy(bool gradients)
  {
    double energy;
    
    IF_OBFF_LOGLVL_MEDIUM
      OBFFLog("\nE N E R G Y\n\n");
    
    if (gradients) {
      ClearGradients();
      energy  = E_Bond<true>();
      energy += E_Angle<true>();
      energy += E_Torsion<true>();
      energy += E_OOP<true>();
      energy += E_VDW<true>();
    } else {
      energy  = E_Bond<false>();
      energy += E_Angle<false>();
      energy += E_Torsion<false>();
      energy += E_OOP<false>();
      energy += E_VDW<false>();
    }

    // The electrostatic term, by default is 0.0
    // You will need to call SetupEletrostatics if you want it
    // energy += E_Electrostatic(gradients);
     
    IF_OBFF_LOGLVL_MEDIUM {
      snprintf(_logbuf, BUFF_SIZE, "\nTOTAL ENERGY = %8.5f %s\n", energy, GetUnit().c_str());
      OBFFLog(_logbuf);
    }

    return energy;
  }
  
  OBFFParameter* OBForceFieldUFF::GetParameterUFF(std::string a, vector<OBFFParameter> &parameter)
  {
    for (unsigned int idx=0; idx < parameter.size(); ++idx) {
      if (a == parameter[idx]._a) {
        return &parameter[idx];
      }
    }
    return NULL;
  }
 
  bool OBForceFieldUFF::ValidateGradients ()
  {
    vector3 numgrad, anagrad, err;
    bool passed = true; // set to false if any component fails
    int coordIdx;
    
    OBFFLog("\nV A L I D A T E   G R A D I E N T S\n\n");
    OBFFLog("ATOM IDX      NUMERICAL GRADIENT           ANALYTICAL GRADIENT        REL. ERROR (%)   \n");
    OBFFLog("----------------------------------------------------------------------------------------\n");
    //     "XX       (000.000, 000.000, 000.000)  (000.000, 000.000, 000.000)  (00.00, 00.00, 00.00)"
   
    FOR_ATOMS_OF_MOL (a, _mol) {
      coordIdx = (a->GetIdx() - 1) * 3;

      // OBFF_ENERGY (i.e., overall)
      numgrad = NumericalDerivative(&*a, OBFF_ENERGY);
      Energy(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "%2d       (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", a->GetIdx(), numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      
      // OBFF_EBOND
      numgrad = NumericalDerivative(&*a, OBFF_EBOND);
      ClearGradients();
      E_Bond(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    bond    (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;
      
      // OBFF_EANGLE
      numgrad = NumericalDerivative(&*a, OBFF_EANGLE);
      ClearGradients();
      E_Angle(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    angle   (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;
      
      // OBFF_ETORSION
      numgrad = NumericalDerivative(&*a, OBFF_ETORSION);
      ClearGradients();
      E_Torsion(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    torsion (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      // 8% tolerance here because some 180 torsions cause numerical instability
      if (err.x() > 8.0 || err.y() > 8.0 || err.z() > 8.0)
        passed = false;
      
      // OBFF_EOOP
      numgrad = NumericalDerivative(&*a, OBFF_EOOP);
      ClearGradients();
      E_OOP(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    oop     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
//      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
//        passed = false;
        
      // OBFF_EVDW
      numgrad = NumericalDerivative(&*a, OBFF_EVDW);
      ClearGradients();
      E_VDW(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    vdw     (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;

      // OBFF_EELECTROSTATIC
      numgrad = NumericalDerivative(&*a, OBFF_EELECTROSTATIC);
      ClearGradients();
      E_Electrostatic(); // compute
      anagrad.Set(_gradientPtr[coordIdx], _gradientPtr[coordIdx+1], _gradientPtr[coordIdx+2]);
      err = ValidateGradientError(numgrad, anagrad);

      snprintf(_logbuf, BUFF_SIZE, "    electro (%7.3f, %7.3f, %7.3f)  (%7.3f, %7.3f, %7.3f)  (%5.2f, %5.2f, %5.2f)\n", numgrad.x(), numgrad.y(), numgrad.z(), 
              anagrad.x(), anagrad.y(), anagrad.z(), err.x(), err.y(), err.z());
      OBFFLog(_logbuf);
      if (err.x() > 5.0 || err.y() > 5.0 || err.z() > 5.0)
        passed = false;
    }
    
    return passed; // did we pass every single component?
  }
  
} // end namespace OpenBabel

//! \file forcefieldUFF.cpp
//! \brief UFF force field
