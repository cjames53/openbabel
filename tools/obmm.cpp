/**********************************************************************
obmm.cpp - openbabel molecular mechanics program

commands:            description:
load <filename>      load a molecule from filename				done
save <filename>      save currently loaded molecule to filename			done
ff <forcefield>      select the force field
forcefields          print the available forcefields

energy               calculate the energy					done
ebond                calculate the bond stretching energy			"
eangle               calculate the angle bending energy				"
estrbnd              calculate the stretch-bending enregy			"
eoop                 calculate the out-of-plane bending energy			"
etorsion             calculate the torsional energy				"
evdw                 calculate the Van der Waals energy				"
eeq                  calculate the electrostatic energy				"

sd <n>               steepest descent energy minimization for n steps
cg <n>               conjugate gradients energy minimization for n steps	todo

addH                 add hydrogens (for smiles, ...)				done

quit                 quit							done

Copyright (C) 2006 Tim Vandermeersch
Some portions Copyright (C) 2006 Geoffrey R. Hutchison
 
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

// used to set import/export for Cygwin DLLs
#ifdef WIN32
#define USING_OBDLL
#endif

#include <openbabel/babelconfig.h>
#include <openbabel/base.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/forcefield.h>
#include <unistd.h>

using namespace std;
using namespace OpenBabel;

// PROTOTYPES /////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//! \brief  Generate rough 3D coordinates for SMILES (or other 0D files)
//          based on bonding network, rings, atom types, etc.
//
int main(int argc,char **argv)
{
  OBForceField* pFF = OBForceField::FindForceField("MMFF94");
  OBMol mol;
  char commandline[100];
  vector<string> vs;

  mol.Clear();

  cout << endl;
  cout << "openbabel                            " << endl;
  cout << "M O L E C U L A R   M E C H A N I C S" << endl;
  cout << "                              program" << endl;
  cout << "                v 0.1                " << endl << endl;

  while (1) {

    cout << "command > ";
    cin.getline(commandline, 100);

    //
    // commands with no parameters
    //
    if (EQn(commandline, "quit", 4) || cin.eof()) {
       cout << "bye." << endl;
       exit(0);
    }

    if (EQn(commandline, "energy", 6)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  total energy = " << pFF->GetEnergy() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }

    if (EQn(commandline, "ebond", 5)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  bond stretching energy = " << pFF->E_Bond() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }

    if (EQn(commandline, "eangle", 6)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  angle bending energy = " << pFF->E_Angle() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }

    if (EQn(commandline, "estrbnd", 7)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  stretch-bending energy = " << pFF->E_StrBnd() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }

    if (EQn(commandline, "eoop", 4)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  out-of-plane bending energy = " << pFF->E_OOP() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }

    if (EQn(commandline, "etorsion", 8)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  torsional energy = " << pFF->E_Torsion() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }

    if (EQn(commandline, "evdw", 4)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  Van der Waals energy = " << pFF->E_VDW() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }
    
    if (EQn(commandline, "eeq", 3)) {
      if (mol.Empty()) {
        cout << "no molecule loaded." << endl;
	continue;
      }
      cout << endl << "  electrostatic energy = " << pFF->E_Electrostatic() << " " << pFF->GetUnit() << endl << endl;
      continue;
    }
    
    if (EQn(commandline, "addH", 4)) {
      int num1, num2;
      num1 = mol.NumAtoms();
      mol.AddHydrogens(false, true);
      num2 = mol.NumAtoms();
      cout << (num2 - num1) << " hydrogens added." << endl;
      continue;
    }


    //
    // commands with parameters
    //
    tokenize(vs, commandline);
    
    // load <filename>
    if (EQn(commandline, "load", 4)) {
      if (vs.size() < 2) {
        cout << "no <filename> specified." << endl;
        continue;
      }
      
      ifstream ifs;
      OBConversion conv;
      OBFormat *format_in = conv.FormatFromExt(vs[1].c_str());
   
      if (!format_in || !conv.SetInFormat(format_in)) {
        cout << "could not detect format." << endl;
	continue;
      }
       
      ifs.open(vs[1].c_str());
      if (!ifs) {
        cout << "could not open '" << vs[1] << "'." <<endl;
        continue;
      }
      
      mol.Clear();
      if (!conv.Read(&mol, &ifs)) {
        cout << "could not read a molecule from '" << vs[1] << "'." <<endl;
        continue;
      }
      
      if (mol.Empty()) {
        cout << "this molecule is empty." <<endl;
        continue;
      }

      if (!pFF->Setup(mol)) {
        cout << "error while initializing the force field for this molecule." <<endl;
        continue;
      }

      cout << "molecule succesfully loaded." << endl;
      cout << "  " << mol.NumAtoms() << " atoms" << endl;
      cout << "  " << mol.NumBonds() << " bonds" << endl;

      ifs.close();
 
      continue;
    }
    
    // save <filename>
    if (EQn(commandline, "save", 4)) {
      if (vs.size() < 2) {
        cout << "no <filename> specified." << endl;
        continue;
      }
      
      ofstream ofs;
      OBConversion conv;
      OBFormat *format_out = conv.FormatFromExt(vs[1].c_str());
   
      if (!format_out || !conv.SetOutFormat(format_out)) {
        cout << "could not detect format." << endl;
	continue;
      }
       
      ofs.open(vs[1].c_str());
      if (!ofs) {
        cout << "could not open '" << vs[1] << "'." <<endl;
        continue;
      }
      
      if (!conv.Write(&mol, &ofs)) {
        cout << "could not read a molecule from '" << vs[1] << "'." <<endl;
        continue;
      }
      
      cout << "molecule succesfully saved." << endl;
      cout << "  " << mol.NumAtoms() << " atoms" << endl;
      cout << "  " << mol.NumBonds() << " bonds" << endl;

      ofs.close();
 
      continue;
    }

    if (EQn(commandline, "sd", 2)) {
      if (vs.size() < 2) {
        cout << "no <n> steps specified." << endl;
        continue;
      }

      pFF->SteepestDescent(atoi(vs[1].c_str()));
      pFF->UpdateCoordinates(mol);

      continue;
    }

    cout << "invalid command." << endl;
  }

  return(1);
}
