/**********************************************************************
Copyright (C) 2003 by Michael Banck <mbanck@gmx.net>
Some portions Copyright (C) 2004 by Chris Morley

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

#include "mol.h"
#include "obconversion.h"
#include <stdlib.h>

using namespace std;
namespace OpenBabel {

class CHTFormat : public OBFormat
{
public:
	//Register this format type ID
	CHTFormat() {OBConversion::RegisterFormat("CHT",this);}

	virtual const char* Description() //required
	{ return
"ChemTool\n \
No comments yet\n \
";
	};

	virtual const char* SpecificationURL(){return
		"http://ruby.chemie.uni-freiburg.de/~martin/chemtool/chemtool.html";}; //optional

	//Flags() can return be any the following combined by | or be omitted if none apply
	// NOTREADABLE  READONEONLY  NOTWRITABLE  WRITEONEONLY 
	virtual unsigned int Flags(){return NOTREADABLE;};

	////////////////////////////////////////////////////
	/// The "API" interface functions
	virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

////////////////////////////////////////////////////
	/// The "Convert" interface functions
	virtual bool WriteChemObject(OBConversion* pConv)
	{
		//Retrieve the target OBMol
		OBBase* pOb = pConv->GetChemObject();
		OBMol* pmol = dynamic_cast<OBMol*> (pOb);
		bool ret=false;
		if(pmol)
			ret=WriteMolecule(pmol,pConv);
		delete pOb; 
		return ret;
	};
};

//Make an instance of the format class
CHTFormat theCHTFormat;

/////////////////////////////////////////////////////////////////

bool CHTFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
{
	OBMol* pmol = dynamic_cast<OBMol*>(pOb);
	if(pmol==NULL) return false;

	//Define some references so we can use the old parameter names
	ostream &ofs = *pConv->GetOutStream();
	OBMol &mol = *pmol;
	const char *dimension = pConv->GetDimension();

	char buffer[BUFF_SIZE];
	int w, h, x, y; 	// to calculate the geometry
	int bondtype;		// type of bond
	int conv_factor = 50;	// please adjust
	int natoms = 0;		// number of additional (non-carbon) atoms
	OBAtom *atom, *atom1, *atom2;
	OBBond *bond;

	ofs << "Chemtool Version 1.4" << endl;

	// get the geometry
	w = 0;
	h = 0;
	vector<OBNodeBase*>::iterator i;
	for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
		x = (int)(atom->GetX()) * conv_factor;
		y = (int)(atom->GetY()) * conv_factor;
		if (x > w) w = x;
		if (y > h) h = y;
		if (atom->GetAtomicNum() != 6) natoms++;
	}
	ofs << "geometry " << w * 1.1 << " " << h * 1.1 << endl;

	// write out bonds
	ofs << "bonds "<< mol.NumBonds() << endl;
	vector<OBEdgeBase*>::iterator j;
	for(bond = mol.BeginBond(j); bond; bond = mol.NextBond(j)) {
		bondtype = 0;
		atom1 = bond->GetBeginAtom();
		atom2 = bond->GetEndAtom();
		if (bond->GetBO() == 2) bondtype = 1;
		if (bond->GetBO() == 3) bondtype = 3;
		// FIXME: use flag-info, too
		sprintf(buffer, "%d\t%d\t%d\t%d\t%1d",
			(int)floor(atom1->GetX() * conv_factor + 0.5),
			(int)floor(atom1->GetY() * conv_factor + 0.5),
			(int)floor(atom2->GetX() * conv_factor + 0.5),
			(int)floor(atom2->GetY() * conv_factor + 0.5),
			bondtype);
		ofs << buffer << endl;
	}

	// start over, write additional atoms
	ofs << "atoms " << natoms << endl;
        for (atom = mol.BeginAtom(i); atom; atom = mol.NextAtom(i)) {
		// Carbon does not need to be treated
		if (atom->GetAtomicNum() != 6) { 
			sprintf(buffer, "%d\t%d\t%s\t%d",
				(int)floor(atom->GetX() * conv_factor + 0.5),
				(int)floor(atom->GetY() * conv_factor + 0.5),
				etab.GetSymbol(atom->GetAtomicNum()),
				-1 // assume centered Text
				);
			ofs << buffer << endl;
		}
	}	

	// We don't have any splines to write
	ofs << "splines 0" << endl;

	return true;
}

} //namespace OpenBabel
