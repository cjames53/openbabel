/**********************************************************************
Copyright (C) 2008 Geoffrey R. Hutchison
Some portions Copyright (C) 2004 by Chris Morley
 
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/babelconfig.h>
#include <openbabel/obmolecformat.h>

#include <vector>
#include <map>

#include <sstream>

using namespace std;
namespace OpenBabel
{

  class PQRFormat : public OBMoleculeFormat
  {
  public:
    //Register this format type ID
    PQRFormat()
    {
      OBConversion::RegisterFormat("pqr",this, "chemical/x-pqr");
    }

    virtual const char* Description() //required
    {
      return
        "PQR format\n"
        "Read Options e.g. -as\n"
        "  s  Output single bonds only\n"
        "  b  Disable bonding entirely\n\n";
    };

    virtual const char* SpecificationURL()
    { return "";};

    virtual const char* GetMIMEType() 
    { return "chemical/x-pqr"; };

    //*** This section identical for most OBMol conversions ***
    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual int SkipObjects(int n, OBConversion* pConv);
    virtual bool ReadMolecule(OBBase* pOb, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pOb, OBConversion* pConv);

  };
  //***

  //Make an instance of the format class
  PQRFormat thePQRFormat;

  /////////////////////////////////////////////////////////////////

  // In atomrecord.cpp (shared with PDB format)
  bool ParseAtomRecord(char *, OBMol &, int);
  double ParseAtomCharge(char *, OBMol &);
  double ParseAtomRadius(char *, OBMol &);
  /////////////////////////////////////////////////////////////////
  int PQRFormat::SkipObjects(int n, OBConversion* pConv)
  {
    if (n == 0)
      ++ n;
    istream &ifs = *pConv->GetInStream();
    char buffer[BUFF_SIZE];
    while (n && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          -- n;
      }
      
    return ifs.good() ? 1 : -1;       
  }
  /////////////////////////////////////////////////////////////////
  bool PQRFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {

    OBMol* pmol = pOb->CastAndClear<OBMol>();
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    istream &ifs = *pConv->GetInStream();
    OBMol &mol = *pmol;
    const char* title = pConv->GetTitle();

    int chainNum = 1;
    char buffer[BUFF_SIZE];
    OBBitVec bs;
    vector<double> charges, radii;
    string line, key, value;

    mol.SetTitle(title);
    mol.SetChainsPerceived(); // It's a PDB-like file, we read all chain/res info.

    mol.BeginModify();
    while (ifs.good() && ifs.getline(buffer,BUFF_SIZE))
      {
        if (EQn(buffer,"ENDMDL",6))
          break;
        if (EQn(buffer,"END",3)) {
          // eat anything until the next ENDMDL
          while (ifs.getline(buffer,BUFF_SIZE) && !EQn(buffer,"ENDMDL",6));
          break;
        }
        if (EQn(buffer,"TER",3)) {
          chainNum++;
          continue;
        }
        if (EQn(buffer,"ATOM",4) || EQn(buffer,"HETATM",6))
          {
            ParseAtomRecord(buffer,mol,chainNum);
            if (EQn(buffer,"ATOM",4))
              bs.SetBitOn(mol.NumAtoms());

            // Read in the partial charge and radius too
            charges.push_back( ParseAtomCharge(buffer, mol) );
            radii.push_back( ParseAtomRadius(buffer, mol) );
            continue;
          }
        }

    if (!mol.NumAtoms()) { // skip the rest of this processing
      mol.EndModify();
      return(false);
    }

    // Use residue definitions to assign bond orders
    resdat.AssignBonds(mol,bs);

    mol.EndModify();

    /*Now assign hetatm bonds based on distance*/
    if (!pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.ConnectTheDots();

    if (!pConv->IsOption("s",OBConversion::INOPTIONS) && !pConv->IsOption("b",OBConversion::INOPTIONS))
      mol.PerceiveBondOrders();

    FOR_ATOMS_OF_MOL(a, mol) {
      // WARNING: Atom index issue here
      a->SetPartialCharge(charges[a->GetIdx() - 1]);

      cerr << " charge : " << charges[a->GetIdx() - 1] << endl;

      if (!a->HasData("Radius")) {
        std::ostringstream s;
        s << radii[ a->GetIdx()-1 ];
        OBPairData *p = new OBPairData;
        p->SetAttribute("Radius");
        p->SetValue( s.str() );
        a->SetData(p);
      }

      cerr << " radius : " << radii[a->GetIdx() - 1] << endl;

    }
    mol.SetPartialChargesPerceived();

    // clean out remaining blank lines
    while(ifs.peek() != EOF && ifs.good() && 
          (ifs.peek() == '\n' || ifs.peek() == '\r'))
      ifs.getline(buffer,BUFF_SIZE);

    return(true);
  }

  double ParseAtomCharge(char *buffer, OBMol &mol)
  // In PQR format, either:
  // Field name, atom number, atom name, residue name, residue number
  //    x y z charge radius
  // OR
  // Field, atom number, atom name, chain id, residue number, X, Y, Z, chg, rad
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 10)
      return atof(vs[8].c_str());
    else if (vs.size() == 11)
      return atof(vs[9].c_str());

    return 0.0;
  }

  double ParseAtomRadius(char *buffer, OBMol &mol)
  {
    vector<string> vs;
    tokenize(vs,buffer);

    OBAtom *atom = mol.GetAtom(mol.NumAtoms());

    if (vs.size() == 10)
      return atof(vs[9].c_str());
    else if (vs.size() == 11)
      return atof(vs[10].c_str());

    return 0.0;
  }

  bool PQRFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    OBMol* pmol = dynamic_cast<OBMol*>(pOb);
    if(pmol==NULL)
      return false;

    //Define some references so we can use the old parameter names
    ostream &ofs = *pConv->GetOutStream();
    OBMol &mol = *pmol;

    unsigned int i;
    char buffer[BUFF_SIZE];
    char type_name[10], padded_name[10];
    char the_res[10];
    char the_chain = ' ';
    const char *element_name;
    int res_num;
    bool het=true;
    int model_num = 0;
    if (!pConv->IsLast() || pConv->GetOutputIndex() > 1)
      { // More than one molecule record
      model_num = pConv->GetOutputIndex(); // MODEL 1-based index
      snprintf(buffer, BUFF_SIZE, "MODEL %8d", model_num);
      ofs << buffer << endl;
      }

    if (strlen(mol.GetTitle()) > 0)
      snprintf(buffer, BUFF_SIZE, "COMPND    %s ",mol.GetTitle());
    else
      snprintf(buffer, BUFF_SIZE, "COMPND    UNNAMED");
    ofs << buffer << endl;

    snprintf(buffer, BUFF_SIZE, "AUTHOR    GENERATED BY OPEN BABEL %s",BABEL_VERSION);
    ofs << buffer << endl;

    // before we write any records, we should check to see if any coord < -1000
    // which will cause errors in the formatting

    double minX, minY, minZ;
    minX = minY = minZ = -999.0f;
    FOR_ATOMS_OF_MOL(a, mol)
      {
        if (a->GetX() < minX)
          minX = a->GetX();
        if (a->GetY() < minY)
          minY = a->GetY();
        if (a->GetZ() < minZ)
          minZ = a->GetZ();
      }
    
    vector3 transV = VZero;
    if (minX < -999.0)
      transV.SetX( -1.0*minX - 900.0 );
    if (minY < -999.0)
      transV.SetY( -1.0*minY - 900.0 );
    if (minZ < -999.0)
      transV.SetZ( -1.0*minZ - 900.0 );

    // if minX, minY, or minZ was never changed, shift will be 0.0f
    // otherwise, move enough so that smallest coord is > -999.0f
    mol.Translate(transV);

    OBAtom *atom;
    OBResidue *res;
    for (i = 1; i <= mol.NumAtoms(); i++)
      {
        atom = mol.GetAtom(i);
        strncpy(type_name, etab.GetSymbol(atom->GetAtomicNum()), sizeof(type_name));
        type_name[sizeof(type_name) - 1] = '\0';

        //two char. elements are on position 13 and 14 one char. start at 14
        if (strlen(type_name) > 1)
          type_name[1] = toupper(type_name[1]);
        else
          {
            char tmp[10];
            strncpy(tmp, type_name, 10);
            snprintf(type_name, sizeof(type_name), " %-3s", tmp);
          }

        if ( (res = atom->GetResidue()) != 0 )
          {
            het = res->IsHetAtom(atom);
            snprintf(the_res,4,"%s",(char*)res->GetName().c_str());
            snprintf(type_name,5,"%s",(char*)res->GetAtomID(atom).c_str());
	    the_chain = res->GetChain();

            //two char. elements are on position 13 and 14 one char. start at 14
            if (strlen(etab.GetSymbol(atom->GetAtomicNum())) == 1)
              {
                if (strlen(type_name) < 4)
                  {
                    char tmp[16];
                    strncpy(tmp, type_name, 16);
                    snprintf(padded_name, sizeof(padded_name), " %-3s", tmp);
                    strncpy(type_name,padded_name,4);
                    type_name[4] = '\0';
                  }
                else
                  {
                    type_name[4] = '\0';
                  }
              }
            res_num = res->GetNum();
          }
        else
          {
            strcpy(the_res,"UNK");
            snprintf(padded_name,sizeof(padded_name), "%s",type_name);
            strncpy(type_name,padded_name,4);
            type_name[4] = '\0';
            res_num = 1;
          }

        element_name = etab.GetSymbol(atom->GetAtomicNum());
        //snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  \n",
        snprintf(buffer, BUFF_SIZE, "%s%5d %-4s %-3s %c%4d    %8.3f%8.3f%8.3f %11.8f%8.3f %2s  \n",
                 het?"HETATM":"ATOM  ",
                 i,
                 type_name,
                 the_res,
		 the_chain,
                 res_num,
                 atom->GetX(),
                 atom->GetY(),
                 atom->GetZ(),
                 atom->GetPartialCharge(),
                 etab.GetVdwRad(atom->GetAtomicNum()),
                 element_name);
        ofs << buffer;
      }

    OBAtom *nbr;
    vector<OBBond*>::iterator k;
    for (i = 1; i <= mol.NumAtoms(); i ++)
      {
        atom = mol.GetAtom(i);
        if (atom->GetValence() == 0)
          continue; // no need to write a CONECT record -- no bonds

        snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
        ofs << buffer;
        // Write out up to 4 real bonds per line PR#1711154
        int currentValence = 0;
        for (nbr = atom->BeginNbrAtom(k);nbr;nbr = atom->NextNbrAtom(k))
          {
            snprintf(buffer, BUFF_SIZE, "%5d", nbr->GetIdx());
            ofs << buffer;
            if (++currentValence % 4 == 0) {
              // Add the trailing space to finish this record
              ofs << "                                       \n";
              // write the start of a new CONECT record
              snprintf(buffer, BUFF_SIZE, "CONECT%5d", i);
              ofs << buffer;              
            }
          }

        // Add trailing spaces
        int remainingValence = atom->GetValence() % 4;
        for (int count = 0; count < (4 - remainingValence); count++) {
          snprintf(buffer, BUFF_SIZE, "     ");
          ofs << buffer;
        }
        ofs << "                                       \n";
      }

    snprintf(buffer, BUFF_SIZE, "MASTER        0    0    0    0    0    0    0    0 ");
    ofs << buffer;
    snprintf(buffer, BUFF_SIZE, "%4d    0 %4d    0\n",mol.NumAtoms(),mol.NumAtoms());
    ofs << buffer;

    ofs << "END\n";
    if (model_num) {
      ofs << "ENDMDL" << endl;
    }

    return(true);
  }


} //namespace OpenBabel
