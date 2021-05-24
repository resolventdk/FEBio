/*This file is part of the FEBio source code and is licensed under the MIT license
listed below.

See Copyright-FEBio.txt for details.

Copyright (c) 2020 University of Utah, The Trustees of Columbia University in 
the City of New York, and others.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.*/



#include "stdafx.h"
#include "FEBioOutputSection.h"
#include <FECore/NodeDataRecord.h>
#include <FECore/FaceDataRecord.h>
#include <FECore/ElementDataRecord.h>
#include <FEBioMech/ObjectDataRecord.h>
#include <FECore/NLConstraintDataRecord.h>
#include <FECore/FEModel.h>
#include <FECore/FEModelData.h>
#include <FECore/FSPath.h>

bool string_to_int_vector(const char* szlist, std::vector<int>& list)
{
	list.clear();
	if (szlist == nullptr) return false;
	if (szlist[0] == 0) return true;

	char* ch = nullptr;
	char* sz = (char*)szlist;
	do
	{
		ch = strchr(sz, ',');
		if (ch) *ch = 0;
		int n0, n1, nn;
		int nread = sscanf(sz, "%d:%d:%d", &n0, &n1, &nn);
		switch (nread)
		{
		case 0: return false;
			break;
		case 1:
			n1 = n0;
			nn = 1;
			break;
		case 2:
			nn = 1;
		}

		for (int i = n0; i <= n1; i += nn) list.push_back(i);

		if (ch) *ch = ',';
		sz = ch + 1;
	}
	while (ch != 0);

	return true;
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::Parse(XMLTag& tag)
{
	++tag;
	do
	{
		if      (tag == "logfile" ) ParseLogfile(tag);
		else if (tag == "plotfile") ParsePlotfile(tag);
		else if (tag == "data"    ) ParseDataSection(tag);
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParseDataSection(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	++tag;
	do
	{
		if (tag == "element_data")
		{
			const char* szdata = tag.AttributeValue("data");

			FELogElemData* pd = fecore_new<FELogElemData>(szdata, &fem);
			if (pd == 0) throw XMLReader::InvalidAttributeValue(tag, "data", szdata);

			vector<int> items;
			tag.value(items);

			FEModelData* data = new FEModelData(&fem, pd, items);
			data->SetName(szdata);
			fem.AddModelData(data);
		}
		else throw XMLReader::InvalidTag(tag);
		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParseLogfile(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// Get the feb file path
	const char* szpath = GetFileReader()->GetFilePath();

	// see if the log file has any attributes
	const char* szlog = tag.AttributeValue("file", true);

	// If the log filename is not a path, but just a filename, make the filename
	// relative to the input file, otherwise, leave the path alone.
	if(szlog)
	{
		if(!FSPath::isPath(szlog))
		{
			char szfile[1024] = {0};
			sprintf(szfile, "%s%s", szpath, szlog);
			GetFEBioImport()->SetLogfileName(szfile);
		}
		else
		{
			GetFEBioImport()->SetLogfileName(szlog);
		}

	}

	if (tag.isleaf()) return;
	++tag;
	do
	{
		// get the (optional) file attribute
		char szfilename[1024] = { 0 };
		const char* szfile = tag.AttributeValue("file", true);
		if (szfile)
		{
			// if we have a path, prepend the path's name
			if (szpath && szpath[0])
			{
				sprintf(szfilename, "%s%s", szpath, szfile);
			}
			else strcpy(szfilename, szfile);
			szfile = szfilename;
		}

		// get other attributes
		const char* szdelim = tag.AttributeValue("delim", true);
		const char* szformat = tag.AttributeValue("format", true);

		bool bcomment = true;
		const char* szcomment = tag.AttributeValue("comments", true);
		if (szcomment != 0)
		{
			if (strcmp(szcomment, "on") == 0) bcomment = true;
			else if (strcmp(szcomment, "off") == 0) bcomment = false;
		}


		if (tag == "node_data")
		{
			NodeDataRecord* prec = new NodeDataRecord(&fem, szfile);

			const char* szdata = tag.AttributeValue("data");
			prec->SetData(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname   != 0) prec->SetName(szname); else prec->SetName(szdata);
			if (szdelim  != 0) prec->SetDelim(szdelim);
			if (szformat != 0) prec->SetFormat(szformat);
			prec->SetComments(bcomment);

			const char* sztmp = "set";
			if (GetFileReader()->GetFileVersion() >= 0x0205) sztmp = "node_set";
			const char* sz = tag.AttributeValue(sztmp, true);
			if (sz)
			{
				FENodeSet* pns = mesh.FindNodeSet(sz);
				if (pns == 0) throw XMLReader::InvalidAttributeValue(tag, sztmp, sz);
				prec->SetNodeSet(pns);
			}
			else
			{
				std::vector<int> items;
				string_to_int_vector(tag.szvalue(), items);
				prec->SetItemList(items);
			}

			GetFEBioImport()->AddDataRecord(prec);
		}
		else if (tag == "face_data")
		{
			FaceDataRecord* prec = new FaceDataRecord(&fem, szfile);

			const char* szdata = tag.AttributeValue("data");
			prec->SetData(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname   != 0) prec->SetName(szname); else prec->SetName(szdata);
			if (szdelim  != 0) prec->SetDelim(szdelim);
			if (szformat != 0) prec->SetFormat(szformat);
			prec->SetComments(bcomment);

			const char* sz = tag.AttributeValue("surface");
			FESurface* surf = mesh.FindSurface(sz);
			if (surf == nullptr) throw XMLReader::InvalidAttributeValue(tag, "surface", sz);

			std::vector<int> items;
			string_to_int_vector(tag.szvalue(), items);

			prec->SetSurface(surf, items);

			GetFEBioImport()->AddDataRecord(prec);
		}
		else if (tag == "element_data")
		{
			ElementDataRecord* prec = new ElementDataRecord(&fem, szfile);

			const char* szdata = tag.AttributeValue("data");
			prec->SetData(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname   != 0) prec->SetName(szname); else prec->SetName(szdata);
			if (szdelim  != 0) prec->SetDelim(szdelim);
			if (szformat != 0) prec->SetFormat(szformat);
			prec->SetComments(bcomment);

			const char* sztmp = "elset";
			if (GetFileReader()->GetFileVersion() >= 0x0205) sztmp = "elem_set";

			const char* sz = tag.AttributeValue(sztmp, true);
			if (sz)
			{
				FEElementSet* pes = mesh.FindElementSet(sz);
				if (pes == 0) throw XMLReader::InvalidAttributeValue(tag, sztmp, sz);
				prec->SetElementSet(pes);
			}
			else
			{
				std::vector<int> items;
				string_to_int_vector(tag.szvalue(), items);
				prec->SetItemList(items);
			}

			GetFEBioImport()->AddDataRecord(prec);
		}
		else if (tag == "rigid_body_data")
		{
			ObjectDataRecord* prec = new ObjectDataRecord(&fem, szfile);

			const char* szdata = tag.AttributeValue("data");
			prec->SetData(szdata);

			const char* szname = tag.AttributeValue("name", true);
			if (szname   != 0) prec->SetName(szname); else prec->SetName(szdata);
			if (szdelim  != 0) prec->SetDelim(szdelim);
			if (szformat != 0) prec->SetFormat(szformat);
			prec->SetComments(bcomment);

			std::vector<int> items;
			string_to_int_vector(tag.szvalue(), items);
			prec->SetItemList(items);

			GetFEBioImport()->AddDataRecord(prec);
		}
        else if (tag == "rigid_connector_data")
        {
            NLConstraintDataRecord* prec = new NLConstraintDataRecord(&fem, szfile);
            
            const char* szdata = tag.AttributeValue("data");
            prec->SetData(szdata);
            
            const char* szname = tag.AttributeValue("name", true);
            if (szname   != 0) prec->SetName(szname); else prec->SetName(szdata);
            if (szdelim  != 0) prec->SetDelim(szdelim);
            if (szformat != 0) prec->SetFormat(szformat);
			prec->SetComments(bcomment);

			std::vector<int> items;
			string_to_int_vector(tag.szvalue(), items);
			prec->SetItemList(items);
            
			GetFEBioImport()->AddDataRecord(prec);
        }
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}

//-----------------------------------------------------------------------------
void FEBioOutputSection::ParsePlotfile(XMLTag &tag)
{
	FEModel& fem = *GetFEModel();

	// get the plot file type. Must be "febio"!
	const char* sz = tag.AttributeValue("type", true);
	if (sz)
	{
		if ((strcmp(sz, "febio") != 0) && (strcmp(sz, "febio2") != 0)) throw XMLReader::InvalidAttributeValue(tag, "type", sz);
	}
	else sz = "febio";
	strcpy(GetFEBioImport()->m_szplot_type, sz);

	// get the optional plot file name
	const char* szplt = tag.AttributeValue("file", true);

	// If the plot filename is not a path, but just a filename, make the filename
	// relative to the input file, otherwise, leave the path alone.
	if(szplt)
	{
		if(!FSPath::isPath(szplt))
		{
			// Get the feb file path
			const char* szpath = GetFileReader()->GetFilePath();

			char szfile[1024] = {0};
			sprintf(szfile, "%s%s", szpath, szplt);
			GetFEBioImport()->SetPlotfileName(szfile);
		}
		else
		{
			GetFEBioImport()->SetPlotfileName(szplt);
		}

	}

	// read and store the plot variables
	if (!tag.isleaf())
	{
		++tag;
		do
		{
			if (tag == "var")
			{
				// get the variable name
				const char* szt = tag.AttributeValue("type");

				// get the item list
				vector<int> item;
				if (tag.isempty() == false) tag.value(item);

                // see if a surface is referenced
                const char* szset = tag.AttributeValue("surface", true);
                if (szset)
                {
                    // make sure this tag does not have any children
                    if (!tag.isleaf()) throw XMLReader::InvalidTag(tag);
                    
                    // see if we can find the facet set
                    FEMesh& m = GetFEModel()->GetMesh();
					FEFacetSet* ps = m.FindFacetSet(szset);
                    
                    // create a surface from the facet set
                    if (ps)
                    {
                        // create a new surface
                        FESurface* psurf = fecore_alloc(FESurface, &fem);
                        fem.GetMesh().AddSurface(psurf);
                        if (GetBuilder()->BuildSurface(*psurf, *ps) == false) throw XMLReader::InvalidTag(tag);

                        // Add the plot variable
                        const std::string& surfName = psurf->GetName();
						GetFEBioImport()->AddPlotVariable(szt, item, surfName.c_str());
                    }
                    else throw XMLReader::InvalidAttributeValue(tag, "set", szset);
                }
                else
                {
                    // Add the plot variable
					GetFEBioImport()->AddPlotVariable(szt, item);
                }
			}
			else if (tag=="compression")
			{
				int ncomp;
				tag.value(ncomp);
				GetFEBioImport()->SetPlotCompression(ncomp);
			}
			++tag;
		}
		while (!tag.isend());
	}
}
