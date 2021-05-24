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
#include "FEBioLoadDataSection.h"
#include <FECore/FEModel.h>
#include <FECore/FELoadCurve.h>
#include <FECore/FEPointFunction.h>

//-----------------------------------------------------------------------------
FEBioLoadDataSection::FEBioLoadDataSection(FEFileImport* pim) : FEFileSection(pim) 
{
	m_redefineCurves = false;
}

//-----------------------------------------------------------------------------
//!  This function reads the load data section from the xml file
//!
void FEBioLoadDataSection::Parse(XMLTag& tag)
{
	FEModel& fem = *GetFEModel();
	
	++tag;
	do
	{
		if (tag == "loadcurve")
		{
			// load curve ID
			int nid;
			tag.AttributeValue("id", nid);

			// default type and extend mode
			FEPointFunction::INTFUNC ntype = FEPointFunction::LINEAR;
			FEPointFunction::EXTMODE nextm = FEPointFunction::CONSTANT;

			// get the (optional) type
			XMLAtt* patt = tag.Attribute("type", true);
			if (patt)
			{
				XMLAtt& type = *patt;
				if      (type == "step"  ) ntype = FEPointFunction::STEP;
				else if (type == "linear") ntype = FEPointFunction::LINEAR;
				else if (type == "smooth") ntype = FEPointFunction::SMOOTH;
                else if (type == "cubic spline") ntype = FEPointFunction::CSPLINE;
                else if (type == "control points") ntype = FEPointFunction::CPOINTS;
                else if (type == "approximation" ) ntype = FEPointFunction::APPROX;
				else throw XMLReader::InvalidAttributeValue(tag, "type", type.cvalue());
			}

			// get the optional extend mode
			patt = tag.Attribute("extend", true);
			if (patt)
			{
				XMLAtt& ext = *patt;
				if      (ext == "constant"     ) nextm = FEPointFunction::CONSTANT;
				else if (ext == "extrapolate"  ) nextm = FEPointFunction::EXTRAPOLATE;
				else if (ext == "repeat"       ) nextm = FEPointFunction::REPEAT;
				else if (ext == "repeat offset") nextm = FEPointFunction::REPEAT_OFFSET;
				else throw XMLReader::InvalidAttributeValue(tag, "extend", ext.cvalue());
			}

			// get the number of load curves
			int nlc = fem.LoadControllers();

			// find or create the load curve
			FELoadCurve* plc = 0;

			// see if this refers to a valid curve
			if (m_redefineCurves)
			{
				if ((nid > 0) && (nid <= nlc)) 
				{
					plc = dynamic_cast<FELoadCurve*>(fem.GetLoadController(nid - 1));
					assert(plc);

					// clear the curve since we're about to read in new data points
					plc->Clear();
				}
			}

			// if the ID does not refer to an existing curve, make sure it defines the next curve
			if (plc == 0)
			{
				// check that the ID is one more than the number of load curves defined
				// This is to make sure that the ID's are in numerical order and no values are skipped.
				if (nid != nlc + 1) throw XMLReader::InvalidAttributeValue(tag, "id");

				// create the loadcurve
				plc = fecore_new<FELoadCurve>("loadcurve", &fem); assert(plc);
				plc->SetID(nid - 1);

				// add the loadcurve
				fem.AddLoadController(plc);
			}

			// set the load curve attributes
			plc->SetInterpolation(ntype);
			plc->SetExtendMode(nextm);

			// read the data points
			double d[2];
			++tag;
			do
			{
				tag.value(d, 2);
				plc->Add(d[0], d[1]);
				++tag;
			}
			while (!tag.isend());
		}
		else throw XMLReader::InvalidTag(tag);

		++tag;
	}
	while (!tag.isend());
}
