// FERestartImport.h: interface for the FERestartImport class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FERESTARTIMPORT_H__A5A88D72_026C_45F5_BECB_5B3C7B3C767C__INCLUDED_)
#define AFX_FERESTARTIMPORT_H__A5A88D72_026C_45F5_BECB_5B3C7B3C767C__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FileImport.h"
#include "XMLReader.h"

//-----------------------------------------------------------------------------
//! Restart input file reader.
class FECORE_EXPORT FERestartImport : public FEFileImport
{
public:
	FERestartImport();
	virtual ~FERestartImport();

	bool Load(FEModel& fem, const char* szfile);

protected:
	bool ParseControlSection (XMLTag& tag);
	bool ParseLoadSection    (XMLTag& tag);

public:
	char		m_szdmp[256];	// user defined restart file name

protected:
	FEModel*	m_pfem;			// point to the FEM
	XMLReader	m_xml;			// the file reader
};

#endif // !defined(AFX_FERESTARTIMPORT_H__A5A88D72_026C_45F5_BECB_5B3C7B3C767C__INCLUDED_)
