#pragma once
#include "FEBioImport.h"

//-----------------------------------------------------------------------------
// Material Section
class FEBIOXML_API FEBioMaterialSection : public FEFileSection
{
public:
	FEBioMaterialSection(FEFileImport* pim) : FEFileSection(pim){}
	void Parse(XMLTag& tag);

protected:
	FEMaterial* CreateMaterial(XMLTag& tag);

protected:
	int	m_nmat;
};
