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
#include "FEContinuousFiberDistributionUC.h"

BEGIN_FECORE_CLASS(FEContinuousFiberDistributionUC, FEUncoupledMaterial)
	// set material properties
	ADD_PROPERTY(m_pFmat, "fibers");
	ADD_PROPERTY(m_pFDD , "distribution");
	ADD_PROPERTY(m_pFint, "scheme");
END_FECORE_CLASS();

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::FEContinuousFiberDistributionUC(FEModel* pfem) : FEUncoupledMaterial(pfem)
{
	m_pFmat = 0;
	m_pFDD = 0;
	m_pFint = 0;
}

//-----------------------------------------------------------------------------
FEContinuousFiberDistributionUC::~FEContinuousFiberDistributionUC() {}

//-----------------------------------------------------------------------------
bool FEContinuousFiberDistributionUC::Init()
{
	// initialize fiber integration scheme
	if (FEUncoupledMaterial::Init() == false) return false;

	return true;
}

//-----------------------------------------------------------------------------
// returns a pointer to a new material point object
FEMaterialPoint* FEContinuousFiberDistributionUC::CreateMaterialPointData() 
{
	return m_pFmat->CreateMaterialPointData();
}

//-----------------------------------------------------------------------------
//! calculate stress at material point
mat3ds FEContinuousFiberDistributionUC::DevStress(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// calculate stress
	mat3ds s; s.zero();

	// get the local coordinate systems
	mat3d QT = GetLocalCS(mp).transpose();

	double IFD = 0.0;

	// obtain an integration point iterator
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// set the fiber direction
			vec3d& n0 = it->m_fiber;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			vec3d n0a = QT*n0;
			double R = m_pFDD->FiberDensity(mp, n0a);

			// integrate the fiber distribution
			IFD += R*it->m_weight;

			// calculate the stress
			double wn = it->m_weight;
			s += m_pFmat->DevFiberStress(pt, n0)*(R*wn);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// prevent division by zero
	if (IFD == 0.0) IFD = 1.0;

	return s / IFD;
}

//-----------------------------------------------------------------------------
//! calculate tangent stiffness at material point
tens4ds FEContinuousFiberDistributionUC::DevTangent(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the local coordinate systems
	mat3d QT = GetLocalCS(mp).transpose();

	// initialize stress tensor
	tens4ds c;
	c.zero();

	double IFD = 0.0;

	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// set fiber direction in global coordinate system
			vec3d& n0e = it->m_fiber;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			vec3d n0a = QT*n0e;
			double R = m_pFDD->FiberDensity(mp, n0a);

			// integrate the fiber distribution
			IFD += R*it->m_weight;

			// calculate the tangent
			c += m_pFmat->DevFiberTangent(mp, n0e)*(R*it->m_weight);
		}
		while (it->Next());
	}
	
	// don't forget to delete the iterator
	delete it;

	// prevent division by zero
    if (IFD == 0.0) IFD = 1.0;

	// we multiply by two to add contribution from other half-sphere
	return c / IFD;
}

//-----------------------------------------------------------------------------
//! calculate deviatoric strain energy density
double FEContinuousFiberDistributionUC::DevStrainEnergyDensity(FEMaterialPoint& mp)
{ 
	FEElasticMaterialPoint& pt = *mp.ExtractData<FEElasticMaterialPoint>();

	// get the local coordinate systems
	mat3d QT = GetLocalCS(mp).transpose();

	double IFD = 0.0;
	double sed = 0.0;
	FEFiberIntegrationSchemeIterator* it = m_pFint->GetIterator(&pt);
	if (it->IsValid())
	{
		do
		{
			// set fiber direction in global coordinate system
			vec3d& n0e = it->m_fiber;

			// rotate to local configuration to evaluate ellipsoidally distributed material coefficients
			vec3d n0a = QT*n0e;
			double R = m_pFDD->FiberDensity(mp, n0a);

			// integrate the fiber distribution
			IFD += R*it->m_weight;

			// calculate the stress
			sed += m_pFmat->DevFiberStrainEnergyDensity(mp, n0e)*(R*it->m_weight);
		}
		while (it->Next());
	}

	// don't forget to delete the iterator
	delete it;

	// prevent division by zero
	if (IFD == 0.0) IFD = 1.0;

	// we multiply by two to add contribution from other half-sphere
	return sed / IFD;
}
