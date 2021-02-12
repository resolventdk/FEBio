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



#include "FEFluid.h"
#include "FEFluidMaterialPoint.h"
#include <FECore/FECoreKernel.h>
#include <FECore/DumpStream.h>

// define the material parameters
BEGIN_FECORE_CLASS(FEFluid, FEFluidMaterial)

	// material parameters
    ADD_PARAMETER(m_k   , FE_RANGE_GREATER_OR_EQUAL(0.0), "k");

END_FECORE_CLASS();

//============================================================================
// FEFluid
//============================================================================

//-----------------------------------------------------------------------------
//! FEFluid constructor

FEFluid::FEFluid(FEModel* pfem) : FEFluidMaterial(pfem)
{ 
    m_k = 0;
}

//-----------------------------------------------------------------------------
//! returns a pointer to a new material point object
FEMaterialPoint* FEFluid::CreateMaterialPointData()
{
	return new FEFluidMaterialPoint();
}

//-----------------------------------------------------------------------------
//! bulk modulus
double FEFluid::BulkModulus(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& vt = *mp.ExtractData<FEFluidMaterialPoint>();
    return -vt.m_Jf*Tangent_Pressure_Strain(mp);
}

//-----------------------------------------------------------------------------
//! elastic pressure
double FEFluid::Pressure(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double e = fp.m_Jf - 1;

    return Pressure(e);
}

//-----------------------------------------------------------------------------
//! elastic pressure from dilatation
double FEFluid::Pressure(const double e, const double T)
{
    return -m_k*e;
}

//-----------------------------------------------------------------------------
//! The stress of a fluid material is the sum of the fluid pressure
//! and the viscous stress.

mat3ds FEFluid::Stress(FEMaterialPoint& mp)
{
	// calculate solid material stress
	mat3ds s = GetViscous()->Stress(mp);
    
    double p = Pressure(mp);
	
	// add fluid pressure
	s.xx() -= p;
	s.yy() -= p;
	s.zz() -= p;
	
	return s;
}

//-----------------------------------------------------------------------------
//! The tangent of stress with respect to strain J of a fluid material is the
//! sum of the tangent of the fluid pressure and that of the viscous stress.

mat3ds FEFluid::Tangent_Strain(FEMaterialPoint& mp)
{
    // get tangent of viscous stress
    mat3ds sJ = GetViscous()->Tangent_Strain(mp);
    
    // add tangent of fluid pressure
    double dp = Tangent_Pressure_Strain(mp);
    sJ.xx() -= dp;
    sJ.yy() -= dp;
    sJ.zz() -= dp;
    
    return sJ;
}

//-----------------------------------------------------------------------------
//! calculate strain energy density (per reference volume)
double FEFluid::StrainEnergyDensity(FEMaterialPoint& mp)
{
    FEFluidMaterialPoint& fp = *mp.ExtractData<FEFluidMaterialPoint>();
    double sed = m_k*pow(fp.m_Jf-1,2)/2;
    return sed;
}

//-----------------------------------------------------------------------------
//! invert pressure-dilatation relation
bool FEFluid::Dilatation(const double T, const double p, const double c, double& e)
{
    e = -p/m_k;
    //for solute cases
    if (c != 0)
    {
        e = -(p-T*c)/m_k;
    }
    
    return true;
}

