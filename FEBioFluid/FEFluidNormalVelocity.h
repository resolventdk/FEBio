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



#pragma once
#include <FECore/FESurfaceLoad.h>
#include <FECore/FESurfaceMap.h>
#include "febiofluid_api.h"

//-----------------------------------------------------------------------------
//! FEFluidNormalVelocity is a fluid surface that has a normal
//! velocity prescribed on it.  This routine simultaneously prescribes a
//! surface load and nodal prescribed velocities
class FEBIOFLUID_API FEFluidNormalVelocity : public FESurfaceLoad
{
public:
    //! constructor
    FEFluidNormalVelocity(FEModel* pfem);
    
    //! Set the surface to apply the load to
    void SetSurface(FESurface* ps) override;
    
    //! calculate traction stiffness (there is none)
    void StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp) override {}
    
    //! calculate load vector
    void LoadVector(FEGlobalVector& R, const FETimeInfo& tp) override;
    
	//! serializatsion
	void Serialize(DumpStream& ar) override;
    
    //! set the velocity
    void Update() override;
    
    //! initialization
    bool Init() override;
    
    //! activate
    void Activate() override;
    
    //! parabolic velocity profile
    bool SetParabolicVelocity();
    
    //! rim pressure
    bool SetRimPressure();

private:
	double NormalVelocity(FESurfaceMaterialPoint& mp);
    
private:
    double			m_velocity;	//!< average velocity
    FESurfaceMap	m_VC;		//!< velocity boundary cards
	bool            m_bpv;      //!< flag for prescribing nodal values
	bool            m_bpar;     //!< flag for parabolic velocity
    bool            m_brim;     //!< flag for setting rim pressure

private:
	vector<double>  m_VN;       //!< nodal scale factors
    vector<vec3d>   m_nu;       //!< nodal normals
    vector<int>     m_rim;      //!< list of nodes on the rim

	FEDofList	m_dofW;
    int         m_dofEF;
    
    DECLARE_FECORE_CLASS();
};
