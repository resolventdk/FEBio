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
#include "FEBackFlowStabilization.h"
#include "FEFluidMaterial.h"
#include "FEBioFluid.h"
#include <FECore/FELinearSystem.h>

//-----------------------------------------------------------------------------
// Parameter block for pressure loads
BEGIN_FECORE_CLASS(FEBackFlowStabilization, FESurfaceLoad)
	ADD_PARAMETER(m_beta, "beta");
END_FECORE_CLASS()

//-----------------------------------------------------------------------------
//! constructor
FEBackFlowStabilization::FEBackFlowStabilization(FEModel* pfem) : FESurfaceLoad(pfem), m_dofW(pfem)
{
    m_beta = 1.0;
    m_rho = 1.0;
    
    // get the degrees of freedom
	m_dofW.AddVariable(FEBioFluid::GetVariableName(FEBioFluid::RELATIVE_FLUID_VELOCITY));

    m_dof.Clear();
    m_dof.AddDofs(m_dofW);
}

//-----------------------------------------------------------------------------
//! initialize
bool FEBackFlowStabilization::Init()
{
    if (FESurfaceLoad::Init() == false) return false;

	// get fluid density from first surface element
    // assuming the entire surface bounds the same fluid
    FESurfaceElement& el = m_psurf->Element(0);
    FEElement* pe = el.m_elem[0];
    if (pe == nullptr) return false;

	// get the material
    FEMaterial* pm = GetFEModel()->GetMaterial(pe->GetMatID());
	FEFluidMaterial* fluid = pm->ExtractProperty<FEFluidMaterial>();
	if (fluid == nullptr) return false;

	// get the density
	m_rho = fluid->ReferentialDensity();
    
    return true;
}

//-----------------------------------------------------------------------------

void FEBackFlowStabilization::Serialize(DumpStream& ar)
{
    FESurfaceLoad::Serialize(ar);
	ar & m_rho;
}

//-----------------------------------------------------------------------------
void FEBackFlowStabilization::StiffnessMatrix(FELinearSystem& LS, const FETimeInfo& tp)
{
	m_psurf->LoadStiffness(LS, m_dofW, m_dofW, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, const FESurfaceDofShape& dof_b, matrix& Kab) {

		FESurfaceElement& el = *mp.SurfaceElement();

		// tangent vectors
		vec3d rt[FEElement::MAX_NODES];
		m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
		vec3d dxr = el.eval_deriv1(rt, mp.m_index);
		vec3d dxs = el.eval_deriv2(rt, mp.m_index);

		vec3d n = dxr ^ dxs;
		double da = n.unit();

		// Fluid velocity
		vec3d v = FluidVelocity(mp, tp.alphaf);

		double vn = v*n;

		if (m_beta*vn < 0) {

			// shape functions and derivatives
			double H_i  = dof_a.shape;

			double H_j  = dof_b.shape;

			mat3d K = dyad(n)*(m_beta*m_rho * 2 * vn*da);

			// calculate stiffness component
			mat3d Kww = K*(H_i * H_j)*tp.alphaf;

			Kab.zero();
			Kab.sub(0, 0, Kww);
		}
	});
}

//-----------------------------------------------------------------------------
vec3d FEBackFlowStabilization::FluidVelocity(FESurfaceMaterialPoint& mp, double alpha)
{
	vec3d vt[FEElement::MAX_NODES];
	FESurfaceElement& el = *mp.SurfaceElement();
	int neln = el.Nodes();
	for (int j = 0; j<neln; ++j) {
		FENode& node = m_psurf->Node(el.m_lnode[j]);
		vt[j] = node.get_vec3d(m_dofW[0], m_dofW[1], m_dofW[2])*alpha + node.get_vec3d_prev(m_dofW[0], m_dofW[1], m_dofW[2])*(1 - alpha);
	}
	return el.eval(vt, mp.m_index);
}

//-----------------------------------------------------------------------------
void FEBackFlowStabilization::LoadVector(FEGlobalVector& R, const FETimeInfo& tp)
{
	m_psurf->LoadVector(R, m_dofW, false, [=](FESurfaceMaterialPoint& mp, const FESurfaceDofShape& dof_a, vector<double>& fa) {

		FESurfaceElement& el = *mp.SurfaceElement();

		// tangent vectors
		vec3d rt[FEElement::MAX_NODES];
		m_psurf->GetNodalCoordinates(el, tp.alphaf, rt);
		vec3d dxr = el.eval_deriv1(rt, mp.m_index);
		vec3d dxs = el.eval_deriv2(rt, mp.m_index);

		// normal and area element
		vec3d n = dxr ^ dxs;
		double da = n.unit();

		// fluid velocity
		vec3d v = FluidVelocity(mp, tp.alphaf);
		double vn = v*n;

		if (m_beta*vn < 0) {

			// force vector (change sign for inflow vs outflow)
			vec3d f = n*(m_beta*m_rho*vn*vn*da);

			double H = dof_a.shape;
			fa[0] = H * f.x;
			fa[1] = H * f.y;
			fa[2] = H * f.z;
		}
	});
}
