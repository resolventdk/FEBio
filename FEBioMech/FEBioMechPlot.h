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
#include <FECore/FEPlotData.h>
#include <FECore/FEElement.h>

//=============================================================================
//                            N O D E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Nodal velocities
//!
class FEPlotNodeVelocity : public FEPlotNodeData
{
public:
	FEPlotNodeVelocity(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal accelerations
//!
class FEPlotNodeAcceleration : public FEPlotNodeData
{
public:
	FEPlotNodeAcceleration(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal reaction forces
class FEPlotNodeReactionForces : public FEPlotNodeData
{
public:
	FEPlotNodeReactionForces(FEModel* pfem) : FEPlotNodeData(pfem, PLT_VEC3F, FMT_NODE){}
	bool Save(FEMesh& m, FEDataStream& a);
};

//=============================================================================
//                         S U R F A C E   D A T A
//=============================================================================

//-----------------------------------------------------------------------------
//! Contact gap
//!
class FEPlotContactGap : public FEPlotSurfaceData
{
public:
    FEPlotContactGap(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Vector gap
//!
class FEPlotVectorGap : public FEPlotSurfaceData
{
public:
    FEPlotVectorGap(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact pressure
//!
class FEPlotContactPressure : public FEPlotSurfaceData
{
public:
    FEPlotContactPressure(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact traction
//!
class FEPlotContactTraction : public FEPlotSurfaceData
{
public:
    FEPlotContactTraction(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal contact gap
//!
class FEPlotNodalContactGap : public FEPlotSurfaceData
{
public:
	FEPlotNodalContactGap(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal vector gap
//!
class FEPlotNodalVectorGap : public FEPlotSurfaceData
{
public:
    FEPlotNodalVectorGap(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_MULT){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal contact pressure
//!
class FEPlotNodalContactPressure : public FEPlotSurfaceData
{
public:
    FEPlotNodalContactPressure(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_MULT){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal contact traction
//!
class FEPlotNodalContactTraction : public FEPlotSurfaceData
{
public:
	FEPlotNodalContactTraction(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_MULT){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Surface traction
//!
class FEPlotSurfaceTraction : public FEPlotSurfaceData
{
public:
    FEPlotSurfaceTraction(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Nodal surface traction
//!
class FEPlotNodalSurfaceTraction : public FEPlotSurfaceData
{
public:
    FEPlotNodalSurfaceTraction(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_MULT){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Stick status
//!
class FEPlotStickStatus : public FEPlotSurfaceData
{
public:
    FEPlotStickStatus(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact force
//!
class FEPlotContactForce : public FEPlotSurfaceData
{
public:
	FEPlotContactForce(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_VEC3F, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact area
//!
class FEPlotContactArea : public FEPlotSurfaceData
{
public:
	FEPlotContactArea(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact penalty parameter
class FEPlotContactPenalty : public FEPlotSurfaceData
{
public:
	FEPlotContactPenalty(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Contact status 
class FEPlotContactStatus : public FEPlotSurfaceData
{
public:
	FEPlotContactStatus(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mortar gap
class FEPlotMortarContactGap : public FEPlotSurfaceData
{
public:
	FEPlotMortarContactGap(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_NODE){}
	bool Save(FESurface& S, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Enclosed volume
//!
class FEPlotEnclosedVolume : public FEPlotSurfaceData
{
private:
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotEnclosedVolume(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){ m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Surface area
//!
class FEPlotSurfaceArea : public FEPlotSurfaceData
{
private:
    FEModel*            m_pfem;
    bool                m_binit;
    vector<FEElement*>  m_elem;
    vector<vec3d>       m_area;
    
public:
    FEPlotSurfaceArea(FEModel* pfem) : FEPlotSurfaceData(pfem, PLT_FLOAT, FMT_REGION){ m_binit = true; }
    bool Save(FESurface& surf, FEDataStream& a);
};

//=============================================================================
//							D O M A I N   D A T A
//=============================================================================
//-----------------------------------------------------------------------------
//! Velocity
class FEPlotElementVelocity : public FEPlotDomainData
{
public:
    FEPlotElementVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Acceleration
class FEPlotElementAcceleration : public FEPlotDomainData
{
public:
    FEPlotElementAcceleration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for G
class FEPlotElementGnorm : public FEPlotDomainData
{
public:
	FEPlotElementGnorm(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element stresses
class FEPlotElementStress : public FEPlotDomainData
{
public:
	FEPlotElementStress(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element uncoupled pressure
class FEPlotElementUncoupledPressure : public FEPlotDomainData
{
public:
    FEPlotElementUncoupledPressure(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for Cauchy stress
class FEPlotElementsnorm : public FEPlotDomainData
{
public:
	FEPlotElementsnorm(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK1 stress
class FEPlotElementPK1norm : public FEPlotDomainData
{
public:
	FEPlotElementPK1norm(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element norm for PK1 stress moment
class FEPlotElementQK1norm : public FEPlotDomainData
{
public:
	FEPlotElementQK1norm(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element micro energy
class FEPlotElementMicroEnergy : public FEPlotDomainData
{
public:
	FEPlotElementMicroEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy density
class FEPlotStrainEnergyDensity : public FEPlotDomainData
{
public:
	FEPlotStrainEnergyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Deviatoric strain energy density
class FEPlotDevStrainEnergyDensity : public FEPlotDomainData
{
public:
	FEPlotDevStrainEnergyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Specific strain energy
class FEPlotSpecificStrainEnergy : public FEPlotDomainData
{
public:
	FEPlotSpecificStrainEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy density
class FEPlotKineticEnergyDensity : public FEPlotDomainData
{
public:
    FEPlotKineticEnergyDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mass density
class FEPlotDensity : public FEPlotDomainData
{
public:
	FEPlotDensity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy
class FEPlotElementStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotElementStrainEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy
class FEPlotElementKineticEnergy : public FEPlotDomainData
{
public:
    FEPlotElementKineticEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass
class FEPlotElementCenterOfMass : public FEPlotDomainData
{
public:
    FEPlotElementCenterOfMass(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum
class FEPlotElementLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotElementLinearMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum
class FEPlotElementAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotElementAngularMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Stress power
class FEPlotElementStressPower : public FEPlotDomainData
{
public:
    FEPlotElementStressPower(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Strain energy at current time
class FEPlotCurrentElementStrainEnergy : public FEPlotDomainData
{
public:
    FEPlotCurrentElementStrainEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Kinetic energy at current time
class FEPlotCurrentElementKineticEnergy : public FEPlotDomainData
{
public:
    FEPlotCurrentElementKineticEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Center of mass at current time
class FEPlotCurrentElementCenterOfMass : public FEPlotDomainData
{
public:
    FEPlotCurrentElementCenterOfMass(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Linear momentum at current time
class FEPlotCurrentElementLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotCurrentElementLinearMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Angular momentum at current time
class FEPlotCurrentElementAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotCurrentElementAngularMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Relative volume
class FEPlotRelativeVolume : public FEPlotDomainData
{
public:
	FEPlotRelativeVolume(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Material fibers
class FEPlotFiberVector : public FEPlotDomainData
{
public:
	FEPlotFiberVector(FEModel* pfem);
	bool SetFilter(const char* szfilter) override;
	bool Save(FEDomain& dom, FEDataStream& a) override;
private:
	std::string		m_matComp;
};

//-----------------------------------------------------------------------------
//! Material axes
class FEPlotMaterialAxes : public FEPlotDomainData
{
public:
	FEPlotMaterialAxes(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3F, FMT_ITEM){}
	bool SetFilter(const char* szfilter) override;
	bool Save(FEDomain& dom, FEDataStream& a) override;
private:
	std::string		m_matComp;
};

//-----------------------------------------------------------------------------
//! fiber stretch
class FEPlotFiberStretch : public FEPlotDomainData
{
public:
	FEPlotFiberStretch(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool SetFilter(const char* szfilter) override;
	bool Save(FEDomain& dom, FEDataStream& a) override;
private:
	std::string		m_matComp;
};

//-----------------------------------------------------------------------------
//! deviatoric fiber stretch
class FEPlotDevFiberStretch : public FEPlotDomainData
{
public:
	FEPlotDevFiberStretch(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool SetFilter(const char* szfilter) override;
	bool Save(FEDomain& dom, FEDataStream& a) override;
private:
	std::string		m_matComp;
};

//-----------------------------------------------------------------------------
//! Shell thicknesses
class FEPlotShellThickness : public FEPlotDomainData
{
public:
	FEPlotShellThickness(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Shell directors
class FEPlotShellDirector : public FEPlotDomainData
{
public:
	FEPlotShellDirector(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Element elasticity tensor
class FEPlotElementElasticity : public FEPlotDomainData
{
public:
	FEPlotElementElasticity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_TENS4FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotDamage : public FEPlotDomainData
{
public:
	FEPlotDamage(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Damage reduction factor
class FEPlotNestedDamage : public FEPlotDomainData
{
public:
    FEPlotNestedDamage(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& m, FEDataStream& a);
    bool SetFilter(int nsol);
protected:
    int			m_nmat;
};

//-----------------------------------------------------------------------------
//! Intact bond fraction (fatigue)
class FEPlotIntactBondFraction : public FEPlotDomainData
{
public:
    FEPlotIntactBondFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Yielded bond fraction (fatigue)
class FEPlotYieldedBondFraction : public FEPlotDomainData
{
public:
    FEPlotYieldedBondFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Octahedral Plastic Strain
class FEPlotOctahedralPlasticStrain : public FEPlotDomainData
{
public:
    FEPlotOctahedralPlasticStrain(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Reactive plasticity heat supply
class FEPlotReactivePlasticityHeatSupply : public FEPlotDomainData
{
public:
    FEPlotReactivePlasticityHeatSupply(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Mixture volume fraction
class FEPlotMixtureVolumeFraction : public FEPlotDomainData
{
public:
	FEPlotMixtureVolumeFraction(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
	bool Save(FEDomain& m, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the element nodal stresses for UT4 domains
class FEPlotUT4NodalStresses : public FEPlotDomainData
{
public:
	FEPlotUT4NodalStresses(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_NODE) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Store shell strains
class FEPlotShellStrain : public FEPlotDomainData
{
public:
	FEPlotShellStrain(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Shell relative volume
class FEPlotShellRelativeVolume : public FEPlotDomainData
{
public:
    FEPlotShellRelativeVolume(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_ITEM){}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRStresses : public FEPlotDomainData
{
public:
	FEPlotSPRStresses(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class the projects stresses from integration points to nodes using
//! SPR (superconvergergent patch recovery)
class FEPlotSPRLinearStresses : public FEPlotDomainData
{
public:
	FEPlotSPRLinearStresses(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! class that projects the principal stresses from integration points to nodes
//! using the SPR projection method
class FEPlotSPRPrincStresses : public FEPlotDomainData
{
public:
	FEPlotSPRPrincStresses(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FD, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body displacement
class FEPlotRigidDisplacement : public FEPlotDomainData
{
public:
	FEPlotRigidDisplacement(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body velocity
class FEPlotRigidVelocity : public FEPlotDomainData
{
public:
	FEPlotRigidVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body acceleration
class FEPlotRigidAcceleration : public FEPlotDomainData
{
public:
	FEPlotRigidAcceleration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body rotation
class FEPlotRigidRotation : public FEPlotDomainData
{
public:
	FEPlotRigidRotation(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body angular velocity
class FEPlotRigidAngularVelocity : public FEPlotDomainData
{
public:
	FEPlotRigidAngularVelocity(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body angular acceleration
class FEPlotRigidAngularAcceleration : public FEPlotDomainData
{
public:
	FEPlotRigidAngularAcceleration(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body kinetic energy
class FEPlotRigidKineticEnergy : public FEPlotDomainData
{
public:
	FEPlotRigidKineticEnergy(FEModel* pfem) : FEPlotDomainData(pfem, PLT_FLOAT, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body linear momentum
class FEPlotRigidLinearMomentum : public FEPlotDomainData
{
public:
    FEPlotRigidLinearMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body angular momentum
class FEPlotRigidAngularMomentum : public FEPlotDomainData
{
public:
    FEPlotRigidAngularMomentum(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid Euler angles
class FEPlotRigidEuler : public FEPlotDomainData
{
public:
	FEPlotRigidEuler(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rotation vector
class FEPlotRigidRotationVector : public FEPlotDomainData
{
public:
	FEPlotRigidRotationVector(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that projects stresses from integration points to the nodes
//! TODO: This only works with tet10 and hex8 -domains
class FEPlotNodalStresses : public FEPlotDomainData
{
public:
	FEPlotNodalStresses(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_MULT){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Lagrange strains
class FEPlotLagrangeStrain : public FEPlotDomainData
{
public:
	FEPlotLagrangeStrain(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_ITEM){}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Lagrange strains
class FEPlotSPRLagrangeStrain : public FEPlotDomainData
{
public:
	FEPlotSPRLagrangeStrain(FEModel* pfem) : FEPlotDomainData(pfem, PLT_MAT3FS, FMT_NODE){}
	bool Save(FEDomain& dom, FEDataStream& a);
};


//-----------------------------------------------------------------------------
//! Rigid body reaction force
class FEPlotRigidReactionForce : public FEPlotDomainData
{
public:
	FEPlotRigidReactionForce(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Rigid body reaction torque
class FEPlotRigidReactionTorque : public FEPlotDomainData
{
public:
    FEPlotRigidReactionTorque(FEModel* pfem) : FEPlotDomainData(pfem, PLT_VEC3F, FMT_REGION) {}
    bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotStressError : public FEPlotDomainData
{
public:
	FEPlotStressError(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the in-situ fiber stretch
class FEPlotFiberTargetStretch : public FEPlotDomainData
{
public:
	FEPlotFiberTargetStretch(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the fiber stretch
class FEPlotPreStrainStretch : public FEPlotDomainData
{
public:
	FEPlotPreStrainStretch(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the error in the fiber stretch
class FEPlotPreStrainStretchError : public FEPlotDomainData
{
public:
	FEPlotPreStrainStretchError(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the pre-strain correction deformation gradient
class FEPlotPreStrainCorrection : public FEPlotDomainData
{
public:
	FEPlotPreStrainCorrection(FEModel* fem) : FEPlotDomainData(fem, PLT_MAT3F, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
//! Class that outputs the pre-strain correction deformation gradient
class FEPlotSPRPreStrainCorrection : public FEPlotDomainData
{
public:
	FEPlotSPRPreStrainCorrection(FEModel* fem) : FEPlotDomainData(fem, PLT_MAT3F, FMT_NODE) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotPreStrainCompatibility : public FEPlotDomainData
{
public:
	FEPlotPreStrainCompatibility(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotDiscreteElementStretch : public FEPlotDomainData
{
public:
	FEPlotDiscreteElementStretch(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotDiscreteElementElongation : public FEPlotDomainData
{
public:
	FEPlotDiscreteElementElongation(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotDiscreteElementPercentElongation : public FEPlotDomainData
{
public:
	FEPlotDiscreteElementPercentElongation(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotDiscreteElementForce : public FEPlotDomainData
{
public:
	FEPlotDiscreteElementForce(FEModel* fem) : FEPlotDomainData(fem, PLT_VEC3F, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotDiscreteElementStrainEnergy : public FEPlotDomainData
{
public:
	FEPlotDiscreteElementStrainEnergy(FEModel* fem) : FEPlotDomainData(fem, PLT_FLOAT, FMT_ITEM) {}
	bool Save(FEDomain& dom, FEDataStream& a);
};

//-----------------------------------------------------------------------------
class FEPlotContinuousDamage : public FEPlotDomainData
{
public:
	FEPlotContinuousDamage(FEModel* fem);
	bool Save(FEDomain& dom, FEDataStream& a) override;
	bool SetFilter(const char* sz) override;
private:
	std::string	m_prop;
	int			m_propIndex;
};

//-----------------------------------------------------------------------------
class FEPlotContinuousDamageBeta : public FEPlotDomainData
{
public:
	FEPlotContinuousDamageBeta(FEModel* fem);
	bool Save(FEDomain& dom, FEDataStream& a) override;
	bool SetFilter(const char* sz) override;
private:
	std::string	m_prop;
	int			m_propIndex;
};

//-----------------------------------------------------------------------------
class FEPlotContinuousDamageGamma : public FEPlotDomainData
{
public:
	FEPlotContinuousDamageGamma(FEModel* fem);
	bool Save(FEDomain& dom, FEDataStream& a) override;
	bool SetFilter(const char* sz) override;
private:
	std::string	m_prop;
	int			m_propIndex;
};
