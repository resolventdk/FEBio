//
//  FEElasticEASShellDomain.hpp
//  FEBioMech
//
//  Created by Gerard Ateshian on 12/6/17.
//  Copyright © 2017 febio.org. All rights reserved.
//

#ifndef FEElasticEASShellDomain_hpp
#define FEElasticEASShellDomain_hpp

#include "FESSIShellDomain.h"
#include "FEElasticDomain.h"
#include "FESolidMaterial.h"

//-----------------------------------------------------------------------------
//! Domain described by 3D shell elements
class FECORE_API FEElasticEASShellDomain : public FESSIShellDomain, public FEElasticDomain
{
public:
    FEElasticEASShellDomain(FEModel* pfem);
    
    //! \todo do I really need this?
    FEElasticEASShellDomain& operator = (FEElasticEASShellDomain& d);
    
    //! Initialize domain
	bool Init() override;
    
    //! Activate the domain
    void Activate() override;
    
    //! Unpack shell element data
    void UnpackLM(FEElement& el, vector<int>& lm) override;
    
    //! get the material (overridden from FEDomain)
    FEMaterial* GetMaterial() override { return m_pMat; }
    
    //! set the material
    void SetMaterial(FEMaterial* pmat) override;
    
public: // overrides from FEElasticDomain
    
    //! calculates the residual
    //    void Residual(FESolver* psolver, vector<double>& R);
    
    //! internal stress forces
    void InternalForces(FEGlobalVector& R) override;
    
    //! Calculates inertial forces for dynamic problems
    void InertialForces(FEGlobalVector& R, vector<double>& F) override;
    
    //! calculate body force
    void BodyForce(FEGlobalVector& R, FEBodyForce& bf) override;
    
    // update stresses
    void Update(const FETimeInfo& tp) override;
    
    //! initialize elements for this domain
    void PreSolveUpdate(const FETimeInfo& timeInfo) override;
    
    //! calculates the global stiffness matrix for this domain
    void StiffnessMatrix(FESolver* psolver) override;
    
    // inertial stiffness
    void MassMatrix(FESolver* psolver, double scale) override;
    
    // body force stiffness
    void BodyForceStiffness  (FESolver* psolver, FEBodyForce& bf) override;
    
    // evaluate strain E and matrix hu and hw
	void EvaluateEh(FEShellElementNew& el, const int n, const vec3d* Gcnt, mat3ds& E, vector<matrix>& hu, vector<matrix>& hw, vector<vec3d>& Nu, vector<vec3d>& Nw);
    
public:
    
    // --- S T I F F N E S S ---
    
    //! calculates the shell element stiffness matrix
    void ElementStiffness(int iel, matrix& ke);
    
    // --- R E S I D U A L ---
    
    //! Calculates the internal stress vector for shell elements
    void ElementInternalForce(FEShellElementNew& el, vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEModel& fem, FEShellElementNew& el, vector<double>& fe);
    
    //! Calculate extenral body forces for shell elements
	void ElementBodyForce(FEBodyForce& BF, FEShellElementNew& el, vector<double>& fe);
    
    //! calculates the solid element mass matrix
	void ElementMassMatrix(FEShellElementNew& el, matrix& ke, double a);
    
    //! calculates the stiffness matrix due to body forces
	void ElementBodyForceStiffness(FEBodyForce& bf, FEShellElementNew& el, matrix& ke);
    
public:
    
    // --- E A S  M E T H O D ---
    
    // Generate the G matrix for the EAS method
	void GenerateGMatrix(FEShellElementNew& el, const int n, const double Jeta, matrix& G);
    
    // Evaluate contravariant components of mat3ds tensor
    void mat3dsCntMat61(const mat3ds s, const vec3d* Gcnt, matrix& S);
    
    // Evaluate contravariant components of tens4ds tensor
    void tens4dsCntMat66(const tens4ds c, const vec3d* Gcnt, matrix& C);
    
    // Evaluate the matrices and vectors relevant to the EAS method
	void EvaluateEAS(FEShellElementNew& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, vector<mat3ds>& S, vector<tens4ds>& c);
    
    // Evaluate the strain using the ANS method
	void CollocationStrainsANS(FEShellElementNew& el, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW, matrix& NS, matrix& NN);
    
	void EvaluateANS(FEShellElementNew& el, const int n, const vec3d* Gcnt, mat3ds& Ec, vector<matrix>& hu, vector<matrix>& hw, vector<double>& E, vector< vector<vec3d>>& HU, vector< vector<vec3d>>& HW);
    
    // Update alpha in EAS method
    void UpdateEAS(vector<double>& ui) override;
    void UpdateIncrementsEAS(vector<double>& ui, const bool binc) override;
    
protected:
    FESolidMaterial*    m_pMat;
    int                 m_nEAS;
};

#endif /* FEElasticEASShellDomain_hpp */
