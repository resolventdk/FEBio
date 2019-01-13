#pragma once
#include <FECore/FESurfaceConstraint.h>
#include "FEContactSurface.h"
#include "FEDeformableSpringDomain.h"

//-----------------------------------------------------------------------------
class FEDiscreteSet;

//-----------------------------------------------------------------------------
class FECORE_API FEDiscreteContactSurface : public FEContactSurface
{
public:
	//! constructor
	FEDiscreteContactSurface(FEModel* fem);

	//! Initialization
	bool Init();
};

//-----------------------------------------------------------------------------
class FECORE_API FEDiscreteContact : public FESurfaceConstraint
{
	struct NODE
	{
		int					nid;	//!< (local) node ID
		FESurfaceElement*	pe;		//!< master surface element
		double				gap;	//!< gap distance
		double				Lm;		//!< Lagrange multiplier
		vec3d				nu;		//!< normal at master projection
		vec3d				q;		//!< projection point
		double				proj[2];	//!< iso-parametric coordinates of projection point
	};

public:
	FEDiscreteContact(FEModel* pfem);

public:
	bool Init() override;
	void Activate() override;

	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
	bool Augment(int naug, const FETimeInfo& tp) override;
	void BuildMatrixProfile(FEGlobalMatrix& M) override;
	void Update(const FETimeInfo& tp);

	void SetDiscreteSet(FEDiscreteSet* pset);

	FESurface* GetSurface() override { return &m_surf; }

protected:
	void ProjectSurface(bool bupseg);
	void ContactNodalForce    (NODE& nodeData, FESurfaceElement& mel, vector<double>& fe);
	void ContactNodalStiffness(NODE& nodeData, FESurfaceElement& mel, matrix& ke);

protected:
	FEDiscreteContactSurface	m_surf;
	vector<NODE>	m_Node;
	double	m_normg0;
	bool	m_bfirst;

protected:
	bool	m_blaugon;	//!< augmentation flag
	double	m_altol;	//!< augmentation tolerance
	double	m_penalty;	//!< penalty parameter
	double	m_gaptol;	//!< gap tolerance
	int		m_naugmin;	//!< minimum number of augmentations
	int		m_naugmax;	//!< maximum number of augmentations
	int		m_nsegup;	//!< number of segment updates (or zero)

	DECLARE_FECORE_CLASS();
};

//-----------------------------------------------------------------------------
class FEDiscreteContact2 : public FESurfaceConstraint
{
	struct NODE
	{
		int		node;				// node index (local ID into discrete domain)
		FESurfaceElement*	pe;		// master element
		double	proj[2];			// natural coordinates of projection
		vec3d	nu;					// normal on master surface
		vec3d	q;					// new position
	};

public:
	FEDiscreteContact2(FEModel* fem);

	bool Init() override;
	void Activate() override;

	void Residual(FEGlobalVector& R, const FETimeInfo& tp) override;
	void StiffnessMatrix(FESolver* psolver, const FETimeInfo& tp) override;
	void BuildMatrixProfile(FEGlobalMatrix& M) override;
	void Update(const FETimeInfo& tp);
	bool Augment(int naug, const FETimeInfo& tp) override { return true; }

	void SetDiscreteDomain(FEDeformableSpringDomain2* dom) { m_dom = dom; }
	FESurface* GetSurface() override { return &m_surf; }

protected:
	void ProjectNodes();

protected:
	FEDiscreteContactSurface	m_surf;
	FEDeformableSpringDomain2*	m_dom;
	vector<NODE>	m_nodeData;
};
