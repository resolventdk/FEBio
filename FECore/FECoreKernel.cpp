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
#include "FECoreKernel.h"
#include "LinearSolver.h"
#include "Timer.h"
#include <stdarg.h>
using namespace std;

// A module defines a context in which features are defined and searched. 
class FECoreKernel::Module
{
public:
	const char*		m_szname;	// name of module
	unsigned int	m_id;		// unqiue ID (starting at one)
	int				m_alloc_id;	// ID of allocator
	vector<int>		m_depMods;	// module dependencies

	void AddDependency(Module& mod)
	{
		AddDependency(mod.m_id);
		AddDependencies(mod.m_depMods);
	}

private:
	void AddDependency(int mid)
	{
		for (size_t i = 0; i < m_depMods.size(); ++i)
		{
			if (m_depMods[i] == mid) return;
		}
		m_depMods.push_back(mid);
	}

	void AddDependencies(const vector<int>& mid)
	{
		for (size_t i = 0; i < mid.size(); ++i)
		{
			AddDependency(mid[i]);
		}
	}
};

//-----------------------------------------------------------------------------
FECoreKernel* FECoreKernel::m_pKernel = 0;

//-----------------------------------------------------------------------------
FECoreKernel& FECoreKernel::GetInstance()
{
	if (m_pKernel == 0) m_pKernel = new FECoreKernel;
	return *m_pKernel;
}

//-----------------------------------------------------------------------------
// This function is used by plugins to make sure that the plugin and the executable
// are using the same kernel
void FECoreKernel::SetInstance(FECoreKernel* pkernel)
{
	m_pKernel = pkernel;
}

//-----------------------------------------------------------------------------
FECoreKernel::FECoreKernel()
{
	m_activeModule = -1;
	m_alloc_id = 0;
	m_next_alloc_id = 1;

	m_default_solver = nullptr;
}

//-----------------------------------------------------------------------------
// Generate a allocator ID
int FECoreKernel::GenerateAllocatorID()
{
	return m_next_alloc_id++;
}

//-----------------------------------------------------------------------------
FECoreFactory* FECoreKernel::SetDefaultSolverType(const char* sztype)
{
	FECoreFactory* fac = FindFactoryClass(FELINEARSOLVER_ID, sztype);
	if (fac) m_default_solver_type = sztype;
	return fac;
}

//-----------------------------------------------------------------------------
void FECoreKernel::SetDefaultSolver(ClassDescriptor* linsolve)
{
	delete m_default_solver;
	m_default_solver = linsolve;

	if (linsolve)
	{
		m_default_solver_type = linsolve->ClassType();
	}
	else
	{
		m_default_solver_type.clear();
	}
}

//-----------------------------------------------------------------------------
//! get the linear solver type
const char* FECoreKernel::GetLinearSolverType() const
{
	return m_default_solver_type.c_str();
}

//-----------------------------------------------------------------------------
LinearSolver* FECoreKernel::CreateDefaultLinearSolver(FEModel* fem)
{
	if (m_default_solver == nullptr)
	{
		const char* sztype = m_default_solver_type.c_str();
		FECoreFactory* fac = FindFactoryClass(FELINEARSOLVER_ID, sztype);
		return (LinearSolver*)fac->Create(fem);
	}
	else
	{
		return (LinearSolver*)Create(FELINEARSOLVER_ID, fem, *m_default_solver);
	}
}

//-----------------------------------------------------------------------------
void FECoreKernel::RegisterFactory(FECoreFactory* ptf)
{
	unsigned int activeID = 0;
	vector<int> moduleDepends;
	if (m_activeModule != -1)
	{
		Module& activeModule = *m_modules[m_activeModule];
		activeID = activeModule.m_id;
		moduleDepends = activeModule.m_depMods;
	}

	// see if the name already exists
	for (int i=0; i<m_Fac.size(); ++i)
	{
		FECoreFactory* pfi = m_Fac[i];

		if ((pfi->GetSuperClassID() == ptf->GetSuperClassID()) && 
			(strcmp(pfi->GetTypeStr(), ptf->GetTypeStr()) == 0))
		{
			// A feature with the same is already registered. 
			// We need to check the module to see if this would create an ambiguity
			unsigned int modId = pfi->GetModuleID();

			// If the same feature is defined in the active module,
			// then this feature will replace the existing one. 
			if ((modId == activeID) && (pfi->GetSpecID() == ptf->GetSpecID()))
			{
#ifdef _DEBUG
				fprintf(stderr, "WARNING: \"%s\" feature is redefined\n", ptf->GetTypeStr());
#endif
				m_Fac[i] = ptf;
				return;
			}
		}
	}

	// it doesn't so add it
	ptf->SetModuleID(activeID);
	ptf->SetAllocatorID(m_alloc_id);
	m_Fac.push_back(ptf);
}

//-----------------------------------------------------------------------------
bool FECoreKernel::UnregisterFactory(FECoreFactory* ptf)
{
	for (vector<FECoreFactory*>::iterator it = m_Fac.begin(); it != m_Fac.end(); ++it)
	{
		FECoreFactory* pfi = *it;
		if (pfi == ptf)
		{
			m_Fac.erase(it);
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------
//! unregister factories from allocator
void FECoreKernel::UnregisterFactories(int alloc_id)
{
	for (vector<FECoreFactory*>::iterator it = m_Fac.begin(); it != m_Fac.end();)
	{
		FECoreFactory* pfi = *it;
		if (pfi->GetAllocatorID() == alloc_id)
		{
			it = m_Fac.erase(it);
		}
		else ++it;
	}
}

//-----------------------------------------------------------------------------
//! set the current allocator ID
void FECoreKernel::SetAllocatorID(int alloc_id)
{
	m_alloc_id = alloc_id;
}

//-----------------------------------------------------------------------------
//! Create an object. An object is created by specifying the super-class id
//! and the type-string. 
void* FECoreKernel::Create(int superClassID, const char* sztype, FEModel* pfem)
{
	if (sztype == 0) return 0;

	unsigned int activeID = 0;
	vector<int> moduleDepends;
	if (m_activeModule != -1)
	{
		Module& activeModule = *m_modules[m_activeModule];
		activeID = activeModule.m_id;
		moduleDepends = activeModule.m_depMods;
	}

	// first check active module
	if ((activeID > 0) || (activeID == 0))
	{
		std::vector<FECoreFactory*>::iterator pf;
		for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
		{
			FECoreFactory* pfac = *pf;
			if (pfac->GetSuperClassID() == superClassID) {

				// see if we can match module first
				unsigned int mid = pfac->GetModuleID();
				if ((mid == activeID) || (mid== 0))
				{
					// see if the type name matches
					if ((strcmp(pfac->GetTypeStr(), sztype) == 0))
					{
						// check the spec (TODO: What is this for?)
						int nspec = pfac->GetSpecID();
						if ((nspec == -1) || (m_nspec <= nspec))
						{
							return pfac->CreateInstance(pfem);
						}
					}
				}
			}
		}
	}

	// check dependencies in order in which they are defined
	std::vector<FECoreFactory*>::iterator pf;
	for (int i = 0; i < moduleDepends.size(); ++i)
	{
		unsigned modId = moduleDepends[i];
		for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
		{
			FECoreFactory* pfac = *pf;
			if (pfac->GetSuperClassID() == superClassID) {

				// see if we can match module first
				unsigned int mid = pfac->GetModuleID();
				if ((mid == 0) || (mid == modId))
				{
					// see if the type name matches
					if ((strcmp(pfac->GetTypeStr(), sztype) == 0))
					{
						// check the spec (TODO: What is this for?)
						int nspec = pfac->GetSpecID();
						if ((nspec == -1) || (m_nspec <= nspec))
						{
							return pfac->CreateInstance(pfem);
						}
					}
				}
			}
		}
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! Create a specific class
void* FECoreKernel::CreateClass(const char* szclassName, FEModel* fem)
{
	std::vector<FECoreFactory*>::iterator pf;
	for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		const char* szfacName = pfac->GetClassName();
		if (szfacName && (strcmp(szfacName, szclassName) == 0))
		{
			return pfac->CreateInstance(fem);
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
//! Create a class from a class descriptor
void* FECoreKernel::Create(int superClassID, FEModel* pfem, const ClassDescriptor& cd)
{
	const ClassDescriptor::ClassVariable* root = cd.Root();
	FECoreBase* pc = (FECoreBase*)Create(superClassID, root->m_type.c_str(), pfem);
	if (pc == nullptr) return nullptr;
	pc->SetParameters(cd);
	return pc;
}

//-----------------------------------------------------------------------------
int FECoreKernel::Count(SUPER_CLASS_ID sid)
{
	int N = 0;
	std::vector<FECoreFactory*>::iterator pf;
	for (pf=m_Fac.begin(); pf!= m_Fac.end(); ++pf)
	{
		FECoreFactory* pfac = *pf;
		if (pfac->GetSuperClassID() == sid) N++;
	}
	return N;
}

//-----------------------------------------------------------------------------
void FECoreKernel::List(SUPER_CLASS_ID sid)
{
  std::vector<FECoreFactory*>::iterator pf;
  for (pf = m_Fac.begin(); pf != m_Fac.end(); ++pf)
    {
      FECoreFactory* pfac = *pf;
      if (pfac->GetSuperClassID() == sid) fprintf(stdout, "%s\n", pfac->GetTypeStr());
    }
}

//-----------------------------------------------------------------------------
int FECoreKernel::FactoryClasses()
{
	return (int) m_Fac.size();
}

//-----------------------------------------------------------------------------
const FECoreFactory* FECoreKernel::GetFactoryClass(int i)
{
	return m_Fac[i];
}

//-----------------------------------------------------------------------------
//! return a factory class
const FECoreFactory* FECoreKernel::GetFactoryClass(int classID, int i)
{
	int n = 0;
	for (int j = 0; j < m_Fac.size(); ++j)
	{
		FECoreFactory* fac = m_Fac[j];
		if (fac->GetSuperClassID() == classID)
		{
			if (i == n) return fac;
			n++;
		}
	}
	return nullptr;
}

//-----------------------------------------------------------------------------
FECoreFactory* FECoreKernel::FindFactoryClass(int classID, const char* sztype)
{
	for (size_t i=0; i<m_Fac.size(); ++i)
	{
		FECoreFactory* fac = m_Fac[i];
		if ((fac->GetSuperClassID() == classID) &&
			(strcmp(fac->GetTypeStr(), sztype) == 0)) return fac;
	}
	return 0;
}

//-----------------------------------------------------------------------------
//! set the active module
bool FECoreKernel::SetActiveModule(const char* szmod)
{
	// See if user want to deactivate modules
	if (szmod == 0)
	{
		m_activeModule = -1;
		return true;
	}

	// see if the module exists or not
	for (size_t i=0; i<m_modules.size(); ++i) 
	{
		Module& mi = *m_modules[i];
		if (strcmp(mi.m_szname, szmod) == 0)
		{
			m_activeModule = (int) i;
			return true;
		}
	}

	// couldn't find it
	m_activeModule = -1;
	return false;
}

//-----------------------------------------------------------------------------
//! count modules
int FECoreKernel::Modules() const
{
	return (int)m_modules.size();
}

//-----------------------------------------------------------------------------
//! create a module
bool FECoreKernel::CreateModule(const char* szmod)
{
	m_activeModule = -1;
	if (szmod == 0) return false;

	// see if this module already exist
	if (SetActiveModule(szmod) == false)
	{
		// The module does not exist, so let's add it.
		unsigned int newID = (unsigned int) m_modules.size() + 1;
		Module* newModule = new Module;
		newModule->m_szname = szmod;
		newModule->m_id = newID;
		m_modules.push_back(newModule);

		// make this the active module
		m_activeModule = (int)m_modules.size() - 1;
	}

	return true;
}

//-----------------------------------------------------------------------------
//! Get a module
const char* FECoreKernel::GetModuleName(int i) const
{
	if ((i<0) || (i >= m_modules.size())) return nullptr;
	return m_modules[i]->m_szname;
}

//! Get a module
const char* FECoreKernel::GetModuleNameFromId(int id) const
{
	for (size_t n = 0; n < m_modules.size(); ++n)
	{
		const Module& mod = *m_modules[n];
		if (mod.m_id == id) return mod.m_szname;
	}
	return 0;
}

//! Get a module's dependencies
vector<int> FECoreKernel::GetModuleDependencies(int i) const
{
	vector<int> md;
	if ((i >= 0) && (i < m_modules.size()))
	{
		md = m_modules[i]->m_depMods;
	}
	return md;
}


//-----------------------------------------------------------------------------
//! set the spec ID. Features with a matching spec ID will be preferred
//! set spec ID to -1 to stop caring
void FECoreKernel::SetSpecID(int nspec)
{
	m_nspec = nspec;
}

//-----------------------------------------------------------------------------
//! set a dependency on a module
bool FECoreKernel::SetModuleDependency(const char* szmodule)
{
	if (m_activeModule == -1) return false;
	Module& activeModule = *m_modules[m_activeModule];

	if (szmodule == 0)
	{
		// clear dependencies
		activeModule.m_depMods.clear();
		return true;
	}

	// find the module
	for (size_t i = 0; i<m_modules.size(); ++i)
	{
		Module& mi = *m_modules[i];
		if (strcmp(mi.m_szname, szmodule) == 0)
		{
			// add the module to the active module's dependency list
			activeModule.AddDependency(mi);

			return true;
		}
	}

	// oh, oh, couldn't find it
	return false;
}

//-----------------------------------------------------------------------------
//! Register a new domain class
void FECoreKernel::RegisterDomain(FEDomainFactory* pf)
{
	m_Dom.push_back(pf); 
}

//-----------------------------------------------------------------------------
FEDomain* FECoreKernel::CreateDomain(const FE_Element_Spec& spec, FEMesh* pm, FEMaterial* pmat)
{
	for (int i=0; i<(int)m_Dom.size(); ++i)
	{
		FEDomain* pdom = m_Dom[i]->CreateDomain(spec, pm, pmat);
		if (pdom != 0) return pdom;
	}
	return 0;
}
