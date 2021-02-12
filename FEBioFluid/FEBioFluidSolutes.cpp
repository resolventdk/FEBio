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
#include "FEBioFluidSolutes.h"
#include <FECore/FECoreKernel.h>
#include "FEFluidSolutesSolver.h"
#include "FEFluidSolutes.h"
#include "FEFluidSolutesDomain3D.h"
#include "FEFluidSolutesDomainFactory.h"
#include "FESoluteBackflowStabilization.h"
#include "FEFluidSolutesNaturalFlux.h"
#include "FEFluidSolutesPressure.h"
#include "FESoluteConvectiveFlow.h"
#include "FEFluidSolutesDomainFactory.h"
#include "FESolutesSolver.h"
#include "FESolutesMaterial.h"
#include "FESolutesDomain.h"
#include "FESolutesDomainFactory.h"
#include <FEBioMix/FESoluteFlux.h>
#include "FEFluidSolutesSolver2.h"
#include "FEFluidSolutesMaterial2.h"
#include "FEFluidSolutesDomain2.h"

//-----------------------------------------------------------------------------
const char* FEBioFluidSolutes::GetVariableName(FEBioFluidSolutes::FLUID_SOLUTES_VARIABLE var)
{
    switch (var)
    {
        case DISPLACEMENT                : return "displacement"               ; break;
        case RELATIVE_FLUID_VELOCITY     : return "relative fluid velocity"    ; break;
        case RELATIVE_FLUID_ACCELERATION : return "relative fluid acceleration"; break;
        case FLUID_DILATATION            : return "fluid dilation"             ; break;
        case FLUID_DILATATION_TDERIV     : return "fluid dilation tderiv"      ; break;
        case FLUID_CONCENTRATION         : return "concentration"              ; break;
        case FLUID_CONCENTRATION_TDERIV  : return "concentration tderiv"       ; break;
    }
    assert(false);
    return nullptr;
}

void FEBioFluidSolutes::InitModule()
{
    FECoreKernel& febio = FECoreKernel::GetInstance();
    
    // register domain
    febio.RegisterDomain(new FEFluidSolutesDomainFactory);
    
    // define the fsi module
    febio.CreateModule("fluid-solutes");
	febio.SetModuleDependency("fluid");
    febio.SetModuleDependency("multiphasic"); // also pulls in solid, biphasic, solutes
    
	// monolithic fluid-solutes solver
    REGISTER_FECORE_CLASS(FEFluidSolutesSolver, "fluid-solutes");
    REGISTER_FECORE_CLASS(FEFluidSolutes, "fluid-solutes");
    REGISTER_FECORE_CLASS(FEFluidSolutesDomain3D, "fluid-solutes-3D");
    
    REGISTER_FECORE_CLASS(FESoluteFlux, "soluteflux");
    REGISTER_FECORE_CLASS(FESoluteBackflowStabilization, "solute backflow stabilization");
    REGISTER_FECORE_CLASS(FEFluidSolutesNaturalFlux, "solute natural flux");
    REGISTER_FECORE_CLASS(FEFluidSolutesPressure, "pressure");
    
    REGISTER_FECORE_CLASS(FESoluteConvectiveFlow, "solute convective flow");

	// solutes solver classes
	febio.RegisterDomain(new FESolutesDomainFactory);
	REGISTER_FECORE_CLASS(FESolutesSolver, "solutes");
	REGISTER_FECORE_CLASS(FESolutesMaterial, "solutes");
	REGISTER_FECORE_CLASS(FESolutesDomain, "solutes-3D");
	febio.SetActiveModule(0);

	febio.CreateModule("fluid-solutes2");
	febio.SetModuleDependency("fluid-solutes");

	// segragated fluid-solutes solver
	REGISTER_FECORE_CLASS(FEFluidSolutesSolver2, "fluid-solutes2");
	REGISTER_FECORE_CLASS(FEFluidSolutesMaterial2, "fluid-solutes2");
	REGISTER_FECORE_CLASS(FEFluidSolutesDomain2, "fluid-solutes2");
    
    febio.SetActiveModule(0);
}
