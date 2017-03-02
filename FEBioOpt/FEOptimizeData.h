#pragma once
#include "FEBioXML/XMLReader.h"
#include "FECore/FEModel.h"
#include "FECore/Logfile.h"
#include "FECore/FECoreTask.h"
#include "FEObjectiveFunction.h"
#include <vector>
#include <string>
using namespace std;

//-----------------------------------------------------------------------------
class FEOptimizeMethod;


//-----------------------------------------------------------------------------
//! This class represents an input parameter. Input parameters define the parameter 
//! space that will be searched for the optimal parameter value.
class FEInputParameter
{
public:
	FEInputParameter(FEModel* fem) : m_fem(fem) { m_min = -1e99; m_max = 1e99; m_scale = 1.0; }
	virtual ~FEInputParameter() {}

	// implement this to initialize the input parameter
	virtual bool Init() { return true; }

	//! return the current value of the input parameter
	virtual double GetValue() = 0;

	//! set the current value of the input parameter.
	//! should return false is the passed value is invalid
	virtual bool SetValue(double newValue) = 0;

	//! get/set the initial value
	double& InitValue() { return m_initVal; }

	//! get/set the min value
	double& MinValue() { return m_min; }

	//! get/set the max value
	double& MaxValue() { return m_max; }

	//! get/set scale factor
	double& ScaleFactor() { return m_scale; }

	//! set the name
	void SetName(const string& name) { m_name = name; }

	//! get the name
	string GetName() { return m_name; }

	//! Get the FEModel
	FEModel* GetFEModel() { return m_fem; }

private:
	string		m_name;			//!< name of input parameter
	double		m_initVal;		//!< initial value
	double		m_min, m_max;	//!< min, max values for parameter
	double		m_scale;		//!< scale factor
	FEModel*	m_fem;			//!< pointer to model data
};

//-----------------------------------------------------------------------------
// Class for choosing model parameters as input parameters
class FEModelParameter : public FEInputParameter
{
public:
	FEModelParameter(FEModel* fem);

	//! Initialize
	bool Init();

public:
	//! return the current value of the input parameter
	double GetValue();

	//! set the current value of the input parameter.
	//! should return false is the passed value is invalid
	bool SetValue(double newValue);

private:
	string	m_name;		//!< variable name
	double*	m_pd;		//!< pointer to variable data
	double	m_val;		//!< value
};

//=============================================================================
#define OPT_MAX_VAR 64
struct OPT_LIN_CONSTRAINT
{
	double	a[OPT_MAX_VAR];
	double	b;
};

//=============================================================================
//! optimization analyses
//! 
class FEOptimizeData
{
public:
	//! constructor
	FEOptimizeData(FEModel& fem);
	~FEOptimizeData(void);

	//! input function
	bool Input(const char* sz);

	//! Initialize data
	bool Init();

	//! solver the problem
	bool Solve();

	//! return the FE Model
	FEModel& GetFEM() { return m_fem; }

	//! solve the FE problem with a new set of parameters
	bool FESolve(const vector<double>& a);

public:
	// return the number of input parameters
	int InputParameters() { return (int)m_Var.size(); }

	//! add an input parameter to optimize
	void AddInputParameter(FEInputParameter* var) { m_Var.push_back(var); }

	//! return an input parameter
	FEInputParameter* GetInputParameter(int n) { return m_Var[n]; }

public:
	//! add a linear constraint
	void AddLinearConstraint(OPT_LIN_CONSTRAINT& con) { m_LinCon.push_back(con); }

	//! return number of constraints
	int Constraints() { return (int) m_LinCon.size(); }

	//! return a linear constraint
	OPT_LIN_CONSTRAINT& Constraint(int i) { return m_LinCon[i]; }

	FEObjectiveFunction& GetObjective() { return *m_obj; }

	void SetObjective(FEObjectiveFunction* obj) { m_obj = obj; }

	void SetSolver(FEOptimizeMethod* po) { m_pSolver = po; }

	bool RunTask();

public:
	int	m_niter;	// nr of minor iterations (i.e. FE solves)

	FECoreTask* m_pTask;	// the task that will solve the FE model

protected:
	FEModel&	m_fem;

	FEObjectiveFunction*	m_obj;		//!< the objective function

	FEOptimizeMethod*	m_pSolver;

	std::vector<FEInputParameter*>	    m_Var;
	std::vector<OPT_LIN_CONSTRAINT>		m_LinCon;
};
