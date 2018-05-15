#include "stdafx.h"
#include "FGMRES_ILUT_Solver.h"

//-----------------------------------------------------------------------------
// We must undef PARDISO since it is defined as a function in mkl_solver.h
#ifdef MKL_ISS
#ifdef PARDISO
#undef PARDISO
#endif
#include "mkl_rci.h"
#include "mkl_blas.h"
#include "mkl_spblas.h"
#endif // MKL_ISS

//=============================================================================
FGMRES_ILUT_Solver::FGMRES_ILUT_Solver() : m_pA(0)
{
	m_maxfill = 1;
	m_fillTol = 1e-16;
	m_maxiter = 0; // use default min(N, 150)
	m_nrestart = 0; // use default = maxiter
	m_print_level = 0;
	m_doResidualTest = true;
	m_tol = 0.0;

	m_checkZeroDiagonal = true;
	m_zeroThreshold = 1e-16;
	m_zeroReplace = 1e-10;

	m_doPreCond = true;
}

//-----------------------------------------------------------------------------
//! Set max nr of iterations
void FGMRES_ILUT_Solver::SetMaxIterations(int n)
{
	m_maxiter = n;
}

//-----------------------------------------------------------------------------
//! Set the nr of non-restarted iterations
void FGMRES_ILUT_Solver::SetNonRestartedIterations(int n)
{
	m_nrestart = n;
}

//-----------------------------------------------------------------------------
// Set the print level
void FGMRES_ILUT_Solver::SetPrintLevel(int n)
{
	m_print_level = n;
}

//-----------------------------------------------------------------------------
// set residual stopping test flag
void FGMRES_ILUT_Solver::DoResidualStoppingTest(bool b)
{
	m_doResidualTest = b;
}

//-----------------------------------------------------------------------------
// set the convergence tolerance for the residual stopping test
void FGMRES_ILUT_Solver::SetResidualTolerance(double tol)
{
	m_tol = tol;
}

//-----------------------------------------------------------------------------
// set the max fill value
void FGMRES_ILUT_Solver::SetMaxFill(int n)
{
	m_maxfill = n;
}

//-----------------------------------------------------------------------------
// Set the fill tolerance
void FGMRES_ILUT_Solver::SetFillTolerance(double fillTol)
{
	m_fillTol = fillTol;
}

//-----------------------------------------------------------------------------
// do the zero diagonal check during preconditioner
void FGMRES_ILUT_Solver::DoZeroDiagonalCheck(bool b)
{
	m_checkZeroDiagonal = b;
}

//-----------------------------------------------------------------------------
// Set the zero diagonal tolerance value
void FGMRES_ILUT_Solver::SetZeroDiagonalTolerance(double tol)
{
	m_zeroThreshold = tol;
}

//-----------------------------------------------------------------------------
// set the zero diagonal replacement value
void FGMRES_ILUT_Solver::SetZeroDiagonalReplacement(double val)
{
	m_zeroReplace = val;
}

//-----------------------------------------------------------------------------
SparseMatrix* FGMRES_ILUT_Solver::CreateSparseMatrix(Matrix_Type ntype)
{
#ifdef MKL_ISS
	if (ntype != REAL_UNSYMMETRIC) return 0;
	m_pA = new CompactUnSymmMatrix(1, true);
	return m_pA;
#else
	return 0;
#endif
}

//-----------------------------------------------------------------------------
void FGMRES_ILUT_Solver::SetSparseMatrix(SparseMatrix* pA)
{
	m_pA = pA;
	assert(m_pA);
	assert(m_pA->Offset() == 1);
	assert(m_pA->isRowBased());
}

//-----------------------------------------------------------------------------
//! Clean up
void FGMRES_ILUT_Solver::Destroy()
{
	m_tmp.clear();
	m_tmp.shrink_to_fit();
}

//-----------------------------------------------------------------------------
bool FGMRES_ILUT_Solver::PreProcess()
{
#ifdef MKL_ISS
	// number of equations
	MKL_INT N = m_pA->Size();

	int M = (N < 150 ? N : 150); // this is the default value of ipar[14]

	if (m_nrestart > 0) M = m_nrestart;
	else if (m_maxiter > 0) M = m_maxiter;

	// allocate temp storage
	m_tmp.resize((N*(2 * M + 1) + (M*(M + 9)) / 2 + 1));

	// set the pre-conditioning flag
	m_doPreCond = true;
	return true;
#else
	return false;
#endif
}

//-----------------------------------------------------------------------------
bool FGMRES_ILUT_Solver::BackSolve(vector<double>& x, vector<double>& b)
{
#ifdef MKL_ISS
	// make sure we have a matrix
	if (m_pA == 0) return false;

	// max number of iterations
	const int MAXITER = 150;

	// number of equations
	MKL_INT N = m_pA->Size();
	double* pa = m_pA->Values();
	int* ia = m_pA->Pointers();
	int* ja = m_pA->Indices();

	// data allocation
	MKL_INT ipar[128] = { 0 };
	double dpar[128] = { 0.0 };
	MKL_INT RCI_request;
	int M = (N < 150 ? N : 150); // this is the default value of ipar[4] and ipar[14]

	int nrestart = M;
	if (m_nrestart > 0) nrestart = m_nrestart;
	else if (m_maxiter > 0) nrestart = m_maxiter;

	int maxIter = M;
	if (m_maxiter > 0) maxIter = m_maxiter;

	vector<double> trvec(N);
	vector<double> b_copy(b);
	vector<double> residual(N);
	double* ptmp = &m_tmp[0];
	MKL_INT ivar = N;

	// initialize the solver
	dfgmres_init(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// parameters affecting the pre-conditioner
	if (m_checkZeroDiagonal)
	{
		ipar[30] = 1;
		dpar[30] = m_zeroThreshold;
		dpar[31] = m_zeroReplace;
	}
	else
	{
		// do this to avoid a warning from the preconditioner
		dpar[30] = m_fillTol;
	}

	// calculate the pre-conditioner
	int ierr;
	const int PCsize = (2 * m_maxfill + 1)*N + m_maxfill*(m_maxfill + 1) + 1;
	static vector<double> bilut; bilut.resize(PCsize, 0.0);
	static vector<int> jbilut; jbilut.resize(PCsize, 0);
	static vector<int> ibilut; ibilut.resize(N + 1, 0);

	if (m_doPreCond)
	{
		dcsrilut(&ivar, pa, ia, ja, &bilut[0], &ibilut[0], &jbilut[0], &m_fillTol, &m_maxfill, ipar, dpar, &ierr);
		if (ierr != 0) { MKL_Free_Buffers(); return false; }
		m_doPreCond = false;
	}

	// Set the desired parameters:
	ipar[ 4] = maxIter;			// max number of iterations
	ipar[14] = nrestart;		// nr of non-restarted iterations
	ipar[ 7] = 1;		// do the stopping test for maximal number of iterations
	ipar[ 8] = (m_doResidualTest ? 1 : 0);	// do residual stopping test
	ipar[ 9] = 0;		// do not request for the user defined stopping test
	ipar[10] = 1;		// do the pre-conditioned version of the FGMRES iterative solver
	ipar[11] = 1;		// do the check of the norm of the next generated vector automatically
	if (m_tol > 0) dpar[0] = m_tol;			// set the relative tolerance

	// Check the correctness and consistency of the newly set parameters
	dfgmres_check(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);
	if (RCI_request != 0) { MKL_Free_Buffers(); return false; }

	// solve the problem
	bool bdone = false;
	bool bconverged = false;
	while (!bdone)
	{
		// compute the solution via FGMRES
		dfgmres(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp);

		switch (RCI_request)
		{
		case 0: // solution converged. 
			bdone = true;
			bconverged = true;
			break;
		case 1:
		{
			// do matrix-vector multiplication
			m_pA->mult_vector(&m_tmp[ipar[21] - 1], &m_tmp[ipar[22] - 1]);

			if (m_print_level == 1)
			{
				fprintf(stderr, "%3d = %lg (%lg)\n", ipar[3], dpar[4], dpar[3]);
			}
		}
		break;
		case 3:	// do the pre-conditioning step
		{
			char cvar1 = 'L';
			char cvar = 'N';
			char cvar2 = 'U';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &bilut[0], &ibilut[0], &jbilut[0], &m_tmp[ipar[21] - 1], &trvec[0]);
			cvar1 = 'U';
			cvar = 'N';
			cvar2 = 'N';
			mkl_dcsrtrsv(&cvar1, &cvar, &cvar2, &ivar, &bilut[0], &ibilut[0], &jbilut[0], &trvec[0], &m_tmp[ipar[22] - 1]);
		}
		break;
		default:	// something went wrong
			bdone = true;
			bconverged = false;
		}
	}

	// get the solution. 
	MKL_INT itercount;
	dfgmres_get(&ivar, &x[0], &b[0], &RCI_request, ipar, dpar, ptmp, &itercount);
	if (m_print_level == 2)
	{
		fprintf(stderr, "%3d = %lg (%lg)\n", ipar[3], dpar[4], dpar[3]);
	}

	MKL_Free_Buffers();
	return bconverged;

#else
	return false;
#endif // MKL_ISS
}