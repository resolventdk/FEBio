#pragma once
#include "SparseMatrix.h"

class FENewtonSolver;

//-----------------------------------------------------------------------------
// This is a class that just mimics a sparse matrix.
// It is only used by the JFNK strategy. 
// The only function it implements is the mult_vector.
class FECORE_API JFNKMatrix : public SparseMatrix
{
public:
	enum MultiplyPolicy {
		ZERO_FREE_DOFS,
		ZERO_PRESCRIBED_DOFS
	};

public:
	JFNKMatrix(FENewtonSolver* pns, SparseMatrix* K = 0);

	//! override multiply with vector (Does not use sparse matrix m_K)
	bool mult_vector(double* x, double* r) override;

	//! set the reference residual
	void SetReferenceResidual(std::vector<double>& R0);

	//! set matrix policy
	void SetPolicy(MultiplyPolicy p);

	//! set the forward difference epsilon
	void SetEpsilon(double eps);

public: // these functions use the actual sparse matrix m_K

	//! set all matrix elements to zero
	void Zero() override { m_K->Zero(); }

	//! Create a sparse matrix from a sparse-matrix profile
	void Create(SparseMatrixProfile& MP) override;

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, std::vector<int>& lm) override { m_K->Assemble(ke, lm); }

	//! assemble a matrix into the sparse matrix
	void Assemble(matrix& ke, std::vector<int>& lmi, std::vector<int>& lmj) override { m_K->Assemble(ke, lmi, lmj); }

	//! check if an entry was allocated
	bool check(int i, int j) override { return m_K->check(i, j); }

	//! set entry to value
	void set(int i, int j, double v) override { m_K->set(i, j, v); }

	//! add value to entry
	void add(int i, int j, double v) override { m_K->add(i, j, v); }

	//! retrieve value
	double get(int i, int j) override { return m_K->get(i, j); }

	//! get the diagonal value
	double diag(int i) override { return m_K->diag(i); }

	//! release memory for storing data
	void Clear() override { m_K->Clear(); }

	// interface to compact matrices
	double* Values() override { return m_K->Values(); }
	int*    Indices() override { return m_K->Indices(); }
	int*    Pointers() override { return m_K->Pointers(); }
	int     Offset() const override { return m_K->Offset(); }

private:
	double			m_eps;		// forward difference epsilon
	SparseMatrix*	m_K;		// the actual sparse matrix (This is only used as a preconditioner and can be null)
	FENewtonSolver*	m_pns;
	vector<double>	m_v, m_R;

	vector<double>	m_R0;

	vector<int>		m_freeNodes, m_prescribedNodes;
	MultiplyPolicy	m_policy;
};
