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
#include "FEMMGSolutionMapper.h"
#include <FECore/FEModel.h>
#include <FECore/FEMeshTopo.h>
#include <FECore/FEMesh.h>
#include <FECore/FEDomain.h>
#include <FECore/FESolidDomain.h>
#include <FECore/FESurface.h>
#include <FECore/log.h>
#include <FECore/FEOctreeSearch.h>
#include <FECore/FENNQuery.h>
#include <FECore/FEMeshAdaptorCriterion.h>
#include <FECore/FEDomainMap.h>
#include "FELeastSquaresInterpolator.h"
#include "FEMeshShapeInterpolator.h"
#include "FEDomainShapeInterpolator.h"
#include <FECore/FECoreKernel.h>

#include <FEBioMech/FEConstPrestrain.h>
#include <FEBioMech/FEMechModel.h>
#include <FEBioMech/FERigidBody.h>

//#ifdef HAS_MMG
#include "mmg/mmg3d/libmmg3d.h"

//
//Implements remeshing via solution mapping to mitigate mesh distortion.
//
//Remeshing here sets both current and reference nodal positions to
//the positions of the new mesh nodes. Then the strain field of the solution
//on the old mesh is mapped from the old to the new mesh and used with a
//prestrain material to account for the displacement not captured by the
//zeroed node displacements. Further, the nodal displacements are also
//mapped and accumulated for plotting purposes.
//
//Solution mapping works by interpolating results from nodes in the old
//mesh to either nodes or integration points in the new mesh:
//1. For integration point variables the solution variables at the nodes
//of the old mesh is obtained by extrapolating and superpositioning
//values from the integration points to the nodes.
//2. Then the location of each point in the new mesh is obtained with
//respect to the old mesh.
//3. The variables are then interpolated from the nodes of the old
//element to the points in the new model.
//
//Mapping of integration point data is implemented in FERefineMesh and activated using the map data flag
//TODO mapping of none nodal fields are diffusive, hence ensure it is only done for elements actually remeshed.


class FEMMGSolutionMapper::MMG
{
public:
	MMG(FEMMGSolutionMapper* mmgRemesh) : m_mmgRemesh(mmgRemesh) {}
	bool build_mmg_mesh(MMG5_pMesh mmgMesg, MMG5_pSol mmgSol, FEMeshTopo& topo);
	bool build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem);

public:
	FEMMGSolutionMapper*	m_mmgRemesh;
	std::vector<double>	m_metric;	// refinement metric
	std::vector<int>	m_nodeSetTag;	// surface tags for node sets
};

//#endif

BEGIN_FECORE_CLASS(FEMMGSolutionMapper, FERefineMesh)
	ADD_PARAMETER(m_maxiter, "max_iters");
	ADD_PARAMETER(m_hmin, "min_element_size");
	ADD_PARAMETER(m_hausd, "hausdorff");
	ADD_PARAMETER(m_hgrad, "gradation");
	ADD_PARAMETER(m_relativeSize, "relative_size");
	ADD_PARAMETER(m_meshCoarsen, "mesh_coarsen");
	ADD_PARAMETER(m_normalizeData, "normalize_data");
	ADD_PROPERTY(m_criterion, "criterion");
	ADD_PROPERTY(m_sfunc, "size_function", 0);
END_FECORE_CLASS();

FEMMGSolutionMapper::FEMMGSolutionMapper(FEModel* fem) : FERefineMesh(fem)
{

	m_bremesh = true;

	m_maxiter = 1;
	m_maxelem = 0;
	m_relativeSize = true;
	m_meshCoarsen = true;
	m_normalizeData = false;

	m_hmin = 0.0;
	m_hausd = 0.01;
	m_hgrad = 1.3;

	m_criterion = nullptr;

	m_bmap_current = true; // use current mesh configuration for mapping
	m_transferMethod = TRANSFER_SHAPE; //TRANSFER_MLQ;
	m_nnc = 8;

	m_sfunc = nullptr;

	m_attachedRid = -1;
	m_atol = 0.01;

//#ifdef HAS_MMG
	mmg = new FEMMGSolutionMapper::MMG(this);
//#endif
}

bool FEMMGSolutionMapper::Init()
{

	FEModel& fem = *GetFEModel();

	FEMesh& mesh = fem.GetMesh();

	if (mesh.IsType(ET_TET4) == false) {
		printf("Mesh is not entirely tets\n");
		return false;
	}

	if (FERefineMesh::Init() == false) {
		printf("Failed in FERefineMesh::Init()\n");
		return false;
	}
	return true;

}

bool FEMMGSolutionMapper::Apply(int iteration)
{

	m_atol = 1.1 * m_hausd;

	FEModel& fem = *GetFEModel();

	FEMesh& mesh = fem.GetMesh();

	string name = "F0_map";

	// get custom data map used to store (pre)strain (deformation gradient)
	FEDomainMap* oldMap = dynamic_cast<FEDomainMap*>(mesh.FindDataMap(name)); // returns a copy of the pointer, so if it is changed it is not stored :(
	assert(oldMap);

	// create backup
	FEDomainMap backupMap = *oldMap;

	if (iteration == 0) {

		// create new map
		FEDomainMap newMap = FEDomainMap(FE_MAT3D, FMT_MATPOINTS); 
		newMap.Create(GetElementSet());
		newMap.SetName(name);

		//compute strain and add to prestrain for prestrain materials
		assert(GenerateGradients(newMap));

		// replace
		*oldMap = newMap;  // avoid pointer nightmare use copy method		

	}

	if (FERefineMesh::Apply(iteration)) {
		return true;
	}
	else {
		// restore
		*oldMap = backupMap;
		return false;
	}



}

bool FEMMGSolutionMapper::RefineMesh()
{

//#ifdef HAS_MMG
	FEModel& fem = *GetFEModel();
	FEMesh& mesh = fem.GetMesh();

	// check rigid body attachments, only 1 allowed
	FEElementSet* elset = GetElementSet();  // get adapted element set
	FENodeList nodeList = elset->GetNodeList();  // get list of nodes in elements adapted
	FEMechModel& mech_fem = dynamic_cast<FEMechModel&>(fem);
	m_attachedRid = -1;  // find attached rigid body
	for (int i = 0; i < nodeList.Size(); i++) {  // check nodes in set
		int rid = nodeList.Node(i)->m_rid;
		if (rid >= 0) { // attached to rigid
			if (m_attachedRid == -1) {  // if none
				m_attachedRid = rid;  // set new
			} else if (rid != m_attachedRid) { // diffferent from existing
				feLogError("Adapted element set cannnot be attached to more than one rigid body.");
			}
		}
	}

	// initialize the MMG mesh
	MMG5_pMesh mmgMesh = NULL;
	MMG5_pSol  mmgSol = NULL;
	MMG3D_Init_mesh(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	// --- build the MMG mesh ---
	FEMeshTopo& topo = *m_topo;
	if (mmg->build_mmg_mesh(mmgMesh, mmgSol, topo) == false) return false;
	
	// set the control parameters
	MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hmin, m_hmin);
	MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hausd, m_hausd);
	MMG3D_Set_dparameter(mmgMesh, mmgSol, MMG3D_DPARAM_hgrad, m_hgrad);

	// save mesh for debugging
	if (MMG3D_saveMshMesh(mmgMesh, mmgSol, "before.msh") == 0)
	{
		feLogError("Failed to save mmgMesh to .msh file.");
	}

	// run the mesher
	if (m_bremesh) {
	
		int ier = MMG3D_mmg3dlib(mmgMesh, mmgSol);

		if (ier == MMG5_STRONGFAILURE) {
			feLogError("MMG was not able to remesh the mesh.");
			return false;
		}
		else if (ier == MMG5_LOWFAILURE)
		{
			feLogError("MMG return low failure error");
		}
	}

	// save mesh for debugging
	if (MMG3D_saveMshMesh(mmgMesh, mmgSol, "after.msh") == 0)
	{
			feLogError("Failed to save mmgMesh to .msh file.");
	}

	// build the new mesh
	bool bret = mmg->build_new_mesh(mmgMesh, mmgSol, fem);

	// Clean up
	MMG3D_Free_all(MMG5_ARG_start,
		MMG5_ARG_ppMesh, &mmgMesh, MMG5_ARG_ppMet, &mmgSol,
		MMG5_ARG_end);

	return bret;

//#else
	return false;
//#endif
}

//#ifdef HAS_MMG

bool FEMMGSolutionMapper::MMG::build_mmg_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEMeshTopo& topo)
{

	FEMesh& mesh = *topo.GetMesh();
	int NN = mesh.Nodes();
	int NE = topo.Elements();
	int NF = topo.SurfaceFaces();

	// get the element set that may be remeshed
	FEElementSet* elset = m_mmgRemesh->GetElementSet();

	// get the criterion
	FEMeshAdaptorCriterion* criterion = m_mmgRemesh->GetCriterion();
	assert(criterion);
	if (criterion == nullptr) return false;
	FEMeshAdaptorSelection elemList = criterion->GetElementSelection(elset);
	if (elemList.size() == 0) {
		printf("empty element list from criterion\n");
		return false;  // nothing to remesh
	}

	// allocate mesh size
	if (MMG3D_Set_meshSize(mmgMesh, NN, NE, 0, NF, 0, 0) != 1)
	{
		assert(false);
		return false;
	}

	// set the vertex coordinates
	// From https://www.mmgtools.org/about-references :
	// Input point references are preserved 
	// (thus, collapse between points with different refs may be refused). 
	// Moreover, new points are created with the ref 0, so if you don’t    (hjs; seems that if the new nodes are located at referenced surfaces they are given surface reference as point reference)
	// use point references, it is better to set it to 0.
	// check for nodes attached to rigid bodies
	FENodeList nodeList = elset->GetNodeList();  // get list of nodes in elements adapted
	int remeshed = 0;
	int not_remeshed = 0;
	for (int i = 0; i < NN; ++i)
	{
		FENode& vi = mesh.Node(i);
		vec3d r = vi.m_rt;
		if (nodeList.GlobalToLocalID(i) != -1) {  // points that might be remeshed --> if 0 their will get random id if not remeshed
			MMG3D_Set_vertex(mmgMesh, r.x, r.y, r.z, 999, i + 1);  
			remeshed++;
		}
		else {
			MMG3D_Set_vertex(mmgMesh, r.x, r.y, r.z, 1000 + i, i + 1);  // set reference to 1000 + node id
			not_remeshed++;
		}
		
	}
	printf("To be remeshed, not to be remeshed: %d, %d\n", remeshed, not_remeshed);

	// set the tetrahedra, ref starts from 1
	int c = 1;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FEDomain& dom = mesh.Domain(i);
		int ne = dom.Elements();
		for (int j = 0; j < ne; ++j, ++c)
		{
			FEElement& e = dom.ElementRef(j);
			int* n = &e.m_node[0];
			MMG3D_Set_tetrahedron(mmgMesh, n[0] + 1, n[1] + 1, n[2] + 1, n[3] + 1, i+1, c);
		}
	}

	// set the facet markers
	vector<int> faceMarker(NF, 0);

	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);
		vector<int> faceIndexList = topo.SurfaceFaceIndexList(surf);
		for (int j = 0; j < faceIndexList.size(); ++j)
		{
			faceMarker[faceIndexList[j]] = faceMark;
		}
		faceMark++;
	}

	// for node sets we are going to create artificial surfaces
	m_nodeSetTag.assign(mesh.NodeSets(), -1);
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		FENodeSet& nset = *mesh.NodeSet(i);
		if (nset.Size() != mesh.Nodes())
		{
			vector<int> nodeTags(mesh.Nodes(), 0);
			for (int j = 0; j < nset.Size(); ++j) nodeTags[nset[j]] = 2;

			// see if this is indeed a surface node set
			for (int j = 0; j < NF; ++j)
			{
				const FEFaceList::FACE& face = topo.SurfaceFace(j);
				const int* fn = face.node;
				if ((nodeTags[fn[0]] != 0) && (nodeTags[fn[1]] != 0) && (nodeTags[fn[2]] != 0))
				{
					nodeTags[fn[0]] = 1;
					nodeTags[fn[1]] = 1;
					nodeTags[fn[2]] = 1;
				}
			}
			int twos = 0;
			for (int j = 0; j < mesh.Nodes(); ++j) if (nodeTags[j] == 2) twos++;

			if (twos == 0)
			{
				for (int j = 0; j < NF; ++j)
				{
					const FEFaceList::FACE& face = topo.SurfaceFace(j);
					const int* fn = face.node;
					if ((nodeTags[fn[0]] == 1) && (nodeTags[fn[1]] == 1) && (nodeTags[fn[2]] == 1))
					{
						if (faceMarker[j] == 0)
						{
							faceMarker[j] = faceMark;
							m_nodeSetTag[i] = faceMark;
						}
						else
						{
							if (m_nodeSetTag[i] == -1)
							{
								m_nodeSetTag[i] = faceMarker[j];
							}
							else if (faceMarker[j] != m_nodeSetTag[i])
							{
								printf("Non-surface nodeset '%s' encountered, returning...\n", nset.GetName().c_str());
								return false;
							}
						}
					}
				}
			}
			faceMark++;
		}
	}

	// create the faces
	for (int i = 0; i < NF; ++i)
	{
		const FEFaceList::FACE& f = topo.SurfaceFace(i);
		const int* n = &f.node[0];
		MMG3D_Set_triangle(mmgMesh, n[0] + 1, n[1] + 1, n[2] + 1, faceMarker[i], i + 1);
	}

	// Now, we build the "solution", i.e. the target element size.
	// If no elements are selected, we set a homogenous remeshing using the element size parameter.
	// set the "solution", i.e. desired element size
	if (MMG3D_Set_solSize(mmgMesh, mmgSol, MMG5_Vertex, NN, MMG5_Scalar) != 1)
	{
		assert(false);
		return false;
	}

	if (m_metric.empty())
	{
		// build the edge length table
		int ET_TET[6][2] = { { 0,1 },{ 1,2 },{ 2,0 },{ 0,3 },{ 1,3 },{ 2,3 } };
		vector<pair<double, int> > edgeLength(NN, pair<double, int>(0.0, 0));
		for (int i = 0; i < NE; ++i)
		{
			FEElement& el = *topo.Element(i);
			for (int j = 0; j < 6; ++j)
			{
				int a = el.m_node[ET_TET[j][0]];
				int b = el.m_node[ET_TET[j][1]];

				vec3d ra = mesh.Node(a).m_rt;
				vec3d rb = mesh.Node(b).m_rt;

				double L = (ra - rb).norm2();

//				if (L > edgeLength[a].first) edgeLength[a].first = L;
//				if (L > edgeLength[b].first) edgeLength[b].first = L;

				edgeLength[a].first += L; edgeLength[a].second++;
				edgeLength[b].first += L; edgeLength[b].second++;
			}
		}

		m_metric.resize(NN, 0.0);
		for (int i = 0; i < NN; ++i)
		{
			if (edgeLength[i].second != 0)
			{
				edgeLength[i].first /= (double)edgeLength[i].second;
				edgeLength[i].first = sqrt(edgeLength[i].first);
			}
			m_metric[i] = edgeLength[i].first;
		}
	}

	if (elset)
	{
		// elements that are not in the element set will be flagged as required.
		FEElementIterator it(&mesh);
		int c = 1;
		for (; it.isValid(); ++it, ++c)
		{
			FEElement& el = *it;
			if (elset->Contains(el) == false)
			{
				MMG3D_Set_requiredTetrahedron(mmgMesh, c);
			}
		}
	}

	// scale factors
	vector<double> nodeScale(NN, 0.0);

	// see if want to normalize the data
	if (m_mmgRemesh->m_normalizeData)
	{
		// Find data range
		double vmin, vmax;
		for (int i = 0; i < (int)elemList.size(); ++i)
		{
			double v = elemList[i].m_elemValue;
			if ((i == 0) || (v < vmin)) vmin = v;
			if ((i == 0) || (v > vmax)) vmax = v;
		}
		if (vmax == vmin) vmax++;

		// normalize data
		for (int i = 0; i < (int)elemList.size(); ++i)
		{
			double v = elemList[i].m_elemValue;
			elemList[i].m_elemValue = (v - vmin) / (vmax - vmin);
		}
	}

	// map to nodal data
	vector<int> tag(NN, 0);
	for (int i = 0; i < (int)elemList.size(); ++i)
	{
		FEElement& el = *mesh.FindElementFromID(elemList[i].m_elementId);
		for (int j = 0; j < el.Nodes(); ++j)
		{
			double s = elemList[i].m_elemValue;
			FEFunction1D* fs = m_mmgRemesh->m_sfunc;
			if (fs)
			{
				s = fs->value(s);
			}
			assert(s > 0.0);
			if (s <= 0.0) return false;
			nodeScale[el.m_node[j]] += s;
			tag[el.m_node[j]]++;
		}
	}
	for (int i = 0; i < NN; ++i)
	{
		if (tag[i] > 0) nodeScale[i] /= (double)tag[i];
		else nodeScale[i] = (m_mmgRemesh->m_relativeSize ? 1.0 : m_metric[i]);
	}

	// adjust for relative scale flag
	if (m_mmgRemesh->m_relativeSize)
	{
		for (int k = 0; k < NN; k++) {
			nodeScale[k] *= m_metric[k];
		}
	}

	// determine new size field
	bool meshCoarsen = m_mmgRemesh->m_meshCoarsen;
	for (int k = 0; k < NN; k++) 
	{
		double s = nodeScale[k];
		if ((meshCoarsen) || (s < m_metric[k])) m_metric[k] = s;
	}

	// set the new metric
	for (int k = 0; k < NN; k++) {
		MMG3D_Set_scalarSol(mmgSol, m_metric[k], k + 1);
	}

	return true;
}

bool FEMMGSolutionMapper::MMG::build_new_mesh(MMG5_pMesh mmgMesh, MMG5_pSol mmgSol, FEModel& fem)
{

	FEMesh& mesh = fem.GetMesh();
	int N0 = mesh.Nodes();

	// get the domain that was remeshed
	const FEElementSet* oldElset = m_mmgRemesh->GetElementSet();
	const FEDomainList& domainList = oldElset->GetDomainList();
	if (domainList.Domains() != 1) {
		printf("Current implementation limited to adapting a single domain only!\n");
		return false;  // 
	}
	FEDomain& remeshedDom = const_cast<FEDomain&>(*domainList.GetDomain(0));

	// get the new mesh sizes
	int nodes, elems, faces;
	MMG3D_Get_meshSize(mmgMesh, &nodes, &elems, NULL, &faces, NULL, NULL);

	// get old node positions
	vector<vec3d> oldNodePos(N0);
	for (int i = 0; i < N0; ++i) oldNodePos[i] = mesh.Node(i).m_rt;

	// copy new nodal positions
	vector<vec3d> newNodePos(nodes);
	vector<int> newNodeRef(nodes);
	for (int i = 0; i < nodes; ++i)
	{
		MMG3D_Get_vertex(mmgMesh, &newNodePos[i].x, &newNodePos[i].y, &newNodePos[i].z,
			&newNodeRef[i], NULL, NULL);
	}

	// copy the metric
	m_metric.resize(nodes, 0.0);
	for (int i = 0; i < nodes; ++i)
	{
		MMG3D_Get_scalarSol(mmgSol, &m_metric[i]);
	}

	int MAX_DOFS = fem.GetDOFS().GetTotalDOFS();
	vector<vec3d> nodeDispl(nodes);
	vector<vec3d> nodeAccel(nodes);
	vector<vector<double> > nodeVal(nodes, vector<double>(MAX_DOFS, 0.0));

	// allocate data mapper
	// In case of contact problems the domains may overlap.
	// So if we do not map/interpolate domains seperately errors occur in the overlapping regions when a point is mapped to the wrong domain.
	// TODO: we still get problems if domain is self overlapping, so be careful.
	FEMeshDataInterpolator* mapper = nullptr;
	switch (m_mmgRemesh->m_transferMethod)
	{
	case TRANSFER_SHAPE: mapper = new FEDomainShapeInterpolator(&remeshedDom, m_mmgRemesh->m_bmap_current, m_mmgRemesh->m_atol); break;
	case TRANSFER_MLQ:
	{
		FELeastSquaresInterpolator* MLQ = new FELeastSquaresInterpolator;
		MLQ->SetNearestNeighborCount(m_mmgRemesh->m_nnc);
		MLQ->SetDimension(m_mmgRemesh->m_nsdim);
		MLQ->SetSourcePoints(oldNodePos); // dont use all nodes
		mapper = MLQ;
	}
	break;
	default:
		assert(false);
		return false;
	}

	// map nodal positions and nodal data
	for (int i = 0; i < nodes; ++i)
	{
		vec3d ri = newNodePos[i];
		int old_nid = newNodeRef[i] - 1000; // TODO use m_offset

		if (old_nid >= 0) {  // existing node --> copy directly from old mesh

			// get the nodal displacement
			nodeDispl[i] = mesh.Node(old_nid).m_rt - mesh.Node(old_nid).m_r0;
			
			// get the nodal acceleration
			nodeAccel[i] = mesh.Node(old_nid).m_at;

			// get dof values
			for (int l = 0; l < MAX_DOFS; ++l)
				nodeVal[i][l] = mesh.Node(old_nid).get(l);

		} else {  // new node --> use interpolation

			if (mapper->SetTargetPoint(ri) == false)
			{
				assert(false);
				delete mapper;
				return false;
			}

			// get the nodal displacement
			nodeDispl[i] = mapper->MapVec3d([&mesh](int sourceNode) {
				return (mesh.Node(sourceNode).m_rt - mesh.Node(sourceNode).m_r0);
				});

			// get the nodal acceleration
			nodeAccel[i] = mapper->MapVec3d([&mesh](int sourceNode) {
				return mesh.Node(sourceNode).m_at;
				});

			// no need to map director, which is unit vector orthogonal to shell mid-surface

			// update values
			for (int l = 0; l < MAX_DOFS; ++l)
			{
				nodeVal[i][l] = mapper->Map([&mesh, l](int sourceNode) {
					return mesh.Node(sourceNode).get(l);
					});
			}

		} // if old_nid >= 0
	} // i in nodes
	delete mapper;

	// reallocate nodes
	mesh.CreateNodes(nodes);

	int zeros = 0;
	int remeshed = 0;
	int not_remeshed = 0;

	// assign dofs to new nodes. hjs TODO check dofs exists
	const int dof_X = fem.GetFEModel()->GetDOFIndex("x");
	const int dof_Y = fem.GetFEModel()->GetDOFIndex("y");
	const int dof_Z = fem.GetFEModel()->GetDOFIndex("z");

	// get attached RB
	FERigidBody* RB;
	quatd rotInv, rot;
	if (m_mmgRemesh->m_attachedRid >= 0) {
		RB = dynamic_cast<FEMechModel&>(fem).GetRigidBody(m_mmgRemesh->m_attachedRid);
		rot = RB->GetRotation();
		rotInv = RB->GetRotation().Inverse();
	}

	for (int i = 0; i < nodes; ++i)
	{
		FENode& node = mesh.Node(i);
		node.SetDOFS(MAX_DOFS);

		if (i == 0) printf("not renumbering nid %d -> %d\n", node.GetID(), i);

		if (newNodeRef[i] == 0)
			zeros++;

		if (newNodeRef[i] < 1000) { // remeshed -> deformation should be reset
			remeshed++;

			node.m_r0 = node.m_rt = newNodePos[i];  // hjs; reset deformation
			node.m_at = nodeAccel[i];

			for (int j = 0; j < node.m_ID.size(); ++j) {
				node.set(j, nodeVal[i][j]);
			}

			// zero displacement dofs, TODO store somewhere else
			node.set(dof_X, 0.0);
			node.set(dof_Y, 0.0);
			node.set(dof_Z, 0.0);


			// If any nodes in remeshed selection is attached to a rigid body
			// do not subtract the rigid body motion, which would then be added twice because of the implementation in FERigidSystem.cpp
			if (m_mmgRemesh->m_attachedRid >= 0) 
			{ 
				node.m_r0 = rotInv * (node.m_rt - RB->m_rt) + RB->m_r0;

				node.set(dof_X, node.m_rt.x - node.m_r0.x);
				node.set(dof_Y, node.m_rt.y - node.m_r0.y);
				node.set(dof_Z, node.m_rt.z - node.m_r0.z);

			}

			if (m_mmgRemesh->m_nsdim == 2) {
				node.m_rt.z = node.m_r0.z;
				node.set(dof_Z, 0.0);
			}

		} else {  // not remeshed -> deformation should NOT be reset
			not_remeshed++;
			node.m_rt = newNodePos[i];
			node.m_r0 = newNodePos[i] - nodeDispl[i]; // TODO use dof_X,Y,Z
			node.m_at = nodeAccel[i];
			if (m_mmgRemesh->m_nsdim == 2) node.m_rt.z = node.m_r0.z;
			for (int j = 0; j < node.m_ID.size(); ++j) {
				node.set(j, nodeVal[i][j]);
			}
		} 

		// update
		node.UpdateValues();
	}
	printf("remeshed, not remeshed, zeros: %d, %d, %d\n", remeshed, not_remeshed, zeros);

	// recreate domains
	int eID = 1;
	for (int i = 0; i < mesh.Domains(); ++i)
	{
		FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(i));

		int nelems = 0;
		for (int j = 0; j < elems; ++j)
		{
			int n[4], gid, breq;
			MMG3D_Get_tetrahedron(mmgMesh, n, n + 1, n + 2, n + 3, &gid, &breq);
			if (gid == i+1) nelems++;
		}

		dom.Create(nelems, dom.GetElementSpec());
		int c = 0;
		for (int j = 0; j < elems; ++j)
		{
			int n[4], gid;
			MMG3D_Get_tetrahedron(mmgMesh, n, n + 1, n + 2, n + 3, &gid, NULL);

			if (gid == i+1)
			{
				FESolidElement& el = dom.Element(c++);
				if (eID == 1) printf("renumbering eid %d -> %d\n", el.GetID(), eID);
				el.SetID(eID++);  // hjs: renumber element ID so that it is consistent with plotfile
				el.m_node[0] = n[0] - 1;
				el.m_node[1] = n[1] - 1;
				el.m_node[2] = n[2] - 1;
				el.m_node[3] = n[3] - 1;
			}
		}

		// re-init domain
		printf("reiniting domain %s (%d)\n", dom.GetName().c_str(), i);
		dom.CreateMaterialPointData();
		dom.Reset();	// NOTE: we need to call this to actually call the Init function on the material points.
		dom.Init();    // hjs: this where Jacobians of ref mesh is initialised and negative determinants are checked!
		dom.Activate();
	}
	mesh.RebuildLUT();

	// recreate element sets
	for (int i = 0; i < mesh.ElementSets(); ++i)
	{
		FEElementSet& eset = mesh.ElementSet(i);

		// get the domain list
		// NOTE: Don't get the reference, since then the same reference
		// is passed to Create below, which causes problems.
		FEDomainList domList = eset.GetDomainList();
		if (domList.IsEmpty()) { throw std::runtime_error("Error in FEMMGSolutionMapper!"); }

		// recreate the element set from the domain list
		eset.Create(domList);
	}

	// get the element set that may be remeshed
	FEElementSet* elset = m_mmgRemesh->GetElementSet();
	FENodeList nodeList = elset->GetNodeList();  // get list of nodes in elements adapted
	for (int i = 0; i < nodes; ++i)
	{
		if (nodeList.GlobalToLocalID(i) != -1) {  // points that might be remeshed --> if 0 their will get random id if not remeshed
			if (newNodeRef[i] >= 1000)
				printf("Error in vertex reference for remeshed domain: ref=%d\n", newNodeRef[i]);
		}
		else {
			if (newNodeRef[i] < 1000)
				printf("Error in vertex reference for not remeshed domain: ref=%d\n", newNodeRef[i]);
		}
	}

	// recreate surfaces
	int faceMark = 1;
	for (int i = 0; i < mesh.Surfaces(); ++i)
	{
		FESurface& surf = mesh.Surface(i);

		// count faces
		int nfaces = 0;
		for (int j = 0; j < faces; ++j)
		{
			int n[3], gid, breq;
			MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, &breq);
			if (gid == faceMark) nfaces++;
		}
		assert(nfaces > 0);
		surf.Create(nfaces);
		int c = 0;
		for (int j = 0; j < faces; ++j)
		{
			int n[3], gid;
			MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, NULL);
			if (gid == faceMark)
			{
				FESurfaceElement& face = surf.Element(c++);
				face.SetType(FE_TRI3G3);
				face.m_node[0] = n[0] - 1;
				face.m_node[1] = n[1] - 1;
				face.m_node[2] = n[2] - 1;
			}
		}
		assert(c == nfaces);
		surf.CreateMaterialPointData();
		surf.Init();

		// also update the facet set if the surface has one
		FEFacetSet* fset = surf.GetFacetSet();
		if (fset)
		{
			fset->Create(surf);
		}

		faceMark++;
	}

	// update nodesets
	for (int i = 0; i < mesh.NodeSets(); ++i)
	{
		int tag = m_nodeSetTag[i];
		FENodeSet& nset = *mesh.NodeSet(i);
		if (nset.Size() != N0)
		{
			vector<int> nodeTags(mesh.Nodes(), 0);
			for (int j = 0; j < faces; ++j)
			{
				int n[3], gid;
				MMG3D_Get_triangle(mmgMesh, n, n + 1, n + 2, &gid, NULL);
				if (gid == tag)
				{
					nodeTags[n[0] - 1] = 1;
					nodeTags[n[1] - 1] = 1;
					nodeTags[n[2] - 1] = 1;
				}
			}

			std::vector<int> nodeList;
			for (int i = 0; i < nodeTags.size(); ++i) if (nodeTags[i] == 1) nodeList.push_back(i);

			if (nodeList.size() > 0)
			{
				nset.Clear();
				nset.Add(nodeList);
			}

			faceMark++;
		}
		else
		{
			// assume this node set is determined by the entire mesh
			FESolidDomain& dom = dynamic_cast<FESolidDomain&>(mesh.Domain(0));
			//			if (nset.Size() == dom.Nodes())
			{
				nset.Clear();
				for (int i = 0; i < dom.Nodes(); ++i) nset.Add(dom.NodeIndex(i));
			}
		}
	}

	return true;
}

//#endif

bool FEMMGSolutionMapper::GenerateGradients(FEDomainMap& map)
{
	const FEElementSet& set = *map.GetElementSet();

	FEMesh& mesh = *set.GetMesh();

	FEDataType dataType = map.DataType();
	if (dataType != FE_MAT3D) return false;

	int storageFormat = map.StorageFormat();
	if (storageFormat != FMT_MATPOINTS) return false;

	// get the domain of the set, works only for a single domain
	const FEDomainList& domainList = set.GetDomainList();
	if (domainList.Domains() != 1) return false;
	FESolidDomain& dom = const_cast<FESolidDomain&>(dynamic_cast<const FESolidDomain&>(*domainList.GetDomain(0)));  
	// developers should make methods const whenever possible, so const_cast is not neccesary
	// TODO check valid domain

	// loop all elements compute deformation gradients in integration points and multiply existing prestrain gradient, then store in map
	int N = set.Elements();
	for (int i = 0; i < N; ++i)
	{
		FESolidElement* pel = dynamic_cast<FESolidElement*>(mesh.FindElementFromID(set[i]));
		if (pel == nullptr) return false;
		FESolidElement& el = *pel;

		int ni = el.GaussPoints();
		for (int j = 0; j < ni; ++j)
		{
			FEPrestrainMaterialPoint* pp = dynamic_cast<FEPrestrainMaterialPoint*>(el.GetMaterialPoint(j));
			assert(pp);
			mat3d F0 = pp->initialPrestrain();  // initial prestrain, do not include correction
			mat3d F; 
			dom.defgrad(el, F, j);  // compute deformation gradient from current mesh displacement relative to reference
			map.setValue(i, j, F*F0);  // product, notice if F0=I, then after remesh 1 you get F*F0=F as intended
		}
	}
	return true;
}
