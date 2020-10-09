#include <stdio.h>

#include <FECore/sdk.h>
#include <FECore/callback.h>
#include <FECore/FEModel.h>
#include <FECore/FECallback.h>
#include <FECore/FESurface.h>
#include <FECore/FENode.h>
#include <FECore/FEElement.h>
#include <FECore/FEPlotData.h>

class MyCallback : public FECallBack
{
public:
    MyCallback(FEModel* pfem) : FECallBack(pfem, CB_ALWAYS) 
    {
    }

    bool Execute(FEModel& fem, int nwhen)
    {
		if (nwhen == CB_SOLVED)
		{
			// Model is solved

			// create the plot variable "contact pressure"
			FEPlotData* ps = fecore_new<FEPlotData>("contact pressure", &fem);
			assert(ps->RegionType() == FE_REGION_SURFACE);
			assert(ps->StorageFormat() == FMT_ITEM);
			int datasize = ps->VarSize(ps->DataType());

			// loop over all surfaces
			FEMesh& m = fem.GetMesh();
			int NS = m.Surfaces();
			for (int i = 0; i < NS; ++i)
			{

				FESurface& S = m.Surface(i);
				int ne = S.Elements();
				int nsize = datasize * ne;
				const std::string name = S.GetName();

				printf("dumping surface %d: (%s)\n", i, name.c_str());

				// get data for this surface, "a" can be accessed like an array
				FEDataStream a; a.reserve(nsize);
				assert(ps->Save(S, a));
				
				// write a legacy VTK file for this surface containing surface data as scalar field
				FILE* pFile;
				int n;
				char fname[100];
				sprintf(fname, "surf%02d_%s.vtk", i, name.c_str());
				pFile = fopen(fname, "w");

				// header
				fprintf(pFile, "# vtk DataFile Version 2.0\n");
				fprintf(pFile, "surface %d with contact pressure\n", i);
				fprintf(pFile, "ASCII\n");
				fprintf(pFile, "DATASET POLYDATA\n");

				// write points / nodes
				fprintf(pFile, "POINTS %i FLOAT\n", S.Nodes());
				for (int n = 0; n < S.Nodes(); n++)
				{
					FENode node = S.Node(n);					
					fprintf(pFile, "%f %f %f\n", node.m_r0.x, node.m_r0.y, node.m_r0.z);					
				}

				// write polygons / surface elements
				int nconn = 0; // count total number entries in the table
				for (int e = 0; e < S.Elements(); e++) 
				{
					FEElement elem = S.Element(e);
					nconn += 1 + elem.Nodes(); // number of element nodes + list of element nodes
				}
				fprintf(pFile, "POLYGONS %i %i\n", ne, nconn);
				// write the table
				for (int e = 0; e < S.Elements(); e++)
				{
					FEElement elem = S.Element(e);
					vector<int> &enodes = elem.m_lnode;
					fprintf(pFile, "%d", (int)enodes.size());
					for(std::vector<int>::iterator it = enodes.begin(); it != enodes.end(); ++it)
						fprintf(pFile, " %d", *it);
					fprintf(pFile, "\n");
				}

				// write cell data
				fprintf(pFile, "CELL_DATA %d\n", (int)a.size());
				fprintf(pFile, "SCALARS contact_pressure float 1\n");
				fprintf(pFile, "LOOKUP_TABLE default\n");
				for (int k = 0; k < (int)a.size(); ++k) {
					fprintf(pFile, "%f\n", a[k]);
				}

				// done writing
				fclose(pFile);

			}
		}
			
		// all done!
        return true;
    }

};

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// feature classes
    REGISTER_FECORE_CLASS(MyCallback, "dump_callback");

}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}