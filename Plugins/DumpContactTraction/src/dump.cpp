#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>

#include <FECore/sdk.h>
#include <FECore/Callback.h>
#include <FECore/FEModel.h>
#include <FECore/FECallBack.h>
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
			FEPlotData* ps = fecore_new<FEPlotData>("contact traction", &fem);
			assert(ps->RegionType() == FE_REGION_SURFACE);
			assert(ps->StorageFormat() == FMT_ITEM);
			int datasize = ps->VarSize(ps->DataType());
			assert(datasize == 3); // ensure 3d vector field

			// get "surfaces"
			FEMesh& m = fem.GetMesh();
			int NS = m.Surfaces();
			
			// first create a list of surface names 
			// because surfaces are duplicated whenever referenced in a surface?/contact?-pair  
			std::vector<std::string> names;
			for (int i = 0; i < NS; ++i)
				names.push_back(m.Surface(i).GetName());

			// find unique names
			std::vector<std::string> unique_names;
			std::sort(names.begin(), names.end());  // need to be sorted before using unique
			std::unique_copy(names.begin(), names.end(), std::back_inserter(unique_names));

			// loop over unique surface names
			for (int i = 0; i < unique_names.size(); ++i)
			{
				const std::string name = unique_names[i];
				printf("found surface: %s\n", name.c_str()); 

				// get reference surface, count elements and size of data
				FESurface* Si = m.FindSurface(name); // what surface will this return???
				int ne = Si->Elements();
				int nsize = datasize * ne;
				
				// allocate and initialize to zero
				double* out = new double[nsize];
				for (int k = 0; k < nsize; ++k) {
					out[k] = 0.0;
				}

                // loop over all surfaces, match name and accumulate
				for (int j = 0; j < NS; ++j)
				{

					FESurface& Sj = m.Surface(j);

					// continue if names does not match
					if (name.compare(Sj.GetName()) != 0)
						continue;

					printf("adding contribution from surface %d\n", j);

					// check that number of elements matches
					assert(ne == Sj.Elements());

					// get data for this surface, "a" can be accessed like an array
					FEDataStream a; a.reserve(nsize);
					assert(ps->Save(Sj, a));
					assert((int)a.size() == nsize);

					// add to output
					for (int k = 0; k < nsize; ++k) {
						out[k] += a[k];
					}

				}

				// write a legacy VTK file for this surface containing surface data as scalar field
				FILE* pFile;
				int n;
				char fname[100];
				sprintf(fname, "surf_%s.vtk", name.c_str());
				pFile = fopen(fname, "w");
				
				// header
				fprintf(pFile, "# vtk DataFile Version 2.0\n");
				fprintf(pFile, "surface %d with contact pressure\n", i);
				fprintf(pFile, "ASCII\n");
				fprintf(pFile, "DATASET POLYDATA\n");

				// write points / nodes
				fprintf(pFile, "POINTS %i FLOAT\n", Si->Nodes());
				for (int n = 0; n < Si->Nodes(); n++)
				{
					FENode node = Si->Node(n);
					fprintf(pFile, "%f %f %f\n", node.m_r0.x, node.m_r0.y, node.m_r0.z);
				}
				
				// write polygons / surface elements
				int nconn = 0; // count total number entries in the table
				for (int e = 0; e < Si->Elements(); e++)
				{
					FEElement elem = Si->Element(e);
					nconn += 1 + elem.Nodes(); // number of element nodes + list of element nodes
				}
				fprintf(pFile, "POLYGONS %i %i\n", ne, nconn);
				// write the table
				for (int e = 0; e < Si->Elements(); e++)
				{
					FEElement elem = Si->Element(e);
					vector<int>& enodes = elem.m_lnode;
					fprintf(pFile, "%d", (int)enodes.size());
					for (std::vector<int>::iterator it = enodes.begin(); it != enodes.end(); ++it)
						fprintf(pFile, " %d", *it);
					fprintf(pFile, "\n");
				}

				//// write cell data
				//fprintf(pFile, "CELL_DATA %d\n", ne);
				//fprintf(pFile, "SCALARS contact_pressure float 1\n");
				//fprintf(pFile, "LOOKUP_TABLE default\n");
				//for (int k = 0; k < nsize; ++k) {
				//	fprintf(pFile, "%f\n", out[k]);
				//}

				// write cell data
				fprintf(pFile, "CELL_DATA %d\n", ne);
				fprintf(pFile, "VECTORS contact_traction float\n");
				for (int e = 0; e < Si->Elements(); ++e) {
					fprintf(pFile, "%f %f %f\n", out[3*e], out[3*e+1], out[3*e+2]);
				}
				fprintf(pFile, "SCALARS contact_pressure float 1\n");
				fprintf(pFile, "LOOKUP_TABLE default\n");
				for (int e = 0; e < Si->Elements(); e++)
				{
					FESurfaceElement& el = Si->Element(e);
					vec3d n = Si->SurfaceNormal(el, 0.5, 0.5);  // use normal in displaced element (natural coord (0.5,0.5) is center)					
					fprintf(pFile, "%f\n", -n.x*out[3*e]-n.y*out[3*e+1]-n.z*out[3*e+2]);
				}

				// done writing
				fclose(pFile);

				// clean up
				delete out;
					
			}

			printf("all done\n");

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
