#include <stdio.h>
#include <vector>

#include <FECore/sdk.h>
#include <FECore/Callback.h>
#include <FECore/FEModel.h>
#include <FECore/FECallBack.h>
#include <FECore/FEDomain.h>

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

			// get "domains"
			FEMesh& mesh = fem.GetMesh();
			int ND = mesh.Domains();
			
			// loop domains get the negative jacobians
            // dump to file
			FILE* pFile;
			pFile = fopen("negjac.dat", "w");
			for (int i = 0; i < ND; ++i)
			{
				FEDomain& dom = mesh.Domain(i);
				string name = dom.GetName();
				for (vector<int>::iterator it = dom.m_NegJacElems.begin(); it != dom.m_NegJacElems.end(); ++it){
					int eid = *it;
					FEElement* el = mesh.Element(eid);
					vec3d pos = mesh.Node(el->m_node[0]).m_rt;
					fprintf(pFile, "%d, %s, %f, %f, %f\n", eid, name.c_str(), pos.x, pos.y, pos.z);
				}
			}
			printf("done dumping negative jacobians to negjac.json\n");

        }



		// all done!
        return true;
    }

};

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// feature classes
    REGISTER_FECORE_CLASS(MyCallback, "negjac_callback");

}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}
