#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>

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
			FEMesh& m = fem.GetMesh();
			int ND = m.Domains();
			
			// loop domains collect the negative jacobians
			vector<int> negJacElems;
			for (int i = 0; i < ND; ++i)
			{
				FEDomain& dom = m.Domain(i);
				for (vector<int>::iterator it = dom.m_NegJacElems.begin(); it != dom.m_NegJacElems.end(); ++it)
					negJacElems.push_back(*it);

			}
			
			// sort
			sort(negJacElems.begin(), negJacElems.end());

			// dump to file
			FILE* pFile;
			pFile = fopen("negjac.dat", "w");
			for (vector<int>::iterator it = negJacElems.begin(); it != negJacElems.end(); ++it)
				fprintf(pFile, "%d\n", *it);

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
