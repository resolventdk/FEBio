#include <FECore/sdk.h>

#include <stdio.h>
#include <FECore/callback.h>
#include <FECore/FECallback.h>

#include <FECore/FEModel.h>
#include <FECore/FEAnalysis.h>
#include <FECore/FETimeStepController.h>
#include <FECore/log.h>

class MyCallback : public FECallBack
{
public:
    MyCallback(FEModel* pfem) : FECallBack(pfem, CB_ALWAYS) 
    {
    }

    bool Execute(FEModel& fem, int nwhen)
    {

        if (nwhen == CB_MINOR_ITERS && fem.GetCurrentStep()->m_timeController) {

            double dt = fem.GetCurrentStep()->m_dt;
            double dtmin = fem.GetCurrentStep()->m_timeController->m_dtmin;

            if (dt < dtmin) {
                feLogError("dt<dtmin (%f<%f)\n", dt, dtmin); 
                throw std::runtime_error("dt<dtmin");
                return false;
            }

        }

        return true;
    }

};

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// feature classes
    REGISTER_FECORE_CLASS(MyCallback, "dtmin_callback");

}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}