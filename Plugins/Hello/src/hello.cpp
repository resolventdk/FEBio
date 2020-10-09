#include <FECore/sdk.h>

#include <stdio.h>
#include <FECore/callback.h>
#include <FECore/FECallback.h>

class MyCallback : public FECallBack
{
public:
    MyCallback(FEModel* pfem) : FECallBack(pfem, CB_ALWAYS) 
    {
    }

    bool Execute(FEModel& fem, int nwhen)
    {
        if (nwhen == CB_INIT)
        {
            // Model was initialized
            printf("hello from plugin!\n");
        }
        else if (nwhen == CB_SOLVED)
        {
            // Model is solved
        }

        return true;
    }

};

FECORE_PLUGIN void PluginInitialize(FECoreKernel& fecore)
{
	FECoreKernel::SetInstance(&fecore);

	// feature classes
    REGISTER_FECORE_CLASS(MyCallback, "hello_callback");

}

FECORE_PLUGIN int GetSDKVersion()
{
	return FE_SDK_VERSION;
}