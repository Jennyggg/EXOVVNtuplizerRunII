#include <string>

#include "TF2.h"

std::string GetApertureParametrisation_version2(const std::string &period, const std::string &arm)
{
	// extracted from config files in:
	//   https://github.com/cms-sw/cmssw/tree/916cb3d20213734a0465240720c8c8c392b92eac/Validation/CTPPS/python/simu_config
	
	if (period == "2016_preTS2")
	{
		if (arm == "arm0") return "3.76296E-05+(([xi]<0.117122)*0.00712775+([xi]>=0.117122)*0.0148651)*([xi]-0.117122)";
		if (arm == "arm1") return "1.85954E-05+(([xi]<0.14324)*0.00475349+([xi]>=0.14324)*0.00629514)*([xi]-0.14324)";
	}

	if (period == "2016_postTS2")
	{
		if (arm == "arm0") return "6.10374E-05+(([xi]<0.113491)*0.00795942+([xi]>=0.113491)*0.01935)*([xi]-0.113491)";
		if (arm == "arm1") return "([xi]-0.110)/130.0";
	}

	if (period == "2017_preTS2")
	{
		if (arm == "arm0") return "-(8.71198E-07*[xangle]-0.000134726)+(([xi]<(0.000264704*[xangle]+0.081951))*-(4.32065E-05*[xangle]-0.0130746)+([xi]>=(0.000264704*[xangle]+0.081951))*-(0.000183472*[xangle]-0.0395241))*([xi]-(0.000264704*[xangle]+0.081951))";
		if (arm == "arm1") return "3.43116E-05+(([xi]<(0.000626936*[xangle]+0.061324))*0.00654394+([xi]>=(0.000626936*[xangle]+0.061324))*-(0.000145164*[xangle]-0.0272919))*([xi]-(0.000626936*[xangle]+0.061324))";
	}

	if (period == "2017_postTS2")
	{
		if (arm == "arm0") return "-(8.92079E-07*[xangle]-0.000150214)+(([xi]<(0.000278622*[xangle]+0.0964383))*-(3.9541e-05*[xangle]-0.0115104)+([xi]>=(0.000278622*[xangle]+0.0964383))*-(0.000108249*[xangle]-0.0249303))*([xi]-(0.000278622*[xangle]+0.0964383))";
		if (arm == "arm1") return "4.56961E-05+(([xi]<(0.00075625*[xangle]+0.0643361))*-(3.01107e-05*[xangle]-0.00985126)+([xi]>=(0.00075625*[xangle]+0.0643361))*-(8.95437e-05*[xangle]-0.0169474))*([xi]-(0.00075625*[xangle]+0.0643361))";
	}

	if (period == "2018")
	{
		if (arm == "arm0") return "-(8.44219E-07*[xangle]-0.000100957)+(([xi]<(0.000247185*[xangle]+0.101599))*-(1.40289E-05*[xangle]-0.00727237)+([xi]>=(0.000247185*[xangle]+0.101599))*-(0.000107811*[xangle]-0.0261867))*([xi]-(0.000247185*[xangle]+0.101599))";
		if (arm == "arm1") return "-(-4.74758E-07*[xangle]+3.0881E-05)+(([xi]<(0.000727859*[xangle]+0.0722653))*-(2.43968E-05*[xangle]-0.0085461)+([xi]>=(0.000727859*[xangle]+0.0722653))*-(7.19216E-05*[xangle]-0.0148267))*([xi]-(0.000727859*[xangle]+0.0722653))";
	}

	throw "GetApertureParametrisation: unknown period or arm";
}
