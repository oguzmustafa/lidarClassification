#ifndef MCC_SOURCE3_CUH
#define MCC_SOURCE3_CUH

namespace mcc
{
	class IUnclassifiedPoints;

	class source3
	{
	public:
		void cells(double ***A, double*** cellsize, double*** res, int Asize, int celS);
	};
}

#endif