#ifndef _BV_H_
#define _BV_H_

#include "MatVec.h"

// A virtual class that defines an interface that needs to be implemented
// by any BV class that is used with the software.
class CBV
{
public:
  virtual REAL computeDistance(const CBV * bv2,
			       const REAL rot[3][3], const REAL trans[3]) = 0;
  virtual void updateBV(const CBV * bv1, const CBV * bv2,
			const REAL rot[3][3], const REAL trans[3]) = 0;

  virtual void undoBV() = 0;

  virtual REAL getVolume() = 0;
};


#endif
