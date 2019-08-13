
#ifndef GBANCROFT_H
#define GBANCROFT_H

#include "../newmat/newmatap.h"

#include "gall/gallobs.h"
#include "gall/gallnav.h"
#include "gutils/gtriple.h"
#include "gutils/gmutex.h"
#include "gset/gsetbase.h"

namespace gnut {   

int gbancroft(const Matrix& BBpass, ColumnVector& pos);

inline double lorentz(const ColumnVector& aa, const ColumnVector& bb);

} // namespace

#endif

