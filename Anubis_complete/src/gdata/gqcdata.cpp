
/* ----------------------------------------------------------------------
  (c) 2018 G-Nut Software s.r.o. (services@gnutsoftware.com)

  Purpose: container for QC data
  Version: $Rev:$

  2018-05-20 /JD: created

-*/

#include <iostream>

#include "gdata/gqcdata.h"

using namespace std;

namespace gnut {

// Constructors
// ---------
t_gqcdata::t_gqcdata()
{
  id_type(t_gdata::QCDATA);
  id_group(t_gdata::GRP_QC);
}
   
// Destructor
// ---------
t_gqcdata::~t_gqcdata()
{
  gtrace("t_gqc::desctruct");
}

} // namespace
