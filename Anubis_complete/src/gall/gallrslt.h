
#ifndef GALLRSLT_H
#define GALLRSLT_H

#include <vector>
#include <set>
#include <iostream>

#include "gutils/gtime.h"

namespace gnut {  

class t_gallrslt
{
 public:
    t_gallrslt();

    struct result
    {
        string type;
        string prn;

        t_gtime beg;
        t_gtime end;

        int index;

        double adj;
        double rms;
        double val;

        bool  operator <(const t_gallrslt::result&) const;
    };

    vector<t_gallrslt::result> v_rslt;

    set<t_gallrslt::result> set_crd;
    set<t_gallrslt::result> set_trp;
    set<t_gallrslt::result> set_clk;
    set<t_gallrslt::result> set_amb;

    void append(const t_gallrslt::result &rslt);
    void append(string type, t_gtime beg, t_gtime end, int idx, 
                string prn, double adj, double rms, double val);

};

} // namespace

#endif // GALLRSLT_H
