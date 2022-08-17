#ifndef ARES_H
#define ARES_H

struct AtomsRes{
      string resName;
      string atomName;
      int chgn;
      int resNum;
      float charge;
      float mass;
   };

struct waterM{
   int sites;
   int nwater;
   vector<int> Oind;
   float Ochg;
   float Hchg;
   float Mchg;
};

#endif
