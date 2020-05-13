#ifndef MYPROJECT_PARTDIST_H
#define MYPROJECT_PARTDIST_H

#include "Particle.h"
#include <vector>
#include <random>

class PartDist {
private:
public:
    PartDist();
    static void setCube( vector<Particle> &v, const double &l);
    void setRand( vector<Particle> &v, const double &l, const double &s);
    static void setSpeed( vector<Particle> &v, const double &MolecMass, const double &T, const double &e);
};

bool checkDistance( Particle &p, vector<Particle> &v, const double &s, const int &i, const double &l);

#endif //MYPROJECT_PARTDIST_H
