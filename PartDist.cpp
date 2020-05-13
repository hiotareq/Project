//
// Created by bobro on 10.05.2020.
//

#include "PartDist.h"
#include <vector>
#include <random>
#include <cmath>

PartDist::PartDist() = default;

void PartDist::setCube(vector<Particle> &v, const double &l){
    int q = cbrt(v.size());
    int c=0;
    for (auto i = 0; i < q; i++){
        for ( auto j = 0; j < q; j++){
            for ( auto k = 0; k < q;k++){
                v[c].setParticle(i * l/(q-1), j * l/(q-1), k * l/(q-1));
                c++;
            }
        }
    }
}

bool checkDistance( Particle &p, vector<Particle> &v, const double &s, const int &i, const double &l){
    int q = v.size();
    for (int j = 0;j < q; j++){
        if ( j == i ) continue;
        else{
            if ( v[i].getDistance(v[j], l)/s < 0.95 ){
                return false;
            }
        }
    }
    return true;
}

//на самом деле, эта функция может сожрать ооооочень много времени, но иначе частицы распределяются неправдоподообным образом
void PartDist::setRand(vector<Particle> &v, const double &l, const double &s) {
    int q = v.size();
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<> dis(0,l);
    for (int i=0;i < q;i++){
        v[i].setParticle(dis(gen), dis(gen), dis(gen));
        while (!checkDistance(v[i],v, s, i, l)){
            v[i].setParticle( dis(gen), dis(gen), dis(gen));
        }
    }
}

void PartDist::setSpeed( vector<Particle> &v, const double &MolecMass, const double &T, const double &e){
    int q = v.size();
    random_device rd;
    mt19937_64 gen(rd());
    normal_distribution<> d(0,(e * T)/MolecMass);
    for (unsigned int i=0;i<q;i++){
        v[i].setvX(d(gen));
        v[i].setvY(d(gen));
        v[i].setvZ(d(gen));
    }
}
