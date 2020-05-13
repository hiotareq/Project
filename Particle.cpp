#include <iostream>
#include <cmath>

#include "Particle.h"


using namespace std;

Particle::Particle(double new_x, double new_y, double new_z) : x(new_x), y(new_y), z(new_z), aX(0),aY(0),aZ(0),vX(0),vY(0), vZ(0)   {}
Particle:: Particle(): Particle(0,0,0){}
Particle::~Particle()=default;

double Particle::getX() const {
	return x;
}
double Particle::getY() const {
	return y;
}
double Particle::getZ() const {
	return z;
}
double Particle::getVX() const {
	return vX;
}
double Particle::getVY() const{
	return vY;
}
double Particle::getVZ() const {
	return vZ;
}
double Particle::getaX() const {
	return aX;
}
double Particle::getaY() const {
	return aY;
}
double Particle::getaZ() const{
	return aZ;
}

void Particle::setX(const double &new_x) {
	x = new_x;
}
void Particle::setY(const double &new_y) {
	y = new_y;
}
void Particle::setZ(const double &new_z) {
	z = new_z;
}
void Particle::addaX(const double &add_aX){
	aX+=add_aX;
}
void Particle::addaY(const double &add_aY) {
	aY+=add_aY;
}
void Particle::addaZ(const double &add_aZ) {
	aZ+=add_aZ;
}
void Particle::setvX(const double &new_vX) {
	vX+=new_vX;
}
void Particle:: setvY(const double &new_vY) {
	vY+=new_vY;
}
void Particle::setvZ(const double &new_vZ) {
	vZ+=new_vZ;
}
void Particle::setaX(const double &newaX) {
    aX = newaX;
}
void Particle::setaY(const double &newaY) {
    aY = newaY;
}
void Particle::setaZ(const double &newaX) {
    aZ = newaX;
}

double Particle::getDistance(const Particle& p, const double &l) const {
    double rx, ry, rz;
    rx = x - p.x;
    ry = y - p.y;
    rz = z - p.z;
    if ( rx > l/2) rx = l -rx;
    else {
        if ( rx < -l/2) rx = l + rx;
    }
    if ( ry > l/2) ry = l -ry;
    else {
        if ( ry < -l/2) ry = l + ry;
    }
    if ( rz > l/2) rz = l -rz;
    else {
        if ( rz < -l/2) rz = l + rz;
    }
    return sqrt( rx*rx + ry*ry + rz*rz);
}
double Particle::getPotential(const Particle& p,  const double &a, const double &b, const double &l) const{
	double r = getDistance(p, l);
	return (4 * b * (pow(a / r, 12) - pow(a / r, 6)));
}
double Particle::getForce(const Particle& p, const double &s, const double &e, const double &l)const {
	double r = getDistance(p, l);
	if ( r == 0) return 0;
	if ( r < 2.5 * s) return (24 * e * (pow( s / r, 12)/( r * r) - pow ( s / r , 6) )/ ( r * r)) ;
	else return 0;
}

void Particle::setParticle(const double &new_x, const double &new_y, const double &new_z) {
    x = new_x;
    y = new_y;
    z = new_z;
}

void Step( vector<Particle> &v,const float &t, const double &s, const double &e, const double &l, const double &mo_mass){
    unsigned int q = v.size();
    for ( unsigned int i =0 ; i < q ; i++){
        v[i].setaX(0);
        v[i].setaY(0);
        v[i].setaZ(0);
        for  (unsigned int j=0; j < q ; j++){
            if ( i == j ) continue;
            else {
                v[i].addaX(v[i].getForce( v[j], s, e, l)*(v[i].getX() - v[j].getX())/mo_mass);
                v[i].addaY(v[i].getForce( v[j], s, e, l)*(v[i].getY() - v[j].getY())/mo_mass);
                v[i].addaZ(v[i].getForce( v[j], s, e, l)*(v[i].getZ() - v[j].getZ())/mo_mass);
            }
        }
    }

    for (unsigned int i = 0 ; i < q; i++){
        v[i].setvX(t * v[i].getaX());
        v[i].setvY(t * v[i].getaY());
        v[i].setvZ(t * v[i].getaZ());
    }

    for (unsigned int i =0 ;i < q; i++){
        v[i].setParticle(v[i].getX() +t*v[i].getVX(), v[i].getY() + t*v[i].getVY(), v[i].getZ() + t*v[i].getVZ());
    }

    for (unsigned int i = 0; i < q ; i++){
        if (v[i].getX() > l ) v[i].setX(v[i].getX() - l);
        if (v[i].getX() < 0 ) v[i].setX(v[i].getX() + l);
        if (v[i].getY() > l ) v[i].setY(v[i].getY() - l);
        if (v[i].getY() < 0 ) v[i].setY(v[i].getY() + l);
        if (v[i].getZ() > l ) v[i].setZ(v[i].getZ() - l);
        if (v[i].getZ() < 0 ) v[i].setZ(v[i].getZ() + l);
    }
}
void impulse_check( vector<Particle> &v, const double &m){
    double i_x = 0, i_y = 0, i_z = 0;
    int q = v.size();
    for ( int i=0; i< q; i++){
        i_x += v[i].getVX();
        i_y += v[i].getVY();
        i_z += v[i].getVZ();
    }
    if ( abs(i_x) > 0.001) {
        if ( abs(i_y) > 0.001) {
            if ( abs(i_z) > 0.001) {
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvX(-i_x/q);
                    v[i].setvY(-i_y/q);
                    v[i].setvZ(-i_z/q);
                }
            }
            else{
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvX(-i_x/q);
                    v[i].setvY(-i_y/q);
                }
            }
        }
        else{
            if ( abs(i_z) > 0.001){
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvX(-i_x/q);
                    v[i].setvZ(-i_z/q);
                }
            }
            else{
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvX(-i_x/q);
                }
            }
        }
    }
    else{
        if ( abs(i_y) > 0.001){
            if ( abs(i_z) > 0.001){
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvY(-i_y/q);
                    v[i].setvZ(-i_z/q);
                }
            }
            else{
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvY(-i_y/q);
                }
            }
        }
        else {
            if ( abs(i_z) > 0.001){
                for ( int i = 0; i < q ; i++ ){
                    v[i].setvZ(-i_z/q);
                }
            }
        }
    }
}

istream& operator>>(istream& is, Particle& p ){
	double x,y,z;
	is>>x>>y>>z;
	p.setX(x);
	p.setY(y);
	p.setZ(z);
	return is;
}

ostream& operator <<(ostream& os, Particle& p) {
	os << p.getX() << " " << p.getY() <<" "<< p.getZ() <<endl;
	return os;
}