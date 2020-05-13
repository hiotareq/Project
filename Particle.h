#ifndef Particle_h
#define Particle_h

#include <iostream>
#include <vector>

using namespace std;

class Particle {
private:

	double x;
	double y;
	double z;
	double vX;
	double vY;
	double vZ;
	double aX;
	double aY;
	double aZ;
public:
	Particle(double new_x, double new_y, double new_z);
	Particle();
    ~Particle();

	double getX() const;
	double getY() const;
	double getZ() const;
    double getVX() const;
	double getVY() const;
	double getVZ() const;
	double getaX() const;
	double getaY() const;
	double getaZ() const;

	void setX(const double &new_x);
	void setY(const double &new_y);
	void setZ(const double &new_z);
	void addaX(const double &new_aX);
	void addaY(const double &new_aY);
	void addaZ(const double &new_aZ);
    void setvX(const double &aX);
	void setvY(const double &aY);
	void setvZ(const double &aZ);
	void setaX(const double &newaX);
    void setaY(const double &newaX);
    void setaZ(const double &newaX);

	double getDistance(const Particle& p, const double &l) const;
	double getPotential(const Particle& p,  const  double &sigma, const  double &e, const double &l) const;
	double getForce(const Particle& p, const double &sigma, const double &e, const double &length_of_cube) const;

	void setParticle( double const &new_x,double const &new_y, double const &new_z );
};

void impulse_check( vector<Particle> &v, const double &mas_part);
void Step( vector<Particle> &v,const float &t, const double &s, const double &e, const double &l, const double &mo_mass);

istream& operator>>(istream& is, Particle& p );
ostream& operator <<(ostream& os, Particle& p);
#endif