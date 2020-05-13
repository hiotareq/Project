#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include "Particle.h"
#include "PartDist.h"
#include <iomanip>

using namespace std;

//26 858 747 789 600 000 000 000 000
//столько атомов аргона в одном кубическом метре
//кубометр делим на 26 858 747 789 600 000 000 000‬
//получаем нужную нам часть объёма, в которой будут "жить" 1000 частиц твёрдого аргона
//плотность взята при нормальных условиях
//одна 26 858 747 789 600 000 000 000‬ -я часть кубического метра
//3.34*10^(-8) (м) - длина стороны рассматриваемого куба (куб для 1000 частиц)



int main() {
    cout << setprecision(3);
	vector<Particle> particles;
	vector<double> Forces;
	vector<double> Potentials;

	unsigned int n;
	double x,sigma = 3.623,V, e = 0.154,T,den;
	double MolecMass=1, mo_mass=1;  //mo_mass - масса одной частицы
	cout<< "Insert number of particles"<<endl;
	cin>>n;
	particles.resize(n);
	Forces.resize((n*n-n)/2);
	Potentials.resize((n*n-n)/2);
	cout << "Insert temperature"<<endl;//температура с графика, T*=kT/e
	cin>>T;
	cout<<"Insert density"<<endl;//плотность тоже с графика
	cin>>den;
	V=(sigma*sigma*sigma*n)/den;
	x=cbrt(V);
    PartDist a;

    //распределение частиц по кубу ( актуально для твёрдого тела )
    //a.setCube(particles, x);

	/* ifstream f("Src.txt");//считываем координаты из файла (актуально для "локальных" исследований)
	if (!(f.is_open()))  
		cout << "ERROR: file 'Src.txt' wasn't open for reading";
	else
	{
		unsigned int i=0;
		while (!f.eof())
		{
			double x0,y0,z0;
			f >> x0>>y0>>z0;
			particles[i].setParticle(x0 ,y0 , z0);
			++i;
		}
	}
	 */

	//задание случайных координат (актуально для газов и (частично) жидкостей)
    a.setRand( particles, x, sigma);

	normal_distribution<> d(0,1);
	//задание скоростей с помощью распределению Максвелла
	//по проекции скорости
	a.setSpeed(particles, MolecMass, T, e);

	//проверка скоростей (импульс системы равен нулю)
	impulse_check(particles, mo_mass);

	//загружаем скорости в файлы -> можем проверить нормальность распределения
	ofstream srcx;
	srcx.open(R"({PROJECTFOLDER}\vx.txt)");
	if ( srcx.is_open()){
	    for ( auto i=0; i < n; i++)
	    srcx << particles[i].getVX()<<endl;
    }

	ofstream srcy;
    srcy.open(R"({PROJECTFOLDER}\vy.txt)");
    if ( srcy.is_open()) {
        for (auto i = 0; i < n; i++)
            srcy << particles[i].getVY()<<endl;
    }

    ofstream srcz;
    srcz.open(R"({PROJECTFOLDER}}\vz.txt)");
    if ( srcz.is_open()){
        for ( auto i=0; i < n; i++)
            srcz << particles[i].getVZ()<<endl;
    }


//основной эксперимент: вычисление сил и потенциалов, 
//вычисление ускорений и скоростей
//задание новых координат
ofstream srcxyz;
    srcxyz.open(R"({PROJECTFOLDER}\1234.xyz)");//файл, для загрузки в визуализатор
for (auto k=0;k<1800;++k){
Step(particles, 0.01, sigma, e, x, mo_mass);
srcxyz<<n<<endl<<k<<endl;                                               //формат записи информации о частицах,
for ( auto i = 0; i < n; i++) srcxyz <<i + 1<<" "<< particles[i];       //"читаемый" визуализатором
}

	for (unsigned int i=0;i<n;i++){
		cout<<particles[i]<<endl;
	}

	return 0;
}