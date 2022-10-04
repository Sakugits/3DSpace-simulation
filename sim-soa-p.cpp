#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <random>
#include <functional>
#include <iomanip>
#include <omp.h>

using namespace std;

//gravedad
const double grav = 6.674e-11;

struct objetos 
{
	vector <double> x;
	vector <double> y;
	vector <double> z;
	vector <double> vx;
	vector <double> vy;
	vector <double> vz;
	vector <double> ax;
	vector <double> ay;
	vector <double> az;
	vector <double> m;
	vector <int> activo;
};


void check_collisions (objetos &obj, int obj1, int obj2){ //calcular colisiones dados 2 objetos, y fusionarlos si es el caso, devolver 1 si se chocan
	double distance;
	
	distance = std::sqrt((obj.x[obj2] - obj.x[obj1]) * (obj.x[obj2] - obj.x[obj1])+ (obj.y[obj2] - obj.y[obj1]) * (obj.y[obj2] - obj.y[obj1]) + (obj.z[obj2] - obj.z[obj1]) * (obj.z[obj2] - obj.z[obj1])); //distancia es raiz cuadrada de la suma de las diferencias elevadas al cuadrado de las posiciones en cada eje

	if (distance < 1 && obj.activo[obj1] && obj.activo[obj2]) //comprobar si colisionan, hay que hacer la norma de las 3 dimensiones y ver si es menor que 1
	{ 
		obj.vx[obj1] = obj.vx[obj1] + obj.vx[obj2]; //actualizar velocidades y masa
		obj.vy[obj1] = obj.vy[obj1] + obj.vy[obj2];
		obj.vz[obj1] = obj.vz[obj1] + obj.vz[obj2];
		obj.m[obj1] = obj.m[obj1] + obj.m[obj2];
        obj.activo[obj2] = 0;
	}
	
}

void calc_grav(objetos &object, int obj1, int obj2) { //calcular fuerza gravitatoria dados 2 objetos y 3 variables en las que depositar las aceleraciones x,y,z

	double mod = grav * object.m[obj1] * object.m[obj2];
	double norma = std::sqrt((object.x[obj2] - object.x[obj1]) * (object.x[obj2] - object.x[obj1]) + (object.y[obj2] - object.y[obj1]) * (object.y[obj2] - object.y[obj1]) + (object.z[obj2] - object.z[obj1]) * (object.z[obj2] - object.z[obj1]));
	double norma3 = norma * norma * norma;

	#pragma omp atomic
	object.ax[obj1] += (( mod * (object.x[obj2] - object.x[obj1]) ) / (norma3)) /object.m[obj1]; //x
	#pragma omp atomic
	object.ay[obj1] += (( mod * (object.y[obj2] - object.y[obj1]) ) / (norma3)) /object.m[obj1]; //y
	#pragma omp atomic
	object.az[obj1] += (( mod * (object.z[obj2] - object.z[obj1]) ) / (norma3)) /object.m[obj1]; //z
	#pragma omp atomic
	object.ax[obj2] += - (( mod * (object.x[obj2] - object.x[obj1]) ) / (norma3)) /object.m[obj2]; //x
	#pragma omp atomic
	object.ay[obj2] += - (( mod * (object.y[obj2] - object.y[obj1]) ) / (norma3)) /object.m[obj2]; //y
	#pragma omp atomic
	object.az[obj2] += - (( mod * (object.z[obj2] - object.z[obj1]) ) / (norma3)) /object.m[obj2]; //z
}

//los argumentos son (por orden) # de objetos, # de iteraciones, semilla, tamaño del cubo, tiempo por iteracion
int main (int argc, char* argv[] ) { //argv es el array de argumentos y el primer argumento es el nombre, argc es el # de args + 1

	//comprobamos y creamos variables para los distintos parametros
	if (argc != 6){// error de número de argumentos
		cerr << "ERROR: el rúmero de argumentos no es 5."<< endl;
		return -1;
	} 

	string num_objects_s; // transforma los parametros de string a int o float(segun el parametro) y comprueba que sean correctos.
	num_objects_s = argv[1];
	size_t index = 0;
	int num_objects = stoi(num_objects_s, &index);
	if (index != num_objects_s.length()){
		cerr << "ERROR: Formato no válido."<< endl;
		return -2;
	}

	if (num_objects<=0 ){
		cerr << "ERROR: el número de objetos debe ser un entero positivo."<< endl;
		return -2;
	}

	string num_iterations_s;
	num_iterations_s = argv[2];
	index = 0;
	int num_iterations = stoi(num_iterations_s, &index);

	if (index != num_iterations_s.length()){
		cerr << "ERROR: Formato no válido."<< endl;
		return -2;
	}

	if (num_iterations<=0 ){
		cerr << "ERROR: el número de interaciones debe ser un entero positivo."<< endl;
		return -2;
	}

	string random_seed_s;
	random_seed_s = argv[3];
	index = 0;
	int random_seed = stoi(random_seed_s, &index);

	if (index != random_seed_s.length()){
		cerr << "ERROR: Formato no válido."<< endl;
		return -2;
	}

	if (random_seed<=0 ){
		cerr << "ERROR: la semilla debe ser un entero positivo."<< endl;
		return -2;
	}

	string size_enclosure_s;
	size_enclosure_s = argv[4];
	index = 0;
	float size_enclosure = stof(size_enclosure_s, &index);

	if (index != size_enclosure_s.length()){
		cerr << "ERROR: Formato no válido."<< endl;
		return -2;
	}

	if (size_enclosure<=0 ){
		cerr << "ERROR: el tamaño del recinto debe ser positivo."<< endl;
		return -2;
	}

	string time_step_s;
	time_step_s =argv[5];
	float time_step = stof(time_step_s, &index);

	if (index != time_step_s.length()){
		cerr << "ERROR: Formato no válido."<< endl;
		return -2;
	}

	if (time_step<=0 ){
		cerr << "ERROR: el tiempo de cada iteración debe ser positivo."<< endl;
		return -2;
	}

	//ahora que los parámetros están bien escribimos el el fichero los necesarios
	ofstream fout("final_config2.txt");
    ofstream iout("init_config2.txt");
    ofstream bout("buffer_soa.txt");
	iout <<fixed << setprecision(3)<< size_enclosure << " " << time_step << " " << num_objects << endl;


	//Calculo del vector / objeto

	objetos objects;
	
	std::mt19937_64 gen(random_seed); //generador de numeros aleatorios
	std::uniform_real_distribution <double> real_rand{0,std::nextafter(size_enclosure,std::numeric_limits<double>::max())}; //posiciones aleatorias, limitadas al volumen del cubo
	double mean_eso = 1e21; //media
	double stddev_eso = 1e15; //desviacion tipica
	std::normal_distribution <> masas( mean_eso, stddev_eso); //distribucion normal usada para generar las masas
    
	//Calculo del objeto

	for (int i = 0; i < num_objects; i++)
	{
		objects.x.push_back(real_rand(gen)); //aleatorio, confinado al volumen del cubo
		objects.y.push_back(real_rand(gen)); //aleatorio
		objects.z.push_back(real_rand(gen)); //aleatorio
		objects.vx.push_back(0); //velocidades nulas
		objects.vy.push_back(0);
		objects.vz.push_back(0);
		objects.ax.push_back(0); //Aceleraciones nulas
		objects.ay.push_back(0);
		objects.az.push_back(0);
		objects.m.push_back(masas(gen)); //aleatorio, distribucion normal
        objects.activo.push_back(1);
        iout << fixed << setprecision(3) << objects.x[i]<< " " << objects.y[i] << " " << objects.z[i] << " " << objects.vx[i] << " " << objects.vy[i] << " " << objects.vz[i] << " " <<objects.m[i] << endl;
		for (int k = 0; k < i; k++) {//el 0 no comprueba, el 1 con el 0, el 2 con el 1 y el 0...
			check_collisions( objects, k, i); //calcular si hay colisiones al comienzo de la simulacion
		}
	}

	//Inicio Bucle de Cálculos
	int num_threads = omp_get_num_procs();

	for (int j = 0; j < num_iterations; j++) //bucle principal
	{
		#pragma omp parallel num_threads(num_threads) 
		{
			#pragma omp for
			for (int i = 0; i < num_objects-1; i++)  //bucle de calcular fuerzas al inicio de cada iteracion
			{
				if (objects.activo[i]) //solo se calculan fuerzas con objetos activos
				{
					for (int k = i+1; k < num_objects; k++) //comprobar fuerzas por cada par de objetos (inacabado, actualmente se compara cada par 2 veces
					{
						if (i != k && objects.activo[k]) { //un objeto no ejerce fuerza sobre si mismo 
							calc_grav(objects, i, k); //calcular fuerza y aceleraciones por objeto
						}
					}
				}
			}
		}
		
		#pragma omp parallel for num_threads(num_threads)
		for (int i = 0; i < num_objects; i++) //bucle de actualizar posiciones
		{
			if (objects.activo[i]) { //solo se actualizan las de objetos activos
				objects.vx[i] = objects.vx[i] + objects.ax[i]  * time_step; //actualizar velocidades
				objects.vy[i] = objects.vy[i] + objects.ay[i]  * time_step;
				objects.vz[i] = objects.vz[i] + objects.az[i]  * time_step;

				objects.ax[i] = 0;
				objects.ay[i] = 0;
				objects.az[i] = 0;

				objects.x[i] = objects.x[i] + objects.vx[i] * time_step; //actualizar posicion

				if (objects.x[i] < 0) { //si se sale del cubo por limite negativo (0)
					objects.x[i] = 0; //volver a meter en el cubo
					objects.vx[i] = -objects.vx[i]; //invertir velocidad en el eje
				}
				else if (objects.x[i] > size_enclosure) { //si se sale del cubo por limite positivo (arista del cubo)
					objects.x[i] = size_enclosure; //meter en el cubo
					objects.vx[i] = -objects.vx[i];
				}

				objects.y[i] = objects.y[i] + objects.vy[i] * time_step;
				
				if (objects.y[i] < 0) {
					objects.y[i] = 0;
					objects.vy[i] = -objects.vy[i];
				}
				else if (objects.y[i] > size_enclosure) {
					objects.y[i] = size_enclosure;
					objects.vy[i] = -objects.vy[i];
				}

			
				objects.z[i] = objects.z[i] + objects.vz[i] * time_step;
				
				if (objects.z[i] < 0) {
					objects.z[i] = 0;
					objects.vz[i] = -objects.vz[i];
				}
				else if (objects.z[i] > size_enclosure) {
					objects.z[i] = size_enclosure;
					objects.vz[i] = -objects.vz[i];
				}
			}
		}
		
		for(int i = 0; i < num_objects; i++){
			for (int k = 0; k < i; k++) { //calcular si hay colisiones al final de cada iteracion
					 check_collisions( objects, k, i); // 0 no se comprueba, 1 compara con 0, 2 con 1 y con 0...
				}
		}

		for (int i = 0; i < num_objects; i++) {
			if (!(objects.activo[i])) {
				//borramos objetos no activos
				objects.x.erase(objects.x.begin()+i);
				objects.y.erase(objects.y.begin()+i);
				objects.z.erase(objects.z.begin()+i);
				objects.vx.erase(objects.vx.begin()+i);
				objects.vy.erase(objects.vy.begin()+i);
				objects.vz.erase(objects.vz.begin()+i);
				objects.ax.erase(objects.ax.begin()+i);
				objects.ay.erase(objects.ay.begin()+i);
				objects.az.erase(objects.az.begin()+i);
				objects.m.erase(objects.m.begin()+i);
				objects.activo.erase(objects.activo.begin()+i);
				num_objects--;
				i--; //habiendo borrado el objeto i, todo el vector de objetos se desplaza y el siguiente objeto tiene el indice i, y no i+1
			}
		}
	}
	// escribimos la primera linea de final config y el resto van a un buffer para evitar  que la primera linea no esté en primer lugar
	for (int i = 0; i < num_objects; i++) {
		if (i == num_objects -1) {
			fout << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl; //aqui el total de objetos no se actualiza si hay añguna colision en la ultima iteracion
		}
		if (objects.activo[i]) {
			//guardamos los datos finales en un buffer
			bout << fixed << setprecision(3) << objects.x[i] << " " << objects.y[i] << " " << objects.z[i] << " " << objects.vx[i] << " " << objects.vy[i] << " " << objects.vz[i] << " " << objects.m[i] << endl;
		}
	}

	// escribimos las lineas del buffer en final config y cerramos las diferentes ofstream que habíamos abierto
	ifstream entrada;
	bout.close();
	iout.close();
	string texto;
	entrada.open("buffer_soa.txt",ios::in);
	
	while(!entrada.eof()){
		getline(entrada,texto);
		if (texto != ""){
			fout << texto<<endl;
		}
	}
	
	fout.close();
	entrada.close();
	return 0;
}