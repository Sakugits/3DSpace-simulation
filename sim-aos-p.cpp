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


struct object {
	double x; 
	double y;
	double z;
	double vx;
	double vy;
	double vz;
	double ax;
	double ay;
	double az;
	double m;
	int activo; //booleano para saber si hace falta actualizar su posicion
};



//Funciones

const double grav = 6.674e-11; //constante gravitacional para el calculo de fuerzas

void calc_grav(object &a, object &b) { //calcular fuerza gravitatoria dados 2 objetos y 3 variables en las que depositar las aceleraciones x,y,z

	double mod = grav * a.m * b.m; //aqui se aplica la formula de la fuerza ejercida por la gravedad
	double norma = std::sqrt((b.x - a.x) * (b.x - a.x) + (b.y - a.y) * (b.y - a.y) + (b.z - a.z) * (b.z - a.z));
	double norma3 = norma * norma * norma;
	//#pragma omp critical 
	//{
	#pragma omp atomic
	a.ax += ( ( mod * (b.x - a.x) ) / (norma3) ) / a.m; //x
	#pragma omp atomic
	a.ay += ( ( mod * (b.y - a.y) ) / (norma3) ) / a.m; //y
	#pragma omp atomic
	a.az += ( ( mod * (b.z - a.z) ) / (norma3) ) / a.m; //z
	
	#pragma omp atomic
	b.ax += - ( ( mod * (b.x - a.x) ) / (norma3) ) / b.m; //x
	#pragma omp atomic
	b.ay += - ( ( mod * (b.y - a.y) ) / (norma3) ) / b.m; //y
	#pragma omp atomic
	b.az += - ( ( mod * (b.z - a.z) ) / (norma3) ) / b.m; //z
	//}
	

}

void check_collisions (object &obj_a, object &obj_b){ //calcular colisiones dados 2 objetos, y fusionarlos si es el caso, devolver 1 si se chocan
	double distance;
	distance = std::sqrt(((obj_b.x - obj_a.x) * (obj_b.x - obj_a.x)) + ((obj_b.y - obj_a.y) * (obj_b.y - obj_a.y)) + ((obj_b.z - obj_a.z) * (obj_b.z - obj_a.z))); //distancia es raiz cuadrada de la suma de las diferencias elevadas al cuadrado de las posiciones en cada eje
	if (distance < 1 && obj_a.activo && obj_b.activo) { //comprobar si colisionan, hay que hacer la norma de las 3 dimensiones y ver si es menor que 1

		//#pragma omp atomic
		obj_a.vx = obj_a.vx + obj_b.vx; //actualizar velocidades y masa
		//#pragma omp atomic
		obj_a.vy = obj_a.vy + obj_b.vy;
		//#pragma omp atomic
		obj_a.vz = obj_a.vz + obj_b.vz;
		//#pragma omp atomic
		obj_a.m = obj_a.m + obj_b.m;

		obj_b.activo = 0; //desactivar objeto absorbido
	}

}

//Funcion principal

//los argumentos son (por orden) # de objetos, # de iteraciones, semilla, tamaño del cubo, tiempo por iteracion

int main (int argc, char* argv[] ) { //argv es el array de argumentos y el primer argumento es el nombre, argc es el # de args + 1

	//Comprobaciones
	
	if (argc != 6){// error de número de argumentos
		cerr << "ERROR: el rúmero de argumentos no es 5."<< endl;
		return -1;
	} 

	string num_objects_s; // transforma los parametros de string a int o float(segun el parametro) y comprueba que sean correctos.
	num_objects_s = argv[1];
	size_t index = 0;
	int num_objects = stoi(num_objects_s, &index);
	if (index != num_objects_s.length())
	{
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
	double size_enclosure = stod(size_enclosure_s, &index);
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
	double time_step = stod(time_step_s, &index);
	if (index != time_step_s.length()){
		cerr << "ERROR: Formato no válido."<< endl;
		return -2;
	}
	if (time_step<=0 ){
		cerr << "ERROR: el tiempo de cada iteración debe ser positivo."<< endl;
		return -2;
	}

	//ahora que los parámetros están bien escribimos el el fichero los necesarios
	ofstream iout("init_config2.txt");
	ofstream fout("final_config2.txt");
	ofstream bout("buffer_aos.txt");
	iout << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl; //linea inicial en el .txt con los parametros usados

	//Calculo del vector / objeto

	vector <object> objetos (num_objects);

	std::mt19937_64 gen(random_seed); //generador de numeros aleatorios
	std::uniform_real_distribution <double> real_rand{0,std::nextafter(size_enclosure,std::numeric_limits<double>::max())}; //posiciones aleatorias, limitadas al volumen del cubo
	double mean_eso = 1e21; //media
	double stddev_eso = 1e15; //desviacion tipica
	std::normal_distribution <> masas( mean_eso, stddev_eso); //distribucion normal usada para generar las masas

	//Calculo del vector / objeto

	for (int i = 0; i < num_objects; i++)
	{
		objetos[i].x = real_rand(gen); //aleatorio, confinado al volumen del cubo
		objetos[i].y = real_rand(gen);
		objetos[i].z = real_rand(gen);
		objetos[i].m = masas(gen); //aleatorio, distribucion normal
		objetos[i].vx = 0; //velocidades nulas
		objetos[i].vy = 0;
		objetos[i].vz = 0;
		objetos[i].ax = 0;
		objetos[i].ay = 0;
		objetos[i].az = 0;
		objetos[i].activo = 1; //todos los objetos activos al inicio
		iout << fixed << setprecision(3) << objetos[i].x << " " << objetos[i].y << " " << objetos[i].z << " " << objetos[i].vx << " " << objetos[i].vy << " " << objetos[i].vz << " " << objetos[i].m << endl; //escribir situacion inicial en init_config
		for (int k = 0; k < i; k++) //el 0 no comprueba, el 1 con el 0, el 2 con el 1 y el 0...
		{
		check_collisions(objetos[k], objetos[i]); //se comprueban colisiones al comienzo de la simulacion
		}
	}

	//Inicio Bucle de Cálculos
	int num_threads = omp_get_num_procs();

	for (int j = 0; j < num_iterations; j++) //bucle principal
	{
		#pragma omp parallel num_threads(num_threads) 
		{  
			#pragma omp for 
			for (int i = 0; i < num_objects-1; i++) //bucle de calcular fuerzas al inicio de cada iteracion
			{
				if (objetos[i].activo) //solo se calculan fuerzas con objetos activos
				{
					for (int k = i+1; k < num_objects; k++) //comprobar fuerzas por cada par de objetos
					{
						//un objeto no ejerce fuerza sobre si mismo
						calc_grav(objetos[i], objetos[k]); //calcular fuerza y aceleraciones por objeto

					}
				}
			}

		} //Barrera OpenMP

		#pragma omp parallel for num_threads(num_threads)
		for (int i = 0; i < num_objects; i++) //bucle de actualizar posiciones
		{
			if (objetos[i].activo) { //solo se actualizan las de objetos activos


				objetos[i].vx = objetos[i].vx + objetos[i].ax * time_step; //actualizar velocidades
				objetos[i].vy = objetos[i].vy + objetos[i].ay * time_step; //actualizar velocidades
				objetos[i].vz = objetos[i].vz + objetos[i].az * time_step; //actualizar velocidades

				
				objetos[i].ax = 0;
				objetos[i].ay = 0;
				objetos[i].az = 0;

				
				objetos[i].x = objetos[i].x + objetos[i].vx * time_step; //actualizar posicion
				
				if (objetos[i].x < 0) { //si se sale del cubo por limite negativo (0)
					objetos[i].x = 0; //volver a meter en el cubo
					objetos[i].vx = -objetos[i].vx; //invertir velocidad en el eje
				}
				else if (objetos[i].x > size_enclosure) { //si se sale del cubo por limite positivo (arista del cubo)
					objetos[i].x = size_enclosure; //meter en el cubo
					objetos[i].vx = -objetos[i].vx;
				}

				objetos[i].y = objetos[i].y + objetos[i].vy * time_step;
				
				if (objetos[i].y < 0) {
					objetos[i].y = 0;
					objetos[i].vy = -objetos[i].vy;
				}
				else if (objetos[i].y > size_enclosure) {
					objetos[i].y = size_enclosure;
					objetos[i].vy = -objetos[i].vy;
				}
			
				objetos[i].z = objetos[i].z + objetos[i].vz * time_step;
				
				if (objetos[i].z < 0) {
					objetos[i].z = 0;
					objetos[i].vz = -objetos[i].vz;
				}
				else if (objetos[i].z > size_enclosure) {
					objetos[i].z = size_enclosure;
					objetos[i].vz = -objetos[i].vz;
				}

				/*for (int k = 0; k < i; k++) { //calcular si hay colisiones al final de cada iteracion
				check_collisions(objetos[k], objetos[i]); // 0 no se comprueba, 1 compara con 0, 2 con 1 y con 0...
				}*/
			}
		}
		
		//#pragma omp parallel for num_threads(num_threads)
		for (int i = 0; i < num_objects; i++) {
			for (int k = 0; k < i; k++) { //calcular si hay colisiones al final de cada iteracion
				check_collisions(objetos[k], objetos[i]); // 0 no se comprueba, 1 compara con 0, 2 con 1 y con 0...
			}
		}
		
		// borramos los valores no activos al final de cada iteracion
		for (int i = 0; i < num_objects; i++) {

			if (!(objetos[i].activo)) {
				objetos.erase(objetos.begin()+i);
				num_objects--;
				i--; //habiendo borrado el objeto i, todo el vector de objetos se desplaza y el siguiente objeto tiene el indice i, y no i+1
			}
			
		}

	}

	// escribimos la primera linea de final config y el resto van a un buffer para evitar  que la primera linea no esté en primer lugar
	for (int i = 0; i < num_objects; i++) { //una vez terminadada la simulacion, imprimimos resultados
		if (i == num_objects -1){
			fout << fixed << setprecision(3) << size_enclosure << " " << time_step << " " << num_objects << endl;
		}
		if (objetos[i].activo) {
			//guardamos los datos finales en un buffer
			bout << fixed << setprecision(3) << objetos[i].x << " " << objetos[i].y << " " << objetos[i].z << " " << objetos[i].vx << " " << objetos[i].vy << " " << objetos[i].vz << " " << objetos[i].m << endl;
		}
	}
	// sacamos los datos del buffer y los metemos en final_config además cerramos los ofstream
	ifstream entrada;
	bout.close();
	iout.close();
	string texto;
	entrada.open("buffer_aos.txt",ios::in);
	
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