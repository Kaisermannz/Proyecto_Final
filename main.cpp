#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

using namespace std;

// Definición de valores del algoritmo
const int Valor_GAP = -2;
const int Valor_Match = 1;
const int Valor_Miss_Match = -1;
const int Valor_Similar_Match = 0;

int valor_optimo = 0;

string secuencia_1 = "";
string secuencia_2 = "";

string alineamento_sec_1 = "";
string alineamento_sec_2 = "";

string nombre_archivo_1 = "";
string nombre_archivo_2 = "";

// Implementacion del algoritmo Needleman-Wunsch
void needlemanWunsch() {
  int n = secuencia_2.length();
  int m = secuencia_1.length();

  vector<vector<int>> matriz(n + 1, vector<int>(m + 1, 0));

  // Inicialización de la matriz con los valores de gaps
  for (int i = 0; i <= n; i++) {
    matriz[i][0] = i * Valor_GAP;
  }
  for (int j = 0; j <= m; j++) {
    matriz[0][j] = j * Valor_GAP;
  }

  // Llenado de la matriz
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      int top = matriz[i - 1][j] + Valor_GAP;
      int left = matriz[i][j - 1] + Valor_GAP;

      int diagonal;
      // Comparación de nucleótidos
      if (secuencia_2[i - 1] == secuencia_1[j - 1]) {
        diagonal = matriz[i - 1][j - 1] + Valor_Match;
      } else if ((secuencia_2[i - 1] == 'G' && secuencia_1[j - 1] == 'C') ||
                 (secuencia_2[i - 1] == 'C' && secuencia_1[j - 1] == 'G')) {
        diagonal =
            matriz[i - 1][j - 1] + Valor_Similar_Match; // Guanina y Citosina
      } else if ((secuencia_2[i - 1] == 'A' && secuencia_1[j - 1] == 'T') ||
                 (secuencia_2[i - 1] == 'T' && secuencia_1[j - 1] == 'A')) {
        diagonal =
            matriz[i - 1][j - 1] + Valor_Similar_Match; // Adenina y Timina
      } else {
        diagonal = matriz[i - 1][j - 1] + Valor_Miss_Match;
      }

      matriz[i][j] = max({top, left, diagonal});
    }
  }

  // el valor optimo es el ultimo de la matriz
  valor_optimo = matriz[n][m];

  int i = n;
  int j = m;

  while (i > 0 || j > 0) {
    // Si hace Match se guarda un caracter en cada secuencia
    if (i > 0 && j > 0 && (secuencia_2[i - 1] == secuencia_1[j - 1])) {
      alineamento_sec_1 = secuencia_1[j - 1] + alineamento_sec_1;
      alineamento_sec_2 = secuencia_2[i - 1] + alineamento_sec_2;
      i--;
      j--;
      // Comparamos el valor de top con un gap si es igual al valor actual
      // significa que es el mayor
    } else if (i > 0 && matriz[i][j] == matriz[i - 1][j] + Valor_GAP) {
      alineamento_sec_1 = "-" + alineamento_sec_1;
      alineamento_sec_2 = secuencia_2[i - 1] + alineamento_sec_2;
      i--;
      // Al ser el que queda es el izquierdo
    } else {
      alineamento_sec_1 = secuencia_1[j - 1] + alineamento_sec_1;
      alineamento_sec_2 = "-" + alineamento_sec_2;
      j--;
    }
  }

  // Se rellenan las secuencias si una es mas grande que la otra
  while (i > 0) {
    alineamento_sec_1 = "-" + alineamento_sec_1;
    alineamento_sec_2 = secuencia_2[i - 1] + alineamento_sec_2;
    i--;
  }

  while (j > 0) {
    alineamento_sec_1 = secuencia_1[j - 1] + alineamento_sec_1;
    alineamento_sec_2 = "-" + alineamento_sec_2;
    j--;
  }
}

// Guarda las las cabeceras de las secuencias que se encuentran en el archivo
// que anteriormente se pidio, luego le pide al usuario que elija una para
// guardarla en la variable de secuencia_1
void procesar_archivos(string nombre_archivo, int secuencia) {
  // Abre el archivo si no lo logra arroja error
  ifstream archivo(nombre_archivo);
  if (!archivo.is_open()) {
    cerr << "No se pudo abrir el archivo." << endl;
  }

  string linea;
  string cabeceras;

  int contador = 1;
  int respuesta;

  // Busca las cabeceras de las secuencias en los archivos identificando ">" y
  // las guardarlos
  while (getline(archivo, linea)) {
    if (linea[0] == '>') {
      cabeceras += to_string(contador) + " : " + linea + "\n";
      contador++;
    }
  }

  archivo.clear();
  archivo.seekg(0);

  // El usuario debe elegir que secuencias alinear
  while (true) {
    cout << cabeceras;
    cout << "Elija que secuencia quiere alinear: ";
    cin >> respuesta;
    if (respuesta <= contador && respuesta >= 1) {
      break;
    } else {
      cout << "Ingrese un valor dentro del rango.";
    }
  }

  // Se busca la secuencia elegida y se guarda
  int contador_2 = 0;
  while (getline(archivo, linea)) {
    if (linea[0] == '>') {
      contador_2++;
    } else if (contador_2 == respuesta) {
      if (secuencia == 1) {
        secuencia_1 += linea;
      } else if (secuencia == 2) {
        secuencia_2 += linea;
      }
    }
  }
  archivo.close();
}

// Funcion que toma el resultado del alineamiento y le da formato para que el
// txt sea mas legible
void guardar_resultado() {
  int ancho = 50;
  // Se guardan las secuencias en el txt dandoles un ancho de 50 caracteres y
  // mostrando cual secuencia es cual
  ofstream archivo("Alinamiento.txt");
  if (!archivo.is_open()) {
    cerr << "Error al abrir el archivo Alinamiento.txt" << endl;
    return;
  }

  long long longitud = max(secuencia_1.length(), secuencia_2.length());

  archivo << "Alineamiento: \n\n";

  for (long long i = 0; i < longitud; i += ancho) {

    string fragmento1 = "";
    string fragmento2 = "";

    // Obtener fragmento de secuencia_1 si está dentro del rango
    if (i < secuencia_1.length()) {
      fragmento1 = alineamento_sec_1.substr(i, ancho);
    }

    // Obtener fragmento de secuencia_2 si está dentro del rango
    if (i < secuencia_2.length()) {
      fragmento2 = alineamento_sec_2.substr(i, ancho);
    }

    archivo << "Secuencia 1: " << fragmento1 << "\n";
    archivo << "Secuencia 2: " << fragmento2 << "\n\n";
  }
  archivo.close();
  cout << endl;

  cout << "Archivo guardado en Alineamiento.txt con exito.";
  cout << endl;
}

void generar_graphviz() {
  int start = 0;
  int end = 0;

  // se le pide al usuario que ingrse un rango para crear la imagen del graphviz
  // por que si no seria demaciado grande
  while (true) {
    cout << "En que intervalo decea generar la imagen (ejemplo 200 300) no "
            "puede ser mayor a 100 ya que perderia mucha definicion y no seria "
            "viable: ";
    cin >> start >> end;
    if (end - start <= 100 && end - start > 0) {
      break;
    } else {
      cout << "Ingrese valores dentro del rango permitido." << endl;
    }
  }

  ofstream dotFile("alineamiento.dot");

  if (!dotFile.is_open()) {
    cerr << "Error al abrir el archivo DOT." << endl;
    return;
  }

  dotFile << "digraph G {\n";

  // Guarda en el graphviz la parte de la secuencia seleccionada indicando si
  // fue match o gap
  for (size_t k = start; k < end && k < alineamento_sec_1.length(); k++) {
    char nucleoA = alineamento_sec_1[k];
    char nucleoB = alineamento_sec_2[k];

    dotFile << "   A" << k << " [label=\"" << nucleoA << "\"];\n";
    dotFile << "   B" << k << " [label=\"" << nucleoB << "\"];\n";

    if (nucleoA != '-' && nucleoB != '-') {
      dotFile << "   A" << k << " -> B" << k << " [label=\"Match\"];\n";
    } else if (nucleoA == '-' || nucleoB == '-') {
      dotFile << "   A" << k << " -> B" << k << " [label=\"Gap\"];\n";
    }
  }

  dotFile << "}\n";
  dotFile.close();
  cout << endl;

  cout << "Archivo DOT guardado como 'alineamiento.dot'." << endl;
  cout << endl;

  system("dot -Tpng alineamiento.dot -o alineamiento.png");
  system("eog alineamiento.png&");
}

int main(int argc, char *argv[]) {

  // elegir_secuencia();

  if (argc < 3) {
    cerr << "Use: " << argv[0] << " -c <archivo1> <archivo2>\n";
    return 1;
  }

  if (string(argv[1]) == "-c" && argc == 4) {
    // Se guardan los 2 archivos para comparar
    nombre_archivo_1 = argv[2];
    nombre_archivo_2 = argv[3];
  } else if (string(argv[1]) == "-w" && argc == 3) {
    //  Se compara dentro del mismo archivo
    nombre_archivo_1 = argv[2];
    nombre_archivo_2 = argv[2];
  } else {
    cout << "Opcion invalida." << endl;
    return 1;
  }

  procesar_archivos(nombre_archivo_1, 1);
  procesar_archivos(nombre_archivo_2, 2);

  needlemanWunsch();

  cout << endl;
  cout << "El alineamiento ya fue realizado." << endl;
  cout << endl;

  while (true) {
    int respuesta;
    cout << "Opciones: " << endl;
    cout
        << "1) Generar una imagen con graphviz de una seccion de las secuencias"
        << endl;
    cout << "2) Guardar resultado del alineamiento en un txt." << endl;
    cout << "3) Calcular el valor optimo del alineamiento." << endl;
    cout << "4) Salir." << endl;
    cin >> respuesta;

    switch (respuesta) {
    case 1: {
      generar_graphviz();
      break;
    }
    case 2: {
      guardar_resultado();
      break;
    }
    case 3: {
      cout << "El valor optimo del alineamiento es : " << valor_optimo << endl;
      break;
    }
    case 4: {
      return 0;
    }
    default:
      cout << "Opción inválida\n";
    }
  }
}
