# Proyecto Final
## Implementación del Algoritmo de Needleman-Wunsch en C++

## Descripcion 
- Este proyecto implementa el algoritmo de Needleman-Wunsch en C++ para realizar el alineamiento global de secuencias. Es comúnmente utilizado en bioinformática para comparar secuencias de ADN, ARN o proteínas y encontrar el mejor alineamiento posible.

------------

## Requisitos previos
- Sistema operativo: Linux/Unix.
- Compilador C++: GCC, Clang, o MSVC con soporte para C++11 o superior.
- Graphviz: Tener instalado [Graphviz](http://graphviz.org/ "Graphviz").

------------
## Características

1. Generación de una representación visual del alineamiento mediante Graphviz.
   
2. Cálculo del puntaje óptimo del alineamiento para evaluar la calidad del resultado.
   
3. Opción para guardar los resultados en un archivo de texto.

------------


## Instalación
#### Clonar el repositorio
`$ git clone https://github.com/Kaisermannz/Proyecto_Final`

`$ cd Proyecto_Final `

#### Compilar el código
`$ g++ main.cpp -o needleman_wunsch`

------------

## Modos de Ejecución

El Programa Cuenta Con Dos Modos De Uso:
1. Modo de Comparación Entre Archivos **(-c)** : Permite comparar secuencias de dos archivos **.fna** distintos.

2. Modo de Comparación Interno **(-w)** : Permite comparar secuencias de un mismo archivo **.fna** .

- `$ ./needleman_wusch -c < ruta/archivo1.txt > < ruta/>archivo2.txt >`

- `$ ./needleman_wusch -w < ruta/archivo1.txt >`

> [!WARNING]
> Solo acepta formatos FASTA .fna

> [!TIP]
> En la carpeta Secuencias hay secuencias para poder probar.

------------

## Ejemplo de uso

**Compilar el código**

`$ g++ main.cpp -o needleman_wunsch`

**Ejecutar con dos secuencias distintas.**

`$ ./needleman -c secuencias/hemoglobina/HBB_hemoglobin_Homo_sapiens_human/ncbi_dataset/data/gene.fna secuencias/hemoglobina/ LOC116926515_Daphnia_magna/ncbi_dataset/data/gene.fna`

**Ejecutar comparando las secuencias de un archivo.**

`$ ./needleman -w secuencias/hemoglobina_mioglobina/HBB_hemoglobin_Homo_sapiens_human/ncbi_dataset/data/gene.fna`

------------


### Autor 
Cristobal Mauricio Araya Saldivia




