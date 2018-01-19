/******************************************************************************/
/* Dinámica molecular velocity-Verlet 1 especie                               */
/* Estructura FCC                                                             */
/* Potencial de Morse                                                         */
/* Versión 1                                                                  */
/*                                                                            */
/* Calcula:                                                                   */
/* g(r)  Función de distribución radial                                       */
/*                                                                            */
/* Parámetros de simulación                                                   */
/* Pi           Valor del número pi                                           */
/* Ncx          Número de celdas en x                                         */
/* Ncy          Número de celdas en y                                         */
/* Ncz          Número de celdas en z                                         */
/* frand()      Generador de números aleatorios                               */
/* maxint       Número máximo para el arreglo de la g(r)                      */
/* Temp         Temperatura reducida                                          */
/* dt           Paso de integración                                           */
/* dr           Intervalo de muestreo                                         */
/* intstep      Número de pasos de integración                                */
/* ncg          Frecuencia de cálculo de la g(r)                              */
/* rprint       Frecuencia de impresión                                       */
/* rprintf      Frecuencia de impresión en el archivo                         */
/* a            Parametro de red                                              */
/* b            Parametro de red                                              */
/* c            Parametro de red                                              */
/* U0           Fondo del pozo de potencial                                   */
/* alpha        Parámetro del potencial de Morse                              */
/* Rij_0        Distancia del mínimo                                          */
/* epsilon      Profundidad del pozo                                          */
/* Npt          Número de partículas                                          */
/* intstepinit  Número para iniciar la estadística                            */
/******************************************************************************/
