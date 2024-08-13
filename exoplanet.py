import requests
from bs4 import BeautifulSoup

# Búsqueda de los datos de Masa y Radio.
def Planet_Data(soup):  
    rows = soup.find_all('tr')

    # Dato de Masa (predeterminado C2 y F5)
    fila_5 = rows [4]
    columna_fila_5 = fila_5.find_all('td')
    Mass = columna_fila_5[1].text.strip()

    # Dato de Radio (predeterminado C2 y F4)
    fila_4 = rows [3]
    columna_fila_4 = fila_4.find_all('td')
    Radium = columna_fila_4[1].text.strip()

    # Datos del sistema volcánico
    VolcanicActivity = soup.find('td', text='Volcanic Activity').find_next('td').text.strip()
    if not VolcanicActivity:
        VolcanicActivity = 0
        print("No se han encontrado datos.") 

    
    return Planet_Data

# Datos Adicionales
def Planet_Data2(RotationPeriod, AxisStability):
    return;

# Datos Atmosféricos
def Atmospheric_Data(Composition, RadiationLevels):
    return;
#Temperatura superficial del planeta en cuestión los datos son obtenido con Blackbody.           
def Temperature(BlackBody1, BlackBody2, BlackBody3, media):
    return ;

# Datos Orbitales
def Orbital(Estimated, OrbitalType, Eccentricity):
    return ;

# Datos Esenciales para la vida.          
def Essentials(Water, O2):
    return ;

# Niveles de Gas.
def Gas_Levels(CO2, CH4):
    return ;

#Datos Estelares
def Stellar(Activity, DistanceFromStar, StarType, InHabitableZone):
    return;


# Unifica los def´s anteriores.
def Datos_Detallados():
    return;






#Llamo a todos los datos obtenidos por los def anteriores,
#los ordena en una tupla almacena en un SQLite
def Datos_Recopilados():
    return


url = 'https://www.exoplanetkyoto.org/exohtml/A_All_Exoplanets.html'
response = requests.get(url)


if response.status_code == 200:
    # Crea un objeto BeautifulSoup
    soup = BeautifulSoup(response.text, 'html.parser')

    # Búsqueda de "table" dentro de la página
    table = soup.find('table')
    # Filas de la tabla
    rows = table.find_all('tr')

# A continuacion se crea una funcion la cual tendrá como propósito poder extraer
#  los Datos Detallados de los ExoPlanetas.

    def PlanetList(link):   # Funcion DatosDetallados con un parámetro "link"
        url_detalles = url_detalles(link)    # Variable que almacena los detalles de los links.
        respuesta_detalles = requests.get(respuesta_detalles)   # Se le da la peticion al HTML con un requests.get y los datos se almacenan dentro de respuesta_detalles
        soup_detalles = BeautifulSoup(respuesta_detalles.text, 'html.parser')  # Se genera una variable soup que almacenara los datos de respuesta_detalles y se trabaja el HTML como parser.
        
        # Se crea una variable que almacena los datos de la tabla.
        tabla_detalles = soup_detalles.find('table')    #soup_detalles.find('value') Busca dentro del HTML el parámetro que necesitemos en este caso un "table"
        datos_detallados = []  #Almacena los datos dentro de una tupla.

        if tabla_detalles:      # Si tabla_detalles es encontrada entonces
            filas_detalles = tabla_detalles.find_all('tr') # Busca dentro del HTML todos los "tr" que son las filas.
            for fila in filas_detalles:          # Para cada fila en filas_detalles, haz lo siguiente:
                celdas = fila.find_all('td')         # Encuentra todas las celdas "td" dentro de la fila.
                if len(celdas)>2:               # Si hay más de 2 celdas, procede.
                    columna_1 = celdas[0].text.strip()      # Obtén el texto de la primera celda y quita los espacios.
                    columna_2 = celdas[1].text.strip()      # Haz lo mismo para la segunda celda.
                    datos_detallados.append((columna_1, columna_2))      # Guarda los datos en la lista.    
        return datos_detallados

    # Imprimir las primeras 10 filas
    for index, row in enumerate(rows[4:5], 1):
        # Encuentra todas las celdas en la fila
        cells = row.find_all('td')
        # Extrae el texto de cada celda
        cell_data = [cell.get_text(strip=True) for cell in cells[:2]]
        # Imprime los datos de la fila
        
        
        # Rearmado de links
        link = cells[1].find('a')['href'] #Genero una variable que analiza los datos de celda 1 y busca href que el faltante del link.
        link = link.lstrip('./')
        base_url =  'https://www.exoplanetkyoto.org/exohtml'
        url_com = base_url + '/' +  link

        print(f"Fila {index}: {cell_data}: {url_com}")
        print(Planet_Data(base_url))



else:
    print('Error al acceder a la página')



                    