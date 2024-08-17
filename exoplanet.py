import requests
from bs4 import BeautifulSoup

class Planet:
    def __init__(self, soup):
        self.soup = soup
        self.mass, self.radium, self.volcanic_activity = self.get_planet_data()
        self.rotation_period, self.axis_stability, self.magnetic_field = self.get_dynamics_data()
        self.composition, self.radiation_levels = self.get_atmospheric_data()
        self.blackbody1, self.blackbody2, self.blackbody3, self.media = self_get_temperatura()
        self.estimated_orbital_period, self.orbital_type, self_eccentricity = self_orbital_data()
        self.tidally_locked_blackbody_T01, self.tidally_locked_blackbody_T02 = self_tidally_locked()
        self.water, self.o2 = self.get_essentials()
        self.co2, self.ch4 = self.get_gas_levels()
        self.stellar_activity, self.diostance_from_star, self.star_type, self.in_habitable_zone = self.get_stellar_data()

    def get_planet_data(self):  
        rows = soup.find_all('tr')

        mass = rows[4].find_all('td')[1].text.strip()
        radium = rows[3].find_all('td')[1].text.strip()
        volcanic_activity = soup.find('td', text='Volcanic Activity').find_next('td').text.strip()
        if not volcanic_activity:
            volcanic_activity = 'None'
            print('Desconocido')

        return  mass, radium, volcanic_activity

    # Datos Adicionales

    def dynamics_data(self):
        rotation_period = soup.find('td', text='Rotation Period').find_next('td').text.strip()
        if rotation_period:
            rotation_period
        else: 
            print('Desconocido')

        axis_stability = soup.find('td', text= 'Axis Stability').find_next('td').text.strip()
        if not axis_stability:
            axis_stability = 'None'
            print('Desconocido')
            
        magnetic_field = soup.find('td', text = 'Magnetic Field').find_next('td').text.strip()
        if not magnetic_field:
            magnetic_field = 'None'
            print('Desconocido')

        return rotation_period, axis_stability, magnetic_field;

    # Datos Atmosféricos
    def atmospheric_data(self):
        composition = soup.find('td',text='Size Class').find_next('td').text.strip()
        if composition == 'Super-Jupiter-size' or 'Jupiter-size':
            composition = 'Gas Planet'
        else:
            composition = 'Rocky-Planet'

        radiation_levels_td = soup.find('td', text = 'Star Radiation at Atmospheric Boundary').find_next('td').text.strip()
        radiation_levels = radiation_levels_td + ' (W/m2)'
        
        return composition, radiation_levels


    #Temperatura superficial del planeta en cuestión los datos son obtenido con Blackbody.           
    def Temperature(BlackBody1, BlackBody2, BlackBody3, media):
        BlackBody1 = soup.find('td', text='Blackbody T 0.1(K)').find_next('td').text.strip()




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
     

else:
    print('Error al acceder a la página')



                    