import requests
import math
from bs4 import BeautifulSoup

class Planet:
    def __init__(self, soup):
        self.soup = soup
        self.mass, self.radium, self.volcanic_activity = self.get_planet_data()
        self.rotation_period, self.axis_stability, self.magnetic_field = self.get_dynamics_data()
        self.composition, self.radiation_levels = self.get_atmospheric_data()
        self.blackbody1, self.blackbody2, self.blackbody3, self.media = self.get_temperature()
        self.observed_orbital_period, self.estimated_orbital_period, self.orbital_type, self.eccentricity = self.get_orbital_data()
        self.tidally_locked_blackbody_T01, self.tidally_locked_blackbody_T02 = self.get_tidally_locked()
        self.water, self.o2 = self.get_essentials()
        self.co2, self.ch4 = self.get_gas_levels()
        self.name, self.star_type, self.distance_from_star, self.star_radiation, self.stellar_activity, self.in_habitable_zone = self.get_stellar_data()

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

    def get_dynamics_data(self):
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
    def get_atmospheric_data(self):
        composition = soup.find('td',text='Size Class').find_next('td').text.strip()
        if composition == 'Super-Jupiter-size' or 'Jupiter-size':
            composition = 'Gas Planet'
        else:
            composition = 'Rocky-Planet'

        radiation_levels_td = soup.find('td', text = 'Star Radiation at Atmospheric Boundary').find_next('td').text.strip()
        radiation_levels = radiation_levels_td + ' (W/m2)'
        
        return composition, radiation_levels


    #Temperatura superficial del planeta en cuestión los datos son obtenido con Blackbody.           
    def get_temperature(self):
        blackBody1 = float(soup.find('td', text='Blackbody T 0.1(K)').find_next('td').text.strip())
        blackBody2 = float(soup.find('td', text = "Blackbody T 0.3(K)").find_next('td').text.strip())
        blackBody3 = float(soup.find('td', text = 'Blackbody T 0.7(K)').find_next('td').text.strip())
        media = (blackBody1 + blackBody2 + blackBody3)/3
        media_s = media + 'K'
        return media_s

    # Datos Orbitales
    def get_orbital_data(self):
        orbital_period = float(soup.find('td', text ='Orbital Period(Observed/Estimated)'.find_next('td')).text.strip())
        if orbital_period == '0.0':
            estimated_orbital_period_days = float(soup.find('td', text = '(Days Obs.)'.find_next('td')).text.strip())
            estimated_orbital_period = estimated_orbital_period_days/365.2425
            estimated_orbital_period = round(estimated_orbital_period,4)
            orbital_period = estimated_orbital_period + ('Estimado')
        elif orbital_period != 0.0:
            new_orbital_period = orbital_period / 365.2425
            orbital_period = new_orbital_period + (' Años')

        eccentricity = soup.find('td', text ='Eccentricity :'.find_next('td')).text.strip()
        if 0.0 <= eccentricity <= 0.1:
            orbital_type = 'Casi Circular'
        elif 0.1 <= eccentricity <= 0.5:
            orbital_type = 'Moderado'
        elif 0.5 <= eccentricity <= 1:
            orbital_type = 'Muy Elíptica'
        elif eccentricity >= 1:
            orbital_type = 'Hiperbólica'

        return orbital_period, eccentricity, orbital_type

    #Datos Estelares
    def get_stellar_data(soup):
        star_name = soup.find('td', text = 'Star Name :	'.find_next('td')).text.strip()
        #Tipo de Estrella.
        star_type_complete = soup.find('td', text = 'Spectral type :'.find_text('td')).text.strip()
        star_type = star_type_complete[0]
        
        distance_from_star_n = round(float(soup.find('td', text='Semi Major Axis / Orbital Distance Calculated	').find_next('td').text.strip()), 4)
        distance_from_star = distance_from_star_n + " UA."

        # Radiacion de Estrella>
        star_radiation_n = round(float(soup.find('td', text = 'Star Radiation at Atmospheric Boundary:').find_next('td').text.strip()),4)
        star_radiation = star_radiation_n + ' (W/m2)'


        if 1 > star_radiation_n:
            stellar_activity = 'Very Low'
        elif 1 < star_radiation_n < 3:
            stellar_activity = 'Low'
        elif 3 < stellar_activity < 10:
            stellar_activity = 'Moderate'
        elif stellar_activity > 10:
            stellar_activity = 'High'

# --------------------------------------------------------------------------------------------------

        ### Calculo de Zona Habitable se buscan las variables necesarias dentro de la tabla:
        # Stellar Radius
        # Stellar Mass
        # star_temperature
        # distance_from_star.

        # Constantes
        sigma = 5.67e-8  # Constante de Stefan-Boltzmann en W/m^2K^4
        R_sun = 6.96e8  # Radio del sol en metros
        L_sun = 3.828e26  # Luminosidad del sol en Watts

        # Variables
        stellar_radius = float(soup.find('td', text = 'Stellar Radius(Rsun) :'.find_next('td')).text.strip())
        stellar_mass = float(soup.find('td', text = 'Stellar Mass(Msun) :'.find_next('td')).text.strip())
        star_temperature = float(soup.find('td', text = 'Temperature :'.find_next('td')).text.strip())

        # Radio de la Estrella a Metros
        stellar_radius_metters = stellar_radius * R_sun

        # Luminosidad Stefan-Boltzmann.
        star_luminocity = 4 * math.pi * (stellar_radius_metters ** 2) * sigma * (star_temperature ** 4)

        # Luminosidad a unidades solares.
        star_luminocity_sun = star_luminocity/L_sun

        # Limites de zona habitable en UA
        inner_habitable_zone = math.sqrt(star_luminocity_sun/1.1)
        outer_habitable_zone = math.sqrt(star_luminocity_sun/0.53)

        if inner_habitable_zone <= distance_from_star_n <= outer_habitable_zone:
            in_habitable_zone = True
        else: 
            in_habitable_zone = False
       
        return star_name, star_type, distance_from_star, star_radiation, stellar_activity, in_habitable_zone

# --------------------------------------------------------------------------------------------------

    def get_tidally_locked(soup):
        return;
    # Datos Esenciales para la vida.          
    def Essentials(Water, O2):
        return ;

    # Niveles de Gas.
    def Gas_Levels(CO2, CH4):
        return ;
    

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



                    