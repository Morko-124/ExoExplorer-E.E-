import requests
import math
from bs4 import BeautifulSoup

class Planet:
    def __init__(self, soup):
        self.soup = soup
        self.j_mass, self.e_mass, self.j_radius, self.e_radium, self.volcanic_activity = self.get_planet_data()
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

        j_mass = rows[4].find_all('td')[1].text.strip()
        j_radius = rows[3].find_all('td')[1].text.strip()
        e_mass = float(soup.find('td', text = '(Rjup)	').find_next('td').text.strip())
        e_radius = float(soup.find('td', text = '(Mjup)	').find_next('td').text.strip()) 
        
        volcanic_activity = soup.find('td', text='Volcanic Activity').find_next('td').text.strip()
        if not volcanic_activity:
            volcanic_activity = 'None'
            print('Desconocido')

        return  j_mass, e_mass, j_radius, e_radius, volcanic_activity

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
    # Datos Estelares
    star_name = soup.find('td', text='Star Name :').find_next('td').text.strip()
    
    # Tipo de Estrella
    star_type_complete = soup.find('td', text='Spectral type :').find_next('td').text.strip()
    star_type = star_type_complete[0]
    
    # Distancia desde la Estrella
    distance_from_star_n = round(float(soup.find('td', text='Semi Major Axis / Orbital Distance Calculated').find_next('td').text.strip()), 4)
    distance_from_star = f"{distance_from_star_n} UA."
    
    # Radiación de la Estrella
    star_radiation_n = round(float(soup.find('td', text='Star Radiation at Atmospheric Boundary:').find_next('td').text.strip()), 4)
    star_radiation = f"{star_radiation_n} (W/m²)"
    
    # Actividad Estelar
    if star_radiation_n <= 1:
        stellar_activity = 'Very Low'
    elif 1 < star_radiation_n <= 3:
        stellar_activity = 'Low'
    elif 3 < star_radiation_n <= 10:
        stellar_activity = 'Moderate'
    else:
        stellar_activity = 'High'


# --------------------------------------------------------------------------------------------------
    # Cálculo de la Zona Habitable
    ### Calculo de Zona Habitable se buscan las variables necesarias dentro de la tabla:
    # Stellar Radius
    # Stellar Mass
    # star_temperature
    # distance_from_star.

    # Constantes
    sigma = 5.67e-8  # Constante de Stefan-Boltzmann en W/m²K⁴
    R_sun = 6.96e8  # Radio del Sol en metros
    L_sun = 3.828e26  # Luminosidad del Sol en Watts

    # Variables
    stellar_radius = float(soup.find('td', text='Stellar Radius (Rsun) :').find_next('td').text.strip())
    star_temperature = float(soup.find('td', text='Temperature :').find_next('td').text.strip())

    # Radio de la Estrella en Metros
    stellar_radius_meters = stellar_radius * R_sun

    # Luminosidad usando la fórmula de Stefan-Boltzmann
    star_luminosity = 4 * math.pi * (stellar_radius_meters ** 2) * sigma * (star_temperature ** 4)

    # Luminosidad en unidades solares
    star_luminosity_sun = star_luminosity / L_sun

    # Límites de la zona habitable en UA
    inner_habitable_zone = math.sqrt(star_luminosity_sun / 1.1)
    outer_habitable_zone = math.sqrt(star_luminosity_sun / 0.53)

    # Determinación de si el exoplaneta está en la zona habitable
    in_habitable_zone = inner_habitable_zone <= distance_from_star_n <= outer_habitable_zone

    return star_name, star_type, distance_from_star, star_radiation, stellar_activity, in_habitable_zone;

# --------------------------------------------------------------------------------------------------
def get_tidally_locked(self, soup):
    # Constantes:
    jupiter_mass = 1
    jupiter_radius = 1
    jupiter_radius_m = 7.1492e7
    AU = 1.496e11  # Unidad Astronómica en metros
    jupiter_density = 1
    saturn_mass = 0.3
    saturn_radius = 0.9
    saturn_density = 0.7
    earth_mass = 0.01
    earth_radius = 0.1
    earth_density = 1.1
    super_earth_mass = 0.02
    super_earth_radius = 0.2
    G = 6.67430e-11  # Constante gravitacional en m^3 kg^−1 s^−2
    M_sun = 1.989e30  # Masa del Sol en kg

    # Datos del planeta
    mass_p = float(self.j_mass)  # Masa del planeta en masas de Júpiter
    radius_p = float(self.j_radius)  # Radio del planeta en radios de Júpiter
    eccentricity = float(self.eccentricity)  # Excentricidad orbital
    periodo_orbital = float(self.orbital_period)  # Período orbital en días
    semi_mayor_axis = float(self.distance_from_star.split()[0])  # Semi-eje mayor en UA
    rotation = float(self.rotation_period.split()[0])  # Período de rotación en horas

    # Datos estelares
    stellar_mass = float(soup.find('td', text='Stellar Mass (Msun) :').find_next('td').text.strip())
    stellar_radius = float(soup.find('td', text='Stellar Radius (Rsun) :').find_next('td').text.strip())
    stellar_temperature = float(soup.find('td', text='Temperature :').find_next('td').text.strip())

    # Conversiones
    j_m_k = 1.898 * 10**27  # Masa de Júpiter en kg
    j_r_m = 7.1482 * 10**7  # Radio de Júpiter en metros

    mass_p_k = mass_p * j_m_k
    radius_p_m = radius_p * j_r_m
    volumen = (4/3) * math.pi * radius_p_m**3
    density = mass_p_k / volumen

    m_star_kg = stellar_mass * M_sun

    # Calcular el tiempo sincrónico
    k = 1  # Valor estándar inicial
    t_sync = ((k * stellar_radius**2 * stellar_mass * periodo_orbital**2))**(1/3)

    # Ajuste del valor de k según la masa, radio y densidad del planeta
    if (mass_p > jupiter_mass and radius_p >= jupiter_radius) or (mass_p > jupiter_mass and density > jupiter_density) or (mass_p == jupiter_mass):
        k = 0.51
    elif (mass_p <= saturn_mass and radius_p > saturn_radius) or (density < saturn_density):
        k = 0.39
    elif (mass_p <= earth_mass and radius_p <= earth_radius and density <= earth_density):
        k = 0.30
    elif (mass_p > earth_mass and mass_p <= super_earth_mass and radius_p > earth_radius and radius_p <= super_earth_radius):
        k = 0.35
    else:
        k = 0.35

    t_sync = ((k * stellar_radius**2 * stellar_mass * periodo_orbital**2))**(1/3)

    if semi_mayor_axis <= t_sync:
        base = "Possible tidal locking"
    elif t_sync <= 10:
        base = "Within the gravitational influence of the star"
    elif t_sync <= 50:
        base = 'Gravitational influence of its star; moderate'
    else:
        base = 'Far from gravitational influence.'

    # Cálculo del factor Q utilizando el calentamiento mareal
    P_diss_min = 1e13  # Potencia disipada mínima estimada en vatios
    P_diss_max = 1e18  # Potencia disipada máxima estimada en vatios
    
    Q_min, Q_max = self.calcular_factor_Q(stellar_mass, radius_p, semi_mayor_axis, eccentricity, P_diss_min, P_diss_max)

    # Determinación de la probabilidad de anclaje por marea
    if Q_min > 1e6:
        probabilidad_anclaje = "Low probability of tidal locking"
    elif Q_max < 1e4:
        probabilidad_anclaje = "High probability of tidal locking"
    else:
        probabilidad_anclaje = "Medium probability of tidal locking"

    return base, probabilidad_anclaje

def calcular_factor_Q(self, M_star, R_planet, a, eccentricity, P_diss_min, P_diss_max):
    """ Calcula el factor de calidad Q usando el calentamiento mareal. """

    # Constantes
    G = 6.67430e-11  # Constante gravitacional en m^3 kg^−1 s^−2
    M_sun = 1.989e30  # Masa del Sol en kg
    R_jup = 7.1492e7  # Radio de Júpiter en metros
    AU = 1.496e11  # Unidad Astronómica en metros

    # Convertir a unidades correctas
    M_star_kg = M_star * M_sun  # Masa de la estrella en kg
    R_planet_m = R_planet * R_jup  # Radio del planeta en metros
    a_m = a * AU  # Semi-eje mayor en metros

    # Calcular el factor Q para el rango de potencia disipada
    Q_min = (63 / 4) * ((G * M_star_kg)**(3/2)) * (R_planet_m**5) * eccentricity / (P_diss_max * a_m**(15/2))
    Q_max = (63 / 4) * ((G * M_star_kg)**(3/2)) * (R_planet_m**5) * eccentricity / (P_diss_min * a_m**(15/2))

    # Interpretación de los resultados
    if Q_min < 10**5:
        q_min_interpretation = "High probability of tidal locking"
    elif 10**5 <= Q_min < 10**7:
        q_min_interpretation = "Medium probability of tidal locking"
    else:
        q_min_interpretation = "Low probability of tidal locking"

    if Q_max < 10**5:
        q_max_interpretation = "High probability of tidal locking"
    elif 10**5 <= Q_max < 10**7:
        q_max_interpretation = "Medium probability of tidal locking"
    else:
        q_max_interpretation = "Low probability of tidal locking"

    # Imprimir interpretaciones para el rango de Q
    print(f"Minimum Q factor: {Q_min} - Interpretation: {q_min_interpretation}")
    print(f"Maximum Q factor:  {Q_max} - Interpretation: {q_max_interpretation}")

    return Q_min, Q_max


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



                    