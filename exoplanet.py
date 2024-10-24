import requests
from bs4 import BeautifulSoup
import json
import math


class Planet:
    def __init__(self, soup):
        self.soup = soup
        self.name = self.get_name_data()
        self.j_mass, self.e_mass, self.j_radius, self.e_radius = self.get_planet_data()
        self.rotation_period, self.axis_stability = self.get_dynamics_data()
        self.sizeClass, self.composition, self.radiation_levels = self.get_atmospheric_data()
        self.blackbody1, self.blackbody2, self.blackbody3, self.media = self.get_temperature()
        self.orbital_period, self.orbital_period_observed, self.orbital_period_estimated, self.eccentricity, self.orbital_type = self.get_orbital_data()
        self.star_name, self.star_type, self.distance_from_star, self.star_radiation_atmospheric_b, self.stellar_activity, self.in_habitable_zone = self.get_stellar_data()
        self.P_diss_min = 1e13  # Potencia disipada mínima estimada
        self.P_diss_max = 1e18  # Potencia disipada máxima estimada
        self.t_sync, self.base, self.probabilidad_anclaje = self.get_tidally_locked(soup)
        self.Q_min, self.q_min_interpretation, self.Q_max, self.q_max_interpretation = self.calcular_factor_Q()
        self.oxygen, self.gravity, self.magnetic_field, self.atmospheric_pressure, self.volcanic_activity, self.has_water = self.get_essentials_data()



# ---------
#   Extracción de datos, se utilizó en ciertos momentos del código este método, pero en algunas extracciones
#   presentaba fallos o datos incongruentes.
# ---------

    def get_data(self):
        tables = self.soup.find_all('table')
        print(f"Found {len(tables)} tables.")
        if len(tables) > 1:
            table = tables[1]
            cells = table.find_all('td')
            print(f"Found {len(cells)} cells in the table.")
        return cells if len(tables) > 1 else []
    

# ---------
# Inicio de extracción de datos.
#
# ---------


    def get_name_data(self):
        tables = self.soup.find_all('table')
        
        # Inicialización
        self.name = None
         
        if len(tables) > 1:
            table = tables[1]
            cells = table.find_all('td')

            if len(cells) > 2:
                c1 = cells[1].text.strip()
                self.name = c1

        return self.name
    

    # Datos planetarios primarios, captura de masa (Jupiters y Tierras) y radio (Jupiters y Tierras.)
    def get_planet_data(self):  
        tables = self.soup.find_all('table')    
        self.j_mass, self.e_mass, self.j_radius, self.e_radius = (None, None, None, None )

    # Posicionamiento en la segunda tabla disponible.
        if len(tables) > 1:
            cells = self.get_data()

            for i, cell in enumerate(cells):
                call_text = cell.get_text(strip=True)
                # Dato de Masa en Jupiters.
                try:
                    if 'Mass   (Mjup / Me):' in call_text:
                        if i+1 < len(cells):
                            adjacent_cell = cells[i+1]
                            self.j_mass = float(adjacent_cell.get_text(strip=True))
                except:
                    self.j_mass = 'Unknown'

                # Dato de Masa en Tierras.
                try:
                    if 'Mearth' in call_text:
                        if i + 1 < len(cells):
                            adjacent_cell = cells [i - 1]
                            self.e_mass = float(adjacent_cell.get_text(strip=True))
                except:
                    self.e_mass = 'Unknown'

                # Dato de radio en Jupiters.
                try:
                    if 'Radius (Rjup / Re):' in call_text:
                        if i + 1 < len(cells):
                            adjacent_cell = cells [i + 1]
                            self.j_radius = float(adjacent_cell.get_text(strip = True))
                except:
                    self.j_radius = 'Unknown'

                # Dato de radio en Tierras.
                try: 
                    if 'Rearth' in call_text:
                        if i + 1 < len(cells):
                            adjacent_cell = cells [i - 1]
                            self.e_radius = float(adjacent_cell.get_text(strip = True))
                except:
                    self.e_radius = 'Unknown'

        return self.j_mass, self.e_mass, self.j_radius, self.e_radius

    def get_dynamics_data(self):
        tables = self.soup.find_all('table')    
        self.rotation_period = None
        self.axis_stability = None

        if len(tables) > 1:
            table = tables [1]
            cells = table.find_all('td')

            for i, cell in enumerate (cells):
                cell_text = cell.get_text(strip=True)

            try:
                if 'Rotation Period' in cell_text:
                    if i + 1 < len(cells):
                        adjacent_cells = cells [i + 1]
                        self.rotation_period = float(adjacent_cells.get_text(strip = True))
            except:
                self.rotation_period = 'Unknown'
            
            try:
                if 'Axis Stability' in cell_text:
                    if i + 1 < len(cells):
                        adjacent_cells = cells [i + 1]
                        self.axis_stability = float(adjacent_cells.get_text(strip = True))
            except:
                self.axis_stability = 'Unknown'

        return self.rotation_period, self.axis_stability;

    # Datos Atmosféricos
    def get_atmospheric_data(self):
        cells = self.get_data()
        self.composition = {}
        self.radiation_levels = None
        self.sizeClass = None

        for i, cell in enumerate(cells):
            call_text = cell.get_text(strip=True)


# ---------------
# Se buscará el Size Class dentro de la tabla dentro del HTLM, la idea del mismo es que nos de una pista de la composición del mismo.
# En otras palabras, usaremos la clasificación del mismo para poder basarnos en la composición y su posible atmósfera, usamos aquí un factor incertidumbre.
# ---------------

            try:
                if 'Size Class' in call_text:
                    if i + 1 < len(cells):
                        adjacent_cells = cells[i + 1]
                        self.sizeClass = adjacent_cells.get_text(strip=True)

                if self.sizeClass == 'Mercury-size':
                    self.composition = {
                        'Composition': 'Primarily rocky with a large metallic core.',
                        'Atmosphere': 'Little or no atmosphere due to low gravity. Possible trace gases: helium, hydrogen.'
                    }
                elif self.sizeClass == 'Mars-size':
                    self.composition = {
                        'Composition': 'Rocky with possible traces of water and minerals.',
                        'Atmosphere': 'Very thin atmospheres mainly composed of carbon dioxide. Possible gases: CO2, N2, Ar.'
                    }
                elif self.sizeClass == 'sub-Earth-size':
                    self.composition = {
                        'Composition': 'Rocky with thin atmospheres or no atmosphere.',
                        'Atmosphere': 'May contain water in the form of surface or subsurface ice. Possible gases: H2, He.'
                    }
                elif self.sizeClass == 'Earth-size':
                    self.composition = {
                        'Composition': 'Rocky with a silicate crust, mantle, and metallic core.',
                        'Atmosphere': 'Mainly composed of iron and nickel. Possible gases: N2, O2, CO2, H2O.'
                    }
                elif self.sizeClass == 'super-Earth-size':
                    self.composition = {
                        'Composition': 'Ranges from completely rocky to those with thick atmospheres of hydrogen and helium.',
                        'Atmosphere': 'Some may have liquid water if located in the habitable zone. Possible gases: H2, He, CO2, CH4.'
                    }
                elif self.sizeClass == 'sub-Neptune-size':
                    self.composition = {
                        'Composition': 'Mix of rocky cores, water, ice, and gas mantles.',
                        'Atmosphere': 'Often have hydrogen and helium atmospheres, and sometimes water vapor or ice clouds. Possible gases: H2, He, CH4, NH3.'
                    }
                elif self.sizeClass == 'Neptune-size':
                    self.composition = {
                        'Composition': 'Rock and ice cores with thick ice mantles and atmospheres rich in hydrogen, helium, and methane.',
                        'Atmosphere': 'Gives them a bluish color. Possible gases: H2, He, CH4, H2O.'
                    }
                elif self.sizeClass == 'sub-Jupiter-size':
                    self.composition = {
                        'Composition': 'Predominantly gaseous with rocky or icy cores.',
                        'Atmosphere': 'Thick hydrogen and helium atmospheres with ammonia and water clouds. Possible gases: H2, He, CH4, NH3.'
                    }
                elif self.sizeClass == 'Jupiter-size':
                    self.composition = {
                        'Composition': 'Large amounts of hydrogen and helium with small rock or ice cores.',
                        'Atmosphere': 'Complex cloud bands and storm systems. Possible gases: H2, He, CH4, NH3, H2O.'
                    }
                elif self.sizeClass == 'Super-Jupiter-size':
                    self.composition = {
                        'Composition': 'Similar to Jupiter-size but even more massive with possible dense cores.',
                        'Atmosphere': 'Denser and more complex atmospheres with intense storms and magnetic activity. Possible gases: H2, He, CH4, NH3, H2O.'
                    }
                else:
                    self.composition = {
                        'Composition': 'Unknown',
                        'Atmosphere': 'Unknown'
                    }
            except Exception as e:
                self.composition = {
                    'Composition': 'Unknown',
                    'Atmosphere': 'Unknown'
                }

            try:
                if 'Star Radiation at Atmospheric Boundary:' in call_text:
                    if i + 1 < len(cells):
                        adjacent_cells = cells[i + 1]
                        self.radiation_levels = float(adjacent_cells.get_text(strip=True))
            except Exception as e:
                self.radiation_levels = 'Unknown'

        return self.sizeClass, self.composition, self.radiation_levels

    def get_temperature(self):
        tables = self.soup.find_all('table')

        # Obtenemos los valores de tabla 1.
        if len(tables) > 1:
            table = tables[1]  # Seleccionamos la segunda tabla (índice 1)
            cells = table.find_all('td')
        else:
            return 'No hay suficientes tablas en el documento.'

        # Inicializadores: 
        self.blackbody1 = None
        self.blackbody2 = None
        self.blackbody3 = None
        self.media = 0

        for i in range(len(cells)):
            cell_text = cells[i].text.replace("\n", "").strip()

            if cell_text == "Blackbody T 0.1(K)":
                self.blackbody1 = float(cells[i + 1].text.strip())  
            elif cell_text == "Blackbody T 0.3(K)":
                self.blackbody2 = float(cells[i + 1].text.strip()) 
            elif cell_text == 'Blackbody T 0.7(K)':
                self.blackbody3 = float(cells[i + 1].text.strip())


        blackbodies = [self.blackbody1, self.blackbody2, self.blackbody3]
        valid_blackbodies = [b for b in blackbodies if b is not None and b != 0]

        if valid_blackbodies:
            self.media = round(sum(valid_blackbodies) / len(valid_blackbodies), 4)
        else:
            self.media = None


        # Retornamos los valores obtenidos
        return self.blackbody1, self.blackbody2, self.blackbody3, self.media
    

    # Datos Orbitales
    def get_orbital_data(self):
        tables = self.soup.find_all('table')
        
        # Inicialización
        self.orbital_period = None
        self.eccentricity = None
        self.orbital_type = None
        self.orbital_period_observed = None
        self.orbital_period_estimated = None

         
        if len(tables) > 1:
            table = tables[1]
            cells = table.find_all('td')

            for i in range(len(cells)):
                cell_text = cells[i].text.replace('\n', '').strip()
                try: 
                    # Buscamos el valor de 'Orbital Period(Yrs)'
                    if cell_text == 'Orbital Period(Yrs)':
                        self.orbital_period = round (float(cells[i + 1].text.strip()), 4)
                except:
                    self.orbital_period = 'Unknown'


                try: 
                    if cell_text == '(Days Obs.)':
                        self.orbital_period_estimated = round(float(cells[i + 1].text.strip()), 4)
                except:
                    self.orbital_period_estimated = 'Unknown'


                try:
                    if 'Orbital Period(Observed/Estimated)' in cell_text:
                        self.orbital_period_observed = round(float(cells[i + 1].text.strip()), 4)
                except:
                    self.orbital_period_observed = 'Unknown'


                try: 
                    if cell_text == 'Eccentricity      :':
                        self.eccentricity = float(cells[i + 1].text.strip())
                except:
                    self.eccentricity = None
                    self.orbital_type = None
                
                if self.eccentricity != None:
                    if 0.0 <= self.eccentricity <= 0.1:
                        self.orbital_type = 'Nearly Circular'
                    elif 0.1 < self.eccentricity <= 0.5:
                        self.orbital_type = 'Moderate'
                    elif 0.5 < self.eccentricity <= 1:
                        self.orbital_type = 'Very Elliptical'
                    elif self.eccentricity > 1:
                        self.orbital_type = 'Hyperbolic'
                    else:
                        self.orbital_type = 'Unknown'
            
        return self.orbital_period, self.orbital_period_observed, self.orbital_period_estimated, self.eccentricity, self.orbital_type

    def get_stellar_data(self):
    # Abro el archivo stellar_activity_mapping.json, donde se clasificó mediante un diccionario todos los valores
#   posibles del tipo de estrella y su subclase.
        try:
            with open('mapping\stellar_activity_mapping.json', 'r') as file:
                stellar_activity_mapping = json.load(file)
        except FileNotFoundError:
            print('Error: The JSON file was not found.')
        except FileNotFoundError:
            print('Error: JSON file is not valid or empty.')


        tables = self.soup.find_all('table')

        self.star_name = None
        self.star_type = None
        self.distance_from_star = None
        self.star_radiation_atmospheric_b = None
        self.stellar_activity = None

        if len(tables) > 1:
            table = tables[1]
            cells = table.find_all('td')

            for i in range(len(cells)):
                cell_text = cells[i].text.replace('\n', '').strip()
                try: 
                    if 'Star Name' in cell_text:
                        self.star_name = (cells[i + 1].text.strip())
                except:
                    self.star_name = 'Unknown'

                try: 
                    if 'Spectral type :' == cell_text:
                        self.star_type = cells[i + 1].text.strip()
                except:
                    self.star_type = 'Unknwon'
                
                if self.star_type:
                    star_type_letter  = self.star_type[0]
                    try: 
                        star_type_number  = int(self.star_type[1])
                    except:
                        star_type_number = str(self.star_type[1])

                    if star_type_letter in stellar_activity_mapping:
                        if isinstance(stellar_activity_mapping[star_type_letter], dict) and star_type_number is not None:
                            self.stellar_activity = stellar_activity_mapping[star_type_letter].get(str(star_type_number), 'Unknown')
                    
                        else:
                            self.stellar_activity = stellar_activity_mapping[star_type_letter]
                
                    else:
                        self.stellar_activity = 'Unknown'

                try: 
                    if 'Atmospheric Boundary' in cell_text:
                        self.star_radiation_atmospheric_b = round(float(cells[i + 1].text.strip()), 2)
                except:
                    self.star_radiation_atmospheric_b = 'Unknown'
                
                
                try:
                    if 'Semi Major Axis / Orbital Distance Calculated' == cell_text or cell_text in 'Semi Major Axis / Orbital Distance Calculated':
                        self.distance_from_star = float(cells[i + 1].text.strip())
                    
                    elif cell_text in ['(AU Est.)']:
                        self.distance_from_star = float(cells[i - 1].text.strip())

                except:
                    if self.distance_from_star in [None]:
                        self.distance_from_star = 'Unknown'
                

        ### Calculo de Zona Habitable se buscan las variables necesarias dentro de la tabla:
        # Stellar Radius
        # Stellar Mass
        # star_temperature
        # distance_from_star.

        try:
            # Constantes
            sigma = 5.67e-8  # Constante de Stefan-Boltzmann en W/m²K⁴
            R_sun = 6.96e8  # Radio del Sol en metros
            L_sun = 3.828e26  # Luminosidad del Sol en Watts
            stellar_radius = ''
            star_temperature = ''
            self.in_habitable_zone = 'N/A'

            # Variables
            for i in range(len(cells)):
                cell_text = cells[i].text.replace('\n', '').strip()

                try:
                    if 'Stellar Radius' in cell_text:
                        stellar_radius = float(cells[i + 1].text.strip())
                except:
                    stellar_radius = None

                try:
                    if 'Temperature' in cell_text:
                        star_temperature = float(cells[i + 1].text.strip())
                except:
                    star_temperature = None

            # Radio de la Estrella en Metros
            stellar_radius_meters = stellar_radius * R_sun

            # Luminosidad usando la fórmula de Stefan-Boltzmann
            star_luminosity = 4 * math.pi * (stellar_radius_meters ** 2) * sigma * (star_temperature ** 4)

            # Luminosidad en unidades solares
            star_luminosity_sun = star_luminosity / L_sun

            # Límites de la zona habitable en UA
            inner_habitable_zone = math.sqrt(star_luminosity_sun / 1.1)
            outer_habitable_zone = math.sqrt(star_luminosity_sun / 0.53)

            #Variables inner y outer.
            self.inner_habitable_zone = inner_habitable_zone
            self.outer_habitable_zone = outer_habitable_zone


            # Determinación de si el exoplaneta está en la zona habitable
            self.in_habitable_zone = inner_habitable_zone <= self.distance_from_star <= outer_habitable_zone

            if self.in_habitable_zone == False:
                if self.distance_from_star <= inner_habitable_zone:
                    self.in_habitable_zone = "Inside the inner border."
                elif self.distance_from_star >= outer_habitable_zone:
                    self.in_habitable_zone = 'Outside the outer border'

        except:
            self.inner_habitable_zone = self.outer_habitable_zone = 0
            self.in_habitable_zone = 'Unknown'


        return self.star_name, self.star_type, self.distance_from_star, self.star_radiation_atmospheric_b, self.stellar_activity, self.in_habitable_zone

    def get_tidally_locked(self, soup):

        tables = self.soup.find_all('table')
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
        mass_p = self.j_mass  # Masa del planeta en masas de Júpiter
        radius_p = self.j_radius  # Radio del planeta en radios de Júpiter
        eccentricity = self.eccentricity  # Excentricidad orbital

        periodo_orbital = self.orbital_period  # Período orbital en días

        if periodo_orbital == 0:
            periodo_orbital = self.orbital_period_observed or self.orbital_period_estimated

        semi_mayor_axis = self.distance_from_star  # Semi-eje mayor en UA
        rotation = self.rotation_period  # Período de rotación en horas

        # Datos estelares
        self.stellar_mass = 'Dato Primario'
        self.stellar_radius = 'Dato Primario'
        self.stellar_temperature = 'Dato Primario'
      
        if len(tables)> 1:
            table = tables [1]
            cells = table.find_all('td')

            for i in range (len(cells)):
                cell_text = cells [i].text.replace('\n', '').strip()

                try:
                    if 'Stellar Mass(Msun)' in cell_text:
                        self.stellar_mass = round(float(cells[i + 1].text.strip()), 4)
                        if self.stellar_mass == 0:
                            cell_text == 'Msun Measured'
                            self.stellar_mass = round(float(cells[i + 1].text.strip()), 4)
                            if self.stellar_mass == 0:
                                self.stellar_mass = None
                except:
                    self.stellar_mass = None
                
                try:
                    if 'Stellar Radius(Rsun) :' in cell_text:
                        self.stellar_radius = round(float(cells[i + 1].text.strip()), 4)
                        if self.stellar_radius == 0:
                            cell_text == 'Msun Measured'
                            self.stellar_radius = round(float(cells[i + 1].text.strip()), 4)
                            if self.stellar_radius == 0:
                                self.stellar_radius = None

                except:
                    self.stellar_radius = None


                try:
                    if 'Temperature   :' in cell_text:
                        self.stellar_temperature = round(float(cells[i + 1].text.strip()), 4)
                        if self.stellar_temperature == 0:
                            cell_text == 'Msun Measured'
                            self.stellar_temperature = round(float(cells[i + 1].text.strip()), 4)
                            if self.stellar_temperature == 0:
                                self.stellar_temperature = None
                except:
                    self.tellar_temperature = None

        # Conversión a unidades 
        j_m_k = 1.898 * 10**27  # Masa de Júpiter en kg
        j_r_m = 7.1482 * 10**7  # Radio de Júpiter en metros

        mass_p_kg = mass_p * j_m_k
        radius_p_m = radius_p * j_r_m
        volumen = (4/3) * math.pi * radius_p_m**3
        density = mass_p_kg / volumen

        # m_star_kg = stellar_mass * M_sun

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

        # Calcular el tiempo sincrónico
        self.t_sync = ((k * self.stellar_radius**2 * self.stellar_mass * periodo_orbital**2))**(1/3)

        if semi_mayor_axis > 1.0:  # Si la distancia del planeta es mayor a 1 UA
            self.probabilidad_anclaje = "Low probability of tidal locking"
            base = 'Far from gravitational influence'
        else:
            if semi_mayor_axis <= self.t_sync:
                base = "Possible tidal locking"
            elif self.t_sync <= 10:
                base = "Within the gravitational influence of the star"
            elif self.t_sync <= 50:
                base = 'Gravitational influence of its star; moderate'
            else:
                base = 'Far from gravitational influence'

            # Cálculo del factor Q utilizando el calentamiento mareal
            self.P_diss_min = 1e13  # Potencia disipada mínima estimada en vatios
            self.P_diss_max = 1e18  # Potencia disipada máxima estimada en vatios

            Q_min, _, Q_max, _ = self.calcular_factor_Q()

            # Determinación de la probabilidad de anclaje por marea
            if Q_min > 1e6:
                self.probabilidad_anclaje = "Low probability of tidal locking"
            elif Q_max < 1e4:
                self.probabilidad_anclaje = "High probability of tidal locking"
            else:
                self.probabilidad_anclaje = "Medium probability of tidal locking"

        return self.t_sync, base, self.probabilidad_anclaje


    def calcular_factor_Q(self):
        """ Calcula el factor de calidad Q usando el calentamiento mareal. """
        P_diss_max = self.P_diss_max
        P_diss_min = self.P_diss_min
        M_star = self.stellar_mass
        R_planet = self.j_radius
        eccentricity = self.eccentricity
        a = self.distance_from_star

        # Constantes
        G = 6.67430e-11  # Constante gravitacional en m^3 kg^−1 s^−2
        M_sun = 1.989e30  # Masa del Sol en kg
        R_jup = 7.1492e7  # Radio de Júpiter en metros
        AU = 1.496e11  # Unidad Astronómica en metros

        # Convertir a unidades correctas
        M_star_kg = M_star * M_sun  # Masa de la estrella en kg
        R_planet_m = R_planet * R_jup  # Radio del planeta en metros
        a_m = a * AU  # Semi-eje mayor en metros

        if a_m > 1 * AU:
            P_diss_min = 1e10  # Disipación mínima para planetas lejanos
            P_diss_max = 1e12  # Disipación máxima para planetas lejanos
        else:
            P_diss_min = 1e13  # Valores más altos para planetas cercanos
            P_diss_max = 1e18

        # Calcular el factor Q para el rango de potencia disipada
        Q_min = (63 / 4) * ((G * M_star_kg)**(3/2)) * (R_planet_m**5) * eccentricity / (P_diss_max * a_m**(15/2))
        Q_max = (63 / 4) * ((G * M_star_kg)**(3/2)) * (R_planet_m**5) * eccentricity / (P_diss_min * a_m**(15/2))

        # Interpretación de los resultados
        if Q_min < 10**5:
            q_min_interpretation = "High energy dissipation; possible significant internal heating"
        elif 10**5 <= Q_min < 10**7:
            q_min_interpretation = "Moderate energy dissipation; moderate geological activity or internal heat"
        else:
            q_min_interpretation = "Low energy dissipation; likely a rigid body with little to no tidal activity"

        if Q_max < 10**5:
            q_max_interpretation = "High energy dissipation; extreme tidal heating conditions"
        elif 10**5 <= Q_max < 10**7:
            q_max_interpretation = "Moderate energy dissipation; potential for mild to moderate geological activity"
        else:
            q_max_interpretation = "Low energy dissipation; stable internal structure, with little to no tide-induced activity"

        return Q_min, q_min_interpretation, Q_max, q_max_interpretation

    def get_essentials_data(self):
        tables = self.soup.find_all('table')
        
        # Inicialización
        self.oxygen = None
        self.gravity = None
        self.atmospheric_pressure = None
        self.magnetic_field = None
        self.volcanic_activity = None
        # self.molar_mass = 'No hay valor'
         
        if len(tables) > 1:
            table = tables[1]
            cells = table.find_all('td')

            for i in range(len(cells)):
                cell_text = cells[i].text.replace('\n', '').strip()
         # --------------------   
                try: 
                    if cell_text == 'Surface Gravity (G Earth)':
                        self.gravity = round (float(cells[i + 1].text.strip()), 4)
                except:
                    self.gravity = 'Unknown'
         # --------------------
                try: 
                    if cell_text == 'Oxygen' or 'Oxygen' in cell_text:
                        self.oxygen = cells[i + 1].text.strip()
                except:
                    self.oxygen = 'Unknown'
         # --------------------
                try:
                    if 'Pressure (bar)' in cell_text:
                        self.atmospheric_pressure = round(float(cells [i + 1].text.strip()), 4)
                except: 
                    self.atmospheric_pressure = 'Unknown'
         # --------------------


         # --------------------
                try:
                    if 'Volcanic Activity' in cell_text:
                        self.volcanic_activity = cells [i + 1].text.strip()
                except:
                    self.volcanic_activity = 'Unknown'

        # --------------------------------------------------------------
    
        if self.atmospheric_pressure in ['Unknown', None]:
            try: 
                if self.sizeClass in ['Mercury-size', 'Sub-Earth-size']:
                    self.molar_mass = (0.002016 + 0.0040026) / 2  #Promedio entre H2 y He
                elif self.sizeClass in ['Mars-size']:
                    self.molar_mass = 0.04401   #CO2
                elif self.sizeClass in ['Earth-size', 'Super-Earth-size']:
                    if 'O2' in self.composition['Atmosphere'] or 'N2' in self.composition['Atmosphere']:
                        self.molar_mass = (0.032 + 0.028014) / 2  # Promedio de N2 y O2
                    else:
                        self.molar_mass = 0.04401 #Atmosfera rica en CO2
                elif self.sizeClass in ['Sub-Neptune-size', 'Neptune-size']:
                    self.molar_mass = (0.002016 + 0.0040026 + 0.01604)/3  #Primedio H2, He y CH4
                elif self.sizeClass in ['Sub-Jupiter-size', 'Jupiter-size', 'Super-Jupiter-size']:
                    self.molar_mass = (0.002016 + 0.0040026) / 2 # Promedio de H2 y He
                else: 
                    self.molar_mass = 0.04401  # CO2 como valor por defecto

                # Constante de las fases ideales en J/(mol-k)
                r = 8.314

                # Estimación de altura en metros.
                height_scale = (r * self.media) / (self.gravity * self.molar_mass)

                # Estimación de presión, la misma puede seguir ajustandose pero se usa un valor estandar.
                base_pressure = 1013.25 #Valor estándar de presión atmosférica en hP
                self.atmospheric_pressure = round(base_pressure * (1/height_scale), 4)

            except:
                self.atmospheric_pressure = 'It is not possible to make a reliable calculation' 


        if self.oxygen in ['Unknown', None]:
            try: 
                if 'O2' in self.composition.get('Atmosphere', ''):
                    self.oxygen = 'Possible oxygen detected directly in atmospheric composition.'
            
                # Ajustar la evaluación según el tamaño del planeta y la zona habitable.
                if self.sizeClass in ['Jupiter-size', 'Sub-Jupiter-size', 'Neptune-size', 'Super-Jupiter-size', 'Sub-Neptune-size']:
                    if self.in_habitable_zone == 'Outside the inner border':
                        if self.media < 150 :
                            self.oxygen = 'Oxygen likely in ice form due to low temperatures and distance'
                        else:
                            self.oxygen = 'Oxygen presence unlikely, characteristics typical of gas giants'
                    elif self.in_habitable_zone == 'Inside the inner border':
                        if self.media > 300:
                            self.oxygen = 'Oxygen unlikely due to high temperature and proximity to the star'
                    else:
                        self.oxygen = 'Oxygen presence in question due to mixed conditions'

                elif self.sizeClass in ['Earth-size', 'Super-Earth-size', 'Sub-Earth-size']:
                    if self.in_habitable_zone:
                        if self.media > 300:
                            self.oxygen = 'Oxygen unlikely due to high temperature'
                        elif self.media < 200:
                            self.oxygen = 'Oxygen possibly retained due to low temperatures'
                        else: 
                            self.oxygen = 'Moderate chance of oxygen presence'
                    else:
                        if self.media < 150:
                            self.oxygen = 'Oxygen may be retained in solid form due to low temperatures outside habitable zone'
                        else: 
                            self.oxygen = 'Oxygen presence unlikely due to distance from habitable zone'
                
                # Ajustar evaluación según el tipo de estrella y la temperatura.
                star_type = self.star_type
                if star_type[0] in ['A', 'B', 'F']:
                    if self.media > 250:
                        self.oxygen = 'Possible oxygen dissociation due to high stellar radiation and temperatures'
                
                # Ajustar la evaluación según los niveles de radiación
                if self.radiation_levels > 10:
                    self.oxygen += ' High radiation levels may reduce the presence of oxygen.'
            except:
                self.oxygen = 'Evaluation failed. Incomplete or inconsistent data.'

        # -------------------------------
        # Obtención de datos para conocer si es probable si el exoplaneta posee Escudo Magnético.
        # Inicialización

        if self.magnetic_field == None:
            try:
                if self.sizeClass in ['Super-Jupiter-size', 'Jupiter-size']:
                    if self.probabilidad_anclaje == 'Low probability of tidal locking':
                        if self.in_habitable_zone == 'Outside the inner border':
                            if self.stellar_activity in ['Low', 'Very Low', 'Moderate']:
                                self.magnetic_field = 'Probably Strong Magnetic Field'
                            elif self.stellar_activity == 'High':
                                self.magnetic_field = 'Probably Eroded Magnetic Field'
                            else:
                                self.magnetic_field = 'Unknown Stellar Activity'
                        elif self.in_habitable_zone == 'Inside the inner border':
                            if self.stellar_activity in ['Low', 'Very Low', 'Moderate']:
                                self.magnetic_field = 'Probably Strong Magnetic Field'
                            elif self.stellar_activity == 'High':
                                self.magnetic_field = 'Probably Eroded Magnetic Field'
                            else:
                                self.magnetic_field = 'Unknown Stellar Activity'
                        elif self.in_habitable_zone == 'Outside the outer border':
                            self.magnetic_field = 'Possibly Weak Magnetic Field due to distance from star'
                        else:
                            self.magnetic_field = 'Unknown Habitable Zone Status'
                    else:  # High probability of tidal locking
                        self.magnetic_field = 'Probably Weak or No Magnetic Field'
                
                elif self.sizeClass in ['Sub-Jupiter-size', 'Neptune-size', 'Sub-Neptune-size']:
                    if self.probabilidad_anclaje == 'Low probability of tidal locking':
                        if self.in_habitable_zone == 'Outside the inner border':
                            if self.stellar_activity in ['Low', 'Very Low', 'Moderate']:
                                self.magnetic_field = 'Possibly Weak Magnetic Field'
                            elif self.stellar_activity == 'High':
                                self.magnetic_field = 'Probably Eroded Magnetic Field'
                            else:
                                self.magnetic_field = 'Unknown Stellar Activity'
                        elif self.in_habitable_zone == 'Inside the inner border':
                            if self.stellar_activity in ['Low', 'Very Low', 'Moderate']:
                                self.magnetic_field = 'Possibly Weak Magnetic Field'
                            elif self.stellar_activity == 'High':
                                self.magnetic_field = 'Probably Eroded Magnetic Field'
                            else:
                                self.magnetic_field = 'Unknown Stellar Activity'
                        elif self.in_habitable_zone == 'Outside the outer border':
                            self.magnetic_field = 'Likely Very Weak Magnetic Field due to distance from star'
                        else:
                            self.magnetic_field = 'Unknown Habitable Zone Status'
                    else:  # High probability of tidal locking
                        self.magnetic_field = 'Probably Weak or No Magnetic Field'
                
                elif self.sizeClass in ['Super-Earth-size', 'Earth-size']:
                    if self.probabilidad_anclaje == 'Low probability of tidal locking':
                        if self.in_habitable_zone == 'Outside the inner border':
                            if self.stellar_activity in ['Low', 'Very Low', 'Moderate']:
                                self.magnetic_field = 'Probably Strong Magnetic Field'
                            elif self.stellar_activity == 'High':
                                self.magnetic_field = 'Probably Eroded Magnetic Field'
                            else:
                                self.magnetic_field = 'Unknown Stellar Activity'
                        elif self.in_habitable_zone == 'Inside the inner border':
                            if self.stellar_activity in ['Low', 'Very Low', 'Moderate']:
                                self.magnetic_field = 'Probably Strong Magnetic Field'
                            elif self.stellar_activity == 'High':
                                self.magnetic_field = 'Probably Eroded Magnetic Field'
                            else:
                                self.magnetic_field = 'Unknown Stellar Activity'
                        elif self.in_habitable_zone == 'Outside the outer border':
                            self.magnetic_field = 'Weak Magnetic Field due to distance from star'
                        else:
                            self.magnetic_field = 'Unknown Habitable Zone Status'
                    else:  # High probability of tidal locking
                        self.magnetic_field = 'Probably Weak or No Magnetic Field'
                
                elif self.sizeClass in ['Sub-Earth-size', 'Mars-size', 'Mercury-size']:
                    if self.probabilidad_anclaje == 'Low probability of tidal locking':
                        if self.in_habitable_zone == 'Inside the inner border':
                            self.magnetic_field = 'Possibly Weak Magnetic Field'
                        else:
                            self.magnetic_field = 'Weak or No Magnetic Field'
                    else:
                        self.magnetic_field = 'Probably Weak or No Magnetic Field'
                
                else:
                    self.magnetic_field = 'Evaluation Inconclusive'

            except Exception as e:
                self.magnetic_field = False


    # -------------------------------
    #   Posibilidad de si posee agua o no. Dado a en que en algunos casos es imposible conocer el dato a ciencia cierta,
    #  Se trabajará haciendo estimaciones con los datos anteriores ya obtenidos, contemplando diferentes escenarios.

        impossible_to_form_water = (self.sizeClass in ['Jupiter-size', 'Sub-Jupiter-size', 'Super-Jupiter-size'] 
                                    or self.radiation_levels > 1000)

        possible_states = []

        # Verificación de datos faltantes
        if self.media in ['Unknown', None] or self.in_habitable_zone in ['Unknown', None] or self.radiation_levels in ['Unknown', None]:
            has_water = 'Data inconclusive'  # Caso de datos faltantes

        else:
            # Agua en estado sólido (hielo)
            if self.media < 273 and self.gravity > 5:
                possible_states.append('Possibility of water in solid state.')

            # Agua en estado líquido
            if (273 <= self.media <= 300 and self.in_habitable_zone == True and 
                self.radiation_levels < 10 and self.stellar_activity in ['Low', 'Moderate'] and 
                self.gravity >= 5 and self.probabilidad_anclaje == 'Low probability of tidal locking'):
                possible_states.append('Possibility of liquid water.')

            # Agua en estado gaseoso (vapor)
            if self.media > 300 or self.radiation_levels > 1:
                possible_states.append('Possibility of water in gaseous state.')

            # Evaluación final
            if impossible_to_form_water:
                has_water = 'Impossible for water to form in any state due to planet type or radiation levels.'
            elif possible_states:
                states_str = ' and '.join(possible_states)
                has_water = f'Possible formation of water in {states_str}'
            else:
                has_water = 'No detectable water in any form based on current data.'

        return self.oxygen, self.gravity, self.magnetic_field, self.atmospheric_pressure, self.volcanic_activity, has_water

    def show_data(self):
            print(f"Masa en Jupiters: {self.j_mass}")
            print(f"Masa en Tierras: {self.e_mass}")
            print(f'Radio en Jupiters: {self.j_radius}')
            print(f"Radio en Tierras: {self.e_radius}")
            print(f"Rotation Period: {self.rotation_period}")
            print(f'Axis Stability: {self.axis_stability}')
            print(f"Size Class: {self.sizeClass}")
            print(f"Composition: {self.composition}")
            print(f'Radiation Levels: {self.radiation_levels}')
            print(f'B 0.1: {self.blackbody1}')
            print(f'B 0.3: {self.blackbody2}')
            print(f'B 0.7: {self.blackbody3}')
            print(f"media: {self.media}")
            print(f'Orbital Period (Yrs): {self.orbital_period}')
            print(f'Orbital Period estimado (Days): {self.orbital_period_estimated}')
            print(f'Periodo orbital observado(Days) : {self.orbital_period_observed}')
            print(f'excentricidad: {self.eccentricity}')
            print(f'Tipo de Orbita {self.orbital_type}')
            print('')
            print('//////////////')
            print('')
            print(f'Nombre de Estrella: {self.star_name}')
            print(f'Tipo de Estrella: {self.star_type}')
            print(f'Distancia de la Estrella: {self.distance_from_star}')
            print(f'Radiacion Stellar a la Atmosfera: {self.star_radiation_atmospheric_b}')
            print(f'Actividad Estelar: {self.stellar_activity}')
            print(f'Zona habitable: {self.in_habitable_zone}')
            print('')
            print('//////////////')
            print('')
            print(f'Radio de la Estrella: {self.stellar_radius}')
            print(f'Masa:{self.stellar_mass}')
            print(f'Temperatura  {self.stellar_temperature}')
            print('')
            print(f'Tiempo de Sincronización: {self.t_sync}')
            print(f'Base: {self.base}')
            print(f'Probabilidad Anclaje {self.probabilidad_anclaje}')
            print(f'{self.Q_min}, {self.q_min_interpretation}')
            print(f'{self.Q_max}, {self.q_max_interpretation}')
            print('')
            print('//////////////')
            print('')
            print(f'Valor nuevo oxigeno: {self.oxygen}')
            print(f'Valor Nuevo gravedad: {self.gravity}')
            print(f'Escudo Magnetico: {self.magnetic_field}')
            print(f'Presion Atmosferica: {self.atmospheric_pressure}')
            print(f'Actividad Volcánica: {self.volcanic_activity}')
            print(f'Posibilidad de agua: {self.has_water}')
            print(f'')
            print(f'')

class Habitability_Score():
    def __init__(self, planet):
        self.planet = planet

    def calculate_score(self):
        score = 0
        p  = self.planet
        m_e_j = 0.00315       #Masa de la tierra en Júpiters.
        r_e_j = 0.0892        # Radio Tierras en Jupiters

        # Score de Masas ----- 10 pts
        if p.j_mass >= 10:      # Super gigante gaseoso
            score += -20

        elif 2 <= p.j_mass < 10:     # Gigante Gaseoso
            score += -15

        elif 0.1 <= p.j_mass < 2:      # Neptuniano
            score += -8

        elif 1.2 <= p.e_mass < 10:      # Supertierras
            score += 5

        elif 0.8 <= p.e_mass < 1.2:      # Tierras
            score += 10

        else:                           # Planeta Enano
            score += 0


        # Score de Radios ---- 10pts
        if p.j_radius >= 1.5:
            score += -10
           
        elif 1 <= p.j_radius < 1.5:
            score += -5

        elif 0.7 <= p.j_radius < 1:      # Neptunianos
            score += -2

        elif 1 <= p.e_radius < 2:      # Supertierras
            score += 5

        elif 0.8 <= p.e_radius < 1.2:      # Tierras
            score += 10

        else:                           # Planeta Enano
            score += -2


        # Score periodo orbital, de primera mano se intentará conocer si el periodo orbital es 0, en caso de que este lo sea
#   se avanzará sobre self.orbital_period (estimated y observed) divido en 365 (años terrestres) ya que estos estan en días.

        if p.orbital_period == 0:
            p.orbital_period = (p.orbital_period_estimated/365) or (p.orbital_period_observed/365)

        # Score habitual: Periodo Orbital --- 8pts.

        if p.orbital_period < 1:
            score += 0              

        elif 1 <= p.orbital_period < 50:
            score += 2  

        elif 50 <= p.orbital_period < 200:
            if p.in_habitable_zone:
                score += 8  
            else:
                score += 4  

        elif 200 <= p.orbital_period < 400:
            if p.in_habitable_zone:
                score += 10  
            else:
                score += 5  

        elif 400 <= p.orbital_period < 1000:
            if p.in_habitable_zone:
                score += 8  
            else:
                score += 2  

        else:
            score += 0


        # Temperatura: 
        if p.media >= 300:
            score += 0
        elif 300 > p.media >= 273:
            score += 15
        elif 273 > p.media > 200:
            score += 10
        elif 200 > p.media >= 150:
            score += 7
        else:
            score += 0

        # Presencia de Agua (15 pts)
        if p.has_water == 'Impossible for water to form in any state':
            score += 0
        elif p.has_water == 'Possible formation of water in ice state':
            score += 3
        elif p.has_water == 'Possible formation of water in vapor state':
            score += 5
        elif p.has_water == 'Possible formation of water in ice and liquid state':
            score += 8
        elif p.has_water == 'Possible formation of water in ice and vapor state':
            score += 10
        elif p.has_water == 'Possible formation of water in liquid and vapor state':
            score += 10
        elif p.has_water == 'Possible formation of water in ice and liquid state, and vapor state':
            score += 15
        elif p.has_water == 'Possible formation of water in ice, liquid, and vapor state':
            score += 15
        elif p.has_water == 'No detectable water in any form':
            score += 0
        else:
            score += 0

        # Presencia de Oxígeno (10 pts)
        if p.oxygen == 'Evaluation failed. Incomplete or inconsistent data':
            score += 0
        elif p.oxygen == 'Evaluation inconclusive without direct atmospheric data':
            score += 0
        elif p.oxygen == 'Oxygen presence unlikely, characteristics typical of gas giants':
            score += 0
        elif p.oxygen == 'Oxygen presence unlikely due to distance from habitable zone':
            score += 1
        elif p.oxygen == 'Oxygen unlikely due to high temperature':
            score += 2
        elif p.oxygen == 'Oxygen unlikely due to high temperature and proximity to the star':
            score += 2
        elif p.oxygen == 'Possible oxygen dissociation due to high stellar radiation and temperatures':
            score += 2
        elif p.oxygen == '"Oxygen may be retained in solid form due to low temperatures outside habitable zone':
            score += 3
        elif p.oxygen == 'High radiation levels may reduce the presence of oxygen':
            score += 3
        elif p.oxygen == 'Oxygen may be present in solid form due to low temperatures':
            score += 3
        elif p.oxygen == 'Oxygen likely in ice form due to low temperatures and distance':
            score += 4
        elif p.oxygen == 'Oxygen presence in question due to mixed conditions':
            score += 5
        elif p.oxygen == 'Oxygen possibly retained due to low temperatures':
            score += 6
        elif p.oxygen == 'Moderate chance of oxygen presence':
            score += 7
        elif p.oxygen == 'Oxygen detected directly in atmospheric composition':
            score += 10
        else:
            score += 0


                       
        # Niveles de Co2 y CH4 (8pts)

        # Excentricidad Orbital
        if p.eccentricity >= 0.5:
            score += 0
        elif 0.5 > p.eccentricity >= 0.3:
            score += 1
        elif 0.3 > p.eccentricity > 0.1:
            score += 3
        else:
            score += 5

        # Actividad Volcánica (5 pts)
        if p.volcanic_activity == 'High':
            score += 0
        elif p.volcanic_activity == 'Moderate':
            score += 5
        elif p.volcanic_activity == 'Low':
            score += 3
        else:
            score += 0



        # Campo Magnético (8 pts)
        if p.magnetic_field == 'Probably Strong Magnetic Field':
            score += 8
        elif p.magnetic_field == 'Possibly Weak Magnetic Field':
            score += 5
        elif p.magnetic_field == 'Probably Eroded Magnetic Field':
            score += 3
        elif p.magnetic_field == 'Probably Weak or No Magnetic Field':
            score += 1
        elif p.magnetic_field == 'Evaluation Inconclusive':
            score += 0
        else: 
            score += 0


        # Periodo de rotación (5 pts)
        if p.distance_from_star <= p.t_sync:
            score += 0
        elif p.probabilidad_anclaje == 'Low probability of tidal locking':
            score += 5
        elif p.probabilidad_anclaje == 'Medium probability of tidal locking':
            score += 3
        else:
            score += 0

        # Estabilidad del eje (5 pts)
        if not p.axis_stability:
            score += 0
        elif p.axis_stability == 'Stable':
            score += 5
        elif p.axis_stability == 'Moderately Stable':
            score += 4
        elif p.axis_stability == 'Precessing':
            score += 3
        elif p.axis_stability == 'Tilting':
            score += 2
        elif p.axis_stability == 'Unstable':
            score += 1
        elif p.distance_from_star <= p.t_sync and p.probabilidad_anclaje == 'High probability of tidal locking':
            score += 0
        else:
            score += 0


        # Composición atmosférica (10pts)----> Detalles.

# ---------------
#   La determinación de la composición atmosférica utilizará los SizeClass como referencia inicial, 
#   pero no será determinante. Un exoplaneta Earth-size, por ejemplo, podría tener una atmósfera rica en O2 y CO2, 
#   pero también es posible que tenga una atmósfera altamente presurizada rica en metano y helio.
#   Composición atmosférica
# ---------------

        if p.composition == 'Mercury-size':
            score += 0  # Posibles gases: Helio (He), Hidrógeno (H2)

        elif p.composition == 'Mars-size':
            score += 5  # Posibles gases: Dióxido de Carbono (CO2), Nitrógeno (N2), Argón (Ar)

        elif p.composition == 'Sub-Earth-size':
            score += 0  # Posibles gases: Hidrógeno (H2), Helio (He).

        elif p.composition == 'Earth-size':
            score += 0  # Posibles gases: Nitrógeno (N2), Oxígeno (O2), Dióxido de carbono (CO2), Vapor de agua (H2O).

        elif p.composition == 'Super-Earth-size':
            score += 5  # Posibles gases: Hidrógeno (H2), Helio (He), Dióxido de carbono (CO2), Metano (CH4).

        elif p.composition == 'Sub-Neptune-size':
            score += 0  # Posibles gases: Hidrógeno (H2), Helio (He), Metano (CH4), Amoníaco (NH3), Vapor de agua (H2O).

        elif p.composition == 'Sub-Jupiter-size':
            score += 0  # Posibles gases: Hidrógeno (H2), Helio (He), Metano (CH4), Amoníaco (NH3).

        elif p.composition == 'Jupiter-size':
            score += 0  # Posibles gases: Hidrógeno (H2), Helio (He), Metano (CH4), Amoníaco (NH3), Vapor de agua (H2O).

        elif p.composition == 'Super-Jupiter-size': 
            score += 0  # Posibles gases: Hidrógeno (H2), Helio (He), Metano (CH4), Amoníaco (NH3), Vapor de agua (H2O).

        else:
            score += 0


        # Radiación ---- 8pts:
# ---------------
#   Nos basamos en la cantidad de radiación que le llega a la atmósfera, mediante los datos otorgados en 
#   self.star_radiation_atmospheric.
# ---------------

        if p.star_radiation_atmospheric_b > 100:
            score += 0
        elif 50 <= p.star_radiation_atmospheric_b <= 99:
            score += 0
        elif 10 <= p.star_radiation_atmospheric_b <= 49:
            score += 2
        elif 1 <= p.star_radiation_atmospheric_b <= 9:
            score += 4
        elif 0.1 <= p.star_radiation_atmospheric_b < 1:
            score += 6
        elif p.star_radiation_atmospheric_b < 0.1:
            score += 8
        else:
            score += 0


        # Actividad Estelar ---- 8pts:
# ---------------
# La actividad estelar es un determinante en el espacio comunal que ocupan los planetas en un sistema,
# esta misma es determinan por los datos que almacenamos en mapping, que de acuerdo con el tipo de estrella,
# podemos determinar o bien, podemos estimar la posible actividad del mismo. 
# ** Con actividad estellar nos referimos a los fenómenos dinámicos que ocurren en la superficie de una estrella anfitriona,
# ** Si una estrella es muy activa, puede emitir grandes cantidades de radiación ultravioleta y particulas cargadas, 
# que pueden erosionar la atmósfera de un planeta.
# ---------------

        if p.stellar_activity == 'Very Low':
            score += 8
        elif p.stellar_activity == 'Low':
            score += 6
        elif p.stellar_activity == 'Low-Moderate':
            score += 5
        elif p.stellar_activity == 'Moderate':
            score += 4
        elif p.stellar_activity == 'Moderate-High':
            score += 1
        elif p.stellar_activity == 'High':
            score += 0

        elif p.stellar_activity == 'Very High':
            score += -2
        else:
            score += 0



        # Distancia de la Estrella y Zona Habitable --- 10pts.

        try:
            # Usar las variables del objeto planet
            if self.planet.inner_habitable_zone <= self.planet.distance_from_star <= self.planet.outer_habitable_zone:
                self.in_habitable_zone = "Optimal zone."
                score += 10  # Zona óptima

            else:
                if self.planet.distance_from_star < self.planet.inner_habitable_zone * 0.5 or self.planet.distance_from_star > self.planet.outer_habitable_zone * 2:
                    self.in_habitable_zone = "Very close or too far from the habitable zone."
                    score += 0  # Muy cerca o muy lejos de la zona habitable.

                elif self.planet.inner_habitable_zone * 0.5 <= self.planet.distance_from_star < self.planet.inner_habitable_zone:
                    self.in_habitable_zone = "Inside the inner border but not optimal."
                    score += 3  # Dentro del rango pero no óptimo.

                elif self.planet.inner_habitable_zone <= self.planet.distance_from_star <= self.planet.outer_habitable_zone:
                    self.in_habitable_zone = "Within the habitable zone but close to the edges."
                    score += 6  # Dentro de la zona pero cercano a los bordes.

                else: 
                    self.in_habitable_zone = 'Unknown'
                    score += 0  # Caso de error.

        except Exception as e:
            # Manejo de excepciones
            self.planet.inner_habitable_zone = self.planet.outer_habitable_zone = 0
            self.in_habitable_zone = 'Unknown'
            score += 0  # Caso de error



 
        #Tipos de Estrellas: 

        type = p.star_type[0]
        number = int(p.star_type[1])

        if type == 'O':
            score += 0

        if type == 'B':
            score += 1 if number <= 4 else 2


        elif p.star_type[0]  == 'A':
            score += 3 if number <= 4 else 4


        elif p.star_type[0]  == 'F':
            score += 5 if number <= 4 else 6


        elif p.star_type[0]  == 'G':
            score += 7 if number <= 4 else 8


        elif p.star_type[0]  == 'K':
            score += 9

        elif p.star_type[0]  == 'M':
            score += 8 if number <= 4 else 10

        else:
            score += 0
         
        return score   
    


#   Se genera una última clase que almacenará los datos dentro de una Tupla,    
#   estos mismos serán utilizados para db_sqlite.py bajo esta clase.

class DataValue:
    def __init__(self, planet, habitability_score):
        
        self.data = (
            planet.name,
            planet.j_mass,
            planet.e_mass,
            planet.j_radius,
            planet.e_radius,
            
            planet.rotation_period,
            planet.axis_stability,
            
            planet.sizeClass,
            planet.composition,
            planet.radiation_levels,
           
            planet.blackbody1,
            planet.blackbody2,
            planet.blackbody3,
            planet.media,

            planet.orbital_period,
            planet.orbital_period_observed,
            planet.orbital_period_estimated,
            planet.eccentricity,
            planet.orbital_type,

            planet.star_name,
            planet.star_type,
            planet.distance_from_star,
            planet.star_radiation_atmospheric_b,
            planet.stellar_activity,
            planet.in_habitable_zone,

            planet.t_sync,
            planet.base,
            planet.probabilidad_anclaje,

            planet.Q_min,
            planet.q_min_interpretation,
            planet.Q_max,
            planet.q_max_interpretation,

            planet.oxygen,
            planet.gravity,
            planet.magnetic_field,
            planet.atmospheric_pressure,
            planet.volcanic_activity,
            planet.has_water,

            habitability_score.calculate_score()
        )

    def get_data(self):
        return self.data

