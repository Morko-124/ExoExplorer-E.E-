import sqlite3
from exoplanet import DataValue, Planet, Habitability_Score
import requests
from bs4 import BeautifulSoup

try:
    connection = sqlite3.connect('Data/ExoExplore_DB.db')
    print('Connection established.')
    cursor = connection.cursor()

    # Corrección en la definición de la tabla
    cursor.execute('''
    CREATE TABLE IF NOT EXISTS ExoPlanets (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                Name TEXT,
                Jupiter_Mass REAL, 
                Earths_Mass REAL, 
                Jupiters_Radius REAL, 
                Earths_Radius REAL,
                Rotation_Period REAL, 
                Axial_Stability REAL,
                Size_Class TEXT, 
                Composition TEXT, 
                Radiation_Levels REAL, 
                Blackbody_T_01_K REAL, 
                Blackbody_T_03_K REAL,  
                Blackbody_T_07_K REAL,		
                Orbital_Period_Yrs REAL, 
                Orbital_Period_Observed_Days REAL, 
                Orbital_Period_Estimated_Days REAL, 
                Eccentricity REAL, 
                Orbital_Type TEXT,
                Star_Name TEXT,
                Star_Type TEXT, 
                Distance_From_Star REAL, 
                Star_Radiation_Atmospheric_Boundary REAL, 
                Stellar_Activity TEXT, 
                Located_in_Habitable_Zone BOOLEAN,  
                Habitable_Zone_Description TEXT,    
                Sync_Time REAL, 
                Star_Influence TEXT, 
                Tidal_Locking_Probability TEXT, 
                Minimum_Quality_Factor_Q REAL, 
                Interpretation_of_Minimum_Quality_Factor_Q TEXT, 
                Maximum_Quality_Factor_Q REAL, 
                Interpretation_of_Maximum_Quality_Factor_Q TEXT, 
                Uncertainty_Factor_Oxygen_Probability TEXT, 
                Gravity REAL, 
                Uncertainty_Factor_Magnetic_Field TEXT, 
                Uncertainty_Factor_Atmospheric_Pressure REAL, 
                Volcanic_Activity REAL, 
                Uncertainty_Factor_Water_Presence_Probability REAL, 
                Habitability_Score REAL
    )''')

    # Función para insertar datos en la base de datos
    def insert_planet_data(data):
        cursor.execute('''
        INSERT INTO ExoPlanets (
            Name, Jupiter_Mass, Earths_Mass, Jupiters_Radius, Earths_Radius, Rotation_Period, Axial_Stability, 
            Size_Class, Composition, Radiation_Levels, Blackbody_T_01_K, Blackbody_T_03_K, Blackbody_T_07_K, 
            Orbital_Period_Yrs, Orbital_Period_Observed_Days, Orbital_Period_Estimated_Days, Eccentricity, 
            Orbital_Type, Star_Name, Star_Type, Distance_From_Star, Star_Radiation_Atmospheric_Boundary, 
            Stellar_Activity, Located_in_Habitable_Zone, Sync_Time, Star_Influence, Tidal_Locking_Probability, 
            Minimum_Quality_Factor_Q, Interpretation_of_Minimum_Quality_Factor_Q, Maximum_Quality_Factor_Q, 
            Interpretation_of_Maximum_Quality_Factor_Q, Uncertainty_Factor_Oxygen_Probability, Gravity, 
            Uncertainty_Factor_Magnetic_Field, Uncertainty_Factor_Atmospheric_Pressure, Volcanic_Activity, 
            Uncertainty_Factor_Water_Presence_Probability, Habitability_Score
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', data)
        connection.commit()
        print('Datos insertados exitosamente.')

    # Prueba con un planeta
    url = 'https://www.exoplanetkyoto.org/exohtml/AB_Aur_b.html'
    response = requests.get(url)

    if response.status_code == 200:
        soup = BeautifulSoup(response.content, 'html.parser')
        planet = Planet(soup)
        habitability_score = Habitability_Score(planet)
        data_value = DataValue(planet, habitability_score)
        
        # Obtener datos y verificar el tamaño
        data = data_value.get_data()
        print(f"Número de elementos en los datos: {len(data)}")
        print(f"Datos: {data}")

        # Insertar los datos
        insert_planet_data(data)
    else:
        print(f"Error al acceder a la página: {response.status_code}")

except Exception as e:
    print(f"Error: {e}")
finally:
    connection.close()
