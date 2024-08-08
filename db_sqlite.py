import sqlite3

try:
    mi_connection=sqlite3.connect('Data/ExoExplore_DB.db')
    print('Conexion Establecida.')

except Exception as ex:
    print(ex)

