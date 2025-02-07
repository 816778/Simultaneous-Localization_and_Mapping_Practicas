import matplotlib.pyplot as plt

def leer_archivo_txt(nombre_archivo):
    datos = {}

    with open(nombre_archivo, 'r', encoding='utf-8') as archivo:
        next(archivo)  # Ignorar la primera línea
        for linea in archivo:
            linea = linea.strip()
            if not linea:
                continue  # Ignorar líneas vacías
            
            partes = linea.split(';')
            clave = int(partes[0])  # Convertimos el primer valor en entero (N_maps)
            valores = [float(x) for x in partes[1:] if x]  # Convertimos los valores en flotantes
            
            datos[clave] = valores  # Guardamos en el diccionario

    return datos

def calcular_medias(datos):
    medias = {}
    for clave, valores in datos.items():
        medias[clave] = sum(valores) / len(valores)  # Calculamos la media de los valores
    return medias


nombre_archivo = "result/elapse_time.txt"  # Cambia esto por el nombre real de tu archivo
datos_leidos = leer_archivo_txt(nombre_archivo)

medias = calcular_medias(datos_leidos)

for key, media in medias.items():
    print(f"N_maps={key}: Media={media:.4f}")

plt.figure(figsize=(8, 5))
plt.plot(medias.keys(), medias.values(), marker='o', linestyle='-', label="Tiempo medio")

# Personalizar la gráfica
plt.xlabel("N_maps")
plt.ylabel("Elapse Time (s)")
plt.title("Elapse Time")
plt.xticks(list(medias.keys()))  # Asegurar que los valores del eje X son correctos
plt.grid(True)
plt.legend()

# Mostrar la gráfica
plt.show()