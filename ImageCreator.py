from PIL import Image
import random

maxvalue = 255

# Getting informations from the prompt
width = int(input("image width: "))
height = int(input("image height: "))

# Crea un oggetto immagine vuoto con una modalità "L" (scala di grigi)
immagine = Image.new("L", (width, height))

for y in range(height):
    for x in range(width):
        pixel = random.randint(0, 255)
        immagine.putpixel((x, y), pixel)

# Salva l'immagine come file PGM
filename = input("PGM file name: ")

with open(filename + ".pgm", "w") as file:
    file.write("P2\n")  # Tipo di immagine
    file.write(f"{width} {height}\n")  # Dimensioni
    file.write(f"{maxvalue}\n")  # Valore massimo dei pixel
immagine.save(filename + "_" + str(width) + "x" + str(height) + ".pgm")

# Chiudi l'oggetto immagine (opzionale ma buona pratica)
immagine.close()