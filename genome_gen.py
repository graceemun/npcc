import struct

file_path = "//Users//kalie_chang//Desktop//genome.bin"

for i in range(7):
    for j in range(7):
        x = i
        y = j
        
for i in range(1024):
    genome = random.randint(0, 1)

data = struct.pack('<3Q', x, y, genome)

# Example binary data as a packed struct
with open(file_path, "wb") as file: #wb = write binary
    file.write(data)

