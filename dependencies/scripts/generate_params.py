from cmath import atan
import os
import yaml
from pathlib import Path

def generate_folder_names(folder_path):
    folder_names = {}
    index = 0
    for folder_name in os.listdir(folder_path):
        if os.path.isdir(os.path.join(folder_path, folder_name)):
            folder_names[index] = folder_name
            index += 1
    return folder_names

# Wygeneruj słownik z nazwami folderów
generated_params = generate_folder_names(Path('DATABASES_HHSUITE'))
print(generated_params)

# Zapisz słownik do pliku YAML
with open("generated_params.yaml", "w") as params_file:
    yaml.dump(generated_params, params_file)

print("Generated params saved to generated_params.yaml")
