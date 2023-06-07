

import yaml

# Python dictionary with the data to write to the YAML file

data = {
    'Lagoon': {
        'coords': '[6, -1.2]',
        'region_type': 'rectangle',
        'frame': 'galactic',
        'region_params': {
            'width': 1.2,
            'height': 1.0,
            'pa': 0
                          },
        'priority': 9,
        'observatory': 'LCO',
        'telescope': 'LVM-160',
        'max_airmass': 1.75,
        'min_shadowheight': 1000.0,
        'exptime': 900,
        'n_exposures': 1,
        'min_exposures': 1,
        'min_moon_dist': 60,
        'max_lunation': 1.0,
        'overhead': 1.1,
        'tiling_strategy': 'lowest_airmass',
        'tile_overlap': 0.02416
                }
        }
        
#file path where I write the YAML file
file_path = '/home/amrita/Desktop/PhD U Chile/U chile/Sem 2 Project/Lagoon.yaml'

# Open the file in write mode
with open(file_path, 'w') as file:
    # Use the yaml.dump() function to convert the dictionary to YAML and write it to the file
    yaml.dump(data, file)

