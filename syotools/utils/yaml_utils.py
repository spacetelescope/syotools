import yaml
import numpy as np
import astropy.units as u

def simplify_data(data):
    print(type(data),data)
    if isinstance(data, str):
        data = str(data)
    elif isinstance(data, u.Quantity):
        data = ["!Astropy", simplify_data(data.value),data.unit]
    elif isinstance(data, (int, np.int32, np.int64)):
        data = int(data)
    elif isinstance(data, (float, np.float32, np.float64)):
        data = float(data)
    elif isinstance(data, (list, tuple, np.ndarray)):
        data = list(data)
        for idx, item in enumerate(data):
            data[idx] = simplify_data(item)
    elif isinstance(data, (dict, set)):
        for item in data:
            data[item] = simplify_data(data[item])
    else:
        print("uncertain")

    return(data)

def complexify_data(data):
    if isinstance(data, list):
        print(data, type(data))
        # for item in data:
        #     print(type(item))
        if data[0] == "!Astropy":
            #print("Astropy", data, type[data[2]])
            data = complexify_data(data[1]) * u.Unit(data[2])
        else:
            for idx,x in enumerate(data):
                data[idx] = complexify_data(x)
    elif isinstance(data, dict):
        for key in data:
            data[key] = complexify_data(data[key])
    return data 


def write_yaml(data, outfilename):
    data=simplify_data(data)
    with open(outfilename, "w") as outfile:
        outfile.write(yaml.dump(data, Dumper=yaml.Dumper))

def read_yaml(infilename):
    with open(infilename, "r") as infile:
        data = yaml.load(infile, Loader=yaml.Loader)
    data = complexify_data(data)

    return data
