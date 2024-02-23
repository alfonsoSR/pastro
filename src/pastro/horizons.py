import numpy as np
from urllib.request import urlopen
from dataclasses import dataclass
from .types import CartesianState


@dataclass
class HorizonsRequest:
    body_id: str
    obj_data: str
    ephem_type: str
    center: str
    ref_plane: str
    start_time: str
    stop_time: str
    step_size: str
    vec_table: str
    vec_labels: str
    csv_format: str
    vec_delta_t: str


def get_ephemeris(config: HorizonsRequest):

    url = ('https://ssd.jpl.nasa.gov/api/horizons.api?format=text&COMMAND=%27'
           f"{config.body_id}%27&OBJ_DATA=%27{config.obj_data}"
           "%27&EPHEM_TYPE=%27"
           f"{config.ephem_type}%27&CENTER=%27{config.center}%27&REF_PLANE=%27"
           f"{config.ref_plane}%27&START_TIME=%27{config.start_time}"
           f"%27&STOP_TIME=%27{config.stop_time}%27&STEP_SIZE=%27"
           f"{config.step_size}%27&VEC_TABLE=%27{config.vec_table}"
           f"%27&VEC_LABELS=%27{config.vec_labels}%27&CSV_FORMAT=%27"
           f"{config.csv_format}%27&VEC_DELTA_T=%27{config.vec_delta_t}%27")

    # Request data to Horizons API
    response = urlopen(url).read().decode('utf-8')

    # Split response into lines
    lines_mod = []
    start = False
    for line in response.splitlines():

        if line == "$$SOE":
            start = True
            continue

        if not start:
            continue
        if line == "$$EOE":
            break
        else:
            lines_mod.append(line.split(","))

    # Transform list into numpy array
    arr = np.array(lines_mod, dtype=object)

    # Delete unnecessary columns
    arr = np.delete(arr, [1, 6, 7, 8, 9], 1).astype(float).T

    # Save cartesian coordinates
    s = CartesianState(arr[2, 0], arr[3, 0], arr[4, 0], 0, 0, 0) * 1e3
    assert (isinstance(s, CartesianState))

    jds = arr[0]
    dts = arr[1]

    return jds, dts, s
